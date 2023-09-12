/*
  This file is part of the BROOMStyx project.

  BROOMStyx is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  BROOMStyx is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with BROOMStyx.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the AUTHORS file
  for the list of copyright holders.
*/

#include "DomainManager.hpp"
#include <chrono>
#include <stdexcept>
#include <omp.h>

#include "AnalysisModel.hpp"
#include "Cell.hpp"
#include "DofManager.hpp"
#include "DomainManager.hpp"
#include "MaterialManager.hpp"
#include "Node.hpp"
#include "NumericsManager.hpp"
#include "SolutionManager.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Numerics/Numerics.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

// Constructor
DomainManager::DomainManager()
{
    _fieldsPerNode = -1;
    _name = "DomainManager";
}

// Destructor
DomainManager::~DomainManager()
{
    for ( auto curNode : _nodeList )
    {
        analysisModel().dofManager().destroyNodalDofsAt( curNode );
        delete curNode;
    }
    
    for ( int dim = 0; dim < 4; dim++ )
    {
        for ( auto curCell : _cellList[ dim ] )
        {
            analysisModel().dofManager().destroyCellDofsAt( curCell );
            delete curCell;
        }
    }
}

// Public methods
// ----------------------------------------------------------------------------
void DomainManager::createPhysicalEntity( int dim, int number, std::string label )
{
    PhysicalEntity newPhysEnt;
    
    newPhysEnt.dimension = dim;
    newPhysEnt.entityNumber = number;
    newPhysEnt.name = label;
    
    _physEnt.push_back( newPhysEnt );
}
// ----------------------------------------------------------------------------
DomainManager::PhysicalEntity DomainManager::givePhysicalEntity( int n )
{
    return _physEnt[ n ];
}
// ----------------------------------------------------------------------------
std::vector<Material*> DomainManager::giveMaterialSetForDomain( int label, int stage )
{
    std::string name = this->givePhysicalEntityNameFor( label );
    try
    {
        auto materialSet = _materialSet[ stage ].at( name );
        return materialSet;
    }
    catch( const std::exception& e )
    {
        throw std::runtime_error( "No Material set defined for '" + name + "' at stage "
            + std::to_string( stage ) + "!\n" );
    }
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfPhysicalNames()
{
    return _physEnt.size();
}
// ----------------------------------------------------------------------------
Numerics* DomainManager::giveNumericsForDomain( int label, int stage )
{
    std::string name = this->givePhysicalEntityNameFor( label );
    Numerics* numerics;

    try
    {
        numerics = _numerics[ stage ].at( name );
    }
    catch( const std::exception& e )
    {
        numerics = nullptr;
    }
    
    return numerics;
}
// ----------------------------------------------------------------------------
std::string DomainManager::givePhysicalEntityNameFor( int physEntNum )
{
    std::string physEntName;
    bool success = false;
    
    for ( int i = 0; i < (int)_physEnt.size(); i++ )
    {
        if ( physEntNum == _physEnt[ i ].entityNumber )
        {
            physEntName = _physEnt[ i ].name;
            success = true;
            break;
        }
    }
    
    if ( !success )
        throw std::runtime_error( "Failed to find name corresponding to physical entity number '" 
                + std::to_string( physEntNum ) + "'!\nSource: " + _name );
    
    return physEntName;
}
// ----------------------------------------------------------------------------
int DomainManager::givePhysicalEntityNumberFor( std::string name )
{
    int physEntNumber;
    bool success = false;
    
    for ( int i = 0; i < (int)_physEnt.size(); i++ )
    {
        if ( name == _physEnt[ i ].name )
        {
            physEntNumber = _physEnt[ i ].entityNumber;
            success = true;
            break;
        }
    }
    
    if ( !success )
        throw std::runtime_error( "Failed to find physical entity number corresponding to '" + name + "'!\nSource: " + _name );
    
    return physEntNumber;
}
// ----------------------------------------------------------------------------
void DomainManager::readDomainAssignmentsFrom( FILE* fp )
{
    // Get number of stages from SolutionManager and initialize maps
    // for numerics and material sets
    int nStages = analysisModel().solutionManager().giveNumberOfStages();

    // Note: Stage numbers start from 1 (and not 0) so we allocate an extra element
    _numerics.resize( nStages + 1 );
    _materialSet.resize( nStages + 1 );

    int nAssign = getIntegerInputFrom( fp, "Failed to read number of domain assignments from input file!", _name );
    
    for ( int i = 0; i < nAssign; i++ )
    {
        // Read stage number
        verifyKeyword( fp, "Stage", _name );
        int stage = getIntegerInputFrom( fp, "Failed reading stage number from input file!", _name );

        // Read domain label
        std::string domainLabel = getStringInputFrom( fp, "Failed reading domain label from input file!", _name );
                
        // Read numerics
        verifyKeyword( fp, "Numerics", _name );
        int numericsLabel = getIntegerInputFrom( fp, "Failed reading numerics label from input file!", _name );
        
        Numerics* numericsPtr = analysisModel().numericsManager().giveNumerics( numericsLabel );
        
        // Create map entry
        std::pair< std::map<std::string, Numerics*>::iterator, bool> entry;
        entry = _numerics[ stage ].insert( std::pair<std::string, Numerics*>( domainLabel, numericsPtr ) );
        if ( !entry.second )
            throw std::runtime_error( "Multiple declaration of numerics for label '" + domainLabel
                + "' detected in input file!\nSource: " + _name );
        
        // Read material set
        int nMat = numericsPtr->requiredNumberOfMaterials();
        if ( nMat > 0 )
        {
            verifyKeyword( fp, "MaterialSet", _name );
        
            std::vector<Material*> matSet;
            matSet.assign( nMat, nullptr );
            for ( int j = 0; j < nMat; j++ )
            {
                int matLabel = getIntegerInputFrom( fp, "Failed to read material label from input file.", _name );
                matSet[ j ] = analysisModel().materialManager().giveMaterial( matLabel );
            }

            // Create map entry
            std::pair<std::map<std::string, std::vector<Material*> >::iterator,bool> tmp;
            tmp = _materialSet[ stage ].insert( std::pair<std::string, std::vector<Material*> >( domainLabel, matSet ) );
            if ( !tmp.second )
                throw std::runtime_error( "Multiple declaration of material sets for label '" + domainLabel
                    + "' detected in input file!\nSource: " + _name );
        }
    }
}
// ----------------------------------------------------------------------------
void DomainManager::setNumberOfStagesTo( int nStage )
{
    _nStage = nStage;

    // Initialize maps for numerics and material sets
    // Note: Stage numbers start from 1 (and not 0) so we allocate an extra element
    _numerics.resize( nStage + 1 );
    _materialSet.resize( nStage + 1 );
}

// Methods involving node access
// ----------------------------------------------------------------------------
void DomainManager::countNodes()
{
    // Determine number of active nodes
    int nActiveNodes = 0;
    for ( auto curNode = _nodeList.begin(); curNode != _nodeList.end(); curNode++ )
        if ( (*curNode)->_isActive )
            (*curNode)->_id = nActiveNodes++;
    
    // Rewrite active node addresses into one contiguous array
    _node.assign( nActiveNodes, nullptr );
    int curCount = 0;
    for ( auto curNode = _nodeList.begin(); curNode != _nodeList.end(); curNode++ )
        if ( (*curNode)->_isActive )
            if ( (*curNode)->_isActive )
                _node[ curCount++ ] = *curNode;
}
// ----------------------------------------------------------------------------
std::set<Cell*> DomainManager::giveCellsAttachedTo( Node* node, int dim )
{
    return node->_attachedCell[ dim ];
}
// ----------------------------------------------------------------------------
RealVector DomainManager::giveCoordinatesOf( Node* node )
{
    return node->_coordinates;
}
// ----------------------------------------------------------------------------
Dof* DomainManager::giveNodalDof( int dofNum, Node* node )
{
    return node->_dof[ dofNum ];
}
// ----------------------------------------------------------------------------
double DomainManager::giveFieldValueAt( Node* node, int fieldNum )
{
    if ( fieldNum == 0 )
        return 0.0;
    else
        return node->_fieldVal( fieldNum - 1 );
}
// ----------------------------------------------------------------------------
int DomainManager::giveIdOf( Node* node )
{
    return node->_id;
}
// ----------------------------------------------------------------------------
Node* DomainManager::giveNode( int nodeNum ) 
{
    return _node[ nodeNum ];
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfNodes()
{
    return _node.size();
}
// ----------------------------------------------------------------------------
void DomainManager::makeNewNodeAt( RealVector& location )
{
    if ( _fieldsPerNode == -1 )
        throw std::runtime_error( "Cannot create new node due to undefined number of fields per node!\nSource: " + _name );

    Node* newNode = new Node();
    
    if ( location.dim() == 2 )
    {
        newNode->_coordinates.init( 3 );
        newNode->_coordinates( 0 ) = location( 0 );
        newNode->_coordinates( 1 ) = location( 1 );
        newNode->_coordinates( 2 ) = 0.0;
    }
    else if ( location.dim() == 3 )
        newNode->_coordinates = location;
    else
        throw std::runtime_error( "Invalid size of vector input for nodal coordinates!\nSource: " + _name );
    
    // Instantiate degrees of freedom for new node
    analysisModel().dofManager().createNodalDofsAt( newNode );
    
    // Initialize nodal fields
    if ( _fieldsPerNode > 0 )
        newNode->_fieldVal.init( _fieldsPerNode );
    
    // Append node to list of nodes
    _nodeList.push_back( newNode );
}
// ----------------------------------------------------------------------------
void DomainManager::performNodalPostProcessing()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf( "    %-40s", "Performing nodal post-processing ..." );
    std::fflush( stdout );
    tic = std::chrono::high_resolution_clock::now();
    
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < (int)_node.size(); i++ )
    {
        // Get field values corresponding to DOF variables
        analysisModel().dofManager().writeConvergedDofValuesTo( _node[i] );
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf( "done (time = %f sec.)\n", tictoc.count() );
}
// ----------------------------------------------------------------------------
void DomainManager::readNumberOfFieldsPerNodeFrom( FILE* fp )
{
    _fieldsPerNode = getIntegerInputFrom( fp, "\nFailed to read number of fields per node from input file", _name );
}
// ----------------------------------------------------------------------------
void DomainManager::setCoordinatesOf( Node* targetNode, const RealVector& coor )
{
    // Need to do the update component wise in order to be able to take advantage of atomic operations
#ifdef _OPENMP
#pragma omp atomic write
#endif
    targetNode->_coordinates(0) = coor(0);
#ifdef _OPENMP
#pragma omp atomic write
#endif
    targetNode->_coordinates(1) = coor(1);
#ifdef _OPENMP
#pragma omp atomic write
#endif
    targetNode->_coordinates(2) = coor(2);

}
// ----------------------------------------------------------------------------
void DomainManager::setFieldValueAt( Node* targetNode, int fieldNum, double val )
{
    targetNode->_fieldVal( fieldNum - 1 ) = val;
}

// Methods for cell access
// ----------------------------------------------------------------------------
void DomainManager::countCells()
{
    int totCount = 0;
    for ( int dim = 0; dim < 4; dim++ )
    {
        _cell[ dim ].assign( _cellList[ dim ].size(), nullptr );
        int vecCount = 0;
        for ( auto it = _cellList[ dim ].begin(); it != _cellList[ dim ].end(); it++ )
        {
            _cell[ dim ][ vecCount++ ] = *it;
            (*it)->_id = totCount++;
        }
    }
}
// ----------------------------------------------------------------------------
void DomainManager::finalizeCellDataAt( const TimeData& time, int stage )
{
    for ( int dim = 0; dim < 4; dim++ )
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < (int)_cell[ dim ].size(); i++ )
        {
            Numerics* numerics = this->giveNumericsFor( _cell[ dim ][ i ], stage );
            if ( numerics )
                numerics->finalizeDataAt( _cell[ dim ][ i ], time );
        }
    }
}
// ----------------------------------------------------------------------------
void DomainManager::findCellAttachments()
{
     std::chrono::time_point<std::chrono::system_clock> tic, toc;
     std::chrono::duration<double> tictoc;

     std::printf( "  %-40s", "Finding cell attachments ..." );
     std::fflush( stdout );
     tic = std::chrono::high_resolution_clock::now();

     for ( int dim = 0; dim < 4; dim ++ )
     {
#ifdef _OPENMP
#pragma omp parallel for
#endif
         for ( int i = 0; i < (int) _cell[ dim ].size(); i++ )
             this->findCellsAttachedTo( _cell[ dim ][ i ] );
     }

     toc = std::chrono::high_resolution_clock::now();
     tictoc = toc - tic;
     std::printf( "done (time = %f sec.)\n", tictoc.count() );
}
// ----------------------------------------------------------------------------
void DomainManager::findCellsAttachedTo( Cell* targetCell )
{
    std::vector<Node*> cellNode = targetCell->_node;
    int nCellNodes = cellNode.size();

    for ( int curDim = 0; curDim < 4; curDim++ )
    {
        for ( int curNode = 0; curNode < nCellNodes; curNode++ )
        {
            std::set<Cell*> candidate = this->giveCellsAttachedTo( cellNode[ curNode ], curDim );

            if ( targetCell->_dim < curDim )
            {
                // All nodes of targetCell must also be nodes of candidate attached cell
                bool candidateIsAttached = true;
                for ( auto it = candidate.begin(); it != candidate.end(); ++it )
                {
                    std::vector<Node*> candNode = (*it)->_node;
                    for ( int j = 0; j < nCellNodes; j++ )
                    {
                        bool nodeIsShared = false;
                        for ( int k = 0; k < (int)candNode.size(); k++ )
                            if ( cellNode[ j ] == candNode[ k ] )
                                nodeIsShared = true;
                        
                        if ( !nodeIsShared )
                            candidateIsAttached = false;
                    }
                    if ( candidateIsAttached )
                    {
                        targetCell->_attachedCell[ curDim ].insert( *it );
                        (*it)->_attachedCell[ targetCell->_dim ].insert( targetCell );
                    }
                }
            }
            else if ( targetCell->_dim > curDim )
            {
                // All nodes of candidate attached cell must also be nodes of targetCell
                bool candidateIsAttached = true;
                for ( auto it = candidate.begin(); it != candidate.end(); ++it )
                {
                    std::vector<Node*> candNode = (*it)->_node;
                    for ( int k = 0; k < (int)candNode.size(); k++ )
                    {
                        bool nodeIsShared = false;
                        for ( int j = 0; j < nCellNodes; j++ )
                            if ( cellNode[ j ] == candNode[ k ] )
                                nodeIsShared = true;
                        
                        if ( !nodeIsShared )
                            candidateIsAttached = false;
                    }
                    if ( candidateIsAttached )
                    {
                        targetCell->_attachedCell[ curDim ].insert( *it );
                        (*it)->_attachedCell[ targetCell->_dim ].insert( targetCell );
                    }
                }
            }
            else
            {
                // ( targetCell->_dim == curDim )
                // targetCell and candidate attached cell must share a face
                int nFaces = analysisModel().meshReader().giveNumberOfFacesForElementType( targetCell->_elType );
                targetCell->_neighbor.assign( nFaces, nullptr );

                for ( int curFace = 0; curFace < nFaces; curFace++ )
                {
                    std::vector<int> faceNodes = analysisModel().meshReader().giveFaceNodeNumbersForElementType( targetCell->_elType, curFace );
                    for ( auto it = candidate.begin(); it != candidate.end(); it++ )
                    {
                        std::vector<Node*> candNode = this->giveNodesOf( *it );

                        bool candidateIsAttached = true;
                        for ( int k = 0; k < (int)faceNodes.size(); k++ )
                        {
                            bool sharesFaceNode = false;
                            for ( int m = 0; m < (int)candNode.size(); m++ )
                                if ( candNode[ m ] == targetCell->_node[ faceNodes[ k ] ] )
                                    sharesFaceNode = true;

                            if ( !sharesFaceNode )
                                candidateIsAttached = false;
                        }

                        if ( candidateIsAttached && (*it) != targetCell )
                        {
                            targetCell->_attachedCell[ curDim ].insert( *it );
                            (*it)->_attachedCell[ curDim ].insert( targetCell );

                            targetCell->_neighbor[ curFace ] = (*it);
                        }
                    }
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
Cell* DomainManager::giveCell( int num, int dim )
{
    return _cell[ dim ][ num ];
}
// ----------------------------------------------------------------------------
Dof* DomainManager::giveCellDof( int dofNum, Cell *targetCell )
{
    return targetCell->_dof[ dofNum ];
}
// ----------------------------------------------------------------------------
int DomainManager::giveElementTypeOf( Cell* targetCell )
{
    return targetCell->_elType;
}
// ----------------------------------------------------------------------------
int DomainManager::giveIdOf( Cell *targetCell )
{
    return targetCell->_id;
}
// ----------------------------------------------------------------------------
int DomainManager::giveLabelOf( Cell *targetCell )
{
    return targetCell->_label;
}
// ----------------------------------------------------------------------------
std::vector<Cell*> DomainManager::giveNeighborsOf( Cell* targetCell )
{
    return targetCell->_neighbor;
}
// ----------------------------------------------------------------------------
std::vector<Node*> DomainManager::giveNodesOf( Cell *targetCell )
{
    return targetCell->_node;
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfCellsWithDimension( int dim )
{
    return (int)_cell[ dim ].size();
}
// ----------------------------------------------------------------------------
int DomainManager::giveNumberOfNodesOf( Cell *targetCell )
{
    return targetCell->_node.size();
}
// ----------------------------------------------------------------------------
Numerics* DomainManager::giveNumericsFor( Cell* targetCell, int stage )
{
    return this->giveNumericsForDomain( targetCell->_label, stage );
}
// ----------------------------------------------------------------------------
void DomainManager::initializeMaterialsAtCells()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf( "  %-40s", "Initializing material data at cells ..." );
    std::fflush( stdout );
    tic = std::chrono::high_resolution_clock::now();

    for ( int dim = 0; dim < 4; dim++ )
    {
#ifdef _OPENMP    
#pragma omp parallel for
#endif
        for ( int i = 0; i < (int)_cell[ dim ].size(); i++ )
        {
            for ( int curStage = 1; curStage <= _nStage; curStage++ )
            {
                Numerics* numerics = this->giveNumericsFor( _cell[ dim ][ i ], curStage );
                if ( numerics )
                    numerics->initializeMaterialsAt( _cell[ dim ][ i ] );
            }
        }
    }

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf( "done (time = %f sec.)\n", tictoc.count() );
}
// ----------------------------------------------------------------------------
void DomainManager::initializeNumericsAtCells()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;
    
    std::printf( "  %-40s", "Initializing numerics at cells ..." );
    std::fflush( stdout );
    tic = std::chrono::high_resolution_clock::now();

    for ( int dim = 0; dim < 4; dim++ )
    {
#ifdef _OPENMP    
#pragma omp parallel for
#endif
        for ( int i = 0; i < (int)_cell[ dim ].size(); i++ )
        {
            for ( int curStage = 1; curStage <= _nStage; curStage++ )
            {
                Numerics* numerics = this->giveNumericsFor( _cell[ dim ][ i ], curStage );
                if ( numerics )
                    numerics->initializeNumericsAt( _cell[ dim ][ i ] );
            }
        }
    }
    
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf( "done (time = %f sec.)\n", tictoc.count() );
}
// ----------------------------------------------------------------------------
Cell* DomainManager::makeNewCell( int elType, int cellLabel, int dim )
{
    // Instantiate new cells
    Cell* newCell = new Cell( elType, cellLabel, dim );

    // Add new cell object to relevant list
    _cellList[ dim ].push_back( newCell );

    // Create DOF objects
    analysisModel().dofManager().createCellDofsAt( newCell );

    return newCell;
}
// ----------------------------------------------------------------------------
void DomainManager::readNumberOfFieldsPerCellFrom( FILE* fp )
{
    _fieldsPerCell = getIntegerInputFrom( fp, "\nFailed to read number of fields per cell in input file!", "DomainManager" );
}
// ----------------------------------------------------------------------------
void DomainManager::removeAllCellConstraints()
{
    for ( int dim = 0; dim < 4; dim++ )
        for ( int i = 0; i < (int)_cell[ dim ].size(); i++ )
            for ( int stage = 1; stage <= _nStage; stage++ )
            {
                Numerics *numerics = analysisModel().domainManager().giveNumericsForDomain( _cell[ dim ][ i ]->_label, stage );
                if ( numerics )
                    numerics->removeConstraintsOn( _cell[ dim ][ i ] );
            }
}
// ----------------------------------------------------------------------------
void DomainManager::reorderNodesOf( Cell* targetCell, std::vector<int>& reordering )
{
    int nNodes = (int)reordering.size();
    std::vector<Node*> reorderedNode( nNodes, nullptr );
    if ( nNodes != (int)targetCell->_node.size() )
        throw std::runtime_error( "ERROR: Number of original and reordered cell nodes don't match!\nSource: " + _name );
    else
        for ( int i = 0; i < nNodes; i++ )
            reorderedNode[ i ] = targetCell->_node[ reordering[ i ] ];

    targetCell->_node = reorderedNode; 

    // Note: Implementation is incomplete in the sense that boundary associations, etc. also have to be
    // updated (which is currently not done).
}
// ----------------------------------------------------------------------------
void DomainManager::reportDetailedStatus()
{
    this->reportStatus();
    
    std::printf( "\n" );
    for ( int i = 0; i < (int)_node.size(); i++ )
    {
        RealVector coor = _node[ i ]->_coordinates;
        std::printf( "    node %d: x = %e, y = %e, z = %e, dofs = ", i, coor( 0 ), coor( 1 ), coor( 2 ) );
        for ( int j = 0; j < (int)_node[ i ]->_dof.size(); j++ )
            std::printf("%d ", analysisModel().dofManager().giveEquationNumberAt( _node[ i ]->_dof[ j ] ) );
        std::printf( "\n" );
    }

    for ( int dim = 0; dim < 4; dim++ )
    {
        for ( int i = 0; i < (int)_cell[ dim ].size(); i++ )
        {
            Cell* curCell = _cell[dim][i];
            std::printf( "   cell %d: dim = %d, nodes = ", curCell->_id, dim );
            std::vector < Node * > node = curCell->_node;
            for ( int j = 0; j < (int)node.size(); j++ )
                std::printf( "%d ", node[j]->_id );

            int nDofs = curCell->_dof.size();
            if ( nDofs > 0 )
            {
                std::printf( ", dofs = " );
                for ( int j = 0; j < nDofs; j++ )
                    std::printf( "%d ", analysisModel().dofManager().giveEquationNumberAt( curCell->_dof[ j ] ) );
            }

            int nNeighbors = curCell->_neighbor.size();
            if ( nNeighbors > 0 )
            {
                std::printf( ", neighbors = " );
                for ( int j = 0; j < (int)curCell->_neighbor.size(); j++ )
                    if ( curCell->_neighbor[ j ] )
                        std::printf( "%d ", curCell->_neighbor[ j ]->_id);
                    else
                        std::printf( "null " );
            }

            std::printf("\n");
        }
    }
}
// ----------------------------------------------------------------------------
void DomainManager::reportStatus()
{
    std::printf( "    Nodes          = %ld\n", _node.size() );
    std::printf( "    3-D cells      = %ld\n", _cell[ 3 ].size() );
    std::printf( "    2-D cells      = %ld\n", _cell[ 2 ].size() );
    std::printf( "    1-D cells      = %ld\n", _cell[ 1 ].size() );
    std::printf( "    0-D cells      = %ld\n", _cell[ 0 ].size() );
}
// ----------------------------------------------------------------------------
void DomainManager::setElementTypeOf( Cell* targetCell, int elemType )
{
    targetCell->_elType = elemType;
}
// ----------------------------------------------------------------------------
void DomainManager::setNodesOf( Cell *targetCell, std::vector<int>& cellNodes )
{
    // Important: this method must not be called from within a loop that is
    // parallelized since std::set is not thread safe!
    
    std::vector<Node*> node( cellNodes.size(), nullptr );
    
    for ( int i = 0; i < (int)cellNodes.size(); i++ )
    {
        node[ i ] = _node[ cellNodes[ i ] ];
        node[ i ]->_attachedCell[ targetCell->_dim ].insert( targetCell );
    }
    
    targetCell->_node = node;
}