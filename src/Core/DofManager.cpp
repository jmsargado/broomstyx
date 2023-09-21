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

#include "DofManager.hpp"
#include <omp.h>
#include <stdexcept>

#include "AnalysisModel.hpp"
#include "Dof.hpp"
#include "Cell.hpp"
#include "Node.hpp"
#include "DomainManager.hpp"
#include "SolutionManager.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

DofManager::DofManager()
{
    _name = "DofManager";
    _multiFreedomConstraint.clear();
    _cellDofInfo[ 0 ].clear();
    _cellDofInfo[ 1 ].clear();
    _cellDofInfo[ 2 ].clear();
    _cellDofInfo[ 3 ].clear();
    _nodalDofInfo.clear();
    _numericsDof.assign( 0, nullptr );
}
// ----------------------------------------------------------------------------
DofManager::~DofManager() {}
// ----------------------------------------------------------------------------
// Public methods
//void DofManager::addToSecondaryVariableAt( Dof* targetDof, double val )
//{
//#ifdef _OPENMP
//#pragma omp atomic
//#endif
//    targetDof->_secVar += val;
//}
// ----------------------------------------------------------------------------
void DofManager::createCellDofsAt( Cell* targetCell )
{
    int dim = targetCell->_dim;
    int nDofs = _cellDofInfo[ dim ].size();
    targetCell->_dof.assign( nDofs, nullptr );
    for ( int i = 0; i < nDofs; i++ )
        targetCell->_dof[ i ] = new Dof( _cellDofInfo[ dim ][ i ].group );
}
//// ----------------------------------------------------------------------------
////void DofManager::createFaceDofsAt( Cell* targetFace )
////{
////    int dofsPerFace = _faceDofInfo.size();
////    targetFace->_dof.assign( dofsPerFace, nullptr );
////    for ( int i = 0; i < dofsPerFace; i++ )
////        targetFace->_dof[ i ] = new Dof( _faceDofInfo[ i ].group );
////}
// ----------------------------------------------------------------------------
void DofManager::createNodalDofsAt( Node* targetNode )
{
    int dofsPerNode = _nodalDofInfo.size();
    targetNode->_dof.assign( dofsPerNode, nullptr );
    for ( int i = 0; i < dofsPerNode; i++ )
        targetNode->_dof[ i ] = new Dof( _nodalDofInfo[ i ].group );
}
//// ----------------------------------------------------------------------------
//Dof* DofManager::createNumericsDofWithGroup( int grp )
//{
//    Dof* newDof = new Dof( grp );
//    _numericsDof.push_back( newDof );
//
//    return newDof;
//}
// ----------------------------------------------------------------------------
void DofManager::destroyCellDofsAt( Cell* targetCell )
{
    int nDofs = targetCell->_dof.size();
    for ( int i = 0; i < nDofs; i++ )
        if ( targetCell->_dof[ i ] )
        {
            delete targetCell->_dof[ i ];
            targetCell->_dof[ i ] = nullptr;
        }
}
//// ----------------------------------------------------------------------------
//void DofManager::destroyFaceDofsAt( Cell* targetFace )
//{
//    for ( int i = 0; i < (int)_faceDofInfo.size(); i++ )
//        if ( targetFace->_dof[ i ] )
//        {
//            delete targetFace->_dof[ i ];
//            targetFace->_dof[ i ] = nullptr;
//        }
//}
// ----------------------------------------------------------------------------
void DofManager::destroyNodalDofsAt( Node* targetNode )
{
    int nDofs = targetNode->_dof.size();
    for ( int i = 0; i < nDofs; i++ )
        if ( targetNode->_dof[ i ] )
        {
            delete targetNode->_dof[ i ];
            targetNode->_dof[ i ] = nullptr;
        }
}
//// ----------------------------------------------------------------------------
//void DofManager::destroyNumericsDof( Dof*& targetDof )
//{
//    delete targetDof;
//    targetDof = nullptr;
//}
// ----------------------------------------------------------------------------
void DofManager::enslave( Dof* targetDof, Dof* masterDof )
{
    targetDof->_isSlave = true;
    targetDof->_masterDof = masterDof;
}
// ----------------------------------------------------------------------------
void DofManager::finalizeDofPrimaryValuesAtStage( int stage )
{
    int nNodes = analysisModel().domainManager().giveNumberOfNodes();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for ( int i = 0; i < nNodes; i++ )
    {
        Node* targetNode = analysisModel().domainManager().giveNode( i );

        for ( int j = 0; j < (int)_nodalDofInfo.size(); j++ )
        {
            Dof *targetDof = analysisModel().domainManager().giveNodalDof(j, targetNode);
            if ( targetDof->_stage == stage )
            {
                if ( targetDof->_isSlave )
                {
                    targetDof->_primVarConverged = targetDof->_masterDof->_primVarCurrent;
                    targetDof->_primVarCurrent = targetDof->_masterDof->_primVarCurrent;
                }
                else
                    targetDof->_primVarConverged = targetDof->_primVarCurrent;
            }
        }
    }

    for ( int dim = 0; dim < 4; dim++ )
    {
        int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for ( int i = 0; i < nCells; i++ )
        {
            Cell *curCell = analysisModel().domainManager().giveCell( i, dim );
            for ( int j = 0; j < (int) _cellDofInfo[ dim ].size(); j++)
            {
                Dof *targetDof = analysisModel().domainManager().giveCellDof( j, curCell );
                if ( targetDof->_stage == stage )
                {
                    if ( targetDof->_isSlave )
                    {
                        targetDof->_primVarConverged = targetDof->_masterDof->_primVarConverged;
                        targetDof->_primVarCurrent = targetDof->_masterDof->_primVarCurrent;
                    }
                    else
                        targetDof->_primVarConverged = targetDof->_primVarCurrent;
                }
            }
        }
    }

    for ( int i = 0; i < (int)_numericsDof.size(); i++ )
    {
        _numericsDof[ i ]->_primVarConverged = _numericsDof[ i ]->_primVarCurrent;
    }
}
// ----------------------------------------------------------------------------
void DofManager::findActiveDofs()
{
    auto nStage = analysisModel().solutionManager().giveNumberOfStages();
//    auto registeredStages = analysisModel().solutionManager().giveRegisteredSolutionStages();

    // First pass: determine number of active dofs at each stage
    _nActiveDof.assign( nStage + 1, 0 );
    _nInactiveDof.assign( nStage + 1, 0 );

//    _nActiveDof.erase(_nActiveDof.begin(), _nActiveDof.end());
//    _nInactiveDof.erase(_nInactiveDof.begin(), _nInactiveDof.end());

//    for ( auto curStage : registeredStages )
//    {
//        _nActiveDof.insert({curStage, 0});
//        _nInactiveDof.insert({curStage, 0});
//    }

    int nNodes = analysisModel().domainManager().giveNumberOfNodes();
    for ( int i = 0; i < nNodes; i++ )
    {
        Node* curNode = analysisModel().domainManager().giveNode( i );
        for ( int j = 0; j < (int)_nodalDofInfo.size(); j++ )
        {
            Dof* curDof = analysisModel().domainManager().giveNodalDof( j, curNode );
            if ( !curDof->_isConstrained && !curDof->_isSlave && curDof->_stage != UNASSIGNED )
                _nActiveDof[ curDof->_stage ]++;
            else if ( curDof->_stage != UNASSIGNED )
                _nInactiveDof[ curDof->_stage ]++;
        }
    }

    for ( int dim = 0; dim < 4; dim++ )
    {
        int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );
        for ( int i = 0; i < nCells; i++ )
        {
            Cell *curCell = analysisModel().domainManager().giveCell( i, dim );
            for ( int j = 0; j < (int) _cellDofInfo[dim].size(); j++ )
            {
                Dof *curDof = analysisModel().domainManager().giveCellDof( j, curCell );
                if ( !curDof->_isConstrained && !curDof->_isSlave && curDof->_stage != UNASSIGNED )
                    _nActiveDof[ curDof->_stage ]++;
                else if ( curDof->_stage != UNASSIGNED )
                    _nInactiveDof[ curDof->_stage ]++;
            }
        }
    }

    for ( int i = 0; i < (int)_numericsDof.size(); i++ )
    {
        Dof* curDof = _numericsDof[i];
        if ( !curDof->_isConstrained && !curDof->_isSlave && curDof->_stage != UNASSIGNED )
            _nActiveDof[ curDof->_stage ]++;
        else if ( curDof->_stage != UNASSIGNED )
            _nInactiveDof[ curDof->_stage ]++;
    }

    _activeDof.assign( nStage, std::vector<Dof*>() );
    _inactiveDof.assign( nStage, std::vector<Dof*>() );

    // 2nd pass: for vectors of Dof pointers
    for ( int stage = 1; stage <= nStage; stage++ )
    {
        _activeDof[ stage ].assign( _nActiveDof[ stage ], nullptr );
        _inactiveDof[ stage ].assign( _nActiveDof[ stage ], nullptr );
    }

    for ( int curStage = 1; curStage <= nStage; curStage++ )
    {
        int curActiveIdx = 0;
        int curInactiveIdx = 0;
        for ( int i = 0; i < nNodes; i++ )
        {
            Node* targetNode = analysisModel().domainManager().giveNode( i );

            for ( int j = 0; j < (int)_nodalDofInfo.size(); j++ )
            {
                Dof* targetDof = analysisModel().domainManager().giveNodalDof( j, targetNode );
                if ( targetDof->_stage == curStage && !targetDof->_isConstrained && !targetDof->_isSlave )
                    _activeDof[ curStage ][ curActiveIdx++ ] = targetDof;
                else if ( targetDof->_stage == curStage )
                    _inactiveDof[ curStage ][ curInactiveIdx++ ] = targetDof;
            }
        }

        for ( int dim = 0; dim < 4; dim++ )
        {
            int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );
            for ( int i = 0; i < nCells; i++ )
            {
                Cell* targetCell = analysisModel().domainManager().giveCell( i, dim );

                for ( int j = 0; j < (int)_cellDofInfo[ dim ].size(); j++ )
                {
                    Dof* targetDof = analysisModel().domainManager().giveCellDof( j, targetCell );
                    if ( targetDof->_stage == curStage && !targetDof->_isConstrained && !targetDof->_isSlave )
                        _activeDof[ curStage ][ curActiveIdx++ ] = targetDof;
                    else if ( targetDof->_stage == curStage )
                        _inactiveDof[ curStage ][ curInactiveIdx++ ] = targetDof;
                }
            }
        }

        for ( int i = 0; i < (int)_numericsDof.size(); i++ )
        {
            Dof* targetDof = _numericsDof[i];
            if ( targetDof->_stage == curStage && !targetDof->_isConstrained && !targetDof->_isSlave )
                _activeDof[ curStage ][ curActiveIdx++ ] = targetDof;
            else if ( targetDof->_stage == curStage )
                _inactiveDof[ curStage ][ curInactiveIdx++ ] = targetDof;
        }
    }
}
//// ----------------------------------------------------------------------------
//std::vector<Dof*> DofManager::giveActiveDofsAtStage( int stg )
//{
//    return _activeDof[ stg ];
//}
//// ----------------------------------------------------------------------------
int DofManager::giveGroupNumberFor( Dof* targetDof )
{
    if ( targetDof->_isSlave )
        targetDof = targetDof->_masterDof;

    return targetDof->_group;
}
// ----------------------------------------------------------------------------
int DofManager::giveIndexForCellDof( const std::string& name )
{
    int index = -1;
    for ( int dim = 0; dim < 4; dim++ )
        for ( int i = 0; i < (int)_cellDofInfo[ dim ].size(); i++ )
        {
            if ( _cellDofInfo[ dim ][ i ].tag == name )
                index = i;
        }
    if ( index < 0 )
        throw std::runtime_error( "ERROR: Cannot give DOF index. Cell DOF '" + name + "' not recognized!\n" );

    return index;
}
//// ----------------------------------------------------------------------------
//int DofManager::giveIndexForFaceDof( const std::string& name )
//{
//    int index = -1;
//    for ( int i = 0; i < (int)_faceDofInfo.size(); i++ )
//    {
//        if ( _faceDofInfo[ i ].tag == name )
//            index = i;
//    }
//    if ( index < 0 )
//        throw std::runtime_error( "ERROR: Cannot give DOF index. Face DOF '" + name + "' not recognized!\n" );
//
//    return index;
//}
// ----------------------------------------------------------------------------
int DofManager::giveIndexForNodalDof( const std::string& name )
{
    int index = -1;
    for ( int i = 0; i < (int)_nodalDofInfo.size(); i++ )
    {
        if ( _nodalDofInfo[ i ].tag == name )
            index = i;
    }
    if ( index < 0 )
        throw std::runtime_error( "ERROR: Nodal DOF '" + name + "' not recognized!\nSource: " + _name );

    return index;
}
// ----------------------------------------------------------------------------
int DofManager::giveEquationNumberAt ( Dof* targetDof )
{
    if ( targetDof->_isSlave )
        targetDof = targetDof->_masterDof;

    return targetDof->_eqNo;
}
//// ----------------------------------------------------------------------------
//int DofManager::giveNumberOfActiveDofsAtStage( int stgNum )
//{
//    return _nActiveDof[ stgNum ];
//}
//// ----------------------------------------------------------------------------
//int DofManager::giveSubsystemNumberFor( Dof* targetDof )
//{
//    if ( targetDof->_isSlave )
//        targetDof = targetDof->_masterDof;
//
//    return targetDof->_subsystem;
//}
//// ----------------------------------------------------------------------------
//double DofManager::giveValueOfConstraintAt( Dof* targetDof, ValueType valType )
//{
//    double val;
//
//    if ( valType == current_value )
//        val = targetDof->_constraintValue;
//    else if ( valType == incremental_value )
//        val = targetDof->_constraintValue - targetDof->_primVarConverged;
//    else
//        val = targetDof->_primVarConverged;
//
//    return val;
//}
// ----------------------------------------------------------------------------
double DofManager::giveValueOfPrimaryVariableAt( Dof* targetDof, ValueType valType )
{
    double val;
    if ( targetDof->_isSlave )
        targetDof = targetDof->_masterDof;

    if ( valType == current_value )
        val = targetDof->_primVarCurrent;
    else if ( valType == incremental_value )
        val = targetDof->_primVarCurrent - targetDof->_primVarConverged;
    else if ( valType == converged_value )
        val = targetDof->_primVarConverged;
    else if ( valType == correction )
        val = targetDof->_primVarCorrection;
    else
        throw std::runtime_error( "ERROR: Cannot request value of replacement correction of primary variable at DOF!" );

    return val;
}
//// ----------------------------------------------------------------------------
//double DofManager::giveValueOfSecondaryVariableAt( Dof* targetDof )
//{
//    return targetDof->_secVar;
//}
//// ----------------------------------------------------------------------------
//double DofManager::giveValueOfResidualAt( Dof* targetDof )
//{
//    if ( targetDof->_isSlave )
//        targetDof = targetDof->_masterDof;
//
//    return targetDof->_residual;
//}
// ----------------------------------------------------------------------------
void DofManager::imposeMultiFreedomConstraints()
{
    for ( int i = 0; i < (int)_multiFreedomConstraint.size(); i++ )
    {
        if ( _multiFreedomConstraint[ i ].type == "NodalDofSlaveConstraint" )
            this->imposeNodalDofSlaveConstraint( _multiFreedomConstraint[ i ] );
        else
            throw std::runtime_error( "ERROR: Unimplemented multi-freedom constraint type!\nSource: " + _name );
    }
}
// ----------------------------------------------------------------------------
void DofManager::putDirichletConstraintOn( Dof* targetDof )
{
    targetDof->_isConstrained = true;
    targetDof->_eqNo = UNASSIGNED;
}
// ----------------------------------------------------------------------------
void DofManager::readCellDofsFrom( FILE* fp )
{
    // Read number of cell DOFs
    int nCellDofs = getIntegerInputFrom( fp, "\nFailed to read number of DOF per cell from input file!", _name );

    // Setup cell DOF info
    for ( int i = 0; i < nCellDofs; i++ )
    {
        DofInfo info;

        // Cell DOF tag
        info.tag = getStringInputFrom( fp, "Failed to read cell DOF tag from input file!", _name );

        // Cell dimension
        info.dim = getIntegerInputFrom( fp, "Failed to read cell dimension from input file!", _name );
        if ( info.dim < 0 || info.dim > 3 )
            throw std::runtime_error("ERROR: Invalid value '" + std::to_string( info.dim ) +
                "' encountered in input file!\nSource: " + _name );

        // Cell DOF group
        verifyKeyword( fp, "DofGroup", _name );
        info.group = getIntegerInputFrom( fp, "Failed to read cell DOF group from input file!", _name );

        _cellDofInfo[ info.dim ].push_back( info );
    }
}
//// ----------------------------------------------------------------------------
//void DofManager::readFaceDofsFrom( FILE* fp )
//{
//    std::string key, src = "DofManager";
//
//    // Read number of face DOFs
//    int nDofsPerFace = getIntegerInputFrom( fp, "\nFailed to read number of DOF per face from input file!", src );
//
//    // Setup cell DOF info
//    _faceDofInfo.assign( nDofsPerFace, DofInfo() );
//
//    for ( int i = 0; i < nDofsPerFace; i++ )
//    {
//        // Face DOF tag
//        _faceDofInfo[ i ].tag = getStringInputFrom( fp, "\nFailed to read face DOF tag from input file!", src );
//
//        // Face DOF group
//        verifyKeyword( fp, key = "DofGroup", src );
//        _faceDofInfo[ i ].group = getIntegerInputFrom( fp, "\nFailed to read face DOF group from input file!", src );
//
//        // Nodal dof primary and secondary field assignment
//        verifyKeyword( fp, key = "FaceField", src );
//        _faceDofInfo[ i ].primField = getIntegerInputFrom( fp, "\nFailed to read primary field assignment for face DOF from input file!", src );
//        _faceDofInfo[ i ].secField = getIntegerInputFrom( fp, "\nFailed reading secondary field assignment for face DOF from input file!", src );
//    }
//
//    if ( nDofsPerFace > 0 )
//        analysisModel().domainManager().mustConstructFaces();
//}
// ----------------------------------------------------------------------------
void DofManager::readMultiFreedomConstraintsFrom( FILE* fp )
{
    std::string src = "DofManager";

    int nMFConstraints = getIntegerInputFrom( fp, "Failed to read number of multi-freedom constraints from input file!", src );
    _multiFreedomConstraint.assign( nMFConstraints, MultiFreedomConstraint() );
    for ( int i = 0; i < nMFConstraints; i++ )
    {
        _multiFreedomConstraint[ i ].type = getStringInputFrom( fp, "Failed to read multi-freedom constraint type from input file!", src );

        if ( _multiFreedomConstraint[ i ].type == "NodalDofSlaveConstraint" )
            this->readNodalDofSlaveConstraintDataFrom( fp, _multiFreedomConstraint[ i ] );
        else
            throw std::runtime_error( "Unrecognized multi-freedom constraint type encountered!\nSource: " + _name );
    }
}
// ----------------------------------------------------------------------------
void DofManager::readNodalDofsFrom( FILE* fp )
{
    // Read number of DOF per node
    int nNodalDofs = getIntegerInputFrom( fp, "\nFailed to read number of DOF per node from input file!", _name );

    // Setup nodal DOF info
    _nodalDofInfo.assign( nNodalDofs, DofInfo() );

    for (  int i = 0; i < nNodalDofs; i++ )
    {
        // Nodal DOF tag
        _nodalDofInfo[ i ].dim = 0;
        _nodalDofInfo[ i ].tag = getStringInputFrom( fp, "\nFailed to read nodal DOF tag from input file!", _name );

        // Nodal DOF group
        verifyKeyword( fp, "DofGroup", _name );
        _nodalDofInfo[ i ].group = getIntegerInputFrom( fp, "\nFailed to read nodal DOF group from input file!", _name );

        // Nodal dof primary and secondary field assignment
        verifyKeyword( fp, "NodalField", _name );
        _nodalDofInfo[ i ].primField = getIntegerInputFrom( fp, "\nFailed to read primary field assignment for\nnodal DOF from input file!", _name );
        _nodalDofInfo[ i ].secField = getIntegerInputFrom( fp, "Failed reading secondary field assignment for\nnodal DOF from input file!", _name );
    }
}
// ----------------------------------------------------------------------------
void DofManager::removeAllDofConstraints()
{
    // Cycle through all nodal DOFs
    int nNodes = analysisModel().domainManager().giveNumberOfNodes();

    for ( int i = 0; i < nNodes; i++ )
    {
        Node* targetNode = analysisModel().domainManager().giveNode( i );

        for ( int j = 0; j < (int)_nodalDofInfo.size(); j++ )
        {
            Dof* targetDof = analysisModel().domainManager().giveNodalDof( j, targetNode );
            targetDof->_isConstrained = false;
        }
    }

    // Cycle through all cell DOFs
    for ( int dim = 0; dim < 4; dim++ )
    {
        int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );
        for ( int i = 0; i < nCells; i++ )
        {
            Cell *targetCell = analysisModel().domainManager().giveCell( i, dim );

            for ( int j = 0; j < (int) _cellDofInfo[ dim ].size(); j++ )
            {
                Dof *targetDof = analysisModel().domainManager().giveCellDof( j, targetCell );
                targetDof->_isConstrained = false;
            }
        }
    }
}
// ----------------------------------------------------------------------------
void DofManager::reportNumberOfActiveDofs()
{
    std::printf( "\n    Stage    Active DOFs" );
    std::printf( "\n    -----------------------" );
    for ( int curStage = 1; curStage <= (int)_nActiveDof.size(); curStage++ )
        std::printf( "\n    %-9d%d\n", curStage, _nActiveDof[ curStage ] );
    std::printf( "\n" );
}
//// ----------------------------------------------------------------------------
//void DofManager::resetDofCurrentPrimaryValues()
//{
//    int nNodes = analysisModel().domainManager().giveNumberOfNodes();
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for (int i = 0; i < nNodes; i++)
//    {
//        Node* targetNode = analysisModel().domainManager().giveNode( i );
//
//        for ( int j = 1; j <= (int)_nodalDofInfo.size(); j++ )
//        {
//            Dof* targetDof = analysisModel().domainManager().giveNodalDof( j, targetNode );
//            targetDof->_primVarCurrent = targetDof->_primVarConverged;
//            targetDof->_secVar = 0.;
//        }
//    }
//
//    int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//    for ( int i = 0; i < nCells; i++ )
//    {
//        Cell* targetCell = analysisModel().domainManager().giveDomainCell( i );
//        for ( int j = 1; j <= (int)_cellDofInfo.size(); j++ )
//        {
//            Dof* targetDof = analysisModel().domainManager().giveCellDof( j, targetCell );
//            targetDof->_primVarCurrent = targetDof->_primVarConverged;
//            targetDof->_secVar = 0.;
//        }
//    }
//}
//// ----------------------------------------------------------------------------
//void DofManager::resetSecondaryVariablesAtStage( int stage )
//{
//    for ( int i = 0; i < (int)_activeDof[ stage ].size(); i++ )
//        _activeDof[ stage ][ i ]->_secVar = 0.;
//
//    for ( int i = 0; i < (int)_inactiveDof[ stage ].size(); i++ )
//        _inactiveDof[ stage ][ i ]->_secVar = 0.;
//}
//// ----------------------------------------------------------------------------
//void DofManager::setConstraintValueAt( Dof* targetDof, double val )
//{
//    targetDof->_constraintValue = val;
//    targetDof->_primVarCurrent = val;
//}
//// ----------------------------------------------------------------------------
//void DofManager::setEquationNumberFor( Dof* targetDof, int eqNo )
//{
//    targetDof->_eqNo = eqNo;
//}
//// ----------------------------------------------------------------------------
//void DofManager::setStageFor( Dof* targetDof, int stgNum )
//{
//    if ( targetDof->_stage != UNASSIGNED && targetDof->_stage != stgNum )
//        throw std::runtime_error( "ERROR: Detected conflict in DOF stage assignment! DOF with previously assigned stage number of '"
//                + std::to_string( targetDof->_stage ) + "' is being reassigned a stage number of '" + std::to_string( stgNum )
//                + std::string( "'.\nSource: DofManager" ) );
//    else
//        targetDof->_stage = stgNum;
//}
//// ----------------------------------------------------------------------------
//void DofManager::setSubsystemFor( Dof* targetDof, int subsysNum )
//{
//    if ( targetDof->_subsystem != UNASSIGNED && targetDof->_subsystem != subsysNum )
//        throw std::runtime_error( "ERROR: Detected conflict in subsystem assignment! DOF with assigned subsystem number of '"
//                + std::to_string( targetDof->_subsystem ) + "' is being reassigned a subsystem number of '" + std::to_string( subsysNum )
//                + std::string( "'.\nSource: DofManager" ) );
//
//    targetDof->_subsystem = subsysNum;
//}
// ----------------------------------------------------------------------------
void DofManager::updatePrimaryVariableAt( Dof*      targetDof
                                        , double    val
                                        , ValueType valType )
{
    if ( valType == current_value )
    {
        targetDof->_primVarCorrection = val - targetDof->_primVarCurrent;
        targetDof->_primVarCurrent = val;
    }
    else if ( valType == incremental_value )
    {
        targetDof->_primVarCorrection = val - targetDof->_primVarCurrent + targetDof->_primVarConverged;
        targetDof->_primVarCurrent = val + targetDof->_primVarConverged;
    }
    else if ( valType == converged_value )
    {
        targetDof->_primVarCorrection = 0.;
        targetDof->_primVarConverged = val;
        targetDof->_primVarCurrent = val;
    }
    else if ( valType == correction )
    {
        targetDof->_primVarCorrection = val;
        targetDof->_primVarCurrent += val;
    }
    else
    {
        targetDof->_primVarCurrent += val - targetDof->_primVarCorrection;
        targetDof->_primVarCorrection = val;
    }
}
//// ----------------------------------------------------------------------------
//void DofManager::updateResidualAt( Dof* targetDof, double val )
//{
//    targetDof->_residual = val;
//}
// ----------------------------------------------------------------------------
void DofManager::writeConvergedDofValuesTo( Node* targetNode )
{
    for ( int i = 0; i < (int)_nodalDofInfo.size(); i++)
    {
        Dof* curDof = targetNode->_dof[i];

        analysisModel().domainManager().setFieldValueAt( targetNode, _nodalDofInfo[ i ].primField, curDof->_primVarConverged );
        analysisModel().domainManager().setFieldValueAt( targetNode, _nodalDofInfo[ i ].secField, curDof->_secVar );
    }
}
// ----------------------------------------------------------------------------
void DofManager::imposeNodalDofSlaveConstraint( MultiFreedomConstraint& mfc )
{
//    int nBCells = analysisModel().domainManager().giveNumberOfBoundaryCells();
    int masterPhysNum = analysisModel().domainManager().givePhysicalEntityNumberFor( mfc.masterTag );
    int slavePhysNum = analysisModel().domainManager().givePhysicalEntityNumberFor( mfc.slaveTag );

    int dim = analysisModel().domainManager().giveDimensionForPhysicalEntity( masterPhysNum );
    int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );

    // Find master node
    Dof* masterDof = nullptr;
    for ( int i = 0; i < nCells; i++ )
    {
        Cell* candCell = analysisModel().domainManager().giveCell( i, dim );
        if ( analysisModel().domainManager().giveLabelOf( candCell ) == masterPhysNum )
        {
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( candCell );
            if ( (int)node.size() != 1 )
                throw std::runtime_error( "ERROR: Detected more than one master node in 'NodalDofSlaveConstraint' assignment!\nSource: DofManager" );

            masterDof = analysisModel().domainManager().giveNodalDof( mfc.masterDofNum, node[ 0 ] );
            break;
        }
    }

    // Find slave nodes
    dim = analysisModel().domainManager().giveDimensionForPhysicalEntity( slavePhysNum );
    nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );

    for ( int i = 0; i < nCells; i++ )
    {
        Cell* candCell = analysisModel().domainManager().giveCell( i, dim );
        if ( analysisModel().domainManager().giveLabelOf( candCell ) == slavePhysNum )
        {
            std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( candCell );
            for ( int j = 0; j < (int)node.size(); j++)
            {
                Dof* slaveDof = analysisModel().domainManager().giveNodalDof( mfc.slaveDofNum, node[ j ] );
                if ( slaveDof != masterDof )
                    this->enslave( slaveDof, masterDof );
            }
        }
    }
}
// ----------------------------------------------------------------------------
void DofManager::readNodalDofSlaveConstraintDataFrom( FILE* fp, MultiFreedomConstraint& mfc )
{
    std::string src = "DofManager", name;
    
    verifyKeyword( fp, "Master", src );
    mfc.masterTag = getStringInputFrom( fp, "Failed to read master node tag from input file!", src );
    name = getStringInputFrom( fp, "Failed to read master DOF name from input file!", src );
    mfc.masterDofNum = this->giveIndexForNodalDof( name );
    verifyKeyword( fp, "Slave", src );
    mfc.slaveTag = getStringInputFrom( fp, "Failed to read slave node tag from input file!", src );
    name = getStringInputFrom( fp, "Failed to read slave DOF name from input file!", src );
    mfc.slaveDofNum = this->giveIndexForNodalDof( name );
}
