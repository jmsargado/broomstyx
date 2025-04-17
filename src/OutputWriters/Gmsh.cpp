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

#include "Gmsh.hpp"

#include <chrono>
#include <stdexcept>
#include "Core/AnalysisModel.hpp"
#include "Core/ObjectFactory.hpp"
#include "Core/DomainManager.hpp"
#include "MeshReaders/MeshReader.hpp"
#include "Numerics/Numerics.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

registerBroomstyxObject( OutputWriter, Gmsh )

// Constructor
Gmsh::Gmsh()
{
    // Initialize counters
    _writeCounter = 0;

    _nElem = 0; // Calculated the first time output has to be written

    _isBinary = false;
    _nNodeData = 0;
    _nodeDataScalars = 0;
    _nodeDataVectors = 0;
    _nodeDataTensors = 0;

    _nElemData = 0;
    _elemDataScalars = 0;
    _elemDataVectors = 0;
    _elemDataTensors = 0;

    _nElemNodeData = 0;
    _elemNodeDataScalars = 0;
    _elemNodeDataVectors = 0;
    _elemNodeDataTensors = 0;

    _name = "Gmsh (OutputWriter)";
}

// Destructor
Gmsh::~Gmsh() = default;

// Public methods
void Gmsh::initialize()
{
    // Make directory
    const int dir_err = system( "mkdir -p Output_Gmsh" );
    if ( dir_err == -1 )
        throw std::runtime_error( "Error creating directory 'Output_Gmsh'!\n" );
}

void Gmsh::readDataFrom( FILE *fp )
{
    std::string str;
    // Output filename
    verifyDeclaration( fp, "FILENAME", _name );
    _outputFilename = getStringInputFrom(fp, "Failed to read Gmsh output filename from input file!", _name );
    
    // Output style
    verifyDeclaration( fp, "OUTPUT_STYLE", _name );
    str = getStringInputFrom( fp, "Failed to read Gmsh output style (binary/ascii) from input file!", _name );
    if ( str == "binary" )
        _isBinary = true;
    else if ( str == "ascii" )
        _isBinary = false;
    else
        throw std::runtime_error("Invalid output style specification '" + str + "' encountered in input file!"
            + std::string("\nValid options: 'binary' or 'ascii'\nSource: ") + _name );
    
    // Node data
    verifyDeclaration( fp, "NODE_DATA", _name );
    _nNodeData = getIntegerInputFrom( fp, "Failed to read number of node data output from input file!", _name );
    _nodeData.assign( _nNodeData, OutputData() );
    
    for ( int i = 0; i < _nNodeData; i++ )
    {
        str = getStringInputFrom( fp, "Failed to read node data type from input file!", _name );
        if ( str == "SCALAR" )
        {
            _nodeData[ i ].dataType = scalar;
            _nodeData[ i ].field.assign( 1, 0 );
        }
        else if ( str == "VECTOR" )
        {
            _nodeData[ i ].dataType = vector;
            _nodeData[ i ].field.assign( 3, 0 );
        }
        else if ( str == "TENSOR" )
        {
            _nodeData[ i ].dataType = tensor;
            _nodeData[ i ].field.assign( 9, 0 );
        }
        else
            throw std::runtime_error( "Invalid point data type '" + str
                + "' encountered in input file!\nSource: " + _name );
        
        // Data name to use in Gmsh
        _nodeData[ i ].name = getStringInputFrom( fp, "Failed reading node data label from input file!", _name );
        
        // Nodal fields comprising components of output data
        for ( int& j : _nodeData[ i ].field )
            j = getIntegerInputFrom(fp, "Failed to read nodal field number for point data from input file.", _name );
    }
    
    // Element data
    verifyDeclaration( fp, "ELEMENT_DATA", _name );
    _nElemData = getIntegerInputFrom( fp, "Failed reading number of cell data output from input file!", _name );
    
    _elemData.assign( _nElemData, OutputData() );
    
    for ( int i = 0; i < _nElemData; i++ )
    {
        str = getStringInputFrom( fp, "Failed reading element data type from input file!", _name );
        if ( str == "SCALAR" )
        {
            _elemData[ i ].dataType = scalar;
            _elemData[ i ].field.assign( 1, 0 );
        }
        else if ( str == "VECTOR" )
        {
            _elemData[ i ].dataType = vector;
            _elemData[ i ].field.assign( 3, 0 );
        }
        else if ( str == "TENSOR" )
        {
            _elemData[ i ].dataType = tensor;
            _elemData[ i ].field.assign( 9, 0 );
        }
        else
            throw std::runtime_error( "Invalid element data type '" + str
                + "' encountered in input file!\nSource: " + _name );
        
        // Data name to use in Paraview
        _elemData[ i ].name = getStringInputFrom( fp, "Failed reading cell data label from input file!", _name );
        
        // Element fields comprising components of output data
        for ( int& j : _elemData[ i ].field )
            j = getIntegerInputFrom( fp, "Failed to read cell field number for element data from input file!", _name );
    }

    // Element node data
    verifyDeclaration( fp, "ELEMENT_NODE_DATA", _name );
    _nElemNodeData = getIntegerInputFrom( fp, "Failed reading number of cell data output from input file!", _name );
    _elemNodeData.assign(_nElemNodeData, OutputData());
    
    for ( int i = 0; i < _nElemNodeData; i++ )
    {
        str = getStringInputFrom( fp, "Failed reading element node data type from input file!", _name );
        if ( str == "SCALAR" )
        {
            _elemNodeData[ i ].dataType = scalar;
            _elemNodeData[ i ].field.assign( 1, 0 );
        }
        else if ( str == "VECTOR" )
        {
            _elemNodeData[ i ].dataType = vector;
            _elemNodeData[ i ].field.assign( 3, 0 );
        }
        else if ( str == "TENSOR" )
        {
            _elemNodeData[ i ].dataType = tensor;
            _elemNodeData[ i ].field.assign( 9, 0 );
        }
        else
            throw std::runtime_error("Invalid element node data type '" + str
                + "' encountered in input file!\nSource: " + _name );
        
        // Data name to use in Gmsh
        _elemNodeData[ i ].name = getStringInputFrom( fp, "Failed reading element node data label from input file!", _name );
        
        // Element fields comprising components of output data
        for ( int& j : _elemNodeData[ i ].field )
            j = getIntegerInputFrom( fp, "Failed to read cell field number for element node data from input file!", _name );
    }
}

void Gmsh::writeOutput( double time )
{
    // build complete string for vtu filename
    std::string mshFilename;
    
    mshFilename = "./Output_Gmsh/" + _outputFilename + "_" + std::to_string(_writeCounter) + ".msh";
    
    // Write .msh file
    std::printf( "\n  %-40s", "Writing results to file ..." );
    
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc{};
    tic = std::chrono::high_resolution_clock::now();
    
    FILE *mshFile = fopen( mshFilename.c_str(), "w" );

    // Mesh format
    std::fprintf( mshFile, "$MeshFormat" );
    std::fprintf( mshFile, "\n2.2 " );
    if ( _isBinary )
        std::fprintf( mshFile, "1 " );
    else
        std::fprintf( mshFile, "0 " );
    
    std::fprintf( mshFile, "%d", (int)sizeof(double) );
    std::fprintf( mshFile, "\n$EndMeshFormat" );

    // Physical Names
    std::fprintf( mshFile, "\n$PhysicalNames" );
    int nPhysNames = analysisModel().domainManager().giveNumberOfPhysicalNames();
    std::fprintf( mshFile, "\n%d\n", nPhysNames );
    for ( int i = 0; i < nPhysNames; i++ )
    {
        DomainManager::PhysicalEntity physEnt = analysisModel().domainManager().givePhysicalEntity( i );
        std::fprintf( mshFile, "%d %d %d %s\n", i + 1, physEnt.dimension, physEnt.entityNumber, physEnt.name.c_str() );
    }
    std::fprintf( mshFile, "$EndPhysicalNames" );

    // Nodes
    std::fprintf( mshFile, "\n$Nodes" );
    int nNodes = analysisModel().domainManager().giveNumberOfNodes();
    std::fprintf( mshFile, "\n%d\n", nNodes );

    for ( int i = 0; i < nNodes; i++ )
    {
        Node* curNode = analysisModel().domainManager().giveNode( i );
        int nodeId = DomainManager::giveIdOf( curNode );
        std::fprintf( mshFile, "%d ", nodeId + 1 );
        
        RealVector coor;
        coor = DomainManager::giveCoordinatesOf( curNode );

        if ( coor.dim() == 2 )
            std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n", coor( 0 ), coor( 1 ), 0.0 );
        else
            std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n", coor( 0 ), coor( 1 ), coor( 2 ) );
    }
    std::fprintf( mshFile, "$EndNodes\n" );
    
    // Elements
    std::fprintf( mshFile, "$Elements\n" );

    if ( _nElem == 0 )
    {
        for ( int dim : { 0, 1, 2, 3 } )
        {
            int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );
            for ( int i = 0; i < nCells; i++ )
            {
                Cell* curCell = analysisModel().domainManager().giveCell( i, dim );
                if ( curCell->hasNumerics() )
                    _nElem++;
            }
        }
    }
    std::fprintf( mshFile, "%d\n", _nElem );

    int curElemCount = 0;
    for ( int dim : { 0, 1, 2, 3 } )
    {
        int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );
        for ( int i = 0; i < nCells; i++ )
        {
            curElemCount++;
            Cell* curCell = analysisModel().domainManager().giveCell( i, dim );
            std::fprintf( mshFile, "%d ", curElemCount );

            int elType = DomainManager::giveElementTypeOf(curCell);
            int label = curCell->label();

            std::fprintf( mshFile, "%d 1 %d", elType, label );

            std::vector<Node*> cellNode = DomainManager::giveNodesOf(curCell);

            // Adjust element nodes for 1-based numbering
            for ( auto & j : cellNode                         )
                std::fprintf( mshFile, " %d", DomainManager::giveIdOf( j ) + 1 );
            std::fprintf( mshFile, "\n" );
        }
    }
    std::fprintf( mshFile, "$EndElements\n" );

    // A. Node Data
    if ( _nNodeData > 0 )
    {
        for ( int i = 0; i < _nNodeData; i++ )
        {
            std::fprintf(mshFile, "$NodeData\n");
            // String tags
            std::fprintf(mshFile, "1\n\"%s\"\n", _nodeData[i].name.c_str());
            
            // Real tags
            std::fprintf(mshFile, "1\n%.15e\n", time);
            
            // Integer tags
            std::fprintf(mshFile, "3\n%d\n", _writeCounter);
            int nComp;
            switch ( _nodeData[i].dataType )
            {
                case scalar:
                    nComp = 1;
                    break;
                case vector:
                    nComp = 3;
                    break;
                default:
                    nComp = 9;
            }
            std::fprintf(mshFile, "%d\n", nComp);
            
            int nNodes = analysisModel().domainManager().giveNumberOfNodes();
            std::fprintf(mshFile, "%d\n", nNodes);
            for (  int j = 0; j < nNodes; j++ )
            {
                std::fprintf(mshFile, "%d ", j+1);
                Node* curNode = analysisModel().domainManager().giveNode(j);
                
                if ( _nodeData[i].dataType == scalar )
                    std::fprintf(mshFile, "%25.15e\n", analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[0]));
                else if ( _nodeData[i].dataType == vector )
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n",
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[0]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[1]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[2]));
                }
                else
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[0]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[1]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[2]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[3]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[4]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[5]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[6]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[7]),
                            analysisModel().domainManager().giveFieldValueAt(curNode, _nodeData[i].field[8]));
                }
            }
            std::fprintf(mshFile, "$EndNodeData\n");
        }
    }
    
    // B. Element Data
    if ( _nElemData > 0 )
    {
        for ( int i = 0; i < _nElemData; i++)
        {
            std::fprintf(mshFile, "$ElementData\n");
            // String tags
            std::fprintf(mshFile, "1\n\"%s\"\n", _elemData[i].name.c_str());
            
            // Real tags
            std::fprintf(mshFile, "1\n%.15e\n", time);
            
            // Integer tags
            std::fprintf(mshFile, "3\n%d\n", _writeCounter);
            int nComp;
            switch ( _elemData[i].dataType )
            {
                case scalar:
                    nComp = 1;
                    break;
                case vector:
                    nComp = 3;
                    break;
                default:
                    nComp = 9;
            }
            std::fprintf(mshFile, "%d\n", nComp);
            
            int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
            std::fprintf(mshFile, "%d\n", nCells);
            for (  int j = 0; j < nCells; j++ )
            {
                std::fprintf(mshFile, "%d ", j+1);
                Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                int label = analysisModel().domainManager().giveLabelOf(curCell);
                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                
                if ( _elemData[i].dataType == scalar )
                    std::fprintf(mshFile, "%25.15e\n", numerics->giveCellFieldValueAt(curCell, _elemData[i].field[0]));
                else if ( _elemData[i].dataType == vector )
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e\n",
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[0]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[1]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[2]));
                }
                else
                {
                    std::fprintf(mshFile, "%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e%25.15e\n",
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[0]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[1]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[2]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[3]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[4]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[5]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[6]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[7]),
                            numerics->giveCellFieldValueAt(curCell, _elemData[i].field[8]));
                }
            }
            std::fprintf(mshFile, "$EndElementData\n");
        }
    }
    
    // C. Element Node Data
    if ( _nElemNodeData > 0 )
    {
        for ( int i = 0; i < _nElemNodeData; i++)
        {
            std::fprintf(mshFile, "$ElementNodeData\n");
            // String tags
            std::fprintf(mshFile, "1\n\"%s\"\n", _elemNodeData[i].name.c_str());
            
            // Real tags
            std::fprintf(mshFile, "1\n%.15e\n", time);
            
            // Integer tags
            std::fprintf(mshFile, "3\n%d\n", _writeCounter);
            int nComp;
            switch ( _elemNodeData[i].dataType )
            {
                case scalar:
                    nComp = 1;
                    break;
                case vector:
                    nComp = 3;
                    break;
                default:
                    nComp = 9;
            }
            std::fprintf(mshFile, "%d\n", nComp);
            
            int nCells = analysisModel().domainManager().giveNumberOfDomainCells();
            std::fprintf(mshFile, "%d\n", nCells);
            
            for (  int j = 0; j < nCells; j++ )
            {
                Cell* curCell = analysisModel().domainManager().giveDomainCell(j);
                int nCellNodes = analysisModel().domainManager().giveNumberOfNodesOf(curCell);
                std::fprintf(mshFile, "%d %d ", j+1, nCellNodes);
                
                int label = analysisModel().domainManager().giveLabelOf(curCell);
                Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain(label);
                
                if ( _elemNodeData[i].dataType == scalar )
                {
                    RealVector cellNodeVals = numerics->giveCellNodeFieldValuesAt(curCell, _elemNodeData[i].field[0]);
                    
                    for (  int k = 0; k < nCellNodes; k++ )
                        std::fprintf(mshFile, "%25.15e", cellNodeVals(k));
                    std::fprintf(mshFile, "\n");
                }
                else if ( _elemNodeData[i].dataType == vector )
                {
                    std::vector<RealVector> cellNodeVals;
                    cellNodeVals.assign(3, RealVector());
                    for (  int p = 0; p < 3; p++ )
                        cellNodeVals[p] = numerics->giveCellNodeFieldValuesAt(curCell, _elemNodeData[i].field[p]);
                    
                    for (  int k = 0; k < nCellNodes; k++ )
                        for ( int p = 0; p < 3; p++ )
                            std::fprintf(mshFile, "%25.15e", cellNodeVals[p](k));
                    std::fprintf(mshFile, "\n");
                }
                else
                {
                    std::vector<RealVector> cellNodeVals;
                    cellNodeVals.assign(9, RealVector());
                    for (  int p = 0; p < 9; p++ )
                        cellNodeVals[p] = numerics->giveCellNodeFieldValuesAt(curCell, _elemNodeData[i].field[p]);
                    
                    for (  int k = 0; k < nCellNodes; k++ )
                        for ( int p = 0; p < 9; p++ )
                            std::fprintf(mshFile, "%25.15e", cellNodeVals[p](k));
                    std::fprintf(mshFile, "\n");
                }
            }
            std::fprintf(mshFile, "$EndElementNodeData\n");
        }
    }
    

    std::fclose(mshFile);
    
    // Increment file count
    _writeCounter += 1;

    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf("done (time = %f sec.)\n", tictoc.count());
}