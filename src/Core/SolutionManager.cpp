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

#include "SolutionManager.hpp"
#include <chrono>
#include <stdexcept>

#include "AnalysisModel.hpp"
#include "Util/Diagnostics.hpp"
#include "ObjectFactory.hpp"
#include "DomainManager.hpp"
#include "LoadStep.hpp"
#include "Node.hpp"
#include "OutputManager.hpp"
#include "Numerics/Numerics.hpp"
#include "User/UserFunction.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

// Constructor
SolutionManager::SolutionManager()
{
    _name = "SolutionManager";
}

// Destructor
SolutionManager::~SolutionManager()
{
    for ( int i = 0; i < (int)_loadStep.size(); i++ )
        if ( _loadStep[ i ] )
            delete _loadStep[ i ];

    for ( int i = 0; i < (int)_userFunction.size(); i++ )
        delete _userFunction[ i ];
}

// Public methods
// ----------------------------------------------------------------------------
void SolutionManager::commenceSolution()
{
    std::chrono::time_point<std::chrono::system_clock> tic, toc;
    std::chrono::duration<double> tictoc;

    // Set stages for nodal and elemental degrees of freedom
    for ( int dim = 0; dim < 4; dim++ )
    {
        tic = std::chrono::high_resolution_clock::now();
        int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );
        for ( int i = 0; i < nCells; i++ )
        {
            Cell *curCell = analysisModel().domainManager().giveCell( i, dim );
            for ( int curStage = 1; curStage <= _nStage; curStage++ )
            {
                Numerics* numerics = analysisModel().domainManager().giveNumericsFor( curCell, curStage );
                if ( numerics )
                    numerics->setDofStagesAt( curCell );
            }
        }
    }
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    diagnostics().addSetupTime( tictoc.count() );

    // Impose initial conditions
    std::printf( "\n  %-40s", "Imposing initial conditions ..." );
    std::fflush( stdout );
    tic = std::chrono::high_resolution_clock::now();
    this->imposeInitialConditions();
    toc = std::chrono::high_resolution_clock::now();
    tictoc = toc - tic;
    std::printf( "done (time = %f sec.)\n", tictoc.count() );
    diagnostics().addSetupTime( tictoc.count() );

    // Output initial (assimed current time is t = 0)
    analysisModel().outputManager().writeOutputQuantities( 0. );

    // Solve load steps
    for ( int i = 0; i < (int)_loadStep.size(); i++ )
    {
        _curLoadStep = _loadStep[ i ];
        _curLoadStep->solveYourself();
    }
}
// ----------------------------------------------------------------------------
LoadStep* SolutionManager::giveCurrentLoadStep()
{
    return _curLoadStep;
}
// ----------------------------------------------------------------------------
int SolutionManager::giveNumberOfStages()
{
    return _nStage;
}
// ----------------------------------------------------------------------------
void SolutionManager::imposeInitialConditions()
{
    for ( int i = 0; i < (int)_initCond.size(); i++ )
    {
        std::string domainLabel = _initCond[ i ].domainLabel();
        std::string condType = _initCond[ i ].conditionType();
        int domainId = analysisModel().domainManager().givePhysicalEntityNumberFor( domainLabel );
        int dim = analysisModel().domainManager().giveDimensionForPhysicalEntity( domainId );

        int nCells = analysisModel().domainManager().giveNumberOfCellsWithDimension( dim );

        if ( condType == "nodalDof" )
        {
            int nNodes = analysisModel().domainManager().giveNumberOfNodes();
            std::vector<bool> nodeIsInitialized( nNodes, false );

            for ( int j = 0; j < nCells; j++ )
            {
                Cell* curCell = analysisModel().domainManager().giveCell( j, dim );
                int cellLabel = analysisModel().domainManager().giveLabelOf( curCell );

                if ( cellLabel == domainId )
                {
                    std::vector<Node*> cellNode = analysisModel().domainManager().giveNodesOf( curCell );
                    for ( auto it = cellNode.begin(); it != cellNode.end(); it++ )
                    {
                        int nodeId = analysisModel().domainManager().giveIdOf( *it );
                        if ( !nodeIsInitialized[ nodeId ] )
                        {
                            RealVector coor = analysisModel().domainManager().giveCoordinatesOf( *it );
                            int dofNum = _initCond[ i ].targetDofNumber();
                            double val = _initCond[ i ].valueAt( coor );

                            Dof* targetDof = analysisModel().domainManager().giveNodalDof( dofNum, *it );
                            analysisModel().dofManager().updatePrimaryVariableAt( targetDof, val, converged_value );
                        }
                    }
                }
            }
        }
        else if ( condType == "CellDof" )
        {
            Numerics* numerics = analysisModel().domainManager().giveNumericsForDomain( domainId );

            for ( int j = 0; j < nCells; j++ )
            {
                Cell* curCell = analysisModel().domainManager().giveDomainCell( j );
                int cellLabel = analysisModel().domainManager().giveLabelOf( curCell );

                if ( cellLabel == domainId )
                {
                    numerics->imposeInitialConditionAt( curCell, _initCond[ i ] );
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
// std::vector<int> SolutionManager::giveRegisteredSolutionStages()
// {
//     std::vector<int> stage;
//     stage.assign( (int)_stage.size(), 0 );
//     int counter = 0;
//     for ( auto i : _stage )
//         stage[ counter++ ] = i;
//
//     return stage;
// }
// ----------------------------------------------------------------------------
UserFunction* SolutionManager::makeNewUserFunction( std::string name )
{
    UserFunction* usrFcn = objectFactory().instantiateUserFunction( name );

    _userFunction.push_back( usrFcn );
    return usrFcn;
}
// ----------------------------------------------------------------------------
void SolutionManager::readInitialConditionsFrom( FILE* fp )
{
    int nInitCond = getIntegerInputFrom( fp, "Failed to read number of initial conditions from input file!", _name );
    _initCond.assign( nInitCond, InitialCondition() );

    for ( int i = 0; i < nInitCond; i++ )
        _initCond[ i ].readDataFrom( fp );
}
// ----------------------------------------------------------------------------
void SolutionManager::readLoadStepsFrom( FILE *fp )
{
    int nLoadSteps = getIntegerInputFrom( fp, "Failed to read number of load steps from input file!", _name );

    _loadStep.assign( nLoadSteps, nullptr );

    for ( int i = 0; i < nLoadSteps; i++ )
    {
        int lsNum = getIntegerInputFrom( fp, "Failed to read load step number from input file!", _name );
        _loadStep[ i ] = new LoadStep( lsNum, _stage.size() );
        _loadStep[ i ]->readDataFrom( fp );
    }
}
// ----------------------------------------------------------------------------
void SolutionManager::readNumberOfStagesFrom( FILE* fp )
{
    _nStage = getIntegerInputFrom( fp, "Failed to read number of solution stages from input file.", _name );
}
