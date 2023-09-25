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

#include "Numerics.hpp"
#include <stdexcept>

#include "Core/AnalysisModel.hpp"
#include "Core/DomainManager.hpp"
#include "Core/SolutionManager.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

// Numerics status
NumericsStatus::NumericsStatus() = default;

NumericsStatus::~NumericsStatus() = default;

// Numerics
Numerics::Numerics()
    : _id( 0 )
    , _stage( 0 )
    , _stageAssigned( false )
    , _dim( 0 )
    , _nDofsPerCell( 0 )
    , _nDofsPerNode( 0 )
    , _nMaterials( 0 )
    , _nNodes( 0 )
    , _nSubsystems( 0 )
    , _nCellFieldOutput( 0 )
{}

Numerics::~Numerics() = default;

//  Public methods
// ----------------------------------------------------------------------------
std::string Numerics::giveName()
{
    return _name;
}
// ----------------------------------------------------------------------------
std::tuple<RealVector,RealVector>
Numerics::giveCellFieldOutputAtEvaluationPointsOf( Cell* targetCell, int fieldNum )
{
    std::string fieldTag;
    try
    {
        fieldTag = _cellFieldOutput.at( fieldNum );
    }
    catch ( std::exception& e )
    {
        fieldTag = "unassigned";
    }
    
    auto result = this->giveFieldOutputAt( targetCell, fieldTag );
    return result;
}
// ----------------------------------------------------------------------------
int Numerics::giveSpatialDimension() const
{
    return _dim;
}
// ----------------------------------------------------------------------------
int Numerics::requiredNumberOfDofsPerCell() const
{
    return _nDofsPerCell;
}
// ----------------------------------------------------------------------------
int Numerics::requiredNumberOfMaterials() const
{
    return _nMaterials;
}
// ----------------------------------------------------------------------------
int Numerics::requiredNumberOfDofsPerNode() const
{
    return _nDofsPerNode;
}
// ----------------------------------------------------------------------------
int Numerics::requiredNumberOfNodes() const
{
    return _nNodes;
}
// ----------------------------------------------------------------------------
void Numerics::setIdTo( int id )
{
    _id = id;
}
// ----------------------------------------------------------------------------
void Numerics::setStageTo( int stage )
{
    if ( !_stageAssigned )
    {
        _stage = stage;
        _stageAssigned = true;
    }
    else
        throw std::runtime_error( "Cannot reassign stage for Numerics # " + std::to_string( _id ) + "!\n" );
}
// ----------------------------------------------------------------------------
void Numerics::imposeConstraintAt( Cell*                    targetCell
                                 , const BoundaryCondition& bndCond
                                 , const TimeData&          time )
{
    // Do nothing.
}
// ----------------------------------------------------------------------------
void Numerics::initializeMaterialsAt( Cell* targetCell ) {}
// ----------------------------------------------------------------------------
bool Numerics::performAdditionalConvergenceCheckAt( Cell* targetCell )
{
    // Does nothing if no special convergence criteria are defined.
    return true;
}
// ----------------------------------------------------------------------------
void Numerics::performPreIterationOperationsAt( int iterNum )
{    
    // Does nothing. Overloading method should be implemented in derived class
    // as needed.
}
// ----------------------------------------------------------------------------
void Numerics::performPrefinalizationCalculationsAt( Cell* targetCell )
{
    // Does nothing. Overloading method should be implemented in derived class
    // as needed.
}
// ----------------------------------------------------------------------------
void Numerics::performPreprocessingAt( Cell* targetCell, std::string tag )
{
    this->error_unimplemented( "performPreprocessingAt(..)" );
}
// ----------------------------------------------------------------------------
void Numerics::printPostIterationMessage()
{
    // Does nothing. Overloading method should be implemented in derived class as needed
}
// ----------------------------------------------------------------------------
void Numerics::readAdditionalDataFrom( FILE* fp ) {}
// ----------------------------------------------------------------------------
void Numerics::readDataFrom( FILE* fp )
{
    std::string str;
    
    // Nodal degrees of freedom to be used for numerics type
    if ( _nDofsPerNode > 0 )
    {
        verifyKeyword( fp, str = "NodalDof", _name );
        _nodalDof.assign( _nDofsPerNode, 0 );
        for ( int i = 0; i < _nDofsPerNode; i++ )
        {
            std::string name = getStringInputFrom( fp, "Failed to read nodal DOF assignment from input file!", _name );
            _nodalDof[ i ] = analysisModel().dofManager().giveIndexForNodalDof( name );
        }
    }
    
    // Cell degrees of freedom to be used for numerics type
    if ( _nDofsPerCell > 0 )
    {
        verifyKeyword( fp, str = "CellDof", _name );
        _cellDof.assign( _nDofsPerCell, 0 );
        for ( int i = 0; i < _nDofsPerCell; i++ )
        {
            std::string name = getStringInputFrom( fp, "Failed to read cell DOF assignment from input file!", _name );
            _cellDof[ i ] = analysisModel().dofManager().giveIndexForCellDof( name );
        }
    }
    
//    // Stages to be used for numerics type
//    _stage.assign( _nStages, 0 );
//    verifyKeyword( fp, str = "Stage", _name );
//    for ( int i = 0; i < _nStages; i++ )
//        _stage[ i ] = getIntegerInputFrom( fp, "Failed to read stage assignment from input file!", _name );
//
//    for ( auto curStage : _stage )
//        analysisModel().solutionManager().registerStage( curStage, _name );
    
    // Subsystem assignments
    if ( _nSubsystems > 0 )
    {
        _subsystem.assign( _nSubsystems, 0 );
        verifyKeyword( fp, str = "Subsystem", _name );
        for ( int i = 0; i < _nSubsystems; i++ )
            _subsystem[ i ] = getIntegerInputFrom( fp, "Failed to read subsystem assignment from input file!", _name );
    }
    
    // B. Cell Field Output
    verifyKeyword( fp, str = "CellFieldOutput", _name );
    _nCellFieldOutput = getIntegerInputFrom( fp, "Failed to read number of cell field output from input file!", _name );
    
    for ( int i = 0; i < _nCellFieldOutput; i++ )
    {
        int fieldNum = getIntegerInputFrom( fp, "Failed to read cell output field number from input file!", _name );
        std::string fieldTag = getStringInputFrom( fp, "Failed to read cell output tag from input file!", _name );
        
        _cellFieldOutput.insert( std::pair<int, std::string>( fieldNum, fieldTag ) );
    }
    
    this->readAdditionalDataFrom( fp );
}
// ----------------------------------------------------------------------------
void Numerics::removeConstraintsOn( Cell* targetCell )
{
    /* This method is only relevant for derived classes in which constraint
       removal is required by the numerics during some part of the solution
       process. In this case, an overriding method must be implemented in
       the derived class.
    */
}

/* -------------------------------------------------------------------------
   The following method implementations are substitute mechanisms for pure 
   virtual functions. The idea is that derived classes should only implement
   methods that are relevant to the specific numerics, however a runtime error
   should result when the some numerics implementation is paired with an incom-
   patible solution method that makes a call to a function not implemented in
   the numerics class.
*/
 
double Numerics::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
{
    this->error_unimplemented( "giveCellFieldValueAt(...)" );
    
    // Return statement just to suppress compilation warnings.
    return 0.;
}
// ----------------------------------------------------------------------------
RealVector Numerics::giveCellNodeFieldValuesAt( Cell* targetCell, int fieldNum )
{
    this->error_unimplemented( "giveCellNodeFieldValuesAt(...)" );
    
    // Return statement just to suppress compilation warnings.
    RealVector dummy;
    return dummy;
}
// ----------------------------------------------------------------------------
std::vector<RealVector> Numerics::giveEvaluationPointsFor( Cell* targetCell )
{
    this->error_unimplemented( "giveEvaluationPointsFor(...)" );
    
    // Return statement just to suppress compilation warnings.
    std::vector<RealVector> dummy;
    return dummy;
}
// ----------------------------------------------------------------------------
std::tuple< RealVector, RealVector >
Numerics::giveFieldOutputAt( Cell* targetCell, const std::string& fieldTag )
{
    this->error_unimplemented( "giveFieldOutputAt(...)" );
    
    // Return statement just to suppress compilation warnings.
    RealVector dummy;
    return std::make_tuple( dummy, dummy );
}
// ----------------------------------------------------------------------------
RealVector Numerics::giveNumericsParameter( const std::string& paramTag )
{
    throw std::runtime_error( "nERROR: Unknown parameter '" + paramTag + "' requested from numerics!\n" );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
Numerics::giveStaticCoefficientMatrixAt( Cell*           targetCell
                                       , int             subsys
                                       , const TimeData& time )
{
    // Return statement just to suppress compilation warnings.
    std::vector<Dof*> dofdummy;
    RealVector dummy;
    return std::make_tuple( dofdummy, dofdummy, dummy );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Numerics::giveStaticLeftHandSideAt( Cell*           targetCell
                                  , int             subsys
                                  , const TimeData& time )
{
    // Return statement just to suppress compilation warnings.
    std::vector<Dof*> dofdummy;
    RealVector dummy;
    return std::make_tuple( dofdummy, dummy );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>
          , std::vector<Dof*>
          , RealVector >
Numerics::giveTransientCoefficientMatrixAt( Cell*           targetCell
                                          , int             subsys
                                          , const TimeData& time )
{
    std::vector<Dof*> dof;
    RealVector val;

    return std::make_tuple( dof, dof, val );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Numerics::giveTransientLeftHandSideAt( Cell*           targetCell
                                     , int             subsys
                                     , const TimeData& time
                                     , ValueType       valType )
{
    std::vector<Dof*> dof;
    RealVector val;

    return std::make_tuple( dof, val );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Numerics::giveStaticRightHandSideAt( Cell*                    targetCell
                                   , int                      subsys
                                   , const BoundaryCondition& bndCond
                                   , const TimeData&          time )
{
    this->error_unimplemented( "giveStatusRightHandSideAt(... bndCond ...)" );
    // Return statement just to suppress compilation warnings.
    std::vector<Dof*> dofdummy;
    RealVector dummy;
    return std::make_tuple( dofdummy, dummy );
}
// ----------------------------------------------------------------------------
std::tuple< std::vector<Dof*>, RealVector >
Numerics::giveStaticRightHandSideAt( Cell*                 targetCell
                                   , int                   subsys
                                   , const FieldCondition& fldCond
                                   , const TimeData&       time )
{
    this->error_unimplemented( "giveStatusRightHandSideAt(... fldCond ...)" );
    // Return statement just to suppress compilation warnings.
    std::vector<Dof*> dofdummy;
    RealVector dummy;
    return std::make_tuple( dofdummy, dummy );
}
// ----------------------------------------------------------------------------
void Numerics::imposeInitialConditionAt( Cell*                   targetCell
                                       , const InitialCondition& initCond )
{
    this->error_unimplemented( "imposeInitialConditionAt(...)" );
}
// ----------------------------------------------------------------------------
void Numerics::performPostprocessingAt( Cell* targetCell, std::string tag )
{
    this->error_unimplemented( "performPostProcessingAt(...)" );
}
// ----------------------------------------------------------------------------
void Numerics::setDofStagesAt( Cell* targetCell )
{
    this->error_unimplemented( "setDofStagesAt(...)" );
}

// Helper methods
// ----------------------------------------------------------------------------
std::vector<Material*> Numerics::giveMaterialSetFor( Cell* targetCell ) const
{
    int label = targetCell->label();
    return analysisModel().domainManager().giveMaterialSetForDomain( label, _stage );
}
// ----------------------------------------------------------------------------
void Numerics::error_unimplemented( const std::string& method )
{
    throw std::runtime_error( "\nError: Call to unimplemented method '"
            + _name + "::" + method + "' encountered!\n" );
}