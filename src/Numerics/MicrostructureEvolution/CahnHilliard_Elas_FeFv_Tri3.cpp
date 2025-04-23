/*
  Copyright (c) 2014 - 2019 University of Bergen
  
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

#include "CahnHilliard_Elas_FeFv_Tri3.hpp"
#include "Core/AnalysisModel.hpp"
#include "Materials/Material.hpp"
#include "Util/linearAlgebra.hpp"
#include "Util/readOperations.hpp"

using namespace broomstyx;

// Constructor for cell numerics status
CahnHilliard_Elas_FeFv_Tri3::CellNumericsStatus::CellNumericsStatus()
    : _area( 0. )
    , _phi( 0. )
    , _phiOld( 0. )
    , _strain( RealVector( 4 ) )
    , _stress( RealVector( 4 ) )
    , _Egy_chem( 0. )
    , _Egy_elas( 0. )
    , _dPsi( RealMatrix( 2,3 ) )
    , _hasConcentrationConstraint( false )
    , _hasConcentrationGradientPrescribedOnFace { false, false, false }
    , _hasPsiGradientPrescribedOnFace{ false, false, false }
    , _cValueOnFace { 0., 0., 0. }
    , _psiValueOnFace { 0., 0., 0. }
    , _hasNotComputedTransmissibilities( true )
    , _transmissibility{ 0., 0., 0. }
    , _materialStatus { nullptr, nullptr, nullptr, nullptr }
{}

// Constructor
CahnHilliard_Elas_FeFv_Tri3::CahnHilliard_Elas_FeFv_Tri3()
{
    _dim = 2;
    _dofPerCell = 2;
    _dofPerNode = 2;

    _nNodes = 3;
    _nMaterials = 4;
    _nStages = 1;
    _nSubsystems = 2;

    _name = "CahnHilliard_Elas_FeFv_Tri3";

    // Lone Gauss point coordinate and weight
    _gpNatCoor = { 1./3., 1./3. };
    _wt = 0.5;

    _M = 0.;
    _kappa = 0.;
    _A = 0.;
    _eps_m = 0.;
    _Nv = 0.;

    // Pre-calculate shape functions and derivatives at Gauss point
    _basisFunctionValues = _basisFunction.giveBasisFunctionsAt( _gpNatCoor );
    std::vector<RealVector> dpsiNat = _basisFunction.giveBasisFunctionDerivativesAt( _gpNatCoor );
    _basisFunctionDerivatives.init( 2,3 );
    for ( int i = 0; i < 3; i++)
    {
        _basisFunctionDerivatives( 0,i ) = dpsiNat[ 0 ]( i );
        _basisFunctionDerivatives( 1,i ) = dpsiNat[ 1 ]( i );
    }
}

// Public methods
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::deleteNumericsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    auto material = this->giveMaterialSetFor( targetCell );

    material[ 0 ]->destroy( cns->_materialStatus[ 0 ] );
    material[ 1 ]->destroy( cns->_materialStatus[ 1 ] );
    material[ 2 ]->destroy( cns->_materialStatus[ 2 ] );
    material[ 3 ]->destroy( cns->_materialStatus[ 3 ] );

    delete cns;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveCellFieldValueAt( Cell* targetCell, int fieldNum )
{
    RealVector val, wt;
    std::string fieldTag;

    try
    {
        fieldTag = _cellFieldOutput.at( fieldNum );
    }
    catch ( std::exception& e )
    {
        fieldTag = "unassigned";
    }

    std::tie( val,wt ) = this->giveFieldOutputAt( targetCell, fieldTag );

    return val( 0 );
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::imposeConstraintAt( Cell*                    targetCell
                                                    , int                      stg
                                                    , const BoundaryCondition& bndCond
                                                    , const TimeData&          time )
{
    // Essential BCs on nodal DOFs
    if ( bndCond.conditionType() == "NodalConstraint" )
    {
        // Retrieve nodes of boundary element
        std::vector< Node* > node = analysisModel().domainManager().giveNodesOf( targetCell );
        int dofNum = analysisModel().dofManager().giveIndexForNodalDof( bndCond.targetDof() );

        for ( auto& curNode : node)
        {
            Dof* targetDof = analysisModel().domainManager().giveNodalDof( dofNum, curNode );

            RealVector coor = analysisModel().domainManager().giveCoordinatesOf( curNode );
            double bcVal = bndCond.valueAt( coor, time );

            DofManager::setConstraintValueAt( targetDof, bcVal );
        }
    }

    // Constraints pertaining to cell DOFS
    if ( bndCond.conditionType() == "CellConstraint" )
    {
        auto cns = this->getNumericsStatusAt( targetCell );

        int dofNum = analysisModel().dofManager().giveIndexForCellDof( bndCond.targetDof() );
        Dof* targetDof = analysisModel().domainManager().giveCellDof( dofNum, targetCell );

        std::vector<RealVector> ep = this->giveEvaluationPointsFor( targetCell );
        double bcVal = bndCond.valueAt( ep[ 0 ], time );
        DofManager::setConstraintValueAt( targetDof, bcVal );

        cns->_hasConcentrationConstraint = true;
    }

    if ( bndCond.conditionType() == "ConcentrationGradient" || bndCond.conditionType() == "PsiGradient" )
    {
        std::vector<Cell*> domCell = analysisModel().domainManager().giveDomainCellsAssociatedWith( targetCell );
        for ( Cell* curDomCell : domCell )
        {
            Numerics* domCellNumerics = analysisModel().domainManager().giveNumericsFor( curDomCell );
            if ( this == domCellNumerics )
            {
                auto cns = this->getNumericsStatusAt( curDomCell );
                std::vector<Node*> bndCellNode = analysisModel().domainManager().giveNodesOf( targetCell );

                // Sanity check on mesh
                if ( bndCellNode.size() != 2 )
                    throw std::runtime_error( "Numerics type '" + _name + "' only accepts boundary cells with exactly two nodes!" );

                std::vector< std::vector<Node*> > face = giveFaceNodesOf( curDomCell );

                // Cycle through faces to determine which one corresponds to current boundary cell
                for ( int i = 0; i < 3; i++ )
                {
                    if ( (bndCellNode[ 0 ] == face[ i ][ 0 ] && bndCellNode[ 1 ] == face[ i ][ 1 ] ) ||
                         (bndCellNode[ 1 ] == face[ i ][ 0 ] && bndCellNode[ 0 ] == face[ i ][ 1 ] ) )
                    {
                        // Set boundary flag for appropriate face
                        if ( bndCond.conditionType() == "ConcentrationGradient" )
                            cns->_hasConcentrationGradientPrescribedOnFace[ i ] = true;
                        else
                            cns->_hasPsiGradientPrescribedOnFace[ i ] = true;
                    }
                }
            }
        }
    }
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::imposeInitialConditionAt( Cell* targetCell, const InitialCondition& initCond )
{
    if ( _cellDof[ 0 ] != initCond.targetDofNumber() )
        throw std::runtime_error( "Error imposing initial condition! Specified DOF in initial condition \
            does not correspond to the cell phase-field DOF.\nSource: " + _name );
    else
    {
        std::vector<RealVector> epCoor = this->giveEvaluationPointsFor( targetCell );
        double initVal = initCond.valueAt( epCoor[ 0 ] );

        Dof* pfDof = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
        DofManager::updatePrimaryVariableAt( pfDof, initVal, converged_value );
    }
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::initializeMaterialsAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );
    std::vector<Material*> material = this->giveMaterialSetFor( targetCell );

    cns->_materialStatus[ 0 ] = material[ 0 ]->createMaterialStatus();
    cns->_materialStatus[ 1 ] = material[ 1 ]->createMaterialStatus();
    cns->_materialStatus[ 2 ] = material[ 2 ]->createMaterialStatus();
    cns->_materialStatus[ 3 ] = material[ 3 ]->createMaterialStatus();
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::initializeNumericsAt( Cell* targetCell )
{
    targetCell->numericsStatus = new CellNumericsStatus();
    auto cns = this->getNumericsStatusAt( targetCell );

    // Pre-calculate det(J) and inv(J);
    RealMatrix Jmat = this->giveJacobianMatrixAt( targetCell );
    double Jdet = Jmat( 0,0 ) * Jmat( 1,1 ) - Jmat( 1,0 ) * Jmat( 0,1 );
    cns->_area = _wt * Jdet;

    // Sanity check
    if ( Jdet <= 0 )
        throw std::runtime_error( "Calculation of negative area detected!\nSource: " + _name );

    // Calculate shape function derivatives in the actual space
    RealMatrix JmatInv = inv( Jmat );
    cns->_dPsi = JmatInv * _basisFunctionDerivatives;
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::readAdditionalDataFrom( FILE* fp )
{
    verifyKeyword( fp, "Mobility", _name );
    _M = getRealInputFrom( fp, "Failed to read mobility from input file!", _name );
    verifyKeyword( fp, "GradientEnergyCoef", _name );
    _kappa = getRealInputFrom( fp, "Failed to read gradient energy coefficient from input file!", _name );
    verifyKeyword( fp, "EnergyBarrier", _name );
    _A = getRealInputFrom( fp, "Failed to read energy barrier from input file!", _name );
    verifyKeyword( fp, "LatticeMisfit", _name );
    _eps_m = getRealInputFrom( fp, "Failed to read lattice misfit from input file!", _name );
    verifyKeyword( fp, "AtomsPerUnitVolume", _name );
    _Nv = getRealInputFrom( fp, "Failed to read number of atoms per unit volume from input file!", _name );
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::removeConstraintsOn( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt(targetCell);
    cns->_hasConcentrationConstraint = false;
    cns->_hasConcentrationGradientPrescribedOnFace[ 0 ] = false;
    cns->_hasConcentrationGradientPrescribedOnFace[ 1 ] = false;
    cns->_hasConcentrationGradientPrescribedOnFace[ 2 ] = false;
    cns->_hasPsiGradientPrescribedOnFace[ 0 ] = false;
    cns->_hasPsiGradientPrescribedOnFace[ 1 ] = false;
    cns->_hasPsiGradientPrescribedOnFace[ 2 ] = false;
    cns->_hasNotComputedTransmissibilities = true;
}
// ----------------------------------------------------------------------------
void CahnHilliard_Elas_FeFv_Tri3::setDofStagesAt( Cell* targetCell )
{
    // A. Nodal DOFs
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    for ( auto& curNode : node )
    {
        Dof *dof_x, *dof_y;

        dof_x = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], curNode );
        dof_y = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], curNode );

        // Set DOF stage numbers
        DofManager::setStageFor( dof_x, _stage[ 0 ] );
        DofManager::setStageFor( dof_y, _stage[ 0 ] );

        // Set DOF subsystem numbers
        DofManager::setSubsystemFor( dof_x, _subsystem[ 0 ] );
        DofManager::setSubsystemFor( dof_y, _subsystem[ 0 ] );
    }

    // B. Cell DOFs
    Dof* dof_c = analysisModel().domainManager().giveCellDof( _cellDof[ 0 ], targetCell );
    Dof* dof_psi = analysisModel().domainManager().giveCellDof( _cellDof[ 1 ], targetCell );

    DofManager::setStageFor( dof_c, _stage[ 0 ] );
    DofManager::setSubsystemFor( dof_c, _subsystem[ 1 ] );
    DofManager::setStageFor( dof_psi, _stage[ 0 ] );
    DofManager::setSubsystemFor( dof_psi, _subsystem[ 1 ] );

}

// Private methods
// ----------------------------------------------------------------------------
CahnHilliard_Elas_FeFv_Tri3::CellNumericsStatus*
CahnHilliard_Elas_FeFv_Tri3::getNumericsStatusAt( Cell* targetCell )
{
    auto cns = dynamic_cast< CellNumericsStatus* >( targetCell->numericsStatus );
    if ( !cns )
        throw std::runtime_error( "Error: Unable to retrieve numerics status at cell!\nSource: " + _name );

    return cns;
}
// ----------------------------------------------------------------------------
RealMatrix CahnHilliard_Elas_FeFv_Tri3::giveBmatAt( Cell* targetCell )
{
    auto cns = this->getNumericsStatusAt( targetCell );

    // Fill entries for Bmat
    RealMatrix& dpsi = cns->_dPsi;
    RealMatrix bmat( { { dpsi( 0,0 ), 0.,          dpsi( 0,1 ), 0.,          dpsi( 0,2 ), 0. },
                       { 0.,          dpsi( 1,0 ), 0.,          dpsi( 1,1 ), 0.,          dpsi( 1,2 ) },
                       { dpsi( 1,0 ), dpsi( 0,0 ), dpsi( 1,1 ), dpsi( 0,1 ), dpsi( 1,2 ), dpsi( 0,2 ) } } );
    return bmat;
}
// ----------------------------------------------------------------------------
std::vector<std::vector<Node*> > CahnHilliard_Elas_FeFv_Tri3::giveFaceNodesOf( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    std::vector<Node*> faceTemplate( 2, nullptr );
    std::vector<std::vector<Node*> > face( 3, faceTemplate );

    face[ 0 ][ 0 ] = node[ 0 ];
    face[ 0 ][ 1 ] = node[ 1 ];
    face[ 1 ][ 0 ] = node[ 1 ];
    face[ 1 ][ 1 ] = node[ 2 ];
    face[ 2 ][ 0 ] = node[ 2 ];
    face[ 2 ][ 1 ] = node[ 0 ];

    return face;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveDistanceToMidpointOf( std::vector<Node*>& face
                                                            , RealVector&         coor )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf( face[ 0 ] );
    coor1 = analysisModel().domainManager().giveCoordinatesOf( face[ 1 ] );
    dx = 0.5 * ( coor0 + coor1 ) - coor;

    return std::sqrt( dx.dot( dx ) );
}
// ----------------------------------------------------------------------------
RealMatrix CahnHilliard_Elas_FeFv_Tri3::giveJacobianMatrixAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );
    RealMatrix coorMat( 3,2 );

    for ( int i = 0; i < 3; i++ )
    {
        RealVector nodeCoor = analysisModel().domainManager().giveCoordinatesOf( node[ i ] );
        coorMat( i,0 ) = nodeCoor( 0 );
        coorMat( i,1 ) = nodeCoor( 1 );
    }

    return _basisFunctionDerivatives * coorMat;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveLengthOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx;
    coor0 = analysisModel().domainManager().giveCoordinatesOf( face[ 0 ] );
    coor1 = analysisModel().domainManager().giveCoordinatesOf( face[ 1 ] );
    dx = coor1 - coor0;

    return std::sqrt( dx.dot( dx ) );
}
// ----------------------------------------------------------------------------
RealVector CahnHilliard_Elas_FeFv_Tri3::giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType )
{
    // Displacements
    RealVector u( 6 );
    u( 0 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 0 ], valType );
    u( 1 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 1 ], valType );
    u( 2 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 2 ], valType );
    u( 3 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 3 ], valType );
    u( 4 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 4 ], valType );
    u( 5 ) = DofManager::giveValueOfPrimaryVariableAt( dof[ 5 ], valType );

    return u;
}
// ----------------------------------------------------------------------------
std::vector<Dof*> CahnHilliard_Elas_FeFv_Tri3::giveNodalDofsAt( Cell* targetCell )
{
    std::vector<Node*> node = analysisModel().domainManager().giveNodesOf( targetCell );

    std::vector<Dof*> dof( 6, nullptr );
    dof[ 0 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 0 ] );
    dof[ 1 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 0 ] );
    dof[ 2 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 1 ] );
    dof[ 3 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 1 ] );
    dof[ 4 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 0 ], node[ 2 ] );
    dof[ 5 ] = analysisModel().domainManager().giveNodalDof( _nodalDof[ 1 ], node[ 2 ] );

    return dof;
}
// ----------------------------------------------------------------------------
RealVector CahnHilliard_Elas_FeFv_Tri3::giveOutwardUnitNormalOf( std::vector<Node*>& face )
{
    RealVector coor0, coor1, dx, normVec;
    coor0 = analysisModel().domainManager().giveCoordinatesOf( face[ 0 ] );
    coor1 = analysisModel().domainManager().giveCoordinatesOf( face[ 1 ] );
    dx = coor1 - coor0;
    double dist = std::sqrt(dx.dot( dx ) );

    normVec = { dx( 1 ) / dist, -dx( 0 ) / dist };

    return normVec;
}
// ----------------------------------------------------------------------------
double CahnHilliard_Elas_FeFv_Tri3::giveTransmissibilityCoefficientAt( std::vector<Node*>& face
        , Cell*               targetCell
        , Cell*               neighborCell )
{
    // A. Target cell
    // Coordinates of cell center
    std::vector<RealVector> ep = this->giveEvaluationPointsFor( targetCell );
    RealVector coor1 = ep[ 0 ];

    // length of face
    double length = giveLengthOf( face );

    // Distance from cell center to midpoint of face
    double d1 = giveDistanceToMidpointOf( face, coor1 );

    // B. Neighbor cell
    // Coordinates of cell center
    ep = this->giveEvaluationPointsFor( neighborCell );
    RealVector coor2 = ep[ 0 ];

    // Distance from cell center to midpoint of face
    double d2 = giveDistanceToMidpointOf( face, coor2 );

    return length / ( d1 + d2 );
}