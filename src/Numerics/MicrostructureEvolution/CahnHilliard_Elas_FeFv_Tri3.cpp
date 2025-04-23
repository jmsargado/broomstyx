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
#include "Util/linearAlgebra.hpp"

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
    , _hasPhsFldConstraint( false )
    , _hasPhsFldGradientPrescribedOnFace { false, false, false }
    , _valueOnFace { 0., 0., 0. }
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