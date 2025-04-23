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
CahnHilliard_Elas_FeFv_Tri3::CahnHilliard_Elas_FeFv_Tri3() = default;


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