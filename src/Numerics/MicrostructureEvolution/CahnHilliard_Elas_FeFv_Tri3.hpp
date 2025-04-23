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

#ifndef CAHNHILLIARD_ELAS_FEFV_TRI3_HPP
#define CAHNHILLIARD_ELAS_FEFV_TRI3_HPP

#include "Numerics/Numerics.hpp"
#include "Core/DofManager.hpp"
#include "Core/NumericsManager.hpp"
#include "BasisFunctions/Triangle_P1.hpp"

namespace broomstyx
{
    class MaterialStatus;
    
    class CahnHilliard_Elas_FeFv_Tri3 final : public Numerics
    {
    public:
        CahnHilliard_Elas_FeFv_Tri3();
        ~CahnHilliard_Elas_FeFv_Tri3() override = default;

    private:

        class CellNumericsStatus : public NumericsStatus
        {
        public:
            CellNumericsStatus();
            ~CellNumericsStatus() override = default;

        private:
            double     _area;
            double     _phi;
            double     _phiOld;
            RealVector _strain;
            RealVector _stress;
            double     _Egy_chem;
            double     _Egy_elas;
            RealMatrix _dPsi;

            bool   _hasPhsFldConstraint;
            bool   _hasPhsFldGradientPrescribedOnFace[ 3 ];
            double _valueOnFace[ 3 ];
            bool   _hasNotComputedTransmissibilities;
            double _transmissibility[ 3 ];

            MaterialStatus* _materialStatus[ 4 ];
        };

        CellNumericsStatus* getNumericsStatusAt( Cell* targetCell );
    };
}

#endif /* CAHNHILLIARD_ELAS_FEFV_TRI3_HPP */