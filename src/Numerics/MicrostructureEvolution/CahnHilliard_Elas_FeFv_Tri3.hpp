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

        void initializeMaterialsAt( Cell* targetCell ) override;
        void initializeNumericsAt( Cell* targetCell ) override;
        void readAdditionalDataFrom( FILE* fp ) override;

    private:
        Triangle_P1 _basisFunction;
        RealVector  _basisFunctionValues;
        RealMatrix  _basisFunctionDerivatives;

        RealVector  _gpNatCoor;
        double      _wt;

        double _M;
        double _kappa;
        double _A;
        double _eps_m;
        double _Nv;

        class CellNumericsStatus : public NumericsStatus
        {
        public:
            CellNumericsStatus();
            ~CellNumericsStatus() override = default;

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
        RealMatrix        giveBmatAt( Cell* targetCell );
        static std::vector< std::vector<Node*> > giveFaceNodesOf( Cell* targetCell );
        static double     giveDistanceToMidpointOf( std::vector<Node*>& face, RealVector& coor);
        RealMatrix        giveJacobianMatrixAt( Cell* targetCell );
        static double     giveLengthOf( std::vector<Node*>& face );
        static RealVector giveLocalDisplacementsAt( std::vector<Dof*>& dof, ValueType valType );
        std::vector<Dof*> giveNodalDofsAt( Cell* targetCell );
        static RealVector giveOutwardUnitNormalOf( std::vector<Node*>& face );
        double            giveTransmissibilityCoefficientAt( std::vector<Node*>& face
                                                           , Cell*               targetCell
                                                           , Cell*               neighborCell );
    };
}

#endif /* CAHNHILLIARD_ELAS_FEFV_TRI3_HPP */