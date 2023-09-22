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

#ifndef QUADRILATERAL_P2_3_HPP
#define QUADRILATERAL_P2_3_HPP

/***********************************************************************
   Basis functions computed according to the following nodal ordering:

                           4 --- 7 --- 3
                           |           |
                           |           |
                           8           6
                           |           |
                           |           |
                           1 --- 5 --- 2

 ***********************************************************************/

#include "ScalarBasisFunction.hpp"

namespace broomstyx
{
    class Quadrilateral_P2_3 final : public ScalarBasisFunction
    {
    public:
        Quadrilateral_P2_3();
        virtual ~Quadrilateral_P2_3();
        
        RealVector giveBasisFunctionsAt( const RealVector& coor ) override;
        std::vector<RealVector> giveBasisFunctionDerivativesAt( const RealVector& coor ) override;
    };
}

#endif /* QUADRILATERAL_P2_3_HPP */