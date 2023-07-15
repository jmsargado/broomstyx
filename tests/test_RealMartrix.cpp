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


#include "Math/RealMatrix.hpp"
#include "Math/floatingPointComparison.hpp"
#include "test_RealMatrix.hpp"
#include "Tester.hpp"

namespace broomstyx
{
    void test_RealMatrix()
    {
        test_RealMatrix_Construction();
        test_RealMatrix_Assignment();
        test_RealMatrix_Transposition();
        test_RealMatrix_SelfAddition();
        test_RealMatrix_SelfSubtraction();
        test_RealMatrix_InPlaceScaling();
        test_RealMatrix_InitAndErase();
    }

    void test_RealMatrix_Construction()
    {
        Tester testObj( "RealMatrix construction" );

        RealMatrix A;
        testObj.checkThat( A.dim1() == 0 );
        testObj.checkThat( A.dim2() == 0 );
        testObj.checkThat( A.ptr() == nullptr );

        RealMatrix B( 3,2 );
        testObj.checkThat( B.dim1() == 3 );
        testObj.checkThat( B.dim2() == 2 );

        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 2; j++ )
                testObj.checkThat( isEqual( B( i,j ), 0. ) );

        RealMatrix C( {{ 1, 2, 3 }, { 4, 5, 6 }} );

        testObj.checkThat( C.dim1() == 2 );
        testObj.checkThat( C.dim2() == 3 );

        testObj.checkThat( isEqual( C( 0,0 ), 1. ) );
        testObj.checkThat( isEqual( C( 0,1 ), 2. ) );
        testObj.checkThat( isEqual( C( 0,2 ), 3. ) );
        testObj.checkThat( isEqual( C( 1,0 ), 4. ) );
        testObj.checkThat( isEqual( C( 1,1 ), 5. ) );
        testObj.checkThat( isEqual( C( 1,2 ), 6. ) );

        testObj.reportResults();
    }

    void test_RealMatrix_Assignment()
    {
        Tester testObj( "RealMatrix assignment" );

        RealMatrix A, B( {{1, 1, 1}, {2, 3, 4}} );
        A = B;
        testObj.checkThat( A.dim1() == 2 );
        testObj.checkThat( A.dim2() == 3 );

        testObj.checkThat( isEqual( A( 0,0 ), 1. ) );
        testObj.checkThat( isEqual( A( 0,1 ), 1. ) );
        testObj.checkThat( isEqual( A( 0,2 ), 1. ) );
        testObj.checkThat( isEqual( A( 1,0 ), 2. ) );
        testObj.checkThat( isEqual( A( 1,1 ), 3. ) );
        testObj.checkThat( isEqual( A( 1,2 ), 4. ) );

        A = {{1,2}, {3,4}, {5,6}};
        testObj.checkThat( A.dim1() == 3 );
        testObj.checkThat( A.dim2() == 2 );

        testObj.checkThat( isEqual( A( 0,0 ), 1. ) );
        testObj.checkThat( isEqual( A( 0,1 ), 2. ) );
        testObj.checkThat( isEqual( A( 1,0 ), 3. ) );
        testObj.checkThat( isEqual( A( 1,1 ), 4. ) );
        testObj.checkThat( isEqual( A( 2,0 ), 5. ) );
        testObj.checkThat( isEqual( A( 2,1 ), 6. ) );

        testObj.reportResults();
    }

    void test_RealMatrix_Transposition()
    {
        Tester testObj( "RealMatrix transposition" );

        RealMatrix A({{1, 2}, {3, 4}, {5, 6}}), B;
        B = A.trp();

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0,0 ), 1. ) );
        testObj.checkThat( isEqual( B( 0,1 ), 3. ) );
        testObj.checkThat( isEqual( B( 0,2 ), 5. ) );
        testObj.checkThat( isEqual( B( 1,0 ), 2. ) );
        testObj.checkThat( isEqual( B( 1,1 ), 4. ) );
        testObj.checkThat( isEqual( B( 1,2 ), 6. ) );

        testObj.reportResults();
    }

    void test_RealMatrix_SelfAddition()
    {
        Tester testObj( "RealMatrix self-addition" );

        RealMatrix A({{1, 2}, {3, 4}, {5, 6}});
        RealMatrix B({{1, 1}, {1, 1}, {1, 1}});

        A += B;
        testObj.checkThat( A.dim1() == 3 );
        testObj.checkThat( A.dim2() == 2 );

        testObj.checkThat( isEqual( A( 0,0 ), 2. ) );
        testObj.checkThat( isEqual( A( 0,1 ), 3. ) );
        testObj.checkThat( isEqual( A( 1,0 ), 4. ) );
        testObj.checkThat( isEqual( A( 1,1 ), 5. ) );
        testObj.checkThat( isEqual( A( 2,0 ), 6. ) );
        testObj.checkThat( isEqual( A( 2,1 ), 7. ) );

        testObj.reportResults();
    }

    void test_RealMatrix_SelfSubtraction()
    {
        Tester testObj( "RealMatrix self-subtraction" );

        RealMatrix A({{1, 2}, {3, 4}, {5, 6}});
        RealMatrix B({{1, 1}, {1, 1}, {1, 1}});

        A -= B;
        testObj.checkThat( A.dim1() == 3 );
        testObj.checkThat( A.dim2() == 2 );

        testObj.checkThat( isEqual( A( 0,0 ), 0. ) );
        testObj.checkThat( isEqual( A( 0,1 ), 1. ) );
        testObj.checkThat( isEqual( A( 1,0 ), 2. ) );
        testObj.checkThat( isEqual( A( 1,1 ), 3. ) );
        testObj.checkThat( isEqual( A( 2,0 ), 4. ) );
        testObj.checkThat( isEqual( A( 2,1 ), 5. ) );

        testObj.reportResults();
    }

    void test_RealMatrix_InPlaceScaling()
    {
        Tester testObj( "RealMatrix in-place scaling" );

        RealMatrix A({{1, 2}, {3, 4}, {5, 6}});
        A *= -0.5;

        testObj.checkThat( A.dim1() == 3 );
        testObj.checkThat( A.dim2() == 2 );

        testObj.checkThat( isEqual( A( 0,0 ), -0.5 ) );
        testObj.checkThat( isEqual( A( 0,1 ), -1.0 ) );
        testObj.checkThat( isEqual( A( 1,0 ), -1.5 ) );
        testObj.checkThat( isEqual( A( 1,1 ), -2.0 ) );
        testObj.checkThat( isEqual( A( 2,0 ), -2.5 ) );
        testObj.checkThat( isEqual( A( 2,1 ), -3.0 ) );

        A /= 0.5;

        testObj.checkThat( A.dim1() == 3 );
        testObj.checkThat( A.dim2() == 2 );

        testObj.checkThat( isEqual( A( 0,0 ), -1. ) );
        testObj.checkThat( isEqual( A( 0,1 ), -2. ) );
        testObj.checkThat( isEqual( A( 1,0 ), -3. ) );
        testObj.checkThat( isEqual( A( 1,1 ), -4. ) );
        testObj.checkThat( isEqual( A( 2,0 ), -5. ) );
        testObj.checkThat( isEqual( A( 2,1 ), -6. ) );

        testObj.reportResults();
    }

    void test_RealMatrix_InitAndErase()
    {
        Tester testObj( "RealMatrix init & erase" );

        RealMatrix A;
        A.init( 3,4 );

        testObj.checkThat( A.dim1() == 3 );
        testObj.checkThat( A.dim2() == 4 );

        for ( int i = 0; i < 3; i++ )
            for ( int j = 0; j < 4; j++ )
                testObj.checkThat( isEqual( A( i,j ), 0. ) );

        A.erase();
        testObj.checkThat( A.dim1() == 0 );
        testObj.checkThat( A.dim2() == 0 );
        testObj.checkThat( A.ptr() == nullptr );

        testObj.reportResults();
    }
}
