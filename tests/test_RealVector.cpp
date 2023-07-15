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


#include "Math/RealVector.hpp"
#include "Math/floatingPointComparison.hpp"
#include "test_RealVector.hpp"
#include "Tester.hpp"

namespace broomstyx
{
    void test_RealVector()
    {
        test_RealVector_Construction();
        test_RealVector_Assignment();
        test_RealVector_SelfAddition();
        test_RealVector_SelfSubtraction();
        test_RealVector_InPlaceScaling();
        test_RealVector_CrossProduct();
        test_RealVector_DotProduct();
        test_RealVector_TensorProduct();
        test_RealVector_InitAndErase();
    }

    void test_RealVector_Construction()
    {
        Tester testObj( "RealVector construction" );

        RealVector A;
        testObj.checkThat( A.dim() == 0 );
        testObj.checkThat( A.ptr() == nullptr );

        RealVector B( 3 );
        testObj.checkThat( B.dim() == 3 );

        for ( int i = 0; i < 3; i++ )
            testObj.checkThat( isEqual( B( i ), 0. ) );

        RealVector C( { 1, 2, 3 } );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 1. ) );
        testObj.checkThat( isEqual( C( 1 ), 2. ) );
        testObj.checkThat( isEqual( C( 2 ), 3. ) );

        testObj.reportResults();
    }

    void test_RealVector_Assignment()
    {
        Tester testObj( "RealVector assignment" );

        RealVector A, B( {2, 3, 4} );
        A = B;
        testObj.checkThat( A.dim() == 3 );

        testObj.checkThat( isEqual( A( 0 ), 2. ) );
        testObj.checkThat( isEqual( A( 1 ), 3. ) );
        testObj.checkThat( isEqual( A( 2 ), 4. ) );

        A = { 10, 11, 12, 13 };
        testObj.checkThat( A.dim() == 4 );

        testObj.checkThat( isEqual( A( 0 ), 10. ) );
        testObj.checkThat( isEqual( A( 1 ), 11. ) );
        testObj.checkThat( isEqual( A( 2 ), 12. ) );
        testObj.checkThat( isEqual( A( 3 ), 13. ) );

        testObj.reportResults();
    }

    void test_RealVector_SelfAddition()
    {
        Tester testObj( "RealVector self-addition" );

        RealVector A( { 1, 2, 3, 4 } );
        RealVector B( { 1, 1, 1, 1 } );

        A += B;
        testObj.checkThat( A.dim() == 4 );

        testObj.checkThat( isEqual( A( 0 ), 2. ) );
        testObj.checkThat( isEqual( A( 1 ), 3. ) );
        testObj.checkThat( isEqual( A( 2 ), 4. ) );
        testObj.checkThat( isEqual( A( 3 ), 5. ) );

        testObj.reportResults();
    }

    void test_RealVector_SelfSubtraction()
    {
        Tester testObj( "RealVector self-subtraction" );

        RealVector A( { 1, 2, 3, 4 } );
        RealVector B( { 1, 1, 1, 1 } );

        A -= B;
        testObj.checkThat( A.dim() == 4 );

        testObj.checkThat( isEqual( A( 0 ), 0. ) );
        testObj.checkThat( isEqual( A( 1 ), 1. ) );
        testObj.checkThat( isEqual( A( 2 ), 2. ) );
        testObj.checkThat( isEqual( A( 3 ), 3. ) );

        testObj.reportResults();
    }

    void test_RealVector_InPlaceScaling()
    {
        Tester testObj( "RealVector in-place scaling" );

        RealVector A( { 1, 2, 3, 4, 5, 6 } );
        A *= -0.5;

        testObj.checkThat( A.dim() == 6 );

        testObj.checkThat( isEqual( A( 0 ), -0.5 ) );
        testObj.checkThat( isEqual( A( 1 ), -1.0 ) );
        testObj.checkThat( isEqual( A( 2 ), -1.5 ) );
        testObj.checkThat( isEqual( A( 3 ), -2.0 ) );
        testObj.checkThat( isEqual( A( 4 ), -2.5 ) );
        testObj.checkThat( isEqual( A( 5 ), -3.0 ) );

        A /= 0.5;

        testObj.checkThat( A.dim() == 6 );

        testObj.checkThat( isEqual( A( 0 ), -1. ) );
        testObj.checkThat( isEqual( A( 1 ), -2. ) );
        testObj.checkThat( isEqual( A( 2 ), -3. ) );
        testObj.checkThat( isEqual( A( 3 ), -4. ) );
        testObj.checkThat( isEqual( A( 4 ), -5. ) );
        testObj.checkThat( isEqual( A( 5 ), -6. ) );

        testObj.reportResults();
    }

    void test_RealVector_CrossProduct()
    {
        Tester testObj( "RealVector cross product" );

        RealVector A( { 1, 2, 3 } );
        RealVector B( { 4, 5, 6 } );

        RealVector C = A.cross( B );

        testObj.checkThat( C.dim() == 3 );
        testObj.checkThat( isEqual( C( 0 ), -3. ) );
        testObj.checkThat( isEqual( C( 1 ), 6. ) );
        testObj.checkThat( isEqual( C( 2 ), -3. ) );

        testObj.reportResults();
    }

    void test_RealVector_DotProduct()
    {
        Tester testObj( "RealVector dot product" );

        RealVector A( { 1, 2, 3, 4 } );
        RealVector B( { 4, 5, 6, 7 } );

        double c = A.dot( B );

        testObj.checkThat( isEqual( c, 60.0 ) );

        testObj.reportResults();
    }

    void test_RealVector_TensorProduct()
    {
        Tester testObj( "RealVector tensor product" );

        RealVector A( { 1, 2, 3, 4 } );
        RealVector B( { 4, 5, 6 } );

        RealMatrix C = A.xMen( B );

        testObj.checkThat( C.dim1() == 4 );
        testObj.checkThat( C.dim2() == 3 );

        testObj.checkThat( isEqual( C( 0,0 ), 4. ) );
        testObj.checkThat( isEqual( C( 0,1 ), 5. ) );
        testObj.checkThat( isEqual( C( 0,2 ), 6. ) );
        testObj.checkThat( isEqual( C( 1,0 ), 8. ) );
        testObj.checkThat( isEqual( C( 1,1 ), 10. ) );
        testObj.checkThat( isEqual( C( 1,2 ), 12. ) );
        testObj.checkThat( isEqual( C( 2,0 ), 12. ) );
        testObj.checkThat( isEqual( C( 2,1 ), 15. ) );
        testObj.checkThat( isEqual( C( 2,2 ), 18. ) );
        testObj.checkThat( isEqual( C( 3,0 ), 16. ) );
        testObj.checkThat( isEqual( C( 3,1 ), 20. ) );
        testObj.checkThat( isEqual( C( 3,2 ), 24. ) );

        testObj.reportResults();
    }

    void test_RealVector_InitAndErase()
    {
        Tester testObj( "RealVector init & erase" );

        RealVector A;
        A.init( 5 );

        testObj.checkThat( A.dim() == 5 );

        for ( int i = 0; i < 3; i++ )
            testObj.checkThat( isEqual( A( i ), 0. ) );

        A.erase();
        testObj.checkThat( A.dim() == 0 );
        testObj.checkThat( A.ptr() == nullptr );

        testObj.reportResults();
    }
}
