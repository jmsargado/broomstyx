#include "test_linearAlgebra.hpp"
#include "Math/linearAlgebra.hpp"
#include "Math/floatingPointComparison.hpp"
#include "Tester.hpp"

namespace broomstyx
{
    void test_linearAlgebra()
    {
        test_matrixAddSubtract();
        test_vectorAddSubtract();
        test_scalarMatrixMultDiv();
        test_scalarVectorMultDiv();
        test_matrixMult();
        test_matVecMult();
        test_vecMatMult();
        test_matrixInverse();
    }

    void test_matrixAddSubtract()
    {
        Tester testObj( "Linear Algebra, matrix add/subtract" );

        RealMatrix A( { { 1, 2, 3 }, { 4, 5, 6 } } );
        RealMatrix B( { { 7, 8, 9 }, { 1, 2, 3 } } );

        RealMatrix D = A + B;

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 8. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 10. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), 12. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 5. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 7. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 9. ) );

        D = ( A + A ) + ( B + B );

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 16. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 20. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), 24. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 10. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 14. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 18. ) );

        D = ( A + B ) + B;

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 15. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 18. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), 21. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 6. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 9. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 12. ) );

        D = A + ( B + B );

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 15. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 18. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), 21. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 6. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 9. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 12. ) );

        D = A - B;

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), -6. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), -6. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), -6. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 3. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 3. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 3. ) );

        D = ( A + B ) - B;

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 1. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 2. ) );
        testObj.checkThat( isEqual( D( 0, 2  ), 3. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 4. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 5. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 6. ) );

        D = A - ( B - B );

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 1. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 2. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), 3. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 4. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 5. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 6. ) );

        D = ( A - A ) - ( B - B );

        testObj.checkThat( D.dim1() == 2 );
        testObj.checkThat( D.dim2() == 3 );

        testObj.checkThat( isEqual( D( 0, 0 ), 0. ) );
        testObj.checkThat( isEqual( D( 0, 1 ), 0. ) );
        testObj.checkThat( isEqual( D( 0, 2 ), 0. ) );
        testObj.checkThat( isEqual( D( 1, 0 ), 0. ) );
        testObj.checkThat( isEqual( D( 1, 1 ), 0. ) );
        testObj.checkThat( isEqual( D( 1, 2 ), 0. ) );

        testObj.reportResults();
    }

    void test_vectorAddSubtract()
    {
        Tester testObj( "Linear algebra, vector add/subtract" );
        RealVector A( { 1, 2, 3 } );
        RealVector B( { 4, 5, 6 } );
        RealVector C;

        C = A + B;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 5. ) );
        testObj.checkThat( isEqual( C( 1 ), 7. ) );
        testObj.checkThat( isEqual( C( 2 ), 9. ) );

        C = ( A + A ) + B;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 6. ) );
        testObj.checkThat( isEqual( C( 1 ), 9. ) );
        testObj.checkThat( isEqual( C( 2 ), 12. ) );

        C = A + ( B + B );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 9. ) );
        testObj.checkThat( isEqual( C( 1 ), 12. ) );
        testObj.checkThat( isEqual( C( 2 ), 15. ) );

        C = ( A + A ) + ( B + B );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 10. ) );
        testObj.checkThat( isEqual( C( 1 ), 14. ) );
        testObj.checkThat( isEqual( C( 2 ), 18. ) );

        C = A - B;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), -3. ) );
        testObj.checkThat( isEqual( C( 1 ), -3. ) );
        testObj.checkThat( isEqual( C( 2 ), -3. ) );

        C = ( A + B ) - B;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 1. ) );
        testObj.checkThat( isEqual( C( 1 ), 2. ) );
        testObj.checkThat( isEqual( C( 2 ), 3. ) );

        C = A - ( B - B );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 1. ) );
        testObj.checkThat( isEqual( C( 1 ), 2. ) );
        testObj.checkThat( isEqual( C( 2 ), 3. ) );

        C = ( A + A ) - ( B - B );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 2. ) );
        testObj.checkThat( isEqual( C( 1 ), 4. ) );
        testObj.checkThat( isEqual( C( 2 ), 6. ) );

        testObj.reportResults();
    }

    void test_scalarMatrixMultDiv()
    {
        Tester testObj( "Linear algebra, scalar-matrix mult/div" );

        RealMatrix A( { { 1, 2, 3 }, { 4, 5, 6 } } );
        RealMatrix B;

        B = 2.0 * A;

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0, 0 ), 2. ) );
        testObj.checkThat( isEqual( B( 0, 1 ), 4. ) );
        testObj.checkThat( isEqual( B( 0, 2 ), 6. ) );
        testObj.checkThat( isEqual( B( 1, 0 ), 8. ) );
        testObj.checkThat( isEqual( B( 1, 1 ), 10. ) );
        testObj.checkThat( isEqual( B( 1, 2 ), 12. ) );

        B = 2.0 * ( A + A );

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0, 0 ), 4. ) );
        testObj.checkThat( isEqual( B( 0, 1 ), 8. ) );
        testObj.checkThat( isEqual( B( 0, 2 ), 12. ) );
        testObj.checkThat( isEqual( B( 1, 0 ), 16. ) );
        testObj.checkThat( isEqual( B( 1, 1 ), 20. ) );
        testObj.checkThat( isEqual( B( 1, 2 ), 24. ) );

        B = A * 2.0;

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0, 0 ), 2. ) );
        testObj.checkThat( isEqual( B( 0, 1 ), 4. ) );
        testObj.checkThat( isEqual( B( 0, 2 ), 6. ) );
        testObj.checkThat( isEqual( B( 1, 0 ), 8. ) );
        testObj.checkThat( isEqual( B( 1, 1 ), 10. ) );
        testObj.checkThat( isEqual( B( 1, 2 ), 12. ) );

        B = ( A + A ) * 2.0;

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0, 0 ), 4. ) );
        testObj.checkThat( isEqual( B( 0, 1 ), 8. ) );
        testObj.checkThat( isEqual( B( 0, 2 ), 12. ) );
        testObj.checkThat( isEqual( B( 1, 0 ), 16. ) );
        testObj.checkThat( isEqual( B( 1, 1 ), 20. ) );
        testObj.checkThat( isEqual( B( 1, 2 ), 24. ) );

        B = A / 0.5;

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0, 0 ), 2. ) );
        testObj.checkThat( isEqual( B( 0, 1 ), 4. ) );
        testObj.checkThat( isEqual( B( 0, 2 ), 6. ) );
        testObj.checkThat( isEqual( B( 1, 0 ), 8. ) );
        testObj.checkThat( isEqual( B( 1, 1 ), 10. ) );
        testObj.checkThat( isEqual( B( 1, 2 ), 12. ) );

        B = ( A + A ) / 0.5;

        testObj.checkThat( B.dim1() == 2 );
        testObj.checkThat( B.dim2() == 3 );

        testObj.checkThat( isEqual( B( 0, 0 ), 4. ) );
        testObj.checkThat( isEqual( B( 0, 1 ), 8. ) );
        testObj.checkThat( isEqual( B( 0, 2 ), 12. ) );
        testObj.checkThat( isEqual( B( 1, 0 ), 16. ) );
        testObj.checkThat( isEqual( B( 1, 1 ), 20. ) );
        testObj.checkThat( isEqual( B( 1, 2 ), 24. ) );

        testObj.reportResults();
    }

    void test_scalarVectorMultDiv()
    {
        Tester testObj( "Linear algebra, scalar-vector mult/div" );

        RealVector A( { 1, 2, 3, 4 } );
        RealVector B;

        B = 2.0 * A;

        testObj.checkThat( B.dim() == 4 );

        testObj.checkThat( isEqual( B( 0 ), 2. ) );
        testObj.checkThat( isEqual( B( 1 ), 4. ) );
        testObj.checkThat( isEqual( B( 2 ), 6. ) );
        testObj.checkThat( isEqual( B( 3 ), 8. ) );

        B = 2.0 * ( A + A );

        testObj.checkThat( B.dim() == 4 );

        testObj.checkThat( isEqual( B( 0 ), 4. ) );
        testObj.checkThat( isEqual( B( 1 ), 8. ) );
        testObj.checkThat( isEqual( B( 2 ), 12. ) );
        testObj.checkThat( isEqual( B( 3 ), 16. ) );

        B = A * 2.0;

        testObj.checkThat( B.dim() == 4 );

        testObj.checkThat( isEqual( B( 0 ), 2. ) );
        testObj.checkThat( isEqual( B( 1 ), 4. ) );
        testObj.checkThat( isEqual( B( 2 ), 6. ) );
        testObj.checkThat( isEqual( B( 3 ), 8. ) );

        B = ( A + A ) * 2.0;

        testObj.checkThat( B.dim() == 4 );

        testObj.checkThat( isEqual( B( 0 ), 4. ) );
        testObj.checkThat( isEqual( B( 1 ), 8. ) );
        testObj.checkThat( isEqual( B( 2 ), 12. ) );
        testObj.checkThat( isEqual( B( 3 ), 16. ) );

        B = A / 0.5;

        testObj.checkThat( B.dim() == 4 );

        testObj.checkThat( isEqual( B( 0 ), 2. ) );
        testObj.checkThat( isEqual( B( 1 ), 4. ) );
        testObj.checkThat( isEqual( B( 2 ), 6. ) );
        testObj.checkThat( isEqual( B( 3 ), 8. ) );

        B = ( A + A ) / 0.5;

        testObj.checkThat( B.dim() == 4 );

        testObj.checkThat( isEqual( B( 0 ), 4. ) );
        testObj.checkThat( isEqual( B( 1 ), 8. ) );
        testObj.checkThat( isEqual( B( 2 ), 12. ) );
        testObj.checkThat( isEqual( B( 3 ), 16. ) );

        testObj.reportResults();
    }

    void test_matrixMult()
    {
        Tester testObj( "Linear algebra, matrix mult" );

        RealMatrix A( { { 1, 2, 3 }, { 4, 5, 6 } } );
        RealMatrix B( { { 6, 3 }, { 5, 2 }, { 4, 1 } } );

        RealMatrix C = A * B;

        testObj.checkThat( C.dim1() == 2 );
        testObj.checkThat( C.dim2() == 2 );

        testObj.checkThat( isEqual( C( 0, 0 ), 28. ) );
        testObj.checkThat( isEqual( C( 0, 1 ), 10. ) );
        testObj.checkThat( isEqual( C( 1, 0 ), 73. ) );
        testObj.checkThat( isEqual( C( 1, 1 ), 28. ) );

        C = trp( A ) * trp( B );

        testObj.checkThat( C.dim1() == 3 );
        testObj.checkThat( C.dim2() == 3 );

        testObj.checkThat( isEqual( C( 0, 0 ), 18. ) );
        testObj.checkThat( isEqual( C( 0, 1 ), 13. ) );
        testObj.checkThat( isEqual( C( 0, 2 ), 8. ) );
        testObj.checkThat( isEqual( C( 1, 0 ), 27. ) );
        testObj.checkThat( isEqual( C( 1, 1 ), 20. ) );
        testObj.checkThat( isEqual( C( 1, 2 ), 13. ) );
        testObj.checkThat( isEqual( C( 2, 0 ), 36. ) );
        testObj.checkThat( isEqual( C( 2, 1 ), 27. ) );
        testObj.checkThat( isEqual( C( 2, 2 ), 18. ) );

        C = A * trp( A );

        testObj.checkThat( C.dim1() == 2 );
        testObj.checkThat( C.dim2() == 2 );

        testObj.checkThat( isEqual( C( 0, 0 ), 14. ) );
        testObj.checkThat( isEqual( C( 0, 1 ), 32. ) );
        testObj.checkThat( isEqual( C( 1, 0 ), 32. ) );
        testObj.checkThat( isEqual( C( 1, 1 ), 77. ) );

        C = trp( A ) * A;

        testObj.checkThat( C.dim1() == 3 );
        testObj.checkThat( C.dim2() == 3 );

        testObj.checkThat( isEqual( C( 0, 0 ), 17. ) );
        testObj.checkThat( isEqual( C( 0, 1 ), 22. ) );
        testObj.checkThat( isEqual( C( 0, 2 ), 27. ) );
        testObj.checkThat( isEqual( C( 1, 0 ), 22. ) );
        testObj.checkThat( isEqual( C( 1, 1 ), 29. ) );
        testObj.checkThat( isEqual( C( 1, 2 ), 36. ) );
        testObj.checkThat( isEqual( C( 2, 0 ), 27. ) );
        testObj.checkThat( isEqual( C( 2, 1 ), 36. ) );
        testObj.checkThat( isEqual( C( 2, 2 ), 45. ) );

        testObj.reportResults();
    }

    void test_matVecMult()
    {
        Tester testObj( "Linear algebra, matrix-vector mult" );

        RealMatrix A( { { 1, 2, 3 }, { 4, 5, 6 } } );
        RealVector B( { -3, -4, 5 } );

        RealVector C = A * B;

        testObj.checkThat( C.dim() == 2 );

        testObj.checkThat( isEqual( C( 0 ), 4. ) );
        testObj.checkThat( isEqual( C( 1 ), -2. ) );

        C = A * ( B + B );

        testObj.checkThat( C.dim() == 2 );

        testObj.checkThat( isEqual( C( 0 ), 8. ) );
        testObj.checkThat( isEqual( C( 1 ), -4. ) );

        C = ( A + A ) * B;

        testObj.checkThat( C.dim() == 2 );

        testObj.checkThat( isEqual( C( 0 ), 8. ) );
        testObj.checkThat( isEqual( C( 1 ), -4. ) );

        C = ( A + A ) * ( B + B );

        testObj.checkThat( C.dim() == 2 );

        testObj.checkThat( isEqual( C( 0 ), 16. ) );
        testObj.checkThat( isEqual( C( 1 ), -8. ) );


        testObj.reportResults();
    }

    void test_vecMatMult()
    {
        Tester testObj( "Linear algebra, vector-matrix mult" );

        RealMatrix A( { { 1, 2, 3 }, { 4, 5, 6 } } );
        RealVector B( { 1, 2 } );

        RealVector C = B * A;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 9. ) );
        testObj.checkThat( isEqual( C( 1 ), 12. ) );
        testObj.checkThat( isEqual( C( 2 ), 15. ) );

        C = B * ( A + A );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 18. ) );
        testObj.checkThat( isEqual( C( 1 ), 24. ) );
        testObj.checkThat( isEqual( C( 2 ), 30. ) );

        C = ( B + B ) * A;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 18. ) );
        testObj.checkThat( isEqual( C( 1 ), 24. ) );
        testObj.checkThat( isEqual( C( 2 ), 30. ) );

        C = ( B + B ) * ( A + A );

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 36. ) );
        testObj.checkThat( isEqual( C( 1 ), 48. ) );
        testObj.checkThat( isEqual( C( 2 ), 60. ) );

        A = { { 5, 4, 3 }, { 2, -5, 6 }, { 4, 7, 8 }, { 1, 9, 2 } };
        B = { 2, 5, 6, -2 };

        C = B * A;

        testObj.checkThat( C.dim() == 3 );

        testObj.checkThat( isEqual( C( 0 ), 42. ) );
        testObj.checkThat( isEqual( C( 1 ), 7. ) );
        testObj.checkThat( isEqual( C( 2 ), 80. ) );

        testObj.reportResults();
    }

    void test_matrixInverse()
    {
        Tester testObj( "Linear algebra, matrix inverse" );

        RealMatrix A( { { 3, 6, 4, 5 }, { 1, 3, 4, 4 }, { 8, 9, 2, 5 }, { 4, 0, 1, -6 } } );
        RealMatrix C = inv( A );

        testObj.checkThat( C.dim1() == 4 );
        testObj.checkThat( C.dim2() == 4 );

        testObj.checkThat( isEqual( C( 0, 0 ), -0.881578947368421 ) );
        testObj.checkThat( isEqual( C( 0, 1 ),  0.697368421052632 ) );
        testObj.checkThat( isEqual( C( 0, 2 ),  0.355263157894737 ) );
        testObj.checkThat( isEqual( C( 0, 3 ),  0.026315789473684 ) );
        testObj.checkThat( isEqual( C( 1, 0 ),  1.118421052631579 ) );
        testObj.checkThat( isEqual( C( 1, 1 ), -0.969298245614035 ) );
        testObj.checkThat( isEqual( C( 1, 2 ), -0.311403508771930 ) );
        testObj.checkThat( isEqual( C( 1, 3 ),  0.026315789473684 ) );
        testObj.checkThat( isEqual( C( 2, 0 ), -0.026315789473684 ) );
        testObj.checkThat( isEqual( C( 2, 1 ),  0.289473684210526 ) );
        testObj.checkThat( isEqual( C( 2, 2 ), -0.078947368421053 ) );
        testObj.checkThat( isEqual( C( 2, 3 ),  0.105263157894737 ) );
        testObj.checkThat( isEqual( C( 3, 0 ), -0.592105263157895 ) );
        testObj.checkThat( isEqual( C( 3, 1 ),  0.513157894736842 ) );
        testObj.checkThat( isEqual( C( 3, 2 ),  0.223684210526316 ) );
        testObj.checkThat( isEqual( C( 3, 3 ), -0.131578947368421 ) );

        testObj.reportResults();
    }
}
