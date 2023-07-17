#include "broomstyx_test.hpp"
#include "test_RealMatrix.hpp"
#include "test_RealVector.hpp"
#include "test_linearAlgebra.hpp"

#include <cstdio>

namespace broomstyx
{
    void performTests()
    {
        std::printf( "\n" );
        std::printf( " --------------------------------------------------------------------------------\n" );
        std::printf( "  Running tests ... \n" );
        std::printf( " --------------------------------------------------------------------------------\n\n" );

        test_RealMatrix();
        test_RealVector();
        test_linearAlgebra();

        std::printf( "\n" );
        std::printf( " --------------------------------------------------------------------------------\n" );
        std::printf( "  Tests finished. \n" );
        std::printf( " --------------------------------------------------------------------------------\n\n" );
    }
}
