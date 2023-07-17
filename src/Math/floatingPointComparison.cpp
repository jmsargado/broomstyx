#ifndef FLOATINGPOINTCOMPARISON_HPP
#define FLOATINGPOINTCOMPARISON_HPP

#include "floatingPointComparison.hpp"
#include <cmath>
#include <cstdio>

namespace broomstyx
{
    bool isEqual( double x, double y )
    {
        const double tol = 1.0e-15;

        if ( std::fabs( x - y ) <= tol * std::max( 1.0, std::max( std::fabs( x ), std::fabs( y ) ) ) )
            return true;
        else
            return false;
    }
}

#endif // FLOATINGPOINTCOMPARISON_HPP
