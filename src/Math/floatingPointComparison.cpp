#ifndef FLOATINGPOINTCOMPARISON_HPP
#define FLOATINGPOINTCOMPARISON_HPP

#include "floatingPointComparison.hpp"
#include <limits>
#include <cmath>

namespace broomstyx
{
    bool isEqual(double x, double y)
    {
        const double tol = std::numeric_limits<double>::epsilon();
        if (std::fabs(x - y) <= tol * std::max(1.0, std::max(std::fabs(x), std::fabs(y))))
            return true;
        else
            return false;
    }
}

#endif // FLOATINGPOINTCOMPARISON_HPP
