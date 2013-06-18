#ifndef BESSEL_H_
#define BESSEL_H_

#include <cmath>

namespace bessel
{
    double j(unsigned l, double x)
    {
        if (l == 0) return std::sin(x) / x;
        if (l == 1) return std::sin(x) / (x * x) - std::cos(x) / x;
        return (2 * l - 1) / x * j(l - 1, x) - j(l - 2, x);
    }

    double n(unsigned l, double x)
    {
        if (l == 0) return std::cos(x) / x;
        if (l == 1) return std::cos(x) / (x * x) - std::sin(x) / x;
        return (2 * l - 1) / x * n (l - 1, x) - n(l - 2, x);
    }
}

#endif
