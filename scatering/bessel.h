#ifndef BESSEL_H_
#define BESSEL_H_

#include <cmath>

namespace bessel
{
    double j(unsigned l, double x)
    {
        unsigned l_top = std::max(3 * (int) x / 2 + 20, (int)(l + 20));

        double s_pone = 0.0;
        double s_curr = 0.1;
        double s_mone;

        double answer;
        double j_0, j_1;
        double denominator;

        for (int i = l_top; i > 0; --i)
        {
            s_mone = (2 * i + 1) / x * s_curr - s_pone;
            s_pone = s_curr;
            s_curr = s_mone;
            if (i == l) answer = s_pone;
            if (i == 1) j_1 = s_pone;
            if (i == 1) j_0 = s_curr;
        }

        denominator = (j_0 - x * j_1) * std::cos(x) + x * j_0 * std::sin(x);

        return answer / denominator;
    }

    double n(unsigned l, double x)
    {
        if (l == 0) return - std::cos(x) / x;
        if (l == 1) return - std::cos(x) / (x * x) - std::sin(x) / x;
        return (2 * l - 1) / x * n (l - 1, x) - n(l - 2, x);
    }
}

#endif
