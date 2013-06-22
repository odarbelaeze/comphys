#ifndef BISECTION_H_
#define BISECTION_H_

#include <cmath>
#include <functional>

namespace root
{
    double pure_bisection (std::function<double(double)> F, double a, double b)
    {
        double mid = a + (b - a) / 2.0;
        double f_mid = F(mid);
        if (f_mid < 1e-15) return mid;
        if (f_mid * F(a) < 0.0) return pure_bisection(F, a, mid);
        else                    return pure_bisection(F, mid, b);
    }

    double bisection(std::function<double(double)> f, double a, double b, int segs = 5)
    {
        // This function looks for a root of the function f
        // in the interval @param a, @param b using the bisection 
        // method first dividing the interval in @param segs
        // and lookigng for the first sign change

        double range = std::abs(b - a);
        double h_seg = range / segs;
        double min = std::min(a, b);
        double max;

        bool sign_change = false;

        for (int i = 0; i < segs; ++i)
        {
            if (f(min) * f(min + h_seg) < 0)
            {
                max = min + h_seg;
                sign_change = true;
            }
            min = min + h_seg;
        }

        if (!sign_change) return std::min(a, b) - 1;
        else              return pure_bisection(f, min, max);

    }


}

#endif
