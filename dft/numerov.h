#ifndef NUEROV_H_
#define NUEROV_H_

#include <cmath>
#include <functional>

namespace numerov
{
    double downward_H (std::function<double(double)> V, double E,
        double r_max, double h = 0.001, bool debug = false)
    {
        // This function solves the radial Schrödinger equation
        // for the Hidrogen atom provided a potential such
        //     - (1 / 2) u''(r) + V(r)u(r) = Eu(r)
        // and a guessed E.

        auto F = [&, V, E](double r) { return 2.0 * (V(r) - E); };

        auto u_func = [&, F, h](double r, double w){
            return std::pow(1.0 - (1.0 / 12.0) * h * h * F(r), - 1) * w;
        };

        auto w_func = [&, F, h](double r, double u){
            return (1.0 - (1.0 / 12.0) * h * h * F(r)) * u;
        };


        // The values of u_curr and u_pone are needed in order to solve the
        // RSE, in this case are casted from the actual solution.
        double u_pone = r_max * std::exp( - r_max);
        double u_curr = (r_max - h) * std::exp( - (r_max - h));

        // The values of w are computed acordingly
        double w_pone = w_func(r_max, u_pone);
        double w_curr = w_func(r_max - h, u_curr);
        double w_mone;

        // Now the iterative process can begin
        double r = r_max - h;

        while(r > 1.0 * h)
        {
            w_mone = 2.0 * w_curr - w_pone + h * h * F(r) * u_curr;
            r = r - h;
            w_pone = w_curr;
            w_curr = w_mone;
            u_curr = u_func(r, w_curr);
            if (debug) std::cout << r << "    " << u_curr << std::endl;
        }


        return u_curr;
    }
}

#endif  