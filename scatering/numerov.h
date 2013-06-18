#ifndef NUMEROV_H_
#define NUMEROV_H_

#include <cmath>
#include <functional>
#include <iostream>

void numerov_rse(std::function<double(double)> V,
    int l, double E, double r_max, double h = 0.001,
    double* r1 = NULL, double* u1 = NULL, 
    double* r2 = NULL, double* u2 = NULL,
    bool debug = true)

{
    double half_lambda = M_PI / std::sqrt(E);
    double r = h;

    // This is a security step for those who wants to use a value
    // for h bellow its optimal value.
    double h_min = 0.001;
    double h_opt = std::max(h_min, h);

    // The value of u for the first stage, it is assumed that there exists
    // a value u_old at the beginin, but that value is zero so...
    double u = std::pow(h_opt, l + 1);

    // The value of w is yet to be computed after the function F is declared.
    double w_old = 0.0;
    double w;
    double w_new;

    auto F = [&, V, l, E](double r){
        return V(r) + l * (l + 1) / (r * r) - E;
    };

    auto u_func = [&, F, h](double r, double w){
        return std::pow(1.0 - (1.0 / 12.0) * h * h * F(r), - 1) * w;
    };

    auto w_func = [&, F, h](double r, double u){
        return (1.0 - (1.0 / 12.0) * h * h * F(r)) * u;
    };

    w = w_func(r, u);

    bool r_max_reached = false;
    bool done = false;

    while (!done)
    {
        if (debug) std::cout << r << "    " << u << std::endl;

        w_new = 2 * w - w_old + h * h * F(r) * u;
        r = r + h;
        w_old = w;
        w = w_new;
        u = u_func(r, w);

        if (r > r_max && !r_max_reached) {
            if (debug) std::cout << "# " << r << "    " << u << std::endl;
            if (r1 != NULL) (*r1) = r;
            if (u1 != NULL) (*u1) = u;
            r_max_reached = true;
        }

        if (r > r_max + half_lambda) {
            if (debug) std::cout << "# " << r << "    " << u << std::endl;
            if (r2 != NULL) (*r2) = r;
            if (u2 != NULL) (*u2) = u;
            done = true;
        }
    }

}

#endif
