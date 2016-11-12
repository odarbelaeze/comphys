#ifndef NUMEROV_H_
#define NUMEROV_H_

#ifndef FACTOR
#define FACTOR 1.0      // hbar^2 / 2m === 1, for generic potentials
#endif

#ifndef ALPHA
#define ALPHA 6.12      // meV^(-1) * sigma ^ -1 (for LJ potential)
#endif

#ifndef EPSILON
#define EPSILON 5.9     // meV (for LJ potential)
#endif

#include <cmath>
#include <functional>
#include <iostream>

namespace numerov
{

    void rse(std::function<double(double)> V,
        int l, double E, double r_max, double h = 0.001,
        double* r1 = NULL, double* u1 = NULL,
        double* r2 = NULL, double* u2 = NULL,
        bool debug = true)

    {
        double half_lambda = std::sqrt(FACTOR) * M_PI / std::sqrt(E);
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
            return (V(r) + FACTOR * l * (l + 1) / (r * r) - E) / FACTOR;
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

    void rse_lj(int l, double E, double r_max, double h = 0.001,
        double* r1 = NULL, double* u1 = NULL,
        double* r2 = NULL, double* u2 = NULL,
        bool debug = true)

    {
        double half_lambda = M_PI / std::sqrt(E * ALPHA);
        double r = 0.55 + h;

        // This is a security step for those who wants to use a value
        // for h bellow its optimal value.
        double h_min = 0.001;
        double h_opt = std::max(h_min, h);

        // The value of u for the first stage, it is assumed that there exists
        // a value u_old at the beginin, but that value is zero so...
        double C = std::sqrt(EPSILON * ALPHA / 25.0);
        double u_old = std::exp( - C * std::pow(r - h, - 5));
        double u     = std::exp( - C * std::pow(r    , - 5));

        // The value of w is yet to be computed after the function F is declared.
        double w_old;
        double w;
        double w_new;

        auto V = [](double r) {
            return EPSILON * (std::pow(1.0 / r, 12)
                               - 2.0 * std::pow(1.0 / r, 6));
        };

        auto F = [&, V, l, E](double r){
            return (V(r) + l * (l + 1) / (r * r * ALPHA) - E) * ALPHA;
        };

        auto u_func = [&, F, h](double r, double w){
            return std::pow(1.0 - (1.0 / 12.0) * h * h * F(r), - 1) * w;
        };

        auto w_func = [&, F, h](double r, double u){
            return (1.0 - (1.0 / 12.0) * h * h * F(r)) * u;
        };

        w_old = w_func(r - h, u_old);
        w     = w_func(r    , u    );

        bool r_max_reached = false;
        bool done = false;

        while (!done)
        {
            if (debug) std::cout << r << "    "
                                 << u << "    "
                                 << V(r) << "    "
                                 << std::endl;

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

}



#endif
