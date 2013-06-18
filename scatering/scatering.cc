#include <cmath>
#include <iomanip>
#include <iostream>

#include "numerov.h"
#include "bessel.h"

double total_cross_section_lj (double E)
{
    double k = std::sqrt(E * ALPHA);
    double r1, r2, u1, u2;
    unsigned l_max = 7;

    double K;
    double tdl;
    double delta_l;

    double sigma_tot = 0.0;

    for (int l = 0; l < l_max; ++l)
    {
        numerov::rse_lj(l, E, 5.0, 0.001, &r1, &u1, &r2, &u2, false);
        K = (r1 * u2) / (r2 * u1);
        tdl = (K * bessel::j(l, k * r1) - bessel::j(l, k * r2)) /
              (K * bessel::n(l, k * r1) - bessel::n(l, k * r2)) ;
        delta_l = std::atan(tdl);
        sigma_tot += (2 * l + 1) * std::pow(std::sin(delta_l), 2);
    }

    return 4 * M_PI / (k * k) * sigma_tot;
}


int main(int argc, char const *argv[])
{
    std::cout << std::setprecision(9) << std::scientific;

    for (double E = 0.2; E < 3.5; E += 0.02)
        std::cout << E << "    " << total_cross_section_lj(E) << std::endl; 

    return 0;
}
