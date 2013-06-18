#include "numerov.h"
#include "bessel.h"

int main(int argc, char const *argv[])
{
    double r1, r2, u1, u2;

    numerov_rse([](double r) { return (r < 1.0)? - 2.0 : 0.0; },
        1, 0.5, 2.0, 0.001, &r1, &u1, &r2, &u2, false);

    std::cout << "# " << r1 << "   " << u1 << std::endl;
    std::cout << "# " << r2 << "   " << u2 << std::endl;

    std::cout << bessel::j(5, 1.5) << std::endl;
    std::cout << bessel::n(5, 1.5) << std::endl;

    return 0;
}
