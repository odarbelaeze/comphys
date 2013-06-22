#include <iostream>

#include "bisection.h"
#include "numerov.h"

int main(int argc, char const *argv[])
{
    for (double E = - 0.55; E < -0.2; E += 0.01)
    {
        double u = numerov::downward_H ([](double r){ return - 1.0 / r; }, E, 25.0);
        std::cout << E << "    " << u << std::endl;
    }

    double root_E = root::bisection([](double E){ 
        return numerov::downward_H ([](double r){ return - 1.0 / r; }, E, 25.0); 
    }, -0.52, - 0.46);

    std::cout << "# " << root_E << std::endl;    

    return 0;
}
