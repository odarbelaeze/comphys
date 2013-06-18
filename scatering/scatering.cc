#include "numerov.h"

int main(int argc, char const *argv[])
{
    numerov_rse([](double r) { return (r < 1.0)? - 2.0 : 0.0; },
        1, 0.5, 2.0);
    return 0;
}
