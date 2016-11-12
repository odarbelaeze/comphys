#include "argon.h"

ArgonSystem::ArgonSystem(int L, double T)
{

    N_ = 4 * L * L * L;
    T_ = T;

    vdarray basis({
        darray({ 0.0, 0.0, 0.0 }),
        darray({ 0.5, 0.5, 0.0 }),
        darray({ 0.0, 0.5, 0.5 }),
        darray({ 0.5, 0.0, 0.5 })
    });

    r_ = vdarray(N_);

    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            for (int k = 0; k < L; ++k)
            {
                for (int b = 0; b < 4; ++b)
                {
                    r_[i * L + j * L + k * 4 + b] = darray({
                        (double) i,
                        (double) j,
                        (double) k
                    }) + basis[b];
                }
            }
        }
    }

    std::mt19937 mte;
    std::normal_distribution<> uniform(0.0, 0.5);

    v_ = vdarray(N_);
    for (int i = 0; i < N_; ++i)
    {
        v_[i] = uniform(mte);
    }

    rescale_speed();
}

ArgonSystem::~ArgonSystem() {}

void ArgonSystem::rescale_speed()
{
    darray v_tot = v_.sum();
    darray v_prom = v_tot * (1.0 / N_);
    v_ = v_ + v_prom;

    std::function<double(darray)> F = [](darray da) -> double {
        return (da * da).sum();
    };

    double lambda = std::sqrt(
        (3.0 * (N_ - 1) * T_) /
        (helpers::apply(v_, F)).sum()
    );

}

