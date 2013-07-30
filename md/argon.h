#ifndef ARGON_H_
#define ARGON_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <valarray>
#include <vector>

typedef std::valarray<double> darray;
typedef std::valarray<darray> vdarray;


class ArgonSystem
{

public:
    ArgonSystem(int L, double T);
    ~ArgonSystem();

private:
    double  T_;

    int     N_;
    vdarray r_;
    vdarray v_;
    vdarray L_;

    void rescale_speed ();

};

namespace helpers
{
    template<typename T_IN, typename T_OUT>
    std::valarray<T_OUT> apply (
        std::valarray<T_IN> array,
        std::function<T_OUT(T_IN)> F
    )
    {
        std::valarray<T_OUT> return_value(array.size());
        for (int i = 0; i < array.size(); ++i)
        {
            return_value[i] = F(array[i]);
        }
        return return_value;
    }
}

#endif
