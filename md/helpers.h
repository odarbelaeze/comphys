#ifndef HELPERS_H_
#define HELPERS_H_ 1

#include <iostream>
#include <valarray>

namespace md
{
    typedef std::valarray<double> darray;
    typedef std::valarray<darray> ddarray;

    template<typename T>
    T v_mean (const std::valarray<T>& array)
    {
        T sum = array[0];
        sum = 0.0;
        for (auto&& item : array)
            sum = sum + item;
        return sum * (1.0 / array.size());
    }
}

template<typename T>
std::ostream& operator<< (std::ostream& os, const std::valarray<T> array)
{
    for (auto&& item : array)
        os << item << " ";
    return os;
}

#endif