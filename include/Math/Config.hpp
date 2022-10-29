#ifndef MATH_CONFIG_HPP
#define MATH_CONFIG_HPP

#if __cplusplus < 201402L
#error "C++14 or higher is required"
#else

// #define USE_GLOBAL_FLOATING_POINT_TYPE

namespace math
{
    typedef double real;
    typedef unsigned int size_type;

} // namespace math

#endif // c++14
#endif // MATH_CONFIG_HPP