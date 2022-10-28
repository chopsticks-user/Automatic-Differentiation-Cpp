#ifndef MATH_CONFIG_HPP
#define MATH_CONFIG_HPP

#if __cplusplus < 201402L
#error "C++14 or higher is required"
#endif // c++14

namespace math
{
    typedef double real;
    typedef int default_integer_type;

} // namespace math

#endif // MATH_CONFIG_HPP