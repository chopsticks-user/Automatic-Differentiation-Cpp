#ifndef MATH_HPP
#define MATH_HPP

#include "Config.hpp"

#if __cplusplus >= 201402L
#include "Calculus.hpp"

#define make_math_function(variable, function) [](auto(variable)) { return (function); }

namespace math
{

} // namespace math
#endif // c++14

#endif // MATH_HPP