#ifndef MATH_CALCULUS_AUTO_DIFF_HPP
#define MATH_CALCULUS_AUTO_DIFF_HPP

#include "Config.hpp"

#include "Algebra/DualNumber.hpp"

#include <iostream>
#include <type_traits>
#include <tuple>

namespace math::calculus
{
    namespace details
    {
        template <typename var_tp, typename func_type>
        auto _single_var_auto_diff(func_type f, var_tp x)
        {
            return (f(algebra::dual_number<var_tp>{x, 1.0}) - f(algebra::dual_number<var_tp>{x}).real).dual;
        }
    } // namespace math::calculus::details

    template <typename var_tp, typename func_type>
    auto auto_diff(func_type f, var_tp x)
    {
        return details::_single_var_auto_diff<
            std::conditional_t<
                std::is_integral<var_tp>::value,
                math::real,
                var_tp>>(f, x); 
    }

#if __cplusplus >= 201703L
    template <int pos, typename func_type, typename... var_tp>
    auto auto_diff(func_type f, var_tp... vars)
    {
        auto var_tuple = std::make_tuple(algebra::dual_number{vars}...);
        std::get<pos>(var_tuple).dual = 1.0f;
        // return (std::apply(f, var_tuple) - algebra::dual_number{std::apply(f, var_tuple).real}).dual;
        return (std::apply(f, var_tuple) - std::apply(f, var_tuple).real).dual;
    }
#else

#endif // c++17

} // namespace math::calculus

#endif // MATH_CALCULUS_AUTO_DIFF_HPP