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
        // can be added into utility
        template <typename func_tp, typename tuple_type, size_t... index>
        auto _pass_tuple_as_function_arguments(func_tp f, tuple_type t, std::index_sequence<index...>)
        {
            return f(std::get<index>(t)...);
        }

        template <typename func_tp, typename tuple_type>
        auto _pass_tuple_as_function_arguments(func_tp f, tuple_type t)
        {
            static constexpr auto size = std::tuple_size<tuple_type>::value;
            return _pass_tuple_as_function_arguments(f, t, std::make_index_sequence<size>{});
        }
#ifdef USE_GLOBAL_FLOATING_POINT_TYPE
        template <typename func_tp>
        auto _single_var_auto_diff(func_tp f, math::real x)
        {
            return (f(algebra::dual_number{x, 1.0}) - f(algebra::dual_number{x}).real).dual;
        }
#else
        template <typename var_tp, typename func_tp>
        auto _single_var_auto_diff(func_tp f, var_tp x)
        {
            return (f(algebra::dual_number<var_tp>{x, 1.0}) - f(algebra::dual_number<var_tp>{x}).real).dual;
        }
#endif // USE_GLOBAL_FLOATING_POINT_TYPE
#ifdef USE_GLOBAL_FLOATING_POINT_TYPE
        // calculate differentiation of single variable functions
        template <typename func_tp>
        auto _auto_diff(func_tp f, math::real x)
        {
            return details::_single_var_auto_diff(f, x);
        }
#else // !defined USE_GLOBAL_FLOATING_POINT_TYPE
       // calculate differentiation of single variable functions
        template <typename var_tp, typename func_tp>
        auto _auto_diff(func_tp f, var_tp x)
        {
            return details::_single_var_auto_diff<
                std::conditional_t<
                    std::is_integral<var_tp>::value,
                    math::real,
                    var_tp>>(f, x);
        }
#endif // USE_GLOBAL_FLOATING_POINT_TYPE
#if __cplusplus >= 201703L
        template <int pos, typename func_tp, typename... var_tp>
        auto _auto_diff(func_tp f, var_tp... vars)
        {
            auto var_tuple = std::make_tuple(algebra::dual_number{vars}...);
            std::get<pos>(var_tuple).dual = 1.0f;
            return (std::apply(f, var_tuple) - std::apply(f, var_tuple).real).dual;
        }
#else // C++14
#ifdef USE_GLOBAL_FLOATING_POINT_TYPE
        template <int pos, typename func_tp, typename... var_type>
        auto _auto_diff(func_tp f, var_type... vars)
        {
            auto var_tuple = std::make_tuple(algebra::dual_number{vars}...);
            std::get<pos>(var_tuple).dual = 1.0f;
            auto temp = details::_pass_tuple_as_function_arguments(f, var_tuple);
            return (temp - temp.real).dual;
        }
#else  // !defined USE_GLOBAL_FLOATING_POINT_TYPE
        template <int pos, typename func_tp, typename... var_tp>
        auto _auto_diff(func_tp f, var_tp... vars)
        {
            auto var_tuple = std::make_tuple(algebra::dual_number<decltype(vars)>{vars}...);
            std::get<pos>(var_tuple).dual = 1.0f;
            auto temp = details::_pass_tuple_as_function_arguments(f, var_tuple);
            return (temp - temp.real).dual;
        }
#endif // USE_GLOBAL_FLOATING_POINT_TYPE
#endif // c++17
#if __cplusplus >= 201703L

#else // c++14

#endif // C++17
    }  // namespace math::calculus::details

} // namespace math::calculus

#endif // MATH_CALCULUS_AUTO_DIFF_HPP