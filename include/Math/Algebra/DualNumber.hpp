#ifndef MATH_ALGEBRA_DUAL_NUMBER_HPP
#define MATH_ALGEBRA_DUAL_NUMBER_HPP

#include "Config.hpp"

#include <iostream>
#include <type_traits>
#include <cmath>

namespace math
{
    namespace algebra
    {
        template <typename value_type = math::real>
        struct dual_number
        {
            static_assert(std::is_floating_point<value_type>::value);

            value_type real;
            value_type dual;

            using type = value_type;

            dual_number() = default;
            dual_number(value_type r, value_type d = value_type{}) : real{r}, dual(d) {}
            dual_number(const dual_number &rhs) = default;
            dual_number(dual_number &&rhs) = default;
            dual_number &operator=(const dual_number &rhs) = default;
            dual_number &operator=(dual_number &&rhs) = default;
            ~dual_number() = default;

            template <typename rhs_value_type>
            dual_number operator+(dual_number<rhs_value_type> rhs)
            {
                return dual_number{real + rhs.real, dual + rhs.dual};
            }

            dual_number operator+(value_type scalar)
            {
                return dual_number{real + scalar, dual};
            }

            friend dual_number operator+(value_type scalar, dual_number num)
            {
                return dual_number{num.real + scalar, num.dual};
            }

            template <typename rhs_value_type>
            dual_number operator-(dual_number<rhs_value_type> rhs)
            {
                return dual_number{real - rhs.real, dual - rhs.dual};
            }

            dual_number operator-(value_type scalar)
            {
                return dual_number{real - scalar, dual};
            }

            template <typename rhs_value_type>
            friend dual_number operator-(value_type scalar, dual_number<rhs_value_type> d_num)
            {
                return dual_number{scalar - d_num.real, -d_num.dual};
            }

            template <typename rhs_value_type>
            dual_number operator*(dual_number<rhs_value_type> rhs)
            {
                return dual_number{real * rhs.real, real * rhs.dual + dual * rhs.real};
            }

            dual_number operator*(value_type rhs)
            {
                return dual_number{real * rhs, dual * rhs};
            }

            friend dual_number operator*(value_type scalar, dual_number d_num)
            {
                return dual_number{d_num.real * scalar, d_num.dual * scalar};
            }

            friend std::ostream &operator<<(std::ostream &os, dual_number d_num)
            {
                return os << "(" << d_num.real << ", " << d_num.dual << ")";
            }
        };
    } // namespace math::algebra

    // x.real != 0
    template <typename var_type>
    algebra::dual_number<var_type> abs(algebra::dual_number<var_type> x)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::abs<dual_number>");
        return algebra::dual_number<var_type>{
            std::abs(x.real),
            x.dual * x.real / std::abs(x.real)};
    }

    // power group

    template <typename var_type>
    algebra::dual_number<var_type> sq(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            x.real * x.real,
            x.dual * 2.0 * x.real};
    }

    template <typename var_type>
    algebra::dual_number<var_type> cb(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            x.real * x.real * x.real,
            x.dual * 3.0 * x.real * x.real};
    }

    template <typename var_type>
    algebra::dual_number<var_type> sqrt(algebra::dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::sqrt<dual_number>");
        auto sqrt_xr = std::sqrt(x.real);
        return algebra::dual_number<var_type>{
            sqrt_xr,
            x.dual * (0.5 / sqrt_xr)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> cbrt(algebra::dual_number<var_type> x)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::cbrt<dual_number>");
        auto cbrt_xr = std::cbrt(x.real);
        return algebra::dual_number<var_type>{
            cbrt_xr,
            x.dual / (3.0 * cbrt_xr * cbrt_xr))};
    }

    // x^x or f(x) ^ f(x)
    template <typename var_type>
    algebra::dual_number<var_type> pow(algebra::dual_number<var_type> x, var_type p)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::pow_x_n<dual_number>");
        auto pow_xr = std::pow(x.real, p);
        return algebra::dual_number<var_type>{
            pow_xr,
            p * x.dual * pow_xr / x.real};
    }

    // exponential and logarithmic group

    template <typename var_type>
    algebra::dual_number<var_type> pow(algebra::dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::pow_x_x<dual_number>");
        auto xr_pow_xr = std::pow(x.real, x.real);
        return algebra::dual_number<var_type>{
            xr_pow_xr,
            x.dual * xr_pow_xr * (1 + std::log(x.real))};
    }

    // pow f(x)^g(x)
    // template <typename var_type>
    // algebra::dual_number<var_type> pow(algebra::dual_number<var_type> x, algebra::dual_number<var_type>g_x)
    // {
    //     auto xr_pow_xr = std::pow(x.real, x.real);
    //     return algebra::dual_number<var_type>{
    //         xr_pow_xr,
    //         x.dual * xr_pow_xr * (1 + std::log(x.real))};
    // }

    template <typename var_type>
    algebra::dual_number<var_type> exp(algebra::dual_number<var_type> x)
    {
        auto exp_xr = std::exp(x.real);
        return algebra::dual_number<var_type>{
            exp_xr,
            x.dual * exp_xr};
    }

    template <typename var_type>
    algebra::dual_number<var_type> exp_n(var_type n, algebra::dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::exp_n_x<dual_number>");
        auto exp_n_xr = std::pow(n, x.real);
        return algebra::dual_number<var_type>{
            exp_n_xr,
            x.dual * std::log(n) * exp_n_xr};
    }

    template <typename var_type>
    algebra::dual_number<var_type> log(algebra::dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::ln<dual_number>");
        return algebra::dual_number<var_type>{
            std::log(x.real),
            x.dual / x.real};
    }

    template <typename var_type>
    algebra::dual_number<var_type> ln(algebra::dual_number<var_type> x)
    {
        return math::log(x);
    }

    // n > 1
    template <typename var_type>
    algebra::dual_number<var_type> log_n(var_type n, algebra::dual_number<var_type> x)
    {
        if (n <= 0.0 || n == 1.0)
            throw std::runtime_error("n <= 0 || n = 1 at math::log_n_x<dual_number>");
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::log_n_x<dual_number>");
        auto ln_n = std::log(n);
        return algebra::dual_number<var_type>{
            std::log(x.real) / ln_n,
            x.dual / (x.real * ln_n)};
    }

    // x.real > 1, log_x_n
    template <typename var_type>
    algebra::dual_number<var_type> log_x_n(algebra::dual_number<var_type> x, var_type n)
    {
        if (n <= 0.0)
            throw std::runtime_error("n <= 0 at math::log_x_n<dual_number>");
        if (x.real <= 0.0 || x.real == 1.0)
            throw std::runtime_error("x.real <= 0 || x.real = 1 at math::log_x_n<dual_number>");
        auto ln_n = std::log(n);
        auto ln_x = std::log(x.real);
        return algebra::dual_number<var_type>{
            ln_n / ln_x,
            x.dual * (-ln_n / (x.real * ln_x * ln_x))};
    }

    // trigonometric group

    template <typename var_type>
    algebra::dual_number<var_type> sin(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::sin(x.real),
            x.dual * (std::cos(x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> cos(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::cos(x.real),
            x.dual * (-std::sin(x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> tan(algebra::dual_number<var_type> x)
    {
        // if (std::cos(x.real) == 0.0)
        //     throw std::runtime_error("x.real == k(pi / 2) at math::tan<dual_number>");
        return algebra::dual_number<var_type>{
            std::tan(x.real),
            x.dual * (2.0 / (1.0 + std::cos(2.0 * x.real)))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> cot(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            1.0 / std::tan(x.real),
            x.dual * (2.0 / (std::cos(2.0 * x.real) - 1))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> sec(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            1.0 / std::cos(x.real),
            x.dual * (std::tan(x.real) / std::cos(x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> csc(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            1.0 / std::sin(x.real),
            x.dual * (-1 / (std::sin(x.real) * std::tan(x.real)))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> asin(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::asin(x.real),
            x.dual / (std::sqrt(1 - x.real * x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> acos(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::acos(x.real),
            -x.dual / (std::sqrt(1.0 - x.real * x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> atan(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::atan(x.real),
            x.dual / (1.0 + x.real * x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> acot(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::atan(1.0 / x.real),
            x.dual / (-1.0 - x.real * x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> asec(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::acos(1.0 / x.real),
            x.dual / (std::abs(x.real) * std::sqrt(x.real * x.real - 1))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> acsc(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::acos(1.0 / x.real),
            -x.dual / (std::abs(x.real) * std::sqrt(x.real * x.real - 1))};
    }

    // hyperbolic group

    template <typename var_type>
    algebra::dual_number<var_type> sinh(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::sinh(x.real),
            x.dual * std::cosh(x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> cosh(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::cosh(x.real),
            x.dual * std::sinh(x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> tanh(algebra::dual_number<var_type> x)
    {
        auto tanh_xr = std::tanh(x.real);
        return algebra::dual_number<var_type>{
            tanh_xr,
            x.dual * (1.0 - tanh_xr * tanh_xr)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> coth(algebra::dual_number<var_type> x)
    {
        auto coth_xr = 1.0 / std::tanh(x.real);
        return algebra::dual_number<var_type>{
            coth_xr,
            x.dual * (1.0 - coth_xr * coth_xr)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> sech(algebra::dual_number<var_type> x)
    {
        auto sech_xr = 1.0 / std::cosh(x.real);
        return algebra::dual_number<var_type>{
            sech_xr,
            x.dual * (-sech_xr * std::tanh(x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> csch(algebra::dual_number<var_type> x)
    {
        auto csch_xr = 1.0 / std::sinh(x.real);
        return algebra::dual_number<var_type>{
            csch_xr,
            x.dual * (-csch_xr / std::tanh(x.real))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> asinh(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::asinh(x.real),
            x.dual * std::sqrt(1.0 + x.real * x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> acosh(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::acosh(x.real),
            x.dual * std::sqrt(x.real * x.real - 1.0)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> atanh(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::atanh(x.real),
            x.dual * (1.0 - x.real * x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> acoth(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::atanh(1.0 / x.real),
            x.dual * (1.0 - x.real * x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> asech(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::acosh(1.0 / x.real),
            -x.dual * x.real * std::sqrt(1.0 - x.real * x.real)};
    }

    template <typename var_type>
    algebra::dual_number<var_type> acsch(algebra::dual_number<var_type> x)
    {
        return algebra::dual_number<var_type>{
            std::acosh(1.0 / x.real),
            -x.dual * std::abs(x.real) * std::sqrt(1.0 + x.real * x.real)};
    }

    // miscellaneous group

} // namespace math

#endif // MATH_ALGEBRA_DUAL_NUMBER_HPP