#ifndef MATH_ALGEBRA_FO_AUTO_DIFF_HPP
#define MATH_ALGEBRA_FO_AUTO_DIFF_HPP

#include "Config.hpp"

#include <iostream>
#include <type_traits>
#include <cmath>

namespace math
{
    namespace calculus::details
    {
#ifdef USE_GLOBAL_FLOATING_POINT_TYPE
        struct _dual_number
        {
            math::real real;
            math::real dual;

            // using type = math::real;

            _dual_number() = default;
            _dual_number(math::real r, math::real d = math::real{}) : real{r}, dual(d) {}
            _dual_number(const _dual_number &rhs) = default;
            _dual_number(_dual_number &&rhs) = default;
            _dual_number &operator=(const _dual_number &rhs) = default;
            _dual_number &operator=(_dual_number &&rhs) = default;
            ~_dual_number() = default;

            _dual_number operator+(_dual_number rhs)
            {
                return _dual_number{real + rhs.real, dual + rhs.dual};
            }

            _dual_number operator+(math::real scalar)
            {
                return _dual_number{real + scalar, dual};
            }

            friend _dual_number operator+(math::real scalar, _dual_number num)
            {
                return _dual_number{num.real + scalar, num.dual};
            }

            _dual_number operator-(_dual_number rhs)
            {
                return _dual_number{real - rhs.real, dual - rhs.dual};
            }

            _dual_number operator-(math::real scalar)
            {
                return _dual_number{real - scalar, dual};
            }

            friend _dual_number operator-(math::real scalar, _dual_number d_num)
            {
                return _dual_number{scalar - d_num.real, -d_num.dual};
            }

            _dual_number operator*(_dual_number rhs)
            {
                return _dual_number{real * rhs.real, real * rhs.dual + dual * rhs.real};
            }

            _dual_number operator*(math::real rhs)
            {
                return _dual_number{real * rhs, dual * rhs};
            }

            friend _dual_number operator*(math::real scalar, _dual_number d_num)
            {
                return _dual_number{d_num.real * scalar, d_num.dual * scalar};
            }

            _dual_number operator/(_dual_number rhs)
            {
                return _dual_number{
                    real / rhs.real,
                    (dual * rhs.real - real * rhs.dual) / (rhs.real * rhs.real)};
            }

            _dual_number operator/(math::real rhs)
            {
                return _dual_number{real / rhs, dual / rhs};
            }

            friend _dual_number operator/(math::real scalar, _dual_number d_num)
            {
                return _dual_number{
                    scalar / d_num.real,
                    -scalar * d_num.dual / (d_num.real * d_num.real)};
            }

            friend std::ostream &operator<<(std::ostream &os, _dual_number d_num)
            {
                return os << "(" << d_num.real << ", " << d_num.dual << ")";
            }
        };
#else  // !defined USE_GLOBAL_FLOATING_POINT_TYPE
        template <typename value_type = math::real>
        struct _dual_number
        {
            static_assert(std::is_floating_point<value_type>::value);

            value_type real;
            value_type dual;

            using type = value_type;

            _dual_number(value_type r, value_type d = value_type{}) : real{r}, dual(d) {}

            _dual_number() = default;
            _dual_number(const _dual_number &rhs) = default;
            _dual_number(_dual_number &&rhs) = default;
            _dual_number &operator=(const _dual_number &rhs) = default;
            _dual_number &operator=(_dual_number &&rhs) = default;
            ~_dual_number() = default;

            template <typename rhs_value_type>
            _dual_number operator+(_dual_number<rhs_value_type> rhs)
            {
                return _dual_number{real + rhs.real, dual + rhs.dual};
            }

            _dual_number operator+(value_type scalar)
            {
                return _dual_number{real + scalar, dual};
            }

            friend _dual_number operator+(value_type scalar, _dual_number num)
            {
                return _dual_number{num.real + scalar, num.dual};
            }

            template <typename rhs_value_type>
            _dual_number operator-(_dual_number<rhs_value_type> rhs)
            {
                return _dual_number{real - rhs.real, dual - rhs.dual};
            }

            _dual_number operator-(value_type scalar)
            {
                return _dual_number{real - scalar, dual};
            }

            template <typename rhs_value_type>
            friend _dual_number operator-(value_type scalar, _dual_number<rhs_value_type> d_num)
            {
                return _dual_number{scalar - d_num.real, -d_num.dual};
            }

            template <typename rhs_value_type>
            _dual_number operator*(_dual_number<rhs_value_type> rhs)
            {
                return _dual_number{real * rhs.real, real * rhs.dual + dual * rhs.real};
            }

            _dual_number operator*(value_type rhs)
            {
                return _dual_number{real * rhs, dual * rhs};
            }

            friend _dual_number operator*(value_type scalar, _dual_number d_num)
            {
                return _dual_number{d_num.real * scalar, d_num.dual * scalar};
            }

            template <typename rhs_value_type>
            _dual_number operator/(_dual_number<rhs_value_type> rhs)
            {
                return _dual_number{
                    real / rhs.real,
                    (dual * rhs.real - real * rhs.dual) / (rhs.real * rhs.real)};
            }

            _dual_number operator/(math::real rhs)
            {
                return _dual_number{real / rhs, dual / rhs};
            }

            friend _dual_number operator/(math::real scalar, _dual_number d_num)
            {
                return _dual_number{
                    scalar / d_num.real,
                    -scalar * d_num.dual / (d_num.real * d_num.real)};
            }

            friend std::ostream &operator<<(std::ostream &os, _dual_number d_num)
            {
                return os << "(" << d_num.real << ", " << d_num.dual << ")";
            }
        };
#endif // USE_GLOBAL_FLOATING_POINT_TYPE
    }  // namespace math::algebra
#ifndef USE_GLOBAL_FLOATING_POINT_TYPE
    // x.real != 0
    template <typename var_type>
    calculus::details::_dual_number<var_type> abs(calculus::details::_dual_number<var_type> x)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::abs<_dual_number>");
        return calculus::details::_dual_number<var_type>{
            std::abs(x.real),
            x.dual * x.real / std::abs(x.real)};
    }

    // power group

    template <typename var_type>
    calculus::details::_dual_number<var_type> sq(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            x.real * x.real,
            x.dual * 2.0 * x.real};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> cb(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            x.real * x.real * x.real,
            x.dual * 3.0 * x.real * x.real};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> sqrt(calculus::details::_dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::sqrt<_dual_number>");
        auto sqrt_xr = std::sqrt(x.real);
        return calculus::details::_dual_number<var_type>{
            sqrt_xr,
            x.dual * (0.5 / sqrt_xr)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> cbrt(calculus::details::_dual_number<var_type> x)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::cbrt<_dual_number>");
        auto cbrt_xr = std::cbrt(x.real);
        return calculus::details::_dual_number<var_type>{
            cbrt_xr,
            x.dual / (3.0 * cbrt_xr * cbrt_xr)};
    }

    // x^x or f(x) ^ f(x)
    template <typename var_type>
    calculus::details::_dual_number<var_type> pow(calculus::details::_dual_number<var_type> x, var_type p)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::pow_x_n<_dual_number>");
        auto pow_xr = std::pow(x.real, p);
        return calculus::details::_dual_number<var_type>{
            pow_xr,
            p * x.dual * pow_xr / x.real};
    }

    // exponential and logarithmic group

    template <typename var_type>
    calculus::details::_dual_number<var_type> pow(calculus::details::_dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::pow_x_x<_dual_number>");
        auto xr_pow_xr = std::pow(x.real, x.real);
        return calculus::details::_dual_number<var_type>{
            xr_pow_xr,
            x.dual * xr_pow_xr * (1 + std::log(x.real))};
    }

    // pow f(x)^g(x)
    // template <typename var_type>
    // calculus::details::_dual_number<var_type> pow(calculus::details::_dual_number<var_type> x, calculus::details::_dual_number<var_type>g_x)
    // {
    //     auto xr_pow_xr = std::pow(x.real, x.real);
    //     return calculus::details::_dual_number<var_type>{
    //         xr_pow_xr,
    //         x.dual * xr_pow_xr * (1 + std::log(x.real))};
    // }

    template <typename var_type>
    calculus::details::_dual_number<var_type> exp(calculus::details::_dual_number<var_type> x)
    {
        auto exp_xr = std::exp(x.real);
        return calculus::details::_dual_number<var_type>{
            exp_xr,
            x.dual * exp_xr};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> exp_n(var_type n, calculus::details::_dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::exp_n_x<_dual_number>");
        auto exp_n_xr = std::pow(n, x.real);
        return calculus::details::_dual_number<var_type>{
            exp_n_xr,
            x.dual * std::log(n) * exp_n_xr};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> log(calculus::details::_dual_number<var_type> x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::ln<_dual_number>");
        return calculus::details::_dual_number<var_type>{
            std::log(x.real),
            x.dual / x.real};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> ln(calculus::details::_dual_number<var_type> x)
    {
        return math::log(x);
    }

    // n > 1
    template <typename var_type>
    calculus::details::_dual_number<var_type> log_n(var_type n, calculus::details::_dual_number<var_type> x)
    {
        if (n <= 0.0 || n == 1.0)
            throw std::runtime_error("n <= 0 || n = 1 at math::log_n_x<_dual_number>");
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::log_n_x<_dual_number>");
        auto ln_n = std::log(n);
        return calculus::details::_dual_number<var_type>{
            std::log(x.real) / ln_n,
            x.dual / (x.real * ln_n)};
    }

    // x.real > 1, log_x_n
    template <typename var_type>
    calculus::details::_dual_number<var_type> log_x_n(calculus::details::_dual_number<var_type> x, var_type n)
    {
        if (n <= 0.0)
            throw std::runtime_error("n <= 0 at math::log_x_n<_dual_number>");
        if (x.real <= 0.0 || x.real == 1.0)
            throw std::runtime_error("x.real <= 0 || x.real = 1 at math::log_x_n<_dual_number>");
        auto ln_n = std::log(n);
        auto ln_x = std::log(x.real);
        return calculus::details::_dual_number<var_type>{
            ln_n / ln_x,
            x.dual * (-ln_n / (x.real * ln_x * ln_x))};
    }

    // trigonometric group

    template <typename var_type>
    calculus::details::_dual_number<var_type> sin(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::sin(x.real),
            x.dual * (std::cos(x.real))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> cos(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::cos(x.real),
            x.dual * (-std::sin(x.real))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> tan(calculus::details::_dual_number<var_type> x)
    {
        // if (std::cos(x.real) == 0.0)
        //     throw std::runtime_error("x.real == k(pi / 2) at math::tan<_dual_number>");
        auto tan_xr = std::tan(x.real);
        return calculus::details::_dual_number<var_type>{
            tan_xr,
            x.dual * (1.0 + tan_xr * tan_xr)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> cot(calculus::details::_dual_number<var_type> x)
    {
        auto cot_xr = 1.0 / std::tan(x.real);
        return calculus::details::_dual_number<var_type>{
            cot_xr,
            x.dual * (-1.0 - cot_xr * cot_xr)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> sec(calculus::details::_dual_number<var_type> x)
    {
        auto cos_xr = std::cos(x.real);
        return calculus::details::_dual_number<var_type>{
            1.0 / cos_xr,
            x.dual * (std::tan(x.real) / cos_xr)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> csc(calculus::details::_dual_number<var_type> x)
    {
        auto sin_xr = std::sin(x.real);
        return calculus::details::_dual_number<var_type>{
            1.0 / sin_xr,
            x.dual * (-1.0 / (sin_xr * std::tan(x.real)))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> asin(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::asin(x.real),
            x.dual / (std::sqrt(1 - x.real * x.real))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> acos(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::acos(x.real),
            -x.dual / (std::sqrt(1.0 - x.real * x.real))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> atan(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::atan(x.real),
            x.dual / (1.0 + x.real * x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> acot(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::atan(1.0 / x.real),
            x.dual / (-1.0 - x.real * x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> asec(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::acos(1.0 / x.real),
            x.dual / (std::abs(x.real) * std::sqrt(x.real * x.real - 1))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> acsc(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::acos(1.0 / x.real),
            -x.dual / (std::abs(x.real) * std::sqrt(x.real * x.real - 1))};
    }

    // hyperbolic group

    template <typename var_type>
    calculus::details::_dual_number<var_type> sinh(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::sinh(x.real),
            x.dual * std::cosh(x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> cosh(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::cosh(x.real),
            x.dual * std::sinh(x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> tanh(calculus::details::_dual_number<var_type> x)
    {
        auto tanh_xr = std::tanh(x.real);
        return calculus::details::_dual_number<var_type>{
            tanh_xr,
            x.dual * (1.0 - tanh_xr * tanh_xr)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> coth(calculus::details::_dual_number<var_type> x)
    {
        auto coth_xr = 1.0 / std::tanh(x.real);
        return calculus::details::_dual_number<var_type>{
            coth_xr,
            x.dual * (1.0 - coth_xr * coth_xr)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> sech(calculus::details::_dual_number<var_type> x)
    {
        auto sech_xr = 1.0 / std::cosh(x.real);
        return calculus::details::_dual_number<var_type>{
            sech_xr,
            x.dual * (-sech_xr * std::tanh(x.real))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> csch(calculus::details::_dual_number<var_type> x)
    {
        auto csch_xr = 1.0 / std::sinh(x.real);
        return calculus::details::_dual_number<var_type>{
            csch_xr,
            x.dual * (-csch_xr / std::tanh(x.real))};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> asinh(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::asinh(x.real),
            x.dual * std::sqrt(1.0 + x.real * x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> acosh(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::acosh(x.real),
            x.dual * std::sqrt(x.real * x.real - 1.0)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> atanh(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::atanh(x.real),
            x.dual * (1.0 - x.real * x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> acoth(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::atanh(1.0 / x.real),
            x.dual * (1.0 - x.real * x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> asech(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::acosh(1.0 / x.real),
            -x.dual * x.real * std::sqrt(1.0 - x.real * x.real)};
    }

    template <typename var_type>
    calculus::details::_dual_number<var_type> acsch(calculus::details::_dual_number<var_type> x)
    {
        return calculus::details::_dual_number<var_type>{
            std::acosh(1.0 / x.real),
            -x.dual * std::abs(x.real) * std::sqrt(1.0 + x.real * x.real)};
    }

    // miscellaneous group

#else // !defined USE_GLOBAL_FLOATING_POINT_TYPE
    // x.real != 0
    calculus::details::_dual_number abs(calculus::details::_dual_number x);

    // power group

    calculus::details::_dual_number sq(calculus::details::_dual_number x);

    calculus::details::_dual_number cb(calculus::details::_dual_number x);

    calculus::details::_dual_number sqrt(calculus::details::_dual_number x);

    calculus::details::_dual_number cbrt(calculus::details::_dual_number x);

    // x^x or f(x) ^ f(x)
    calculus::details::_dual_number pow(calculus::details::_dual_number x, math::real p);

    // exponential and logarithmic group

    calculus::details::_dual_number pow(calculus::details::_dual_number x);

    // pow f(x)^g(x)
    // calculus::details::_dual_number pow(calculus::details::_dual_number x, algebra::dual_numberg_x)

    calculus::details::_dual_number exp(calculus::details::_dual_number x);

    calculus::details::_dual_number exp_n(math::real n, calculus::details::_dual_number x);

    calculus::details::_dual_number log(calculus::details::_dual_number x);

    calculus::details::_dual_number ln(calculus::details::_dual_number x);

    // n > 1
    calculus::details::_dual_number log_n(math::real n, calculus::details::_dual_number x);

    // x.real > 1, log_x_n
    calculus::details::_dual_number log_x_n(calculus::details::_dual_number x, math::real n);

    // trigonometric group

    calculus::details::_dual_number sin(calculus::details::_dual_number x);

    calculus::details::_dual_number cos(calculus::details::_dual_number x);

    calculus::details::_dual_number tan(calculus::details::_dual_number x);

    calculus::details::_dual_number cot(calculus::details::_dual_number x);

    calculus::details::_dual_number sec(calculus::details::_dual_number x);

    calculus::details::_dual_number csc(calculus::details::_dual_number x);

    calculus::details::_dual_number asin(calculus::details::_dual_number x);

    calculus::details::_dual_number acos(calculus::details::_dual_number x);

    calculus::details::_dual_number atan(calculus::details::_dual_number x);

    calculus::details::_dual_number acot(calculus::details::_dual_number x);

    calculus::details::_dual_number asec(calculus::details::_dual_number x);

    calculus::details::_dual_number acsc(calculus::details::_dual_number x);

    // hyperbolic group

    calculus::details::_dual_number sinh(calculus::details::_dual_number x);

    calculus::details::_dual_number cosh(calculus::details::_dual_number x);

    calculus::details::_dual_number tanh(calculus::details::_dual_number x);

    calculus::details::_dual_number coth(calculus::details::_dual_number x);

    calculus::details::_dual_number sech(calculus::details::_dual_number x);

    calculus::details::_dual_number csch(calculus::details::_dual_number x);

    calculus::details::_dual_number asinh(calculus::details::_dual_number x);

    calculus::details::_dual_number acosh(calculus::details::_dual_number x);

    calculus::details::_dual_number atanh(calculus::details::_dual_number x);

    calculus::details::_dual_number acoth(calculus::details::_dual_number x);

    calculus::details::_dual_number asech(calculus::details::_dual_number x);

    calculus::details::_dual_number acsch(calculus::details::_dual_number x);

    // miscellaneous group

#endif // !defined USE_GLOBAL_FLOATING_POINT_TYPE
} // namespace math

#endif // MATH_ALGEBRA_FO_AUTO_DIFF_HPP