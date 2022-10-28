#include "DualNumber.hpp"

#ifdef USE_GLOBAL_FLOATING_POINT_TYPE
namespace math
{
    // x.real != 0
    algebra::dual_number abs(algebra::dual_number x)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::abs<dual_number>");
        return algebra::dual_number{
            std::abs(x.real),
            x.dual * x.real / std::abs(x.real)};
    }

    // power group

    algebra::dual_number sq(algebra::dual_number x)
    {
        return algebra::dual_number{
            x.real * x.real,
            x.dual * 2.0 * x.real};
    }

    algebra::dual_number cb(algebra::dual_number x)
    {
        return algebra::dual_number{
            x.real * x.real * x.real,
            x.dual * 3.0 * x.real * x.real};
    }

    algebra::dual_number sqrt(algebra::dual_number x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::sqrt<dual_number>");
        auto sqrt_xr = std::sqrt(x.real);
        return algebra::dual_number{
            sqrt_xr,
            x.dual * (0.5 / sqrt_xr)};
    }

    algebra::dual_number cbrt(algebra::dual_number x)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::cbrt<dual_number>");
        auto cbrt_xr = std::cbrt(x.real);
        return algebra::dual_number{
            cbrt_xr,
            x.dual / (3.0 * cbrt_xr * cbrt_xr)};
    }

    // x^x or f(x) ^ f(x)
    algebra::dual_number pow(algebra::dual_number x, math::real p)
    {
        if (x.real == 0.0)
            throw std::runtime_error("x.real = 0 at math::pow_x_n<dual_number>");
        auto pow_xr = std::pow(x.real, p);
        return algebra::dual_number{
            pow_xr,
            p * x.dual * pow_xr / x.real};
    }

    // exponential and logarithmic group

    algebra::dual_number pow(algebra::dual_number x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::pow_x_x<dual_number>");
        auto xr_pow_xr = std::pow(x.real, x.real);
        return algebra::dual_number{
            xr_pow_xr,
            x.dual * xr_pow_xr * (1 + std::log(x.real))};
    }

    // pow f(x)^g(x)
    // algebra::dual_number pow(algebra::dual_number x, algebra::dual_numberg_x)
    // {
    //     auto xr_pow_xr = std::pow(x.real, x.real);
    //     return algebra::dual_number{
    //         xr_pow_xr,
    //         x.dual * xr_pow_xr * (1 + std::log(x.real))};
    // }

    algebra::dual_number exp(algebra::dual_number x)
    {
        auto exp_xr = std::exp(x.real);
        return algebra::dual_number{
            exp_xr,
            x.dual * exp_xr};
    }

    algebra::dual_number exp_n(math::real n, algebra::dual_number x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::exp_n_x<dual_number>");
        auto exp_n_xr = std::pow(n, x.real);
        return algebra::dual_number{
            exp_n_xr,
            x.dual * std::log(n) * exp_n_xr};
    }

    algebra::dual_number log(algebra::dual_number x)
    {
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::ln<dual_number>");
        return algebra::dual_number{
            std::log(x.real),
            x.dual / x.real};
    }

    algebra::dual_number ln(algebra::dual_number x)
    {
        return math::log(x);
    }

    // n > 1
    algebra::dual_number log_n(math::real n, algebra::dual_number x)
    {
        if (n <= 0.0 || n == 1.0)
            throw std::runtime_error("n <= 0 || n = 1 at math::log_n_x<dual_number>");
        if (x.real <= 0.0)
            throw std::runtime_error("x.real <= 0 at math::log_n_x<dual_number>");
        auto ln_n = std::log(n);
        return algebra::dual_number{
            std::log(x.real) / ln_n,
            x.dual / (x.real * ln_n)};
    }

    // x.real > 1, log_x_n
    algebra::dual_number log_x_n(algebra::dual_number x, math::real n)
    {
        if (n <= 0.0)
            throw std::runtime_error("n <= 0 at math::log_x_n<dual_number>");
        if (x.real <= 0.0 || x.real == 1.0)
            throw std::runtime_error("x.real <= 0 || x.real = 1 at math::log_x_n<dual_number>");
        auto ln_n = std::log(n);
        auto ln_x = std::log(x.real);
        return algebra::dual_number{
            ln_n / ln_x,
            x.dual * (-ln_n / (x.real * ln_x * ln_x))};
    }

    // trigonometric group

    algebra::dual_number sin(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::sin(x.real),
            x.dual * (std::cos(x.real))};
    }

    algebra::dual_number cos(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::cos(x.real),
            x.dual * (-std::sin(x.real))};
    }

    algebra::dual_number tan(algebra::dual_number x)
    {
        // if (std::cos(x.real) == 0.0)
        //     throw std::runtime_error("x.real == k(pi / 2) at math::tan<dual_number>");
        auto tan_xr = std::tan(x.real);
        return algebra::dual_number{
            tan_xr,
            x.dual * (1.0 + tan_xr * tan_xr)};
    }

    algebra::dual_number cot(algebra::dual_number x)
    {
        auto cot_xr = 1.0 / std::tan(x.real);
        return algebra::dual_number{
            cot_xr,
            x.dual * (-1.0 - cot_xr * cot_xr)};
    }

    algebra::dual_number sec(algebra::dual_number x)
    {
        auto cos_xr = std::cos(x.real);
        return algebra::dual_number{
            1.0 / cos_xr,
            x.dual * (std::tan(x.real) / cos_xr)};
    }

    algebra::dual_number csc(algebra::dual_number x)
    {
        auto sin_xr = std::sin(x.real);
        return algebra::dual_number{
            1.0 / sin_xr,
            x.dual * (-1.0 / (sin_xr * std::tan(x.real)))};
    }

    algebra::dual_number asin(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::asin(x.real),
            x.dual / (std::sqrt(1 - x.real * x.real))};
    }

    algebra::dual_number acos(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::acos(x.real),
            -x.dual / (std::sqrt(1.0 - x.real * x.real))};
    }

    algebra::dual_number atan(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::atan(x.real),
            x.dual / (1.0 + x.real * x.real)};
    }

    algebra::dual_number acot(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::atan(1.0 / x.real),
            x.dual / (-1.0 - x.real * x.real)};
    }

    algebra::dual_number asec(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::acos(1.0 / x.real),
            x.dual / (std::abs(x.real) * std::sqrt(x.real * x.real - 1))};
    }

    algebra::dual_number acsc(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::acos(1.0 / x.real),
            -x.dual / (std::abs(x.real) * std::sqrt(x.real * x.real - 1))};
    }

    // hyperbolic group

    algebra::dual_number sinh(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::sinh(x.real),
            x.dual * std::cosh(x.real)};
    }

    algebra::dual_number cosh(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::cosh(x.real),
            x.dual * std::sinh(x.real)};
    }

    algebra::dual_number tanh(algebra::dual_number x)
    {
        auto tanh_xr = std::tanh(x.real);
        return algebra::dual_number{
            tanh_xr,
            x.dual * (1.0 - tanh_xr * tanh_xr)};
    }

    algebra::dual_number coth(algebra::dual_number x)
    {
        auto coth_xr = 1.0 / std::tanh(x.real);
        return algebra::dual_number{
            coth_xr,
            x.dual * (1.0 - coth_xr * coth_xr)};
    }

    algebra::dual_number sech(algebra::dual_number x)
    {
        auto sech_xr = 1.0 / std::cosh(x.real);
        return algebra::dual_number{
            sech_xr,
            x.dual * (-sech_xr * std::tanh(x.real))};
    }

    algebra::dual_number csch(algebra::dual_number x)
    {
        auto csch_xr = 1.0 / std::sinh(x.real);
        return algebra::dual_number{
            csch_xr,
            x.dual * (-csch_xr / std::tanh(x.real))};
    }

    algebra::dual_number asinh(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::asinh(x.real),
            x.dual * std::sqrt(1.0 + x.real * x.real)};
    }

    algebra::dual_number acosh(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::acosh(x.real),
            x.dual * std::sqrt(x.real * x.real - 1.0)};
    }

    algebra::dual_number atanh(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::atanh(x.real),
            x.dual * (1.0 - x.real * x.real)};
    }

    algebra::dual_number acoth(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::atanh(1.0 / x.real),
            x.dual * (1.0 - x.real * x.real)};
    }

    algebra::dual_number asech(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::acosh(1.0 / x.real),
            -x.dual * x.real * std::sqrt(1.0 - x.real * x.real)};
    }

    algebra::dual_number acsch(algebra::dual_number x)
    {
        return algebra::dual_number{
            std::acosh(1.0 / x.real),
            -x.dual * std::abs(x.real) * std::sqrt(1.0 + x.real * x.real)};
    }

    // miscellaneous group

}; // namespace math

#endif