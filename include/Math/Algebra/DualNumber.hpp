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
                return dual_number{scalar - d_num.real, d_num.dual};
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

    // power group

    template <typename var_type>
    algebra::dual_number<var_type> sqrt(algebra::dual_number<var_type> x)
    {
        auto sqrt_xr = std::sqrt(x.real);
        return algebra::dual_number<var_type>{
            sqrt_xr,
            x.dual * (0.5 / sqrt_xr)};
    }

    // x only
    template <typename var_type>
    algebra::dual_number<var_type> cbrt(algebra::dual_number<var_type> x)
    {
        auto cbrt_xr = std::cbrt(x.real);
        return algebra::dual_number<var_type>{
            cbrt_xr,
            x.dual * (1.0 / (3 * cbrt_xr * cbrt_xr))};
    }

    template <typename var_type>
    algebra::dual_number<var_type> pow(algebra::dual_number<var_type> x, var_type p)
    {
        auto pow_xr = std::pow(x.real, p);
        return algebra::dual_number<var_type>{
            pow_xr,
            p * x.dual * pow_xr / x.real};
    }

    // exponential and logarithmic group

    // x only
    template <typename var_type>
    algebra::dual_number<var_type> pow(algebra::dual_number<var_type> x)
    {
        auto xr_pow_xr = std::pow(x.real, x.real);
        return algebra::dual_number<var_type>{
            xr_pow_xr,
            x.dual * xr_pow_xr * (1 + std::log(x.real))};
    }

    // template <typename var_type>
    // algebra::dual_number<var_type> exp(algebra::dual_number<var_type> x, var_type p)
    // {
    //     auto pow_xr = std::pow(x.real, p);
    //     return algebra::dual_number<var_type>{pow_xr, p * x.dual * pow_xr / x.real};
    // }

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

} // namespace math

#endif // MATH_ALGEBRA_DUAL_NUMBER_HPP