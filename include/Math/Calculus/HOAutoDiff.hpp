#ifndef MATH_CALCULUS_HO_AUTO_DIFF_HPP
#define MATH_CALCULUS_HO_AUTO_DIFF_HPP

#include "Config.hpp"

#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <type_traits>

namespace math
{
    namespace calculus::details
    {
#ifdef USE_GLOBAL_FLOATING_POINT_TYPE

#else  // !defined USE_GLOBAL_FLOATING_POINT_TYPE
        template <typename value_type = math::real, size_type highest_order = 1>
        class _high_order_dual_number
        {
            static_assert(std::is_floating_point<value_type>::value);
            typedef std::array<value_type, highest_order> vl_type;
            typedef _high_order_dual_number<value_type, highest_order> same_type;

        public:
            _high_order_dual_number() = default;
            explicit _high_order_dual_number(value_type value) : _value_list{}
            {
                _value_list[0] = value;
                _value_list[1] = 1.0;
            }

            value_type derivative(size_type order) const noexcept
            {
                size_type fact = 1.0;
                for (size_type i = 2.0; i <= order; ++i)
                    fact *= i;
                return fact * _value_list[order];
            }

            same_type operator+(const same_type &rhs) const noexcept
            {
                same_type result{};

                // size_type index = 0;
                // std::for_each(result._value_list.begin(), result._value_list.end(),
                //               [&](value_type &n)
                //               {n=this->_value_list[index]+rhs._value_list[index];++index; });

                result._value_list[0] = _value_list[0] + rhs._value_list[0];
                for (size_type i = 1; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] + rhs._value_list[i];
                return result;
            }

            same_type operator-(const same_type &rhs) const noexcept
            {
                same_type result{};
                result._value_list[0] = _value_list[0] - rhs._value_list[0];
                for (size_type i = 1; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] - rhs._value_list[i];
                return result;
            }

            same_type operator*(const same_type &rhs) const noexcept
            {
                same_type result{};
                result._value_list[0] = _value_list[0] * rhs._value_list[0];
                for (size_type i = 1; i < highest_order; ++i)
                    for (size_type j = 0; j <= i; ++j)
                        result._value_list[i] += _value_list[j] * rhs._value_list[i - j];
                return result;
            }

            same_type operator*(value_type scalar) const noexcept
            {
                same_type result{};
                result._value_list[0] = scalar * _value_list[0];
                for (size_type i = 1; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] * scalar;
                return result;
            }

            friend same_type operator*(value_type scalar, const same_type &rhs) noexcept
            {
                same_type result{};
                result._value_list[0] = scalar * rhs._value_list[0];
                for (size_type i = 1; i < highest_order; ++i)
                    result._value_list[i] = rhs._value_list[i] * scalar;
                return result;
            }

            // zi must be calculated sequetially
            same_type operator/(const same_type &rhs) const noexcept
            {
                same_type result{_value_list[0] / rhs._value_list[0]};
                for (size_type i = 1; i < highest_order; ++i)
                {
                    value_type temp = 0.0;
                    for (size_type j = 1; j <= i; ++j)
                        temp += rhs._value_list[j] * result._value_list[i - j];
                    result._value_list[i] = (_value_list[i] - temp) / rhs._value_list[0];
                }
                return result;
            }

            same_type operator/(value_type scalar) const noexcept
            {
                same_type result{_value_list[0] / scalar};
                for (size_type i = 1; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] / scalar;
                return result;
            }

            friend same_type operator/(value_type scalar, const same_type &rhs) noexcept
            {
                same_type result{scalar};
                result._value_list[1] = 0.0;
                return result / rhs;
            }

            // power group

            // exponential and logarithmic group

            // trigonometric group

            friend same_type sin(const same_type &x)
            {
                value_type sin_result = std::sin(x._value_list[0]);
                same_type result{sin_result};
                if (highest_order == 1)
                    return result;
                value_type cos_result = std::cos(x._value_list[0]);
                value_type fact = 1.0;
                for (size_type i = 1; i < highest_order; ++i, fact *= i)
                {
                    switch (i % 4)
                    {
                    case 0:
                        result._value_list[i] = sin_result;
                        break;
                    case 1:
                        result._value_list[i] = cos_result;
                        break;
                    case 2:
                        result._value_list[i] = -sin_result;
                        break;
                    case 3:
                        result._value_list[i] = -cos_result;
                        break;
                    }
                    result._value_list[i] /= fact;
                }
                return result;
            }

            friend same_type cos(const same_type &x)
            {
                value_type cos_result = std::cos(x._value_list[0]);
                same_type result{cos_result};
                if (highest_order == 1)
                    return result;
                value_type sin_result = std::sin(x._value_list[0]);
                value_type fact = 1.0;
                for (size_type i = 1; i < highest_order; ++i, fact *= i)
                {
                    switch (i % 4)
                    {
                    case 0:
                        result._value_list[i] = cos_result;
                        break;
                    case 1:
                        result._value_list[i] = -sin_result;
                        break;
                    case 2:
                        result._value_list[i] = -cos_result;
                        break;
                    case 3:
                        result._value_list[i] = sin_result;
                        break;
                    }
                    result._value_list[i] /= fact;
                }
                return result;
            }

            friend same_type tan(const same_type &x)
            {
                return sin(x) / cos(x);
            }

            friend same_type cot(const same_type &x)
            {
                return cos(x) / sin(x);
            }

            friend same_type sec(const same_type &x)
            {
                return 1.0 / cos(x);
            }

            friend same_type csc(const same_type &x)
            {
                return 1.0 / sin(x);
            }

        private:
            std::array<value_type, highest_order> _value_list;
        };
#endif // USE_GLOBAL_FLOATING_POINT_TYPE
    }  // namespace math::calculus::details

    // friend auto math::sin(const same_type &x)
    // {

    // }

} // namespace math

#endif // MATH_CALCULUS_HO_AUTO_DIFF_HPP