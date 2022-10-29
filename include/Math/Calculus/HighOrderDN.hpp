#ifndef MATH_CALCULUS_HIGH_ORDER_DN_HPP
#define MATH_CALCULUS_HIGH_ORDER_DN_HPP

#include "Config.hpp"

#include <vector>
#include <array>
#include <algorithm>

namespace math::calculus
{
    namespace details
    {
#ifdef USE_GLOBAL_FLOATING_POINT_TYPE

#else  // !defined USE_GLOBAL_FLOATING_POINT_TYPE
        template <typename value_type = math::real, size_type highest_order = 1>
        class _high_order_dual_number
        {
            typedef std::array<value_type, highest_order> vl_type;
            typedef _high_order_dual_number<value_type, highest_order> same_type;

        public:
            _high_order_dual_number() = default;
            _high_order_dual_number(value_type value) : _value_list{}
            {
                _value_list[0] = value;
                _value_list[1] = 1.0;
            }

            value_type derivative(size_type order) const noexcept
            {
                value_type fact = 1.0;
                for (size_type i = 2; i <= order; ++i)
                    fact *= i;
                return _value_list[order] * fact;
            }

            same_type operator+(const same_type &rhs) const noexcept
            {
                same_type result{};

                // size_type index = 0;
                // std::for_each(result._value_list.begin(), result._value_list.end(),
                //               [&](value_type &n)
                //               {n=this->_value_list[index]+rhs._value_list[index];++index; });

                for (size_type i = 0; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] + rhs._value_list[i];
                return result;
            }

            same_type operator-(const same_type &rhs) const noexcept
            {
                same_type result{};
                for (size_type i = 0; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] - rhs._value_list[i];
                return result;
            }

            same_type operator*(const same_type &rhs) const noexcept
            {
                same_type result{};
                for (size_type i = 0; i < highest_order; ++i)
                    for (size_type j = 0; j <= i; ++j)
                        result._value_list[i] += _value_list[j] * rhs._value_list[i - j];
                return result;
            }

            same_type operator*(value_type scalar) const noexcept
            {
                same_type result{};
                for (size_type i = 0; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] * scalar;
                return result;
            }

            friend same_type operator*(value_type scalar, const same_type &rhs) noexcept
            {
                same_type result{};
                for (size_type i = 0; i < highest_order; ++i)
                    result._value_list[i] = rhs._value_list[i] * scalar;
                return result;
            }

            // zi must be calculated sequetially
            same_type operator/(const same_type &rhs) const noexcept
            {
                same_type result{};
                for (size_type i = 0; i < highest_order; ++i)
                {
                    for (size_type j = 0; j <= i; ++j)
                        result._value_list[i] += rhs._value_list[j] * result._value_list[i - j];
                    result._value_list[i] = _value_list[i] - result._value_list[i];
                    result._value_list[i] /= rhs._value_list[0];
                    if (result._value_list[i] == 0)
                        return result;
                }
                return result;
            }

            same_type operator/(value_type scalar) const noexcept
            {
                same_type result{};
                for (size_type i = 0; i < highest_order; ++i)
                    result._value_list[i] = _value_list[i] / scalar;
                return result;
            }

            friend same_type operator/(value_type scalar, const same_type &rhs) noexcept
            {
                same_type result{scalar / rhs._value_list[0]};
                value_type temp = result._value_list[0];
                for (size_type i = 1; i < highest_order; ++i)
                {
                    temp *= i * -1.0 / rhs._value_list[0];
                    result._value_list[i] = temp;
                }
                return result;
            }

        private:
            std::array<value_type, highest_order> _value_list;
        };
#endif // USE_GLOBAL_FLOATING_POINT_TYPE
    }  // namespace math::calculus::details
} // namespace math::calculus

#endif // MATH_CALCULUS_HIGH_ORDER_DN_HPP