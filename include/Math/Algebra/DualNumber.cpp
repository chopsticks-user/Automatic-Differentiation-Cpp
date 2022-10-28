#ifdef AAA

#include "DualNumber.hpp"

namespace math
{
    template struct algebra::dual_number<float>;
    template struct algebra::dual_number<double>;
    template struct algebra::dual_number<long double>;

    template algebra::dual_number<float> sin(algebra::dual_number<float> x);
    template algebra::dual_number<double> sin(algebra::dual_number<double> x);
    template algebra::dual_number<long double> sin(algebra::dual_number<long double> x);

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

}; // namespace math

#endif