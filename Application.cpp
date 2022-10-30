#include <Math.hpp>

#include <iostream>
#include <typeinfo>

using namespace math;
using calculus::details::_high_order_dual_number;

int main([[maybe_unused]] int argc,
         [[maybe_unused]] char **argv)
{
    _high_order_dual_number<double, 6> x{0.2758};
    auto f = sec(x);
    std::cout << f.derivative(0) << '\n';
    std::cout << f.derivative(1) << '\n';
    std::cout << f.derivative(2) << '\n';
    std::cout << f.derivative(3) << '\n';

    return 0;
}