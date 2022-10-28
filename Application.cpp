#include <Math.hpp>

#include <iostream>
#include <typeinfo>

using namespace math;

auto func = [](auto x, auto y, auto z)
{ return x * 2.1f * z * z + 0.5 * y; };

auto sfunc = [](auto x)
{ return math::sin(math::sqrt(math::pow(x, 3.0))); };

int main([[maybe_unused]] int argc,
         [[maybe_unused]] char **argv)
{
    auto x = 3.0;
    // auto y = 3.0;
    // auto z = 1.0;

    // std::cout << calculus::auto_diff<0>(func, x, y, z) << '\n';
    // std::cout << calculus::auto_diff<1>(func, x, y, z) << '\n';
    // std::cout << calculus::auto_diff<2>(func, x, y, z) << '\n';

    std::cout << calculus::auto_diff(sfunc, x) << '\n';

    return 0;
}