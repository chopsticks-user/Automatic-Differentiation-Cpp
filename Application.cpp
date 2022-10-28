#include <Math.hpp>

#include <iostream>
#include <typeinfo>

using namespace math;

auto func = [](auto x)
{
    return x * x;
};

auto sfunc = [](auto x, auto y)
{
    return x * x / math::sin(x);
};

using calculus::details::_auto_diff;

int main([[maybe_unused]] int argc,
         [[maybe_unused]] char **argv)
{
    double x = -5;
    double y = 3.0;

    std::cout << _auto_diff(func, x) << '\n';
    std::cout << _auto_diff<0>(sfunc, x, y) << '\n';

    return 0;
}