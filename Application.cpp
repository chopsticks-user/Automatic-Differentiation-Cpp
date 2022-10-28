#include <Math.hpp>

#include <iostream>
#include <typeinfo>

using namespace math;

auto func = [](auto x)
{
    return math::tan(x);
};

auto sfunc = [](auto x)
{
    return math::sqrt(x);
};

using calculus::details::_auto_diff;

int main([[maybe_unused]] int argc,
         [[maybe_unused]] char **argv)
{
    double x = 1.1;

    std::cout << _auto_diff(func, x) << '\n';
    std::cout << _auto_diff(sfunc, x) << '\n';

    return 0;
}