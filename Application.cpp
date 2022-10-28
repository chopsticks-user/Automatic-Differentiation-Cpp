#include <Math.hpp>

#include <iostream>
#include <typeinfo>

using namespace math;

auto func = [](auto x)
{ return math::log(x); };

auto sfunc = [](auto x)
{ return 1 - 2 * x - x * x; };

// f(x) = a - g(x)
// auto sfunc = [](auto x)
// { return math::sqrt(1.0 - math::pow(x, -3.0)); };

int main([[maybe_unused]] int argc,
         [[maybe_unused]] char **argv)
{
    auto x = 3.0;

    std::cout << calculus::auto_diff(func, x) << '\n';
    std::cout << calculus::auto_diff(sfunc, x) << '\n';

    return 0;
}