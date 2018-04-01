#include <chislaki/transcendental/utility.hpp>

#include <cmath>
#include <iostream>

using namespace chislaki;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Bad arguments count! Needed epsilon!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    double epsilon = std::stod(std::string(argv[1]));
    auto f = [](auto&& x) { return std::log(x + 2) - x * x; };

    auto x1 = fixed_point_iteration(-1.0, -1.0, 0.0, epsilon, f);
    auto x2 = fixed_point_iteration(0.5, 0.5, 1.5, epsilon, f);

    std::cout << "x1 = " << x1 << ", f(x1) = " << f(x1) << std::endl;
    std::cout << "x2 = " << x2 << ", f(x2) = " << f(x2) << std::endl;
}
