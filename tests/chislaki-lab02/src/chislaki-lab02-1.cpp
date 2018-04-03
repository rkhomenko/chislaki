#include <chislaki/transcendental/utility.hpp>

#include <cmath>
#include <iostream>

using namespace chislaki;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Bad arguments count! Needed method type and epsilon!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    double epsilon = std::stod(std::string(argv[2]));
    auto f = [](auto&& x) { return std::log(x + 2) - x * x; };

    if (std::string(argv[1]) == "i") {
        auto x1 = fixed_point_iteration_solver(-1.0, -1.0, 0.0, epsilon, f);
        auto x2 = fixed_point_iteration_solver(0.5, 0.5, 1.5, epsilon, f);

        std::cout << "************************ Fixed point iteration "
                     "***********************"
                  << std::endl;
        std::cout << "x1 = " << x1 << ", f(x1) = " << f(x1) << std::endl;
        std::cout << "x2 = " << x2 << ", f(x2) = " << f(x2) << std::endl;
    } else if (std::string(argv[1]) == "n") {
        auto x1 = newton_solver(-1.0, epsilon, f);
        auto x2 = newton_solver(0.5, epsilon, f);

        std::cout << "**************************** Newton method "
                     "***************************"
                  << std::endl;
        std::cout << "x1 = " << x1 << ", f(x1) = " << f(x1) << std::endl;
        std::cout << "x2 = " << x2 << ", f(x2) = " << f(x2) << std::endl;
    } else {
        std::cout << "Bad method name!" << std::endl;
    }
}
