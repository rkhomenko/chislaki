#include <chislaki/transcendental/utility.hpp>

#include <array>
#include <cmath>
#include <functional>
#include <iostream>

using namespace chislaki;

template <class T>
matrix<T> fixed_point_iteration(matrix<T> x0, T epsilon) {
    T a = 3;

    auto phi1 = [&a](auto&& x) {
        // return std::sqrt(a * a - std::pow(x(1) - a / 2, 2)) + a / 2;
        return (x(0) * x(0) + x(1) * x(1)) / a - x(1) - a + 0.5;
    };

    auto phi2 = [&a](auto&& x) { return a * a * a / (x(0) * x(0) + a * a); };

    auto x_k = x0;
    auto x_k_1 = x0;

    while (true) {
        x_k_1(0) = phi1(x_k);
        x_k_1(1) = phi2(x_k);
        if (norm_inf(x_k_1 - x_k) < epsilon) {
            break;
        }
        x_k = x_k_1;
    }

    return x_k_1;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Bad arguments count! Needed method type and epsilon!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    double epsilon = std::stod(std::string(argv[2]));
    double a = 3;

    using functor = std::function<double(const matrix<double>&)>;

    functor f1 = [&a](const matrix<double>& x) {
        return (x(0) * x(0) + a * a) * x(1) - a * a * a;
    };

    functor f2 = [&a](const matrix<double>& x) {
        return std::pow(x(0) - a / 2, 2) + std::pow(x(1) - a / 2, 2) - a * a;
    };

    auto func_array = std::array<functor, 2>{f1, f2};

    auto x01 = make_column<double>(2, std::vector<double>{-1, 2});
    auto x02 = make_column<double>(2, std::vector<double>{4.5, 1});

    if (std::string(argv[1]) == "i") {
        auto x1 = fixed_point_iteration_solver(x01, epsilon, func_array);
        auto x2 = fixed_point_iteration_solver(x02, epsilon, func_array);

        std::cout << "x1 = \n"
                  << x1 << "f1(x1) = " << f1(x1) << std::endl
                  << "f2(x1) = " << f2(x1) << std::endl;
        std::cout << "x2 = \n"
                  << x2 << "f1(x2) = " << f1(x2) << std::endl
                  << "f2(x2) = " << f2(x2) << std::endl;
    } else if (std::string(argv[1]) == "n") {
    } else {
        std::cout << "Bad method name!" << std::endl;
    }
}
