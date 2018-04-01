#include <chislaki/linalg/utility.hpp>

#include <cstdlib>
#include <string>

using namespace chislaki;

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Bad arguments count! Needed method type and epsilon!"
                  << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (std::string(argv[1]) == "i") {
        long double epsilon = std::stold(std::string(argv[2]));
        size_type n;
        std::cin >> n;

        auto matr = matrix<long double>(n);
        auto b = make_column<long double>(n);
        std::cin >> matr >> b;

        auto x = fixed_point_iteration(matr, b, epsilon);
        std::cout << "Checking Ax = b:\n" << matr * x << "b:\n" << b;
    } else if (std::string(argv[1]) == "z") {
        long double epsilon = std::stold(std::string(argv[2]));
        size_type n;
        std::cin >> n;

        auto matr = matrix<long double>(n);
        auto b = make_column<long double>(n);
        std::cin >> matr >> b;

        auto x = gauss_seidel(matr, b, epsilon);
        std::cout << "Checking Ax = b:\n" << matr * x << "b:\n" << b;
    } else {
        std::cout << "Bad method type!" << std::endl;
    }
}
