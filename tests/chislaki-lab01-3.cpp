#include <chislaki/linalg/utility.hpp>

#include <cstdlib>
#include <string>

using namespace chislaki;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Bad arguments count! Epsilon value needed!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    long double epsilon = std::stold(std::string(argv[1]));
    size_type n;
    std::cin >> n;

    auto matr = matrix<long double>(n);
    auto b = make_column<long double>(n);
    std::cin >> matr >> b;

    auto x = fixed_point_iteration(matr, b, epsilon);
    std::cout << "Checking Ax = b:\n" << matr * x << b;
}
