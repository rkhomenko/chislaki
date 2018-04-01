#include <chislaki/linalg/utility.hpp>

#include <cstdlib>
#include <string>

using namespace chislaki;

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Bad arguments count! Needed epsilon!" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    long double epsilon = std::stold(std::string(argv[1]));
    size_type n;
    std::cin >> n;

    auto matr = matrix<long double>(n);
    std::cin >> matr;

    jacobi_eigenvalues<long double> je(matr, epsilon);

    for (index_type i = 0; i < je.count(); i++) {
        auto lambda = je.get_value(i);
        auto x_i = je.get_vector(i);
        std::cout << "l_" << i + 1 << " = " << lambda << std::endl
                  << "x_" << i + 1 << ":\n"
                  << x_i << "Ax\n"
                  << matr * x_i << "lx\n"
                  << lambda * x_i << "***********************************\n";
    }
}
