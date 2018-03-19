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

    qr_eigenvalues qre(matr, epsilon);

    for (index_type i = 0; i < matr.rows(); i++) {
        std::cout << "l_" << i << " = " << qre.get_value(i) << std::endl;
    }
}
