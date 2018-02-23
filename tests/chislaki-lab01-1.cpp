#include <chislaki/linalg/utility.hpp>

using namespace chislaki;

int main() {
    size_type n;
    std::cin >> n;

    auto A = matrix<long double>(n);
    auto b = make_column<long double>(n);

    for (size_type i = 0; i < n; i++) {
        for (size_type j = 0; j < n; j++) {
            std::cin >> A(i, j);
        }
    }

    for (size_type i = 0; i < n; i++) {
        std::cin >> b(i);
    }

    auto[L, U, P, swap_count] = LUP_decomposition(A);
    auto x = LUP_solver(A, b);

    std::cout << "************************** LU decomposition "
                 "**************************\n"
              << "L:\n"
              << L << "U:\n"
              << U << "P:\n"
              << P << "Checking: A = LUP\n"
              << L * U * P
              << "************************* Solve linear system "
                 "*************************\n"
                 "x:\n"
              << x << "Checking: A * x\n"
              << A * x
              << "***************************** Determinant "
                 "*****************************\n"
              << determinant(A) << std::endl
              << "****************************** Iversion "
                 "******************************\n"
              << "inverse(A):\n"
              << inverse(A) << "Checking: A * inverse(A)\n"
              << A * inverse(A);
}
