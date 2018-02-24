#include <chislaki/linalg/decompositions.hpp>

#include <eigen3/Eigen/LU>

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

    lup_decomposition lup(A);

    const auto& L = lup.matrix_l();
    const auto& U = lup.matrix_u();
    const auto& P = lup.matrix_p();
    auto x = lup.solve(b);

    std::cout << "******************************** A, b "
                 "********************************\n"
              << "A:\n"
              << A << "b:\n"
              << b
              << "************************** LU decomposition "
                 "**************************\n"
              << "L:\n"
              << L << "U:\n"
              << U << "P:\n"
              << P << "Checking: A = P^(T)LU\n"
              << transpose(P) * L * U
              << "************************* Solve linear system "
                 "*************************\n"
                 "x:\n"
              << x << "Checking: A * x\n"
              << A * x
              << "***************************** Determinant "
                 "*****************************\n"
              << lup.determinant() << std::endl
              << "****************************** Iversion "
                 "******************************\n"
              << "inverse(A):\n"
              << lup.inverse() << "Checking: A * inverse(A)\n"
              << A * lup.inverse();
}
