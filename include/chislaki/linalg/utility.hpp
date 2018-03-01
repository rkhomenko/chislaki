#ifndef CHISLAKI_LINALG_UTILITY_HPP_
#define CHISLAKI_LINALG_UTILITY_HPP_

#include <chislaki/linalg/matrix.hpp>

#include <cmath>

namespace chislaki {

// **********************************************************************
// ******************** Tridiagonal matrix algorithm ********************
// ************************* (Thomas algorithm) *************************
// **********************************************************************

template <class T>
matrix<T> thomas_algorithm(const matrix<T>& matr, const matrix<T>& d) {
    if (matr.rows() != d.rows()) {
        throw bad_size{};
    }

    auto a = make_row<T>(d.rows());
    auto b = make_row<T>(d.rows());

    // a(0) = -c_1 / b_1, b(0) = d_1 / b_1
    a(0) = -matr(0, 1) / matr(0, 0);
    b(0) = d(0) / matr(0, 0);

    for (index_type i = 1; i < matr.rows() - 1; i++) {
        const auto& a_i = matr(i, i - 1);
        const auto& b_i = matr(i, i);
        const auto& c_i = matr(i, i + 1);

        a(i) = -c_i / (b_i + a_i * a(i - 1));
        b(i) = (d(i) - a_i * b(i - 1)) / (b_i + a_i * a(i - 1));
    }

    auto n = a.columns();
    a(n - 1) = 0;
    b(n - 1) = (d(n - 1) - matr(n - 1, n - 2) * b(n - 2)) /
               (matr(n - 1, n - 1) + matr(n - 1, n - 2) * a(n - 2));

    auto x = make_column<T>(d.rows());
    n = x.rows();

    x(n - 1) = b(n - 1);
    for (index_type i = n - 1; i != 0; i--) {
        // x_n = a_n * x_n+1 + b_n
        x(i - 1) = a(i - 1) * x(i) + b(i - 1);
    }

    return x;
}

// **********************************************************************
// ************************ Fixed point iteration ************************
// **********************************************************************

template <class T>
matrix<T> fixed_point_iteration(const matrix<T>& matr,
                                const matrix<T>& b,
                                T epsilon) {
    if (matr.rows() != b.rows()) {
        throw bad_size{};
    }

    auto betta = make_column<T>(b.rows());
    auto alpha = matrix<T>(matr.rows(), matr.columns());

    for (index_type i = 0; i < b.rows(); i++) {
        betta(i) = b(i) / matr(i, i);
    }

    for (index_type i = 0; i < matr.rows(); i++) {
        for (index_type j = 0; j < matr.columns(); j++) {
            if (i != j) {
                alpha(i, j) = -matr(i, j) / matr(i, i);
            } else {
                alpha(i, j) = 0;
            }
        }
    }

    if (norm_inf(alpha) < 1) {
        auto iter_count = (std::log10(epsilon) - std::log10(norm_inf(betta)) +
                           std::log10(1 - norm_inf(alpha))) /
                              std::log10(norm_inf(alpha)) -
                          1;

        std::cout << "Expected iteration count: " << std::ceil(iter_count)
                  << std::endl;
    } else {
        std::cout << "Cannot expect iteration count! ||alpha|| = 1!"
                  << std::endl;
    }

    auto x_k_1 = betta;
    auto x_k = make_column<T>(betta.rows());
    T norm = 0;
    size_type count = 1;
    while ((norm = norm_inf(x_k - x_k_1)) && norm > epsilon) {
        auto tmp = x_k;
        x_k = betta + alpha * x_k_1;
        x_k_1 = tmp;
        count++;
        std::cout << "||x_k - x_k-1|| = " << norm << std::endl;
        std::cout << "x_k:\n" << x_k;
    }

    std::cout << "Iteration count: " << count << std::endl;

    return x_k;
}

}  // namespace chislaki

#endif  // CHISLAKI_LINALG_UTILITY_HPP_
