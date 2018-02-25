#ifndef CHISLAKI_LINALG_UTILITY_HPP_
#define CHISLAKI_LINALG_UTILITY_HPP_

#include <chislaki/linalg/matrix.hpp>

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

}  // namespace chislaki

#endif  // CHISLAKI_LINALG_UTILITY_HPP_
