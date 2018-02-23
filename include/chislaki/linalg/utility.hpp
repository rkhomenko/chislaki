#ifndef CHISLAKI_LINALG_UTILITY_HPP_
#define CHISLAKI_LINALG_UTILITY_HPP_

#include <chislaki/linalg/decompositions.hpp>

namespace chislaki {

// **********************************************************************
// ***************************** LUP solver *****************************
// **********************************************************************

template <class T>
matrix<T> LUP_solver(const matrix<T>& A, const matrix<T>& b) {
    using index_type = typename matrix<T>::index_type;

    if (A.rows() != b.rows()) {
        throw bad_size{};
    }

    auto[L, U, P] = LUP_decomposition(A);

    // Ax = b, A = LU =>
    // LUx = b, Lz = Pb, Ux = z
    auto z = matrix<T>::column(A.columns());

    auto b1 = P * b;
    z(0) = b1(0);
    for (index_type i = 1; i < A.rows(); i++) {
        T sum = 0;
        for (index_type j = 0; j < i; j++) {
            sum += L(i, j) * z(j);
        }
        z(i) = b1(i) - sum;
    }

    // Ux = z
    auto x = matrix<T>::column(A.columns());
    const auto n = x.rows();
    auto U1 = U;
    x(n - 1) = z(n - 1) / U1(n - 1, n - 1);
    for (index_type i = n - 1; i != 0; i--) {
        T sum = 0;
        for (index_type j = (i - 1) + 1; j < n; j++) {
            sum += U1(i - 1, j) * x(j);
        }
        x(i - 1) = (z(i - 1) - sum) / U1(i - 1, i - 1);
    }

    return x;
}
}  // namespace chislaki

#endif  // CHISLAKI_LINALG_UTILITY_HPP_
