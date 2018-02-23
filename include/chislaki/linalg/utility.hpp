#ifndef CHISLAKI_LINALG_UTILITY_HPP_
#define CHISLAKI_LINALG_UTILITY_HPP_

#include <chislaki/linalg/decompositions.hpp>

namespace chislaki {

// **********************************************************************
// ***************************** LUP solver *****************************
// **********************************************************************

template <class T>
matrix<T> LUP_solver(const matrix<T>& A, const matrix<T>& b) {
    if (A.rows() != b.rows()) {
        throw bad_size{};
    }

    auto P = make_eye<T>(A.rows());
    auto L = P;
    auto U = P;
    std::tie(L, U, P, std::ignore) = LUP_decomposition(A);

    // Ax = b, A = LU =>
    // LUx = b, Lz = Pb, Ux = z
    auto z = make_column<T>(A.columns());

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
    auto x = make_column<T>(A.columns());
    const auto n = x.rows();
    x(n - 1) = z(n - 1) / U(n - 1, n - 1);
    for (index_type i = n - 1; i != 0; i--) {
        T sum = 0;
        for (index_type j = (i - 1) + 1; j < n; j++) {
            sum += U(i - 1, j) * x(j);
        }
        x(i - 1) = (z(i - 1) - sum) / U(i - 1, i - 1);
    }

    return x;
}

// ***********************************************************************
// ***************************** Determinant *****************************
// ***********************************************************************

template <class T>
T determinant(const matrix<T>& matr) {
    if (matr.rows() != matr.columns()) {
        throw bad_size{};
    }

    auto U = make_eye<T>(matr.columns());
    auto P = U;
    size_type swap_count = 0;
    std::tie(std::ignore, U, std::ignore, swap_count) = LUP_decomposition(matr);

    // PA = LUP, det(A) = det(L)det(U)det(P)
    // det(L) = 1
    // det(P) = (-1)^swap_count

    auto det = [](auto&& m) {
        T result = 1;
        for (index_type i = 0; i < m.rows(); i++) {
            result *= m(i, i);
        }
        return result;
    };

    return det(U) * (swap_count % 2 == 0 ? 1 : -1);
}

// ***********************************************************************
// ****************************** Inversion ******************************
// ***********************************************************************

template <class T>
matrix<T> inverse(const matrix<T>& matr) {
    if (determinant(matr) == 0) {
        throw degenerate{};
    }

    auto result = make_eye<T>(matr.rows());

    for (index_type j = 0; j < result.columns(); j++) {
        auto b = make_column(j, result);
        result.set_column(j, LUP_solver(matr, b));
    }

    return result;
}

}  // namespace chislaki

#endif  // CHISLAKI_LINALG_UTILITY_HPP_
