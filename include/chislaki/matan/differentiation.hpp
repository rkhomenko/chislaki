#ifndef CHISLAKI_MATAN_DIFFERENTIATION_HPP_
#define CHISLAKI_MATAN_DIFFERENTIATION_HPP_

#include <chislaki/linalg/matrix.hpp>

#include <array>

namespace chislaki {

template <class T, class F>
T diff(T x, T delta, F&& f) {
    return (f(x + delta) - f(x)) / delta;
}

template <class T, class F>
T diff(const matrix<T>& x, T delta, index_type i, F&& f) {
    auto delta_x = make_column<T>(x.rows());
    delta_x(i) = delta;

    return (f(x + delta_x) - f(x)) / delta;
}

template <class T, class F, size_type N = 2>
matrix<T> jacobi_matrix(const matrix<T>& x, T delta, const std::array<F, N>& f) {
    auto result = matrix<T>(N, x.rows());

    for (index_type i = 0; i < N; i++) {
        for (index_type j = 0; j < result.columns(); j++) {
            result(i, j) = diff(x, delta, j, f[i]);
        }
    }

    return result;
}

}  // namespace chislaki

#endif  // CHISLAKI_MATHAN_DIFFERENTIATION_HPP_
