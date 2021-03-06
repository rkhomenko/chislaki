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
matrix<T> jacobi_matrix(const matrix<T>& x,
                        T delta,
                        const std::array<F, N>& f) {
    auto result = matrix<T>(N, x.rows());

    for (index_type i = 0; i < N; i++) {
        for (index_type j = 0; j < result.columns(); j++) {
            result(i, j) = diff(x, delta, j, f[i]);
        }
    }

    return result;
}

template <class T>
T diff1(T x0, const matrix<T>& x, const matrix<T>& y) {
    index_type i = 0;
    for (index_type j = 0; j < x.rows() - 1; j++) {
        if (x(j) <= x0 && x0 <= x(j + 1)) {
            i = j;
            break;
        }
    }

    auto y_ = [&y](auto&& i1, auto&& i2) { return y(i1) - y(i2); };
    auto x_ = [&x](auto&& i1, auto&& i2) { return x(i1) - x(i2); };

    return y_(i + 1, i) / x_(i + 1, i) +
           (y_(i + 2, i + 1) / x_(i + 2, i + 1) - y_(i + 1, i) / x_(i + 1, i)) /
               x_(i + 2, i) * (2 * x0 - x(i) - x(i + 1));
}

template <class T>
T diff2(T x0, const matrix<T>& x, const matrix<T>& y) {
    index_type i = 0;
    for (index_type j = 0; j < x.rows() - 1; j++) {
        if (x(j) <= x0 && x0 <= x(j + 1)) {
            i = j;
            break;
        }
    }

    auto y_ = [&y](auto&& i1, auto&& i2) { return y(i1) - y(i2); };
    auto x_ = [&x](auto&& i1, auto&& i2) { return x(i1) - x(i2); };

    return 2 *
           (y_(i + 2, i + 1) / x_(i + 2, i + 1) - y_(i + 1, i) / x_(i + 1, i)) /
           x_(i + 2, i);
}

}  // namespace chislaki

#endif  // CHISLAKI_MATHAN_DIFFERENTIATION_HPP_
