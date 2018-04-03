#ifndef CHISLAKI_TRANSCENDENTAL_UTILITY_HPP_
#define CHISLAKI_TRANSCENDENTAL_UTILITY_HPP_

#include <chislaki/linalg/decompositions.hpp>
#include <chislaki/matan/differentiation.hpp>
#include <chislaki/optimize/minimize.hpp>

namespace chislaki {

// **********************************************************************
// ************************ Fixed point iteration ***********************
// **********************************************************************

template <class T, class F>
T fixed_point_iteration_solver(T x0, T a, T b, T epsilon, F&& f) {
    auto x_max_derivative = minimize(a, 0.1, a, b, epsilon, [&f](auto&& x) {
        return -std::abs(diff(x, 1e-10, f));
    });

    auto max_derivative = std::abs(diff(x_max_derivative, 1e-10, f));

    auto lambda = [&max_derivative, &f](auto&& x) {
        return std::copysign(1.0 / max_derivative, diff(x, 1e-10, f));
    };

    auto phi = [&f, &lambda](auto&& x) { return x - lambda(x) * f(x); };

    auto x_q = minimize(a, 0.1, a, b, epsilon, [&phi](auto&& x) {
        return -std::abs(diff(x, 1e-10, phi));
    });

    auto q = std::abs(diff(x_q, 1e-10, phi));

    auto x_k = x0;
    auto x_k_1 = std::numeric_limits<T>::max();

    while (true) {
        x_k_1 = x_k - lambda(x_k) * f(x_k);
        if (q / (1 - q) * std::abs(x_k_1 - x_k) < epsilon) {
            break;
        }
        x_k = x_k_1;
    }

    return x_k_1;
}

// **********************************************************************
// **************************** Newton method ***************************
// **********************************************************************

template <class T, class F>
T newton_solver(T x0, T epsilon, F&& f) {
    auto x_k = x0;
    auto x_k_1 = x0 + 1e-2;
    auto x_k_2 = x0;

    while (true) {
        x_k_2 = x_k_1 - f(x_k_1) * (x_k_1 - x_k) / (f(x_k_1) - f(x_k));
        if (std::abs(x_k_2 - x_k_1) < epsilon) {
            break;
        }
        x_k = x_k_1;
        x_k_1 = x_k_2;
    }

    return x_k_2;
}

// **********************************************************************
// ************************ Fixed point iteration ***********************
// **************************** system solver ***************************
// **********************************************************************

template <class T, class F, size_type N>
matrix<T> fixed_point_iteration_solver(const matrix<T>& x0,
                                       T epsilon,
                                       std::array<F, N> f) {
    auto J = jacobi_matrix(x0, epsilon, f);

    lup_decomposition<T> lup(J);
    auto Lambda = lup.inverse();

    auto phi_v = [&Lambda, &f](auto&& x) {
        auto f_v = make_column<T>(N);
        for (index_type i = 0; i < N; i++) {
            f_v(i) = f[i](x);
        }

        return x - Lambda * f_v;
    };

    auto phi = std::array<F, N>{};
    for (index_type i = 0; i < N; i++) {
        phi[i] = [i, &phi_v](auto&& x) { return phi_v(x)(i); };
    }

    auto x_k = x0;
    auto x_k_1 = x0;
    while (true) {
        x_k_1 = phi_v(x_k);
        if (norm_inf(x_k_1 - x_k) < epsilon) {
            break;
        }
        x_k = x_k_1;
    }

    return x_k_1;
}

// **********************************************************************
// *************************** Newton method ****************************
// *************************** system solver ****************************
// **********************************************************************

template <class T, class F, size_type N>
matrix<T> newton_solver(const matrix<T>& x0,
                                       T epsilon,
                                       std::array<F, N> f) {
    auto J = jacobi_matrix(x0, epsilon, f);

    lup_decomposition<T> lup(J);
    auto Jinv = lup.inverse();

    auto f_v = [&f](auto&& x) {
        auto result = make_column<T>(N);
        for (index_type i = 0; i < N; i++) {
            result(i) = f[i](x);
        }

        return result;
    };

    auto x_k = x0;
    auto x_k_1 = x0;
    while (true) {
        x_k_1 = x_k - Jinv * f_v(x_k);
        if (norm_inf(x_k_1 - x_k) < epsilon) {
            break;
        }
        x_k = x_k_1;
    }

    return x_k_1;
}

}  // namespace chislaki

#endif  // CHISLAKI_TRANSCENDENTAL_UTILITY_HPP_
