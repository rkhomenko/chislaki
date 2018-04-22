#include <chislaki/matan/interpolation.hpp>

#include <cmath>

using namespace chislaki;

template <class F>
matrix<double> apply_functor(const matrix<double>& x, F&& f) {
    auto y = make_column<double>(x.rows());

    for (index_type i = 0; i < x.rows(); i++) {
        y(i) = f(x(i));
    }

    return y;
}

matrix<double> chebyshev_nodes(size_type n) {
    auto pi = 4 * std::atan(1.0);
    auto result = make_column<double>(n);

    for (index_type i = 1; i <= n; i++) {
        result(i - 1) = std::cos((2 * i - 1) * pi / (2 * n));
    }

    return result;
}

int main() {
    auto pi = 4 * std::atan(1.0);
    auto x0 = pi / 4;
    auto x1 = make_column<double>(
        4, std::vector<double>{0, pi / 6, pi / 3, pi / 2});
    auto x2 = make_column<double>(
        4, std::vector<double>{0, pi / 6, 5 * pi / 12, pi / 2});
    auto x_che = chebyshev_nodes(20);

    auto y1 = apply_functor(x1, [](auto&& x) { return std::cos(x); });
    auto y2 = apply_functor(x2, [](auto&& x) { return std::cos(x); });
    auto y_che = apply_functor(x_che, [](auto&& x) { return std::cos(x); });

    auto l1 = lagrange_polinom(x1, y1);
    auto l2 = lagrange_polinom(x2, y2);
    auto l_che = lagrange_polinom(x_che, y_che);
    auto n1 = newton_polynom(x1, y1);
    auto n2 = newton_polynom(x2, y2);
    auto n_che = newton_polynom(x_che, y_che);

    std::cout << "************************** Lagrange polynom "
                 "**************************\n"
              << "y1:\n"
              << y1 << "y1 by lagrange:\n"
              << apply_functor(x1, l1) << "y2:\n"
              << y2 << "y2 by lagrange:\n"
              << apply_functor(x2, l2) << "y1 chebyshev:\n"
              << apply_functor(x1, l_che) << "y2 chebyshev:\n"
              << apply_functor(x2, l_che)
              << "f(x*) - L1(x*) = " << std::cos(x0) - l1(x0) << std::endl
              << "f(x*) - L2(x*) = " << std::cos(x0) - l2(x0) << std::endl
              << "f(x*) - LChe(x*) = " << std::cos(x0) - l_che(x0) << std::endl
              << "*************************** Newton polynom "
                 "***************************\n"
              << "y1 by newton:\n"
              << apply_functor(x1, n1) << "y2 by newton:\n"
              << apply_functor(x2, n2) << "y1 chebyshev:\n"
              << apply_functor(x1, n_che) << "y2 chebyshev:\n"
              << apply_functor(x2, n_che)
              << "f(x*) - N1(x*) = " << std::cos(x0) - n1(x0) << std::endl
              << "f(x*) - N2(x*) = " << std::cos(x0) - n2(x0) << std::endl
              << "f(x*) - NChe(x*) = " << std::cos(x0) - n_che(x0) << std::endl;
}
