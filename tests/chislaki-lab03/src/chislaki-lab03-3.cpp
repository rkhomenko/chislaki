#include <chislaki/matan/interpolation.hpp>

#include <iostream>

using namespace chislaki;

template <class F>
double err_calc(const matrix<double>& x, const matrix<double>& y, F&& f) {
    double err = 0;
    for (index_type i = 0; i < x.rows(); i++) {
        err += std::pow(f(x(i)) - y(i), 2);
    }
    return err;
}

int main() {
    size_type n;

    std::cin >> n;

    auto x = make_column<double>(n);
    auto y = x;

    for (index_type i = 0; i < x.rows(); i++) {
        std::cin >> x(i);
    }

    for (index_type i = 0; i < y.rows(); i++) {
        std::cin >> y(i);
    }

    auto lls = linear_least_squares(x, y);
    auto sls = square_least_squares(x, y);

    std::cout << "linear coeffs:\n"
              << lls.get_coeffs() << "err: " << err_calc(x, y, lls) << std::endl
              << "square coeffs:\n"
              << sls.get_coeffs() << "err: " << err_calc(x, y, sls)
              << std::endl;
}
