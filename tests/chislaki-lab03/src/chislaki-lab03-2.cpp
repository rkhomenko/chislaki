#include <chislaki/matan/interpolation.hpp>

#include <iostream>

using namespace chislaki;

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

    auto s = cubic_spline(x, y);

    std::cout << s(1.5) << std::endl;
}

