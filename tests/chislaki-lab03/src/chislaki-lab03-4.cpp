#include <chislaki/matan/differentiation.hpp>

#include <iostream>

using namespace chislaki;

int main() {
    size_type n;
    long double x0{};

    std::cin >> n >> x0;

    auto x = make_column<long double>(n);
    auto y = x;

    for (index_type i = 0; i < x.rows(); i++) {
        std::cin >> x(i);
    }

    for (index_type i = 0; i < y.rows(); i++) {
        std::cin >> y(i);
    }

    std::cout << "y'(x0) = " << diff1(x0, x, y) << std::endl
              << "y''(x0) = " << diff2(x0, x, y) << std::endl;
}
