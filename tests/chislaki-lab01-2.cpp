#include <chislaki/linalg/utility.hpp>

using namespace chislaki;

int main() {
    size_type n;

    std::cin >> n;

    auto A = matrix<long double>(n);
    auto d = make_column<long double>(n);

    std::cin >> A(0, 0) >> A(0, 1);
    for (size_type i = 1; i < n - 1; i++) {
        std::cin >> A(i, i - 1) >> A(i, i) >> A(i, i + 1);
    }
    std::cin >> A(n - 1, n - 2) >> A(n - 1, n - 1);

    for (size_type i = 0; i < n; i++) {
        std::cin >> d(i);
    }

    auto x = thomas_algorithm(A, d);

    std::cout << "******************************** A, b "
                 "********************************\n"
              << "A:\n"
              << A << "d:\n"
              << d
              << "************************** Thomas algorithm "
                 "**************************\n"
              << "x:\n"
              << x << "Checking: Ax = d\n"
              << A * x;
}
