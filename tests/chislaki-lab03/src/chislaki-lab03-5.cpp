#include <chislaki/matan/intеgration.hpp>

#include <cmath>

using namespace chislaki;

template <class T>
T runge(T Ih, T Ikh, index_type k) {
    return Ih + (Ih - Ikh) / (k * k - 1);
}

int main() {
    auto f = [](auto&& x) { return x / std::pow(3 * x + 4, 2); };

    auto R1 = integrate_riemann(f, 0.0, 4.0, 0.5);
    auto R2 = integrate_riemann(f, 0.0, 4.0, 1.0);
    auto T1 = integrate_trapezoidal(f, 0.0, 4.0, 0.5);
    auto T2 = integrate_trapezoidal(f, 0.0, 4.0, 1.0);
    auto S1 = integrate_simpson(f, 0.0, 4.0, 0.5);
    auto S2 = integrate_simpson(f, 0.0, 4.0, 1.0);

    std::cout << "Riemann: " << R1 << ", " << R2
              << ", err = " << runge(R1, R2, 2) << std::endl;
    std::cout << "Trapezoidal: " << T1 << ", " << T2
              << ", err = " << runge(T1, T2, 2) << std::endl;
    std::cout << "Simpson: " << S1 << ", " << S2
              << ", err = " << runge(S1, S2, 2) << std::endl;
}
