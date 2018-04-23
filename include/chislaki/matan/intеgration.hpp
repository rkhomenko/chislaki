#ifndef CHISLAKI_MATAN_INTEGRATION_HPP_
#define CHISLAKI_MATAN_INTEGRATION_HPP_

#include <chislaki/linalg/matrix.hpp>

namespace chislaki {

template <class T, class F>
T integrate_riemann(F&& f, T a, T b, T h) {
    T result = 0;
    index_type i = 0;
    while (a + i * h < b) {
        result += h * f((h * i + a + h * (i + 1) + a) / 2.);
        i++;
    }

    return result;
}

template <class T, class F>
T integrate_trapezoidal(F&& f, T a, T b, T h) {
    T result = 0;
    index_type i = 0;
    while (a + i * h < b) {
        result += (f(a + i * h) + f(a + (i + 1) * h)) * h;
        i++;
    }

    return 0.5 * result;
}

template <class T, class F>
T integrate_simpson(F&& f, T a, T b, T h) {
    T I2 = {};
    T I4 = f(a + h);

    for (index_type i = 2; i < static_cast<index_type>((b - a) / h); i += 2) {
        I4 += f(a + (i + 1) * h);
        I2 += f(a + i * h);
    }

    return h / 3 * (f(a) + f(b) + 4 * I4 + 2 * I2);
}

}  // namespace chislaki

#endif  // CHISLAKI_MATAN_INTEGRATION_HPP_
