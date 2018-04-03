#ifndef CHISLAKI_OPTIMIZE_MINIMIZE_HPP_
#define CHISLAKI_OPTIMIZE_MINIMIZE_HPP_

#include <chislaki/matan/differentiation.hpp>

#include <algorithm>
#include <limits>
#include <set>

namespace chislaki {

template <class T, class F>
T minimize(T x0, T step, T a, T b, T epsilon, F&& f) {
    auto x_k = x0;
    auto x_k_1 = std::numeric_limits<T>::max();

    T result = 0;
    while (true) {
        x_k_1 = x_k - step * diff(x_k, 1e-10, f);
        if (!(a < x_k_1 && x_k_1 < b)) {
            result = (f(a) < f(b)) ? a : b;
            break;
        }
        if (std::abs(x_k_1 - x_k) < epsilon) {
            result = x_k_1;
            break;
        }
        x_k = x_k_1;
    }

    return result;
}

}  // namespace chislaki

#endif  // CHISLAKI_OPTIMIZE_MINIMIZE_HPP_
