#ifndef CHISLAKI_MATAN_DIFFERENTIATION_HPP_
#define CHISLAKI_MATAN_DIFFERENTIATION_HPP_

namespace chislaki {

template <class T, class F>
T diff(T x, T delta, F&& f) {
    return (f(x + delta) - f(x)) / delta;
}

}  // namespace chislaki

#endif  // CHISLAKI_MATHAN_DIFFERENTIATION_HPP_
