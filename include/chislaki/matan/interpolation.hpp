#ifndef CHISLAKI_MATAN_INTERPOLATION_HPP_
#define CHISLAKI_MATAN_INTERPOLATION_HPP_

#include <chislaki/linalg/decompositions.hpp>
#include <chislaki/linalg/matrix.hpp>
#include <chislaki/linalg/utility.hpp>

#include <cmath>

namespace chislaki {

// **********************************************************************
// ************************** Lagrange polynom **************************
// **********************************************************************

template <class T>
class lagrange_polinom {
public:
    lagrange_polinom(const matrix<T>& x, const matrix<T>& y) : x_{x}, y_{y} {}

    T operator()(T x) {
        T result = 0;
        for (index_type i = 0; i < x_.rows(); i++) {
            result += y_(i) * l(x, i);
        }
        return result;
    }

private:
    T l(T x, index_type i) {
        T result = 1;
        for (index_type j = 0; j < i; j++) {
            result *= (x - x_(j)) / (x_(i) - x_(j));
        }
        for (index_type j = i + 1; j < x_.rows(); j++) {
            result *= (x - x_(j)) / (x_(i) - x_(j));
        }

        return result;
    }

    matrix<T> x_;
    matrix<T> y_;
};

// **********************************************************************
// *************************** Newton polynom ***************************
// **********************************************************************

template <class T>
class newton_polynom {
public:
    newton_polynom(const matrix<T>& x, const matrix<T>& y) : x_{x}, y_{y} {
        calc_diff();
    }

    T operator()(T x) {
        T result = y_(0);

        for (index_type i = 1; i < x_.rows(); i++) {
            T prod = 1;
            for (index_type j = 0; j < i; j++) {
                prod *= (x - x_(j));
            }
            result += prod * diff_[i - 1][0];
        }

        return result;
    }

private:
    void calc_diff() {
        diff_.push_back({});
        for (index_type i = 0; i < x_.rows() - 1; i++) {
            diff_[0].push_back((y_(i + 1) - y_(i)) / (x_(i + 1) - x_(i)));
        }

        for (index_type i = 1; i < x_.rows() - 1; i++) {
            diff_.push_back({});
            for (index_type j = 0; j < x_.rows() - i - 1; j++) {
                diff_[i].push_back((diff_[i - 1][j + 1] - diff_[i - 1][j]) /
                                   (x_(j + i + 1) - x_(j)));
            }
        }
    }

    matrix<T> x_;
    matrix<T> y_;
    std::vector<std::vector<T>> diff_;
};

// **********************************************************************
// **************************** Cubic spline ****************************
// **********************************************************************

template <class T>
class cubic_spline {
public:
    cubic_spline(const matrix<T>& x, const matrix<T>& y) : x_{x}, y_{y} {
        calculate_coeffs();
    }

    T operator()(T x) {
        index_type index = 0;
        for (index_type i = 1; i < x_.rows(); i++) {
            if (x_(i - 1) <= x && x <= x_(i)) {
                index = i;
                break;
            }
        }

        return coeffs_[index - 1][0] +
               coeffs_[index - 1][1] * (x - x_(index - 1)) +
               coeffs_[index - 1][2] * std::pow(x - x_(index - 1), 2) +
               coeffs_[index - 1][3] * std::pow(x - x_(index - 1), 3);
    }

private:
    void calculate_coeffs() {
        auto matr = matrix<T>(x_.rows() - 2);
        auto h = [this](auto&& i) { return x_(i) - x_(i - 1); };

        auto b = make_column<T>(x_.rows() - 2);
        for (index_type i = 0; i < b.rows(); i++) {
            b(i) = 3 * ((y_(i + 2) - y_(i + 1)) / h(i + 2) -
                        (y_(i + 1) - y_(i)) / h(i + 1));
        }

        auto n = matr.rows();

        matr(0, 0) = 2 * (h(1) + h(2));
        matr(0, 1) = h(2);

        for (index_type i = 1; i < matr.rows() - 1; i++) {
            matr(i, i - 1) = h(i);
            matr(i, i) = 2 * (h(i) + h(i + 1));
            matr(i, i + 1) = h(i + 1);
        }

        matr(n - 1, n - 2) = h(n);
        matr(n - 1, n - 1) = 2 * (h(n + 1) + h(n));

        auto c = thomas_algorithm(matr, b);

        coeffs_.resize(n + 1);

        // a_i
        for (index_type i = 0; i < coeffs_.size(); i++) {
            coeffs_[i].resize(4);
            coeffs_[i][0] = y_(i);
        }

        auto c_i = [&c](auto&& i) {
            if (i == 0) {
                return T{0};
            }
            return c(i - 1);
        };

        // b_i
        for (index_type i = 0; i < coeffs_.size() - 1; i++) {
            coeffs_[i][1] = (y_(i + 1) - y_(i)) / h(i + 1) -
                            1.0 / 3 * h(i + 1) * (c_i(i + 1) + 2 * c_i(i));
        }

        // c_i
        for (index_type i = 0; i < coeffs_.size(); i++) {
            coeffs_[i][2] = c_i(i);
        }

        // d_i
        for (index_type i = 0; i < coeffs_.size() - 1; i++) {
            coeffs_[i][3] =
                (coeffs_[i + 1][2] - coeffs_[i][2]) / (3 * h(i + 1));
        }

        coeffs_[n][1] =
            (y_(n) - y_(n - 1)) / h(n) - 1.0 / 3 * h(n) * coeffs_[n][2];
        coeffs_[n][3] = -coeffs_[n][2] / (3 * h(n));

        for (const auto& v : coeffs_) {
            for (const auto& item : v) {
                std::cout << item << " ";
            }
            std::cout << std::endl;
        }
    }

    matrix<T> x_;
    matrix<T> y_;
    std::vector<std::vector<T>> coeffs_;
};

// **********************************************************************
// ************************ Linear least squares ************************
// **********************************************************************

template <class T>
class linear_least_squares {
public:
    linear_least_squares(const matrix<T>& x, const matrix<T>& y)
        : x_{x}, y_{y} {
        calculate_coeffs();
    }

    T operator()(T x) { return coeffs_(0) + coeffs_(1) * x; }

    const matrix<T>& get_coeffs() const { return coeffs_; }

private:
    void calculate_coeffs() {
        auto x_sum = [&](auto&& n) {
            T result = 0;
            for (index_type i = 0; i < x_.rows(); i++) {
                result += std::pow(x_(i), n);
            }
            return result;
        };

        auto x_y_sum = [&](auto&& n) {
            T result = 0;
            for (index_type i = 0; i < x_.rows(); i++) {
                result += y_(i) * std::pow(x_(i), n);
            }
            return result;
        };

        auto matr = matrix<T>(2);
        for (index_type i = 0; i < matr.rows(); i++) {
            for (index_type j = 0; j < matr.columns(); j++) {
                matr(i, j) = x_sum(i + j);
            }
        }

        auto b = make_column<T>(2);
        b(0) = x_y_sum(0);
        b(1) = x_y_sum(1);

        coeffs_ = lup_decomposition(matr).solve(b);
    }

    matrix<T> x_;
    matrix<T> y_;
    matrix<T> coeffs_;
};

// **********************************************************************
// ************************ Square least squares ************************
// **********************************************************************

template <class T>
class square_least_squares {
public:
    square_least_squares(const matrix<T>& x, const matrix<T>& y)
        : x_{x}, y_{y} {
        calculate_coeffs();
    }

    T operator()(T x) {
        return coeffs_(0) + coeffs_(1) * x + coeffs_(2) * x * x;
    }

    const matrix<T>& get_coeffs() const { return coeffs_; }

private:
    void calculate_coeffs() {
        auto x_sum = [&](auto&& n) {
            T result = 0;
            for (index_type i = 0; i < x_.rows(); i++) {
                result += std::pow(x_(i), n);
            }
            return result;
        };

        auto x_y_sum = [&](auto&& n) {
            T result = 0;
            for (index_type i = 0; i < x_.rows(); i++) {
                result += y_(i) * std::pow(x_(i), n);
            }
            return result;
        };

        auto matr = matrix<T>(3);
        for (index_type i = 0; i < matr.rows(); i++) {
            for (index_type j = 0; j < matr.columns(); j++) {
                matr(i, j) = x_sum(i + j);
            }
        }

        auto b = make_column<T>(3);
        for (index_type i = 0; i < b.rows(); i++) {
            b(i) = x_y_sum(i);
        }

        coeffs_ = lup_decomposition(matr).solve(b);
    }

    matrix<T> x_;
    matrix<T> y_;
    matrix<T> coeffs_;
};

}  // namespace chislaki

#endif  // CHISLAKI_MATAN_INTERPOLATION_HPP_
