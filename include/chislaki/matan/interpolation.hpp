#ifndef CHISLAKI_MATAN_INTERPOLATION_HPP_
#define CHISLAKI_MATAN_INTERPOLATION_HPP_

#include <chislaki/linalg/matrix.hpp>

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
        for (index_type i = 0; i < diff_.size(); i++) {
            for (index_type j = 0; j < diff_[i].size(); j++) {
                std::cout << diff_[i][j] << " ";
            }
            std::cout << std::endl;
        }
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

}  // namespace chislaki

#endif  // CHISLAKI_MATAN_INTERPOLATION_HPP_
