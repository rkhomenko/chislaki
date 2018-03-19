#ifndef CHISLAKI_LINALG_DECOMPOSITIONS_HPP_
#define CHISLAKI_LINALG_DECOMPOSITIONS_HPP_

#include <chislaki/linalg/matrix.hpp>

#include <cmath>
#include <tuple>

namespace chislaki {

// Exceptions
class not_computed {};

// ***********************************************************************
// ************************** LUP decomposition **************************
// ***********************************************************************

template <class T>
class lup_decomposition {
public:
    using value_type = matrix<T>;
    using const_reference = const value_type&;

    // **********************************************************************
    // **************************** Constructors ****************************
    // **********************************************************************

    lup_decomposition() : computed_{false}, permutations_count_{0} {}

    lup_decomposition(const_reference matr) : lup_decomposition() {
        compute(matr);
        computed_ = true;
    }

    // **********************************************************************
    // ******************** LU decomposition calcualtion ********************
    // **********************************************************************

    void compute(const_reference matr) {
        if (matr.rows() != matr.columns()) {
            throw bad_size{};
        }

        u_ = matr;
        l_ = make_eye<T>(u_.rows());
        p_ = l_;

        auto max = [](auto&& matr, auto&& col_index, auto&& row_current) {
            index_type max_index = 0;
            T max_value = 0;
            for (index_type i = row_current; i < matr.rows(); i++) {
                auto current_value = std::abs(matr(i, col_index));
                if (current_value > max_value) {
                    max_value = current_value;
                    max_index = i;
                }
            }
            return max_index;
        };

        for (index_type j = 0; j < u_.columns(); j++) {
            auto max_index = max(u_, j, j);
            auto need_swap = j != max_index;

            if (need_swap) {
                u_.swap_rows(j, max_index);
                p_.swap_rows(j, max_index);
                permutations_count_++;
            }

            for (index_type i = j + 1; i < u_.rows(); i++) {
                l_(i, j) = u_(i, j) / u_(j, j);

                for (index_type k = j; k < u_.columns(); k++) {
                    u_(i, k) -= l_(i, j) * u_(j, k);
                }
            }

            if (need_swap) {
                for (index_type k = 0; k < j; k++) {
                    std::swap(l_(j, k), l_(max_index, k));
                }
            }
        }
    }

    // ***********************************************************************
    // *********************** L, U, P matrices access ***********************
    // ***********************************************************************

    inline const_reference matrix_l() const noexcept { return l_; }

    inline const_reference matrix_u() const noexcept { return u_; }

    inline const_reference matrix_p() const noexcept { return p_; }

    inline size_type permutations_count() const noexcept {
        return permutations_count_;
    }

    // **********************************************************************
    // ***************************** Algorithms *****************************
    // **********************************************************************

    value_type solve(const_reference b) const {
        if (!is_computed()) {
            throw not_computed{};
        }

        const auto& l = matrix_l();
        const auto& u = matrix_u();
        const auto& p = matrix_p();

        if (u.rows() != b.rows()) {
            throw bad_size{};
        }

        // Ax = b, PA = LU =>
        // PAx = LUx = Pb, Lz = Pb, Ux = z

        auto z = make_column<T>(u.columns());
        auto b1 = p * b;

        z(0) = b1(0);
        z(0) = b1(0);
        for (index_type i = 1; i < u.rows(); i++) {
            T sum = 0;
            for (index_type j = 0; j < i; j++) {
                sum += l(i, j) * z(j);
            }
            z(i) = b1(i) - sum;
        }

        // Ux = z
        auto x = make_column<T>(u.columns());
        const auto n = x.rows();
        x(n - 1) = z(n - 1) / u(n - 1, n - 1);
        for (index_type i = n - 1; i != 0; i--) {
            T sum = 0;
            for (index_type j = (i - 1) + 1; j < n; j++) {
                sum += u(i - 1, j) * x(j);
            }
            x(i - 1) = (z(i - 1) - sum) / u(i - 1, i - 1);
        }

        return x;
    }

    T determinant() const {
        if (!is_computed()) {
            throw not_computed{};
        }

        const auto& u = matrix_u();
        T result = 1;
        for (index_type i = 0; i < u.rows(); i++) {
            result *= u(i, i);
        }
        return result * (permutations_count() % 2 == 0 ? 1 : -1);
    }

    value_type inverse() const {
        if (!is_computed()) {
            throw not_computed{};
        }

        if (determinant() == 0) {
            throw degenerate{};
        }

        auto result = make_eye<T>(matrix_u().rows());

        for (index_type j = 0; j < result.columns(); j++) {
            auto b = make_column(j, result);
            result.set_column(j, solve(b));
        }

        return result;
    }

private:
    // **********************************************************************
    // ********************** Private member functions **********************
    // **********************************************************************

    inline bool is_computed() const noexcept { return computed_; }

    // **********************************************************************
    // ************************** Member variables **************************
    // **********************************************************************

    bool computed_;
    matrix<T> l_;
    matrix<T> u_;
    matrix<T> p_;
    size_type permutations_count_;
};  // class lup_decomposition

// **********************************************************************
// ************************** QR decomposition **************************
// **********************************************************************

template <class T>
class qr_decomposition {
public:
    using value_type = matrix<T>;
    using const_reference = const value_type&;

    // **********************************************************************
    // **************************** Constructors ****************************
    // **********************************************************************

    qr_decomposition() : computed_{false} {}

    qr_decomposition(const_reference matr) : qr_decomposition() {
        compute(matr);
    }

    // **********************************************************************
    // ******************** QR decomposition calcualtion ********************
    // **********************************************************************

    void compute(const_reference matr) {
        auto A = matr;
        std::vector<matrix<T>> h_i;

        for (index_type i = 0; i < matr.rows() - 1; i++) {
            auto H = householder(A, i);
            A = H * A;
            h_i.push_back(H);
        }
        r_ = A;

        q_ = make_eye<T>(matr.rows());
        for (auto&& H : h_i) {
            q_ = q_ * H;
        }

        computed_ = true;
    }

    // **********************************************************************
    // ************************* Q, R matrix access *************************
    // **********************************************************************

    const_reference matrix_q() const {
        if (!is_computed()) {
            throw not_computed{};
        }

        return q_;
    }

    const_reference matrix_r() const {
        if (!is_computed()) {
            throw not_computed{};
        }
        return r_;
    }

    // **********************************************************************
    // ************************* Householder matrix *************************
    // **********************************************************************

    static matrix<T> householder(const matrix<T>& matr, index_type index) {
        auto sign = [](T value) {
            return static_cast<int>((T(0) < value) - (value < T(0)));
        };

        auto calculate_v = [&sign](auto&& matr, index_type index) {
            auto v = make_column<T>(matr.rows());

            for (index_type i = 0; i < index; i++) {
                v(i) = 0;
            }

            T sum = 0;
            for (index_type i = index; i < matr.rows(); i++) {
                sum += std::pow(matr(i, index), 2);
            }

            v(index) =
                matr(index, index) + sign(matr(index, index)) * std::sqrt(sum);

            for (index_type i = index + 1; i < matr.rows(); i++) {
                v(i) = matr(i, index);
            }

            return v;
        };

        auto v = calculate_v(matr, index);

        auto E = make_eye<T>(v.rows());
        auto value = (v * transpose(v)) / (transpose(v) * v)(0);

        return E - static_cast<T>(2) * value * E;
    }

private:
    // **********************************************************************
    // ********************** Private member functions **********************
    // **********************************************************************

    inline bool is_computed() const noexcept { return computed_; }

    // **********************************************************************
    // ************************** Member variables **************************
    // **********************************************************************

    bool computed_;
    matrix<T> q_;
    matrix<T> r_;
};  // class qr_decomposition

}  // namespace chislaki

#endif  // CHISLAKI_LINALG_DECOMPOSITIONS_HPP_
