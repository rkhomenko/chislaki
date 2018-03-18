#ifndef CHISLAKI_LINALG_UTILITY_HPP_
#define CHISLAKI_LINALG_UTILITY_HPP_

#include <chislaki/linalg/matrix.hpp>

#include <cmath>

namespace chislaki {

// **********************************************************************
// ******************** Tridiagonal matrix algorithm ********************
// ************************* (Thomas algorithm) *************************
// **********************************************************************

template <class T>
matrix<T> thomas_algorithm(const matrix<T>& matr, const matrix<T>& d) {
    if (matr.rows() != d.rows()) {
        throw bad_size{};
    }

    auto a = make_row<T>(d.rows());
    auto b = make_row<T>(d.rows());

    // a(0) = -c_1 / b_1, b(0) = d_1 / b_1
    a(0) = -matr(0, 1) / matr(0, 0);
    b(0) = d(0) / matr(0, 0);

    for (index_type i = 1; i < matr.rows() - 1; i++) {
        const auto& a_i = matr(i, i - 1);
        const auto& b_i = matr(i, i);
        const auto& c_i = matr(i, i + 1);

        a(i) = -c_i / (b_i + a_i * a(i - 1));
        b(i) = (d(i) - a_i * b(i - 1)) / (b_i + a_i * a(i - 1));
    }

    auto n = a.columns();
    a(n - 1) = 0;
    b(n - 1) = (d(n - 1) - matr(n - 1, n - 2) * b(n - 2)) /
               (matr(n - 1, n - 1) + matr(n - 1, n - 2) * a(n - 2));

    auto x = make_column<T>(d.rows());
    n = x.rows();

    x(n - 1) = b(n - 1);
    for (index_type i = n - 1; i != 0; i--) {
        // x_n = a_n * x_n+1 + b_n
        x(i - 1) = a(i - 1) * x(i) + b(i - 1);
    }

    return x;
}

// **********************************************************************
// ************************ Jacobi transformation ***********************
// **********************************************************************

template <class T>
std::pair<matrix<T>, matrix<T>> jacobi_transform(const matrix<T>& matr,
                                                 const matrix<T>& b) {
    if (matr.rows() != b.rows()) {
        throw bad_size{};
    }

    auto betta = make_column<T>(b.rows());
    auto alpha = matrix<T>(matr.rows(), matr.columns());

    for (index_type i = 0; i < b.rows(); i++) {
        betta(i) = b(i) / matr(i, i);
    }

    for (index_type i = 0; i < matr.rows(); i++) {
        for (index_type j = 0; j < matr.columns(); j++) {
            if (i != j) {
                alpha(i, j) = -matr(i, j) / matr(i, i);
            } else {
                alpha(i, j) = 0;
            }
        }
    }

    return std::make_pair(alpha, betta);
}

// **********************************************************************
// ************************ Fixed point iteration ***********************
// **********************************************************************

template <class T>
matrix<T> fixed_point_iteration(const matrix<T>& matr,
                                const matrix<T>& b,
                                T epsilon) {
    auto[alpha, betta] = jacobi_transform(matr, b);

    if (norm_inf(alpha) < 1) {
        auto iter_count = (std::log10(epsilon) - std::log10(norm_inf(betta)) +
                           std::log10(1 - norm_inf(alpha))) /
                              std::log10(norm_inf(alpha)) -
                          1;

        std::cout << "Expected iteration count: " << std::ceil(iter_count)
                  << std::endl;
    } else {
        std::cout << "Cannot expect iteration count! ||alpha|| = 1!"
                  << std::endl;
    }

    auto x_k_1 = betta;
    auto x_k = make_column<T>(betta.rows());
    T norm = 0;
    size_type count = 1;
    while ((norm = norm_inf(x_k - x_k_1)) && norm > epsilon) {
        auto tmp = x_k;
        x_k = betta + alpha * x_k_1;
        x_k_1 = tmp;
        count++;
        std::cout << "||x_k - x_k-1|| = " << norm_inf(x_k - x_k_1) << std::endl;
        std::cout << "x_k:\n" << x_k;
    }

    std::cout << "Iteration count: " << count << std::endl;

    return x_k;
}

// **********************************************************************
// ************************* Gauss-Seidel method ************************
// **********************************************************************

template <class T>
matrix<T> gauss_seidel(const matrix<T>& matr, const matrix<T>& b, T epsilon) {
    auto[alpha, betta] = jacobi_transform(matr, b);

    auto x_k_1 = betta;
    auto x_k = make_column<T>(betta.rows());
    T norm = 0;
    size_type count = 0;
    while ((norm = norm_inf(x_k - x_k_1)) && norm > epsilon) {
        auto tmp = x_k;
        for (size_type i = 0; i < x_k.rows(); i++) {
            x_k(i) = betta(i);
            for (size_type k = 0; k < i; k++) {
                x_k(i) += alpha(i, k) * x_k(k);
            }
            for (size_type k = i; k < x_k_1.rows(); k++) {
                x_k(i) += alpha(i, k) * x_k_1(k);
            }
        }
        x_k_1 = tmp;
        std::cout << "||x_k - x_k-1|| = " << norm_inf(x_k - x_k_1) << std::endl;
        std::cout << "x_k:\n" << x_k;
        count++;
    }

    std::cout << "Iteration count: " << count << std::endl;

    return x_k;
}

// **********************************************************************
// *************************** Max for matrix ***************************
// **********************************************************************

template <class T>
std::tuple<T, index_type, index_type> max_nondiag(const matrix<T>& matr) {
    T max = 0;
    index_type max_i = 0;
    index_type max_j = 0;
    for (index_type i = 0; i < matr.rows(); i++) {
        for (index_type j = 0; j < matr.columns(); j++) {
            T value = std::abs(matr(i, j));
            if (value > max && i != j) {
                max = value;
                max_i = i;
                max_j = j;
            }
        }
    }

    return std::make_tuple(max, max_i, max_j);
}

// **********************************************************************
// ********************* Jacobi Eigenvalue Algorithm ********************
// **********************************************************************

template <class T>
class jacobi_eigenvalues {
public:
    using value_type = matrix<T>;
    using const_reference = const value_type&;

    // **********************************************************************
    // **************************** Constructors ****************************
    // **********************************************************************

    jacobi_eigenvalues() : computed_{false}, epsilon_{1e-6} {}

    jacobi_eigenvalues(const_reference matr, T epsilon) : jacobi_eigenvalues() {
        compute(matr, epsilon);
        computed_ = true;
    }

    // **********************************************************************
    // *********************** Eigenvalues calculation **********************
    // **********************************************************************

    void compute(const_reference matr, T epsilon) {
        epsilon_ = epsilon;

        auto off = [](auto&& matr) {
            T sum = 0;

            for (index_type i = 0; i < matr.rows(); i++) {
                for (index_type j = i + 1; j < matr.columns(); j++) {
                    T value = std::abs(matr(i, j));
                    sum += value * value;
                }
            }

            return std::sqrt(sum);
        };

        auto generate_u = [](auto&& matr) {
            auto[m, i, j] = max_nondiag(matr);
            auto u = make_eye<T>(matr.rows());

            T phi = (matr(i, i) - matr(j, j) == 0)
                        ? std::atan(1.0)
                        : 0.5 * std::atan(2 * matr(i, j) /
                                          (matr(i, i) - matr(j, j)));

            u(i, i) = std::cos(phi);
            u(j, i) = std::sin(phi);
            u(i, j) = -std::sin(phi);
            u(j, j) = std::cos(phi);

            return u;
        };

        a_ = matr;
        u_ = make_eye<T>(matr.rows());

        auto u_i = std::vector<matrix<T>>();
        while (off(a_) > get_epsilon()) {
            u_ = generate_u(a_);
            a_ = transpose(u_) * a_ * u_;
            u_i.push_back(u_);
        }

        u_ = make_eye<T>(matr.rows());
        for (index_type i = 0; i < u_i.size(); i++) {
            u_ = u_ * u_i[i];
        }
    }

    // **********************************************************************
    // ****************** Eigenvalues, eigenvectors access ******************
    // **********************************************************************

    inline size_type count() const noexcept { return a_.rows(); }

    T get_value(index_type i) const {
        if (i >= count()) {
            throw bad_index{};
        }
        return a_(i, i);
    }

    matrix<T> get_vector(index_type i) const {
        if (i >= count()) {
            throw bad_index{};
        }
        return u_.column(i);
    }

private:
    // **********************************************************************
    // ********************** Private member functions **********************
    // **********************************************************************

    inline bool is_computed() const noexcept { return computed_; }
    inline T get_epsilon() const noexcept { return epsilon_; }

    // **********************************************************************
    // ************************** Member variables **************************
    // **********************************************************************

    bool computed_;
    T epsilon_;
    matrix<T> u_;
    matrix<T> a_;
};  // class jacobi_eigenvalues

}  // namespace chislaki

#endif  // CHISLAKI_LINALG_UTILITY_HPP_
