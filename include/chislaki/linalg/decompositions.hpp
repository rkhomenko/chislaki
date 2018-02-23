#ifndef CHISLAKI_LINALG_DECOMPOSITION_HPP_
#define CHISLAKI_LINALG_DECOMPOSITION_HPP_

#include <chislaki/linalg/matrix.hpp>

#include <tuple>

namespace chislaki {
template <class T>
std::tuple<matrix<T>, matrix<T>, matrix<T>> LUP_decomposition(
    const matrix<T>& matr) {
    using index_type = typename matrix<T>::index_type;

    if (matr.rows() != matr.columns()) {
        throw bad_size{};
    }

    auto U = matr;
    auto L = matrix<T>::eye(U.rows());
    auto P = L;

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

    for (index_type j = 0; j < U.columns(); j++) {
        auto max_index = max(U, j, j);

        if (j != max_index) {
            U.swap_rows(j, max_index);
            P.swap_rows(j, max_index);
        }

        for (index_type i = j + 1; i < U.rows(); i++) {
            L(i, j) = U(i, j) / U(j, j);

            for (index_type k = j; k < U.columns(); k++) {
                U(i, k) -= L(i, j) * U(j, k);
            }
        }
    }

    return std::make_tuple(L, U, P);
}
}  // namespace chislaki

#endif  // CHISLAKI_LINALG_DECOMPOSITION_HPP_
