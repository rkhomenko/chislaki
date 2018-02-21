#ifndef CHISLAKI_MATRIX_HPP_
#define CHISLAKI_MATRIX_HPP_

#include <iostream>
#include <numeric>
#include <vector>

namespace chislaki {

class bad_index {};
class bad_size {};

template <class T>
class matrix {
public:
    using index_type = std::size_t;
    using value_type = T;
    using reference = T&;
    using const_reference = const T&;
    using iterator = T*;
    using const_iterator = const T*;

    matrix(index_type rows, index_type columns)
        : rows_{rows}, columns_{columns}, data_(rows * columns) {
        if (rows == 0 || columns == 0) {
            throw bad_size{};
        }
    }

    matrix(index_type rows, index_type columns, const std::vector<T>& data)
        : rows_{rows}, columns_{columns}, data_{data} {}

    matrix(const matrix& matr) = default;

    matrix(matrix&& matr) : rows_{matr.rows_}, columns_{matr.columns_} {
        data_ = matr.data_;
        matr.data_ = std::vector<T>();
    }

    index_type columns() const noexcept { return columns_; }

    index_type rows() const noexcept { return rows_; }

    const_iterator data() const noexcept { return data_.data(); }

    iterator data() noexcept { return data_.data(); }

    const_iterator begin() const noexcept { return &data_[0]; }

    iterator begin() noexcept { return &data_[0]; }

    const_iterator end() const noexcept {
        return &data_[rows() * columns() - 1];
    }

    iterator end() noexcept { return &data_[rows() * columns() - 1]; }

    const_reference operator()(index_type row, index_type column) const
        noexcept {
        return data_[row * columns() + column];
    }

    const_reference at(index_type row, index_type column) const {
        if (!(row < rows() && column < columns())) {
            throw bad_index{};
        }
        return data_[row * columns() + column];
    }

    reference operator()(index_type row, index_type column) noexcept {
        return data_[row * columns() + column];
    }

    reference at(index_type row, index_type column) {
        if (!(row << rows() && column < rows())) {
            throw bad_index{};
        }
        return data_[row * columns() + column];
    }

    T row_max_value(index_type row) const { return row_max(row).second; }

    T row_max_index(index_type row) const { return row_max(row).first; }

    std::pair<index_type, T> row_max(index_type row) const {
        return row_max_min(
            std::numeric_limits<T>::min(), row,
            [this](T left, T right) { return left > right; });
    }

    T col_max_value(index_type column) const {
        return column_max(column).second;
    }

    T col_max_index(index_type column) const {
        return column_max(column).first;
    }

    std::pair<index_type, T> column_max(index_type index) const {
        return column_max_min(std::numeric_limits<T>::min(), index,
                              [this](T left, T right) { return left > right; });
    }

    T row_min_value(index_type row) const { return row_min(row).second; }

    T row_min_index(index_type row) const { return row_min(row).first; }

    std::pair<index_type, T> row_min(index_type row) const {
        return row_max_min(std::numeric_limits<T>::max(), row,
                           [this](T left, T right) { return left < right; });
    }

    T col_min_value(index_type column) const {
        return column_min(column).second;
    }

    T col_min_index(index_type column) const {
        return column_min(column).first;
    }

    std::pair<index_type, T> column_min(index_type column) const {
        return column_max_min(std::numeric_limits<T>::max(), column,
                              [this](T left, T right) { return left < right; });
    }

    void swap_rows(std::size_t row1, std::size_t row2) {
        // TODO: Implement!
    }

private:
    template <class P>
    std::pair<std::size_t, T> column_max_min(value_type start_value,
                                             index_type index,
                                             P&& compare) const {
        T current_value = start_value;
        std::size_t current_index = 0;
        for (std::size_t i = 0; i < rows(); i++) {
            if (compare(at(i, index), current_value)) {
                current_value = at(i, index);
                current_index = i;
            }
        }
        return std::make_pair(current_index, current_value);
    }

    template <class P>
    std::pair<std::size_t, T> row_max_min(value_type start_value,
                                          index_type index,
                                          P&& compare) const {
        T current_value = start_value;
        std::size_t current_index = 0;
        for (std::size_t j = 0; j < columns(); j++) {
            if (compare(at(index, j), current_value)) {
                current_value = at(index, j);
                current_index = j;
            }
        }
        return std::make_pair(current_index, current_value);
    }

    index_type rows_;
    index_type columns_;
    std::vector<T> data_;
};  // class matrix

template <class T>
std::ostream& operator<<(std::ostream& os, const matrix<T>& matr) {
    for (std::size_t i = 0; i < matr.rows(); i++) {
        for (std::size_t j = 0; j < matr.columns(); j++) {
            os << matr(i, j) << " ";
        }
        os << std::endl;
    }

    return os;
}

template <class T>
std::istream& operator>>(std::istream& is, matrix<T>& matr) {
    for (std::size_t i = 0; i < matr.rows(); i++) {
        for (std::size_t j = 0; j < matr.columns(); j++) {
            is >> matr(i, j);
        }
    }

    return is;
}

}  // namespace chislaki

#endif  // CHISLAKI_MATRIX_HPP_