#ifndef CHISLAKI_MATRIX_HPP_
#define CHISLAKI_MATRIX_HPP_

#include <cstring>
#include <memory>
#include <numeric>

namespace chislaki {

class bad_index {};
class bad_size {};

template <class T>
class matrix {
public:
    /*
     * Constructors
     */

    matrix(std::size_t rows, std::size_t columns)
        : rows_{rows}, columns_{columns} {
        if (rows == 0 || columns == 0) {
            throw bad_size{};
        }
        data_ = std::shared_ptr<T>{new T[rows * columns],
                                   std::default_delete<T[]>()};
    }

    matrix(std::size_t rows, std::size_t columns, const T* data)
        : matrix(rows, columns) {
        std::memcpy(data_.get(), data, rows * columns);
    }

    matrix(const matrix& matr) = default;

    matrix(matrix&& matr) : rows_{matr.rows_}, columns_{matr.columns_} {
        data_ = nullptr;
        std::swap(data_, matr.data_);
    }

    /*
     * Size
     */

    std::size_t colums() const noexcept { return columns_; }

    std::size_t rows() const noexcept { return rows_; }

    /*
     * Indexing
     */

    const T& operator()(std::size_t row, std::size_t column) const noexcept {
        return data_[row * rows() + column];
    }

    const T& at(std::size_t row, std::size_t column) const {
        if (!(row < rows() && column < colums())) {
            throw bad_index{};
        }
        return data_[row * rows() + column];
    }

    T& operator()(std::size_t row, std::size_t column) noexcept {
        return data_[row * rows() + column];
    }

    T& at(std::size_t row, std::size_t column) {
        if (!(row << rows() && column < rows())) {
            throw bad_index{};
        }
        return data_[row * rows() + column];
    }

    /*
     * Algorithms
     */

    T row_max_value(std::size_t row) const { return row_max(row).first; }

    T row_max_index(std::size_t row) const { return row_max(row).second; }

    std::pair<std::size_t, T> row_max(std::size_t row) const {
        return column_row_max_min(
            std::numeric_limits<T>::min(), row, colums(),
            [this](std::size_t index, std::size_t current_index) {
                return at(index, current_index);
            },
            [this](T left, T right) { return left > right; });
    }

    T col_max_value(std::size_t column) const { return row_max(column).first; }

    T col_max_index(std::size_t column) const {
        return column_max(column).second;
    }

    std::pair<std::size_t, T> column_max(std::size_t column) const {
        return column_row_max_min(
            std::numeric_limits<T>::min(), column, rows(),
            [this](std::size_t index, std::size_t current_index) {
                return at(current_index, index);
            },
            [this](T left, T right) { return left > right; });
    }

    T row_min_value(std::size_t row) const { return row_min(row).first; }

    T row_min_index(std::size_t row) const { return row_min(row).second; }

    std::pair<std::size_t, T> row_min(std::size_t row) const {
        return column_row_max_min(
            std::numeric_limits<T>::max(), row, colums(),
            [this](std::size_t index, std::size_t current_index) {
                return at(index, current_index);
            },
            [this](T left, T right) { return left < right; });
    }

    T col_min_value(std::size_t column) const { return row_min(column).first; }

    T col_min_index(std::size_t column) const {
        return column_min(column).second;
    }

    std::pair<std::size_t, T> column_min(std::size_t column) const {
        return column_row_max_min(
            std::numeric_limits<T>::max(), column, rows(),
            [this](std::size_t index, std::size_t current_index) {
                return at(current_index, index);
            },
            [this](T left, T right) { return left < right; });
    }

private:
    template <class F1, class F2>
    std::pair<std::size_t, T> column_row_max_min(T start_value,
                                                 std::size_t index,
                                                 std::size_t max_index,
                                                 F1&& get,
                                                 F2&& compare) {
        T current_value = start_value;
        std::size_t current_index = 0;
        for (std::size_t i = 0; i < max_index; i++) {
            if (compare(get(index, i), current_value)) {
                current_value = get(i);
                current_index = i;
            }
        }
        return std::make_pair(current_index, current_value);
    }

    std::size_t rows_;
    std::size_t columns_;
    std::shared_ptr<T> data_;
};  // class matrix
}  // namespace chislaki

#endif  // CHISLAKI_MATRIX_HPP_
