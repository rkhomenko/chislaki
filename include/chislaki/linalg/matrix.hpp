#ifndef CHISLAKI_LINALG_MATRIX_HPP_
#define CHISLAKI_LINALG_MATRIX_HPP_

#include <iostream>
#include <limits>
#include <vector>

namespace chislaki {

// exceptions
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

    // **********************************************************************
    // ************************** Friend functions **************************
    // **********************************************************************

    template <class Type>
    friend matrix<Type> operator*(const matrix<Type>& matr,
                                  Type alpha) noexcept;
    template <class Type>
    friend matrix<Type> operator*(Type alpha,
                                  const matrix<Type>& matr) noexcept;

    // **********************************************************************
    // **************************** Constructors ****************************
    // **********************************************************************

    matrix(index_type n) : matrix(n, n) {}

    matrix(index_type rows, index_type columns)
        : rows_{rows}, columns_{columns}, data_(rows * columns) {
        if (rows == 0 || columns == 0) {
            throw bad_size{};
        }
    }

    matrix(index_type n, const std::vector<T>& data) : matrix(n, n, data) {}

    matrix(index_type rows, index_type columns, const std::vector<T>& data)
        : rows_{rows}, columns_{columns}, data_{data} {}

    matrix(const matrix& matr) = default;

    matrix(matrix&& matr) : rows_{matr.rows_}, columns_{matr.columns_} {
        data_ = matr.data_;
        matr.data_ = std::vector<T>();
    }

    // ***********************************************************************
    // ***************************** Matrix size *****************************
    // ***********************************************************************

    inline index_type columns() const noexcept { return columns_; }

    inline index_type rows() const noexcept { return rows_; }

    // **********************************************************************
    // ****************************** Raw data ******************************
    // **********************************************************************

    inline const_iterator data() const noexcept { return data_.data(); }

    inline iterator data() noexcept { return data_.data(); }

    // ***********************************************************************
    // ****************************** Iterators ******************************
    // ***********************************************************************

    inline const_iterator begin() const noexcept { return data(); }

    inline iterator begin() noexcept { return data(); }

    inline const_iterator end() const noexcept {
        return &data()[rows() * columns() - 1];
    }

    inline iterator end() noexcept { return &data()[rows() * columns() - 1]; }

    // **********************************************************************
    // ****************************** Indexing ******************************
    // **********************************************************************

    inline const_reference operator()(index_type row, index_type column) const
        noexcept {
        return data_[row * columns() + column];
    }

    const_reference at(index_type row, index_type column) const {
        if (!(row < rows() && column < columns())) {
            throw bad_index{};
        }
        return data_[row * columns() + column];
    }

    inline reference operator()(index_type row, index_type column) noexcept {
        return data_[row * columns() + column];
    }

    reference at(index_type row, index_type column) {
        if (!(row < rows() && column < columns())) {
            throw bad_index{};
        }
        return data_[row * columns() + column];
    }

    // **********************************************************************
    // ***************************** Algorithms *****************************
    // **********************************************************************

    // **********************************************************************
    // *********************** Row and column maximum ***********************
    // **********************************************************************

    T row_max_value(index_type row) const { return row_max(row).second; }

    T row_max_index(index_type row) const { return row_max(row).first; }

    std::pair<index_type, T> row_max(index_type row) const {
        return row_max_min(std::numeric_limits<T>::min(), row,
                           [](T left, T right) { return left > right; });
    }

    T row_max_abs_value(index_type row) const {
        return row_max_abs(row).second;
    }

    T row_max_abs_index(index_type row) const { return row_max_abs(row).first; }

    std::pair<index_type, T> row_max_abs(index_type row) const {
        return row_max_min(0, row, [](T left, T right) {
            return std::abs(left) > std::abs(right);
        });
    }

    T col_max_value(index_type column) const {
        return column_max(column).second;
    }

    T col_max_index(index_type column) const {
        return column_max(column).first;
    }

    std::pair<index_type, T> column_max(index_type index) const {
        return column_max_min(std::numeric_limits<T>::min(), index,
                              [](T left, T right) { return left > right; });
    }

    T col_max_abs_value(index_type column) const {
        return column_max_abs(column).second;
    }

    T col_max_abs_index(index_type column) const {
        return column_max_abs(column).first;
    }

    std::pair<index_type, T> column_max_abs(index_type index) const {
        return column_max_min(0, index, [](T left, T right) {
            return std::abs(left) > std::abs(right);
        });
    }

    // **********************************************************************
    // *********************** Row and column minimum ***********************
    // **********************************************************************

    T row_min_value(index_type row) const { return row_min(row).second; }

    T row_min_index(index_type row) const { return row_min(row).first; }

    std::pair<index_type, T> row_min(index_type row) const {
        return row_max_min(std::numeric_limits<T>::max(), row,
                           [](T left, T right) { return left < right; });
    }

    T row_min_abs_value(index_type row) const {
        return row_min_abs(row).second;
    }

    T row_abs_index(index_type row) const { return row_min_abs(row).first; }

    std::pair<index_type, T> row_min_abs(index_type row) const {
        return row_max_min(
            std::numeric_limits<T>::max(), row,
            [](T left, T right) { return std::abs(left) < (right); });
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

    T col_min_abs_value(index_type column) const {
        return column_min_abs(column).second;
    }

    T col_min_abs_index(index_type column) const {
        return column_min_abs(column).first;
    }

    std::pair<index_type, T> column_min_abs(index_type column) const {
        return column_max_min(
            std::numeric_limits<T>::max(), column,
            [](T left, T right) { return std::abs(left) < std::abs(right); });
    }

    // ***********************************************************************
    // ************************ Rows and columns swap ************************
    // ***********************************************************************

    void swap_rows(index_type row1_index, index_type row2_index) {
        std::vector<T> row1(columns());

        std::copy(data_.begin() + row1_index * columns(),
                  data_.begin() + row1_index * columns() + columns(),
                  row1.begin());

        std::copy(data_.begin() + row2_index * columns(),
                  data_.begin() + row2_index * columns() + columns(),
                  data_.begin() + row1_index * columns());

        std::copy(row1.begin(), row1.end(),
                  data_.begin() + row2_index * columns());
    }

    void swap_columns(index_type col1_index, index_type col2_index) {
        for (index_type i = 0; i < rows(); i++) {
            std::swap(at(i, col1_index), at(i, col2_index));
        }
    }

    // **********************************************************************
    // ************************ Arithmetic operators ************************
    // **********************************************************************

    matrix operator+(const matrix& matr) const {
        return plus_minus_operator(matr,
                                   [](T lhs, T rhs) { return lhs + rhs; });
    }

    matrix operator-(const matrix& matr) const {
        return plus_minus_operator(matr,
                                   [](T lhs, T rhs) { return lhs - rhs; });
    }

    matrix operator-() const { return -1 * (*this); }

    // ***********************************************************************
    // *********************** Static member functions ***********************
    // ***********************************************************************

    static matrix eye(index_type size) {
        auto result = matrix(size);
        for (index_type i = 0; i < result.rows(); i++) {
            result(i, i) = 1;
        }
        return result;
    }

    static inline matrix row(index_type size) { return matrix(1, size); }

    static inline matrix row(index_type size, const std::vector<T>& data) {
        return matrix(1, size, data);
    }

    static inline matrix column(index_type size) { return matrix(size, 1); }

    static inline matrix column(index_type size, const std::vector<T>& data) {
        return matrix(size, 1, data);
    }

private:
    template <class P>
    std::pair<index_type, T> column_max_min(value_type start_value,
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
    std::pair<index_type, T> row_max_min(value_type start_value,
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

    template <class F>
    matrix plus_minus_operator(const matrix& matr, F&& functor) const {
        matrix result{rows(), columns()};

        if (rows() != matr.rows() && columns() != matr.columns()) {
            throw bad_size{};
        }

        for (index_type i = 0; i < rows() * columns(); i++) {
            result.data_[i] = functor(data_[i], matr.data_[i]);
        }

        return result;
    }

    // **********************************************************************
    // ************************** Member variables **************************
    // **********************************************************************

    index_type rows_;
    index_type columns_;
    std::vector<T> data_;
};  // class matrix

// ***********************************************************************
// *********************** Matrix friend functions ***********************
// ***********************************************************************

template <class T>
matrix<T> operator*(const matrix<T>& matr, T alpha) noexcept {
    matrix<T> result{matr.rows(), matr.columns()};

    using index_type = typename matrix<T>::index_type;
    for (index_type i = 0; i < matr.rows() * matr.columns(); i++) {
        result.data_[i] = alpha * matr.data_[i];
    }

    return result;
}

template <class T>
matrix<T> operator*(T alpha, const matrix<T>& matr) noexcept {
    return matr * alpha;
}

// **********************************************************************
// **************************** IO operators ****************************
// **********************************************************************

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

#endif  // CHISLAKI_LINALG_MATRIX_HPP_
