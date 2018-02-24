#ifndef CHISLAKI_LINALG_MATRIX_HPP_
#define CHISLAKI_LINALG_MATRIX_HPP_

#include <iomanip>
#include <iostream>
#include <vector>

namespace chislaki {

// exceptions
class bad_index {};
class bad_size {};
class degenerate {};

// global types
using size_type = std::size_t;
using index_type = std::size_t;

template <class T>
class matrix {
public:
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

    matrix(size_type n) : matrix(n, n) {}

    matrix(size_type rows, size_type columns)
        : rows_{rows}, columns_{columns}, data_(rows * columns) {
        if (rows == 0 || columns == 0) {
            throw bad_size{};
        }
    }

    matrix(size_type n, const std::vector<T>& data) : matrix(n, n, data) {}

    matrix(size_type rows, size_type columns, const std::vector<T>& data)
        : rows_{rows}, columns_{columns}, data_{data} {}

    matrix(const matrix& matr) = default;

    matrix& operator=(const matrix& matr) = default;

    matrix(matrix&& matr) : rows_{matr.rows_}, columns_{matr.columns_} {
        data_ = matr.data_;
        matr.data_ = std::vector<T>();
    }

    // ***********************************************************************
    // ***************************** Matrix size *****************************
    // ***********************************************************************

    inline size_type columns() const noexcept { return columns_; }

    inline size_type rows() const noexcept { return rows_; }

    // ***********************************************************************
    // ************************ Row and column access ************************
    // ***********************************************************************

    matrix row(index_type index) const noexcept {
        auto result = matrix(1, columns());
        for (index_type j = 0; j < columns(); j++) {
            result(0, j) = (*this)(index, j);
        }
        return result;
    }

    matrix column(index_type index) const noexcept {
        auto result = matrix(rows(), 1);
        for (index_type i = 0; i < rows(); i++) {
            result(i, 0) = (*this)(i, index);
        }
        return result;
    }

    void set_row(index_type index, const matrix& matr) noexcept {
        for (index_type j = 0; j < columns(); j++) {
            (*this)(index, j) = matr(0, j);
        }
    }

    void set_column(index_type index, const matrix& matr) noexcept {
        for (index_type i = 0; i < rows(); i++) {
            (*this)(i, index) = matr(i, 0);
        }
    }

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

    const_reference operator()(index_type index) const noexcept {
        return (rows() == 1) ? (*this)(0, index) : (*this)(index, 0);
    }

    reference operator()(index_type index) noexcept {
        return (rows() == 1) ? (*this)(0, index) : (*this)(index, 0);
    }

    const_reference at(index_type index) const {
        if (rows() != 1 && columns() != 1) {
            throw bad_size{};
        }
        return (rows() == 1) ? at(0, index) : at(index, 0);
    }

    reference at(index_type index) {
        if (rows() != 1 && columns() != 1) {
            throw bad_size{};
        }
        return (rows() == 1) ? at(0, index) : at(index, 0);
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

    matrix operator*(const matrix& right) const {
        const auto& left = *this;

        if (left.columns() != right.rows()) {
            throw bad_size{};
        }

        auto result = matrix(left.rows(), right.columns());

        for (index_type i = 0; i < left.rows(); i++) {
            for (index_type j = 0; j < right.columns(); j++) {
                for (index_type k = 0; k < left.columns(); k++) {
                    result(i, j) += left(i, k) * right(k, j);
                }
            }
        }

        return result;
    }

    matrix operator-() const { return -1 * (*this); }

    // ***********************************************************************
    // *********************** Static member functions ***********************
    // ***********************************************************************

    static matrix make_eye(size_type size) {
        auto result = matrix(size);
        for (index_type i = 0; i < result.rows(); i++) {
            result(i, i) = 1;
        }
        return result;
    }

    static inline matrix make_row(size_type size) { return matrix(1, size); }

    static inline matrix make_row(size_type size, const std::vector<T>& data) {
        return matrix(1, size, data);
    }

    static inline matrix make_column(size_type size) { return matrix(size, 1); }

    static inline matrix make_column(size_type size,
                                     const std::vector<T>& data) {
        return matrix(size, 1, data);
    }

private:
    // **********************************************************************
    // ********************** Private member functions **********************
    // **********************************************************************

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

    size_type rows_;
    size_type columns_;
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
            os << std::setw(15) << std::setprecision(8) << matr(i, j) << " ";
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

// ***********************************************************************
// *********************** Special matrix creation ***********************
// ***********************************************************************

template <class T>
matrix<T> make_eye(size_type size) {
    return matrix<T>::make_eye(size);
}

template <class T>
matrix<T> make_row(size_type size) {
    return matrix<T>::make_row(size);
}

template <class T>
matrix<T> make_row(size_type size, const std::vector<T>& data) {
    return matrix<T>::make_row(size, data);
}

template <class T>
matrix<T> make_row(index_type index, const matrix<T>& matr) {
    return matr.row(index);
}

template <class T>
matrix<T> make_column(size_type size) {
    return matrix<T>::make_column(size);
}

template <class T>
matrix<T> make_column(size_type size, const std::vector<T>& data) {
    return matrix<T>::make_column(size, data);
}

template <class T>
matrix<T> make_column(index_type index, const matrix<T>& matr) {
    return matr.column(index);
}

// **********************************************************************
// ************************ Matrix transposition ************************
// **********************************************************************

template <class T>
matrix<T> transpose(const matrix<T>& matr) noexcept {
    auto result = matrix<T>(matr.columns(), matr.rows());
    for (index_type j = 0; j < matr.columns(); j++) {
        for (index_type i = 0; i < matr.rows(); i++) {
            result(j, i) = matr(i, j);
        }
    }
    return result;
}

}  // namespace chislaki

#endif  // CHISLAKI_LINALG_MATRIX_HPP_
