#include "Matrix.h"
#include <algorithm>

using namespace std;

// constructors
template <typename T>
Matrix<T>::Matrix(unsigned _rows_num, unsigned _cols_num, const T _initial[]){
    static_assert(std::is_same<T, double>::value ||
                  std::is_same<T, float>::value ||
                  std::is_same<T, int>::value ||
                  std::is_same<T, char>::value, "Martix elements must be numbers.");

    this->f_matrix_table.resize(_rows_num);
    // moving elements to new matrix
    int k = 0;
    for (size_t i=0; i<_rows_num; i++) {
        this->f_matrix_table[i].resize(_cols_num);
        for (size_t j=0; j<_cols_num; j++) {
            this->f_matrix_table[i][j] = _initial[k];
            k++;
        }
    }
    f_num_rows = _rows_num;
    f_num_columns = _cols_num;
}

template <typename T>
Matrix<T>::Matrix(): f_num_rows(0), f_num_columns(0){};

// copy constructor
template <typename T>
Matrix<T>::Matrix(const Matrix<T>& other) {
    this->f_matrix_table = other.f_matrix_table;
    f_num_rows = other.f_num_rows;
    f_num_columns = other.f_num_columns;
}

// destructor
template <typename T>
Matrix<T>::~Matrix() = default;

// get dimensional
template <typename T>
unsigned Matrix<T>::m_rows_size() const noexcept {
    return this->f_num_rows;
}

template <typename T>
unsigned Matrix<T>::m_cols_size() const noexcept  {
    return this->f_num_columns;
}

template <typename T>
unsigned Matrix<T>::m_size() const noexcept  {
    return this->f_num_columns * this->f_num_rows;
}

//excess to elements
template<typename T>
T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) noexcept {
    if (f_num_columns <= col || f_num_rows <= row) {
        throw out_of_range("Out of range exception.\nIndexing starting at (0,0)."); }
    return this->f_matrix_table[row][col];
}

template<typename T>
vector<T> Matrix<T>::column(const unsigned& n) const noexcept {
    if (this->f_num_columns <= n) {
        throw out_of_range("Out of range exception.\nIndexing starting at (0,0)."); }

    vector<T> column_vect(this->f_num_columns, 0.0);
    for (unsigned i = 0; i < this->f_num_columns; i++) {
        column_vect[i] = this->f_matrix_table[i][n];
    }
    return column_vect;
}

template<typename T>
vector<T> Matrix<T>::row(const unsigned& n) const noexcept {
    if (this->f_num_rows <= n) {
        throw out_of_range("Out of range exception.\nIndexing starting at (0,0)."); }

    vector<T> row_vect(this->f_num_rows, 0.0);
    for (unsigned i = 0; i < this->f_num_rows; i++) {
        row_vect[i] = this->f_matrix_table[n][i];
    }
    return row_vect;
}

template<typename T>
Matrix<T> Matrix<T>::m_get_submatrix(size_t from_row, size_t to_row, size_t from_column, size_t to_colunm) const noexcept {
    if (f_num_columns < to_colunm || f_num_rows < to_row ||
        to_colunm < 0 || to_row < 0 || from_column < 0 || from_row < 0 ||
        from_row > to_row || from_column > to_colunm) {
        throw out_of_range("Out of range exception.\nWrong sub-matrix sizes."); }

    T empty_list[(to_colunm-from_column)*(to_row-from_row)] = {  };
    Matrix<T> result(to_row-from_row, to_colunm-from_column, empty_list);

    for (unsigned i=from_row; i<to_row; i++) {
        for (unsigned j=from_column; j<to_colunm; j++) {
            result(i - from_row,j - from_column) = this->f_matrix_table[i][j];
        }
    }

    return result;
}

// Access the individual elements (const)
template<typename T>
const T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) const noexcept {
    if (f_num_columns <= col || f_num_rows <= row) {
        throw out_of_range("Out of range exception.\nIndexing starting at (0,0)."); }
    return this->f_matrix_table[row][col];
}

// matrix operators
template <typename T>
void Matrix<T>::m_append_columns(size_t n, const T initial_value) {
    for (size_t i=0; i<this->f_num_rows; i++) {
        this->f_matrix_table[i].insert(this->f_matrix_table[i].end(), n, initial_value);
    }
    this->f_num_columns += n;
}

template <typename T>
void Matrix<T>::m_append_rows(size_t n, const T initial_value) {
    this->f_matrix_table.resize(this->m_rows_size() + n);
    this->f_num_rows += n;
    for (size_t i=0; i<n; i++) {
        this->f_matrix_table[this->m_rows_size()-n+i]=vector<T>(this->m_cols_size(), initial_value);
    }
}

template <typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& other) noexcept {
    if ( &other == this) {return *this;}

    this->f_matrix_table.assign(other.f_matrix_table.begin(), other.f_matrix_table.end());
    f_num_rows = other.f_num_rows;
    f_num_columns = other.f_num_columns;

    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }
    T empty_list[other.f_num_rows*other.f_num_columns] = {  };
    Matrix<T> result(other.f_num_rows, other.f_num_columns, empty_list);

    for (unsigned i=0; i<other.f_num_rows; i++) {
      for (unsigned j=0; j<other.f_num_columns; j++) {
          result(i,j) = this->f_matrix_table[i][j] + other(i,j);
      }
    }

  return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }
    T empty_list[other.f_num_rows*other.f_num_columns] = {  };
    Matrix<T> result(other.f_num_rows, other.f_num_columns, empty_list);

    for (unsigned i=0; i<other.f_num_rows; i++) {
        for (unsigned j=0; j<other.f_num_columns; j++) {
            result(i,j) = this->f_matrix_table[i][j] - other(i,j);
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }
    T empty_list[other.f_num_rows*other.f_num_columns] = {  };
    Matrix<T> result(other.f_num_rows, other.f_num_columns, empty_list);

    for (unsigned i=0; i<other.f_num_rows; i++) {
        for (unsigned j=0; j<other.f_num_columns; j++) {
            T sum = 0;
            for (unsigned k = 0; k < other.f_num_rows; k++) {
                sum += this->f_matrix_table[i][k] * other(k, j);
            }
            result(i, j) += sum;
        }
    }

    return result;
}

template<typename T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& other) {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }
    for (unsigned i=0; i<other.f_num_rows; i++) {
        for (unsigned j=0; j<other.f_num_columns; j++) {
            this->f_matrix_table[i][j] += other(i,j);
        }
    }

    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& other) {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }
    for (unsigned i=0; i<other.f_num_rows; i++) {
        for (unsigned j=0; j<other.f_num_columns; j++) {
            this->f_matrix_table[i][j] -= other(i,j);
        }
    }

    return *this;
}

template<typename T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& other) {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }
    Matrix result = (*this) * other;
    (*this) = result;
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::transpose() {
    T empty_list[this->f_num_rows*this->f_num_columns] = {  };
    Matrix<T> result(this->f_num_columns, this->f_num_rows, empty_list);

    for (unsigned i=0; i<this->f_num_rows; i++) {
        for (unsigned j=0; j<this->f_num_columns; j++) {
            result(j,i) = this->f_matrix_table[i][j];
        }
    }

    std::swap(this->f_num_rows, this->f_num_columns);
    return result;
}

// some algorithms
template<typename T>
Matrix<T>  Matrix<T>::strassen_multiplication(const Matrix<T>& other) const {
    if (f_num_columns != other.f_num_columns || f_num_rows != other.f_num_rows) {
        throw MatrixDimensionalException();
    }

    //base case
    if (this->f_num_rows <= 2) { return (*this)*other; }

    Matrix<T> first_matrix = *(this);
    Matrix<T> second_matrix = other;

    // append matrixs to [2^n * 2^n]
    {
        int size_diff = first_matrix.f_num_columns - first_matrix.f_num_rows;
        if (size_diff < 0) {
            first_matrix.m_append_columns(abs(size_diff), 0);
            second_matrix.m_append_columns(abs(size_diff), 0);
        }else if (size_diff > 0) {
            first_matrix.m_append_rows(size_diff, 0);
            second_matrix.m_append_rows(size_diff, 0);
        }

        size_diff = 0;
        while ((first_matrix.m_cols_size()+size_diff)&(first_matrix.m_cols_size()-1+size_diff) != 0) { size_diff++; }
        first_matrix.m_append_columns(size_diff, 0);
        second_matrix.m_append_columns(size_diff, 0);
        first_matrix.m_append_rows(size_diff, 0);
        second_matrix.m_append_rows(size_diff, 0);
    }

    Matrix<T> a_1_1 = first_matrix.m_get_submatrix(0,                            first_matrix.m_rows_size()/2,
                                                   0,                            first_matrix.m_cols_size()/2);
    Matrix<T> a_2_1 = first_matrix.m_get_submatrix(first_matrix.m_rows_size()/2, first_matrix.m_rows_size(),
                                                   0,                            first_matrix.m_cols_size()/2);
    Matrix<T> a_1_2 = first_matrix.m_get_submatrix(0,                            first_matrix.m_rows_size()/2,
                                                   first_matrix.m_cols_size()/2, first_matrix.m_cols_size());
    Matrix<T> a_2_2 = first_matrix.m_get_submatrix(first_matrix.m_rows_size()/2, first_matrix.m_rows_size(),
                                                   first_matrix.m_cols_size()/2, first_matrix.m_cols_size());

    Matrix<T> b_1_1 = second_matrix.m_get_submatrix(0,                             second_matrix.m_rows_size()/2,
                                                    0,                             second_matrix.m_cols_size()/2);
    Matrix<T> b_2_1 = second_matrix.m_get_submatrix(second_matrix.m_rows_size()/2, second_matrix.m_rows_size(),
                                                    0,                             second_matrix.m_cols_size()/2);
    Matrix<T> b_1_2 = second_matrix.m_get_submatrix(0,                             second_matrix.m_rows_size()/2,
                                                    second_matrix.m_cols_size()/2, second_matrix.m_cols_size());
    Matrix<T> b_2_2 = second_matrix.m_get_submatrix(second_matrix.m_rows_size()/2, second_matrix.m_rows_size(),
                                                    second_matrix.m_cols_size()/2, second_matrix.m_cols_size());

    Matrix<T> m_1 = (a_1_1 + a_2_2).strassen_multiplication(b_1_1 + b_2_2);
    Matrix<T> m_2 = (a_2_1 + a_2_2).strassen_multiplication(b_1_1);
    Matrix<T> m_3 = a_1_1.strassen_multiplication(b_1_2 - b_2_2);
    Matrix<T> m_4 = a_2_2.strassen_multiplication(b_2_1 - b_1_1);
    Matrix<T> m_5 = (a_1_1 + a_1_2).strassen_multiplication(b_2_2);
    Matrix<T> m_6 = (a_2_1 - a_1_1).strassen_multiplication(b_1_1 + b_1_2);
    Matrix<T> m_7 = (a_1_2 - a_2_2).strassen_multiplication(b_2_1 + b_2_2);

    Matrix<T> c_1_1 = m_1 + m_4 - m_5 + m_7;
    Matrix<T> c_1_2 = m_3 + m_5;
    Matrix<T> c_2_1 = m_2 + m_4;
    Matrix<T> c_2_2 = m_1 - m_2 + m_3 + m_6;

    // gather result(C) matrix
    Matrix<T> result = c_1_1;
    result.m_append_columns(result.m_cols_size());
    result.m_append_rows(result.m_rows_size());
        // c_2_2
    for (int i=c_2_2.m_rows_size(); i < result.m_rows_size(); i++) {
        for (int j=c_2_2.m_cols_size(); j < result.m_cols_size(); j++) {
            result.f_matrix_table[i][j] = c_2_2.f_matrix_table[i-c_2_2.m_rows_size()][j-c_2_2.m_cols_size()];
        }
    }
        // c_1_2
    for (int i=0; i < c_1_2.m_rows_size(); i++) {
        for (int j=c_1_2.m_cols_size(); j < result.m_cols_size(); j++) {
            result(i, j) = c_1_2(i, j-c_1_2.m_cols_size());
        }
    }
        // c_2_1
    for (int i=c_2_1.m_rows_size(); i < result.m_rows_size(); i++) {
        for (int j=0; j < c_2_1.m_cols_size(); j++) {
            result(i, j) = c_1_2(i-c_2_1.m_rows_size(), j);
        }
    }


    return result;
}

// vector operators
template<typename T>
std::vector<T> Matrix<T>::operator*(const std::vector<T>& vect) const {
    std::vector<T> result(this->f_num_rows, 0.0);

    for (unsigned i=0; i<this->f_num_rows; i++) {
        for (unsigned j = 0; j < this->f_num_columns; j++) {
            result[i] = this->f_matrix_table[i][j] * vect[j];
        }
    }
    return result;
}

template<typename T>
std::vector<T> Matrix<T>::diag_vec() const noexcept{
    std::vector<T> result(std::min(this->f_num_rows, this->f_num_columns), 0.0);

    for (int i=0; i < std::min(this->f_num_rows, this->f_num_columns); i++) {
        result[i] = this->f_matrix_table[i][i];
    }
    return result;
}

// scalar operators
template<typename T>
Matrix<T> Matrix<T>::operator+(const T& scalar) const noexcept {
    T empty_list[this->f_num_rows*this->f_num_columns] = {  };
    Matrix<T> result(this->f_num_rows, this->f_num_columns, empty_list);

    for (unsigned i=0; i<this->f_num_rows; i++) {
        for (unsigned j=0; j<this->f_num_columns; j++) {
            result(i,j) = this->f_matrix_table[i][j] + scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const T& scalar) const noexcept {
    T empty_list[this->f_num_rows*this->f_num_columns] = {  };
    Matrix<T> result(this->f_num_rows, this->f_num_columns, empty_list);

    for (unsigned i=0; i<this->f_num_rows; i++) {
        for (unsigned j=0; j<this->f_num_columns; j++) {
            result(i,j) = this->f_matrix_table[i][j] - scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) const noexcept {
    T empty_list[this->f_num_rows*this->f_num_columns] = {  };
    Matrix<T> result(this->f_num_rows, this->f_num_columns, empty_list);

    for (unsigned i=0; i<this->f_num_rows; i++) {
        for (unsigned j=0; j<this->f_num_columns; j++) {
            result(i,j) = this->f_matrix_table[i][j] * scalar;
        }
    }

    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator/(const T& scalar) const noexcept {
    T empty_list[this->f_num_rows*this->f_num_columns] = {  };
    Matrix<T> result(this->f_num_rows, this->f_num_columns, empty_list);

    for (unsigned i=0; i<this->f_num_rows; i++) {
        for (unsigned j=0; j<this->f_num_columns; j++) {
            result(i,j) = this->f_matrix_table[i][j] / scalar;
        }
    }

    return result;
}

// Overloading vector printing
template<typename T>
inline std::ostream& operator<< (std::ostream& out, const std::vector<T>& vector){
    out << "( ";
    std::copy(vector.begin(), vector.end(),
              std::ostream_iterator<T>(out, ", "));
    out << ")" << std::endl;
    return out;
}
