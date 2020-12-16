#ifndef MATRIX_MATRIX_H
#define MATRIX_MATRIX_H

#include <vector>
#include <array>
#include <exception>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <string>

class MatrixDimensionalException : public std::exception {
    const char * what () const noexcept override
    {
        return "Matrix has different dimensions. Operation is impossible.";
    }
};

class MatrixPropertiesException : public std::exception {
    const char * what () const noexcept override
    {
        return "Matrix doesn`t have necessary properties. Operation is impossible.";
    }
};
template <typename T>
class Matrix {
public:
    Matrix(unsigned _rows_num, unsigned _cols_num, const T[]);
    Matrix();
    Matrix(const Matrix<T>&);
    virtual ~Matrix();

    // Comparison
    bool       operator==(const Matrix<T>&) const noexcept;
//    bool       operator>(const Matrix<T>&)  const noexcept;
//    bool       operator<(const Matrix<T>&)  const noexcept;

    // Mathematical operations
    Matrix<T>& operator=(const Matrix<T>&) noexcept;
    Matrix<T>  operator+(const Matrix<T>&)   const ;
    Matrix<T>& operator+=(const Matrix<T>&)        ;
    Matrix<T>  operator-(const Matrix<T>&)   const ;
    Matrix<T>& operator-=(const Matrix<T>&)        ;
    Matrix<T>  operator*(const Matrix<T>&)   const ;
    Matrix<T>& operator*=(const Matrix<T>&)        ;

    // Some matrix algorithms
    Matrix<T>  strassen_multiplication(const Matrix<T>&)   const ;
    Matrix<T>  cholesky_decomposition()                          ;

    // Matrix operations & properties
    Matrix<T>  transpose();
    bool       is_symmetric();
    void       m_append_rows(size_t n, T initial_value=0);
    void       m_append_columns(size_t n, T initial_value=0);

    //Matrix-scalar operations
    Matrix<T>  operator+(const T&)  const noexcept;
    Matrix<T>  operator-(const T&)  const noexcept;
    Matrix<T>  operator*(const T&)  const noexcept;
    Matrix<T>  operator/(const T&)  const noexcept;

    // Matrix-vector operations
    std::vector<T> operator*(const std::vector<T>&)  const;
    std::vector<T> diag_vec()  const noexcept;

    // Access the individual elements
    T& operator()(const unsigned& row, const unsigned& col)             noexcept;
    const T& operator()(const unsigned& row, const unsigned& col) const noexcept;

    // Access the row and column sizes
    unsigned m_rows_size() const noexcept;
    unsigned m_cols_size() const noexcept;
    unsigned m_size()      const noexcept;

    // Access the row and column
    std::vector<T> row(const unsigned& n)       const noexcept;
    std::vector<T> column(const unsigned& n)    const noexcept;

    // Sub-matrix access
    Matrix<T> m_get_submatrix(size_t from_row, size_t to_row, size_t from_column, size_t to_column) const noexcept;

    // Printing matrix
    friend std::ostream& operator<< (std::ostream& out, const Matrix<T>& matrix){
        size_t max_elem_size = 1;
        for (auto i =matrix.f_matrix_table.begin(); i != matrix.f_matrix_table.end(); i++){
            std::for_each((*i).begin(), (*i).end(),
                          [&](T elem){
                              if (std::to_string(elem).size() > max_elem_size) {
                                  max_elem_size = std::to_string(elem).size();
                              }
                          });
        }
        max_elem_size += max_elem_size % 2 == 1 ? 1 : 0;

        for (auto line = matrix.f_matrix_table.begin(); line != matrix.f_matrix_table.end(); line++){
            out << "| ";
            for (auto elem = (*line).begin(); elem != (*line).end(); elem++){
                int element_size = std::to_string(*elem).size();
                out << std::setw((max_elem_size-element_size)/2) << "" << *elem << std::setw((max_elem_size-element_size)/2) << " ";
            }
            out << "|";
//            if (line != matrix.f_matrix_table.end() - 1) { out << "\n"; }
            out << "\n";
        }
        return out;
}
    // vector printing in .cpp

private:
    std::vector<std::vector<T>> f_matrix_table;
    unsigned                    f_num_rows    ;
    unsigned                    f_num_columns ;

};

template<typename T>
std::vector<double> solve_system_by_cholesky_decomposition(Matrix<T>&, std::vector<T>&);


#endif //MATRIX_MATRIX_H
