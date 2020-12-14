#include <iostream>
#include <iterator>
#include "Matrix.cpp"
typedef int T;
int main() {
    T matrix[]  = {1, 2, 3,
                   4, 6, 7,
                   7, 1, 5,
                   8, 9, 10};
    T matrix_2[]  = {9, 4, 3,
                     4, 6, 7,
                     1, 1, 5,
                     8, 9, 10};
    // initializer test
    Matrix<T> test_matrix(4, 3, matrix);
    Matrix<T> test_matrix_2(4, 3, matrix_2);
    Matrix<T> coppied_matrix = test_matrix;
    Matrix<T> default_matrix;

    //copy, *=, +=, -= tests
    coppied_matrix = test_matrix;
    coppied_matrix *= coppied_matrix;
    coppied_matrix += coppied_matrix;
    coppied_matrix -= test_matrix;

    // access tests
    vector<T> test_row = coppied_matrix.row(1);
    vector<T> test_column = coppied_matrix.column(1);
    std::copy(test_row.begin(), test_row.end(),
              std::ostream_iterator<T>(std::cout, ", "));
    std::copy(test_column.begin(), test_column.end(),
              std::ostream_iterator<T>(std::cout, ", "));
    std::cout << coppied_matrix(1, 2) << std::endl;
    std::cout << "Matrix: " << std::endl;
    std::cout << coppied_matrix << std::endl;
    std::cout << "Sub-matrix: " << std::endl;
    std::cout << coppied_matrix.m_get_submatrix(1, 4, 0, coppied_matrix.m_cols_size() ) << std::endl;

    // matrix operators
    std::cout << "Matrix test: " << std::endl;
    std::cout << test_matrix << std::endl;
    std::cout << "Appending matrix test: " << std::endl;
    test_matrix.append_columns(2);
    test_matrix.append_rows(3, 1);
    std::cout << test_matrix << "\n" << std::endl;

    // -, *, +, << tests
    std::cout << coppied_matrix + coppied_matrix - coppied_matrix * coppied_matrix<< std::endl;

    // transpose test

    std::cout << (test_matrix_2) << std::endl;
    std::cout << (test_matrix_2).transpose()<< std::endl;

    // vector func tests
    vector<T> diagonal= coppied_matrix.diag_vec();
    diagonal = coppied_matrix * diagonal;
    std::copy(diagonal.begin(), diagonal.end(),
              std::ostream_iterator<T>(std::cout, ", "));

    // scalar tests
    std::cout << "\n" << test_matrix + 4 - 10 / 4 * 4 << std::endl;

    // if dimensional is right
    std::cout << (int)(test_matrix.m_size() == test_matrix.m_rows_size() * test_matrix.m_cols_size()) << endl;




    return 0;
}