#include <iostream>
#include <iterator>
// #include "Matrix.cpp"
// Чомко Василь К-21
/*
 * 1) Клас (загальний опис)
        Клас являє собою абстрактний тип даних, що визначається користувачем і являє собою модель реального об’єкта у вигляді даних
        та функцій для роботи з ними.

        Схему об'єкта і його реалізацію розділені, що дозволяє описати, що може робити об'єкт, не
        конкретизуючи, як це робити. У цьому полягає принцип абстракції даних. Принципова схема об'єкта не
        повинна залежати від її конкретного наповнення. Необхідно лише вказати, що об'єкт повинний містити певну
        операцію. Якщо алгоритм знадобиться змінити, ми модифікуємо його реалізацію, але загальна структура об'єкта
        від цього не постраждає.
        Успадкування - створення об'єкта, що успадковує усі властивості свого попередника із додаванням в нього нових можливостей.
        Об'єкт може одержувати різноманітну інформацію. Зрозуміло, для її обробки можна заздалегідь передбачити
        його реакцію, написавши відповідні функції. Однак в об’єктно-орієнтованих мовах існує можливість
        перевантажувати операції і функції, що дозволяє об'єкту самому конкретизувати їхній зміст у ході виконання
        програми.
        Об'єкт — це фізична сутність, що виникає при виконанні програми, тобто сукупність комірок пам'яті, що зберігають дані і
        код. Для того щоб створити об'єкт, необхідна його схема: які дані він містить, які функції обробляють ці дані, як
        організований доступ до цих даних і функцій, що називаються членами об'єкта. У мові С++ схемою об'єкта
        називається клас.

        Для оголошення класу в мові С++ призначене ключове слово class. Усе, що розташовано між фігурними
        дужками, що слідують за цим ключовим словом, являє собою тіло класу, що складається з його членів. Варто
        пам'ятати, що клас — це логічна схема об'єкта. Отже, виділення пам'яті відбудеться лише тоді, коли об'єкт буде
        визначений.

        Загальний вид оголошення класу:
        class ім'я{
         члени класу;
         private:
         члени класу;
         protected:
         члени класу;
         public:
         члени класу;
         } <список об'єктів>;

        Визначати об'єкти відразу після оголошення класів не обов’язково — це можна зробити в придатному місці
        програми. Як показане вище, тіло класу розділяється на три розділи, позначені наступними ключовими словами
        (специфікаторами доступу).

        public // Відкритий розділ
        private // Закритий розділ
        protected // Захищений розділ

        Ключове слово public позначає розділ, у якому розміщаються функції і дані, доступні з будь-якої точки
        програми, — відкриті члени. Ключове слово private оголошує розділ, у якому розташовані закриті члени
        класу, що є доступними лише функціям-членам самого класу, а також дружнім функціям і класам. За замовчуванням
        функції і дані, оголошені в private розділі.

 * 2) Лямбда-функції
        Лямбда-функція — це опис функціональної можливості, яку можна визначити в операторі і виразі.
        Таким чином, лямбда-функцію можна використовуватияк функцію, що підставляється.
        Лямбда-функції забезпечують більш інтуїтивний підхід до визначення функціонального
        поводження алгоритмів з бібліотеки STL. Крім того, лямбда-функції працюють швидше функціональних об'єктів.
        Лямбда-функції завжди передує так називаний ініціатор лямбда-функції (lambda introducer):
        квадратні дужки, усередині яких можна визначити так називане захоплення для доступу до нестатичних зовнішніх об'єктів у лямбда-функції.

        Якщо доступ до зовнішніх даних не потрібний, квадратні дужки залишаються порожніми.
        У лямбда-функціях можна використовувати статичні об'єкти.
        Між ініціатором лямбда-функції і тілом лямбда-функції можна вказати параметри, ключове слово mutable,
        специфікацію винятку, специфікатори атрибутів і тип значення, що повертається.
        Усе це не обов’язково, але якщо один з перерахованих пунктів присутній,
        то круглі дужки для параметрів стають обов'язковими.
        Лямбда-функції не можуть бути шаблонними. У них завжди необхідно указувати всі типи.
        Лямбда-функція може повертати який-небудь об'єкт.
        Якщо тип об'єкта, що повертається, не зазначений, він виводиться з його значення.
        В ініціаторі лямбда-функції (квадратних дужках перед лямбда-функцією) можна задати список захоплення (capture)
        для доступу до даних із зовнішньої області видимості, що не передаються як аргументи.
        Символи [=] означають, що зовнішня область видимості передається в лямбдафункцію за значенням.
        Символи [&] означають, що зовнішня область видимості передається в лямбдафункцію по посиланню.
        Для кожного об'єкта в лямбда-функції необхідно вказати режим доступу до нього: за значенням чи за посиланням.
        Для того щоб змішати передачу за посиланням і за значенням, лямбда-функцію можна оголосити з ключовим словом mutable.
        У цьому випадку об'єкти передаються за значенням, але у функції-об'єкті, визначеній в лямбда-функції,
        передане значення можна змінити.
        Типом лямбда-функції є анонімна функція-об'єкт (чи функтор), що є унікальною для кожного лямбда-
        виразу.
        Для оголошення об'єктів такого типу необхідні шаблони чи ключове слово auto.
        Якщо необхідний тип, можна використовувати ключове слово decltype, що потрібно при передачі лямбда-функції
        як хеш-функції чи критерію сортування для асоціативних чи неупорядкованих контейнерів.
 */
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
    Matrix<T>  cholesky_decomposition()   const ;

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
    std::vector<T> solve_by_cholesky_decomposition(const Matrix<T>&, std::vector<T>&);

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

#include "Matrix.h"
#include <algorithm>
#include <cmath>

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

// matrix comperisson
template<typename T>
bool Matrix<T>::operator==(const Matrix<T>& other) const noexcept {
    for (int i=0; i < this->m_rows_size(); i++) {
        for (int j=0; j < this->m_rows_size(); j++) {
            if (this->f_matrix_table[i][j] != other(i, j)) { return false;}
        }
    }
    return true;
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

template<typename T>
bool Matrix<T>::is_symmetric() {
    Matrix<T> transposed = this->transpose();
    if (this->operator==(transposed)) {
        return true;
    }
    return false;
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

template<typename T>
Matrix<T> Matrix<T>::cholesky_decomposition() const {
    if (!this->is_symmetric()) { throw MatrixPropertiesException(); }

    T empty_list[this->f_num_rows*this->f_num_columns] = {  };
    Matrix<T> L(this->f_num_rows, this->f_num_columns, empty_list);

    // get diagonal elements
    for (int i=0;i<this->m_rows_size() && i<this->m_cols_size(); i++) {
        for (int j=0;j<=i; j++) {
            int sum = 0;

            if (i==j) {
                for (int k=0; k<i;k++) { sum += L(i, k)*L(i, k); }
                L(i, i) = sqrt(this->f_matrix_table[i][i] - sum);
            }
            else {
                for (int k=0; k<j;k++) { sum += L(i, k)*L(j, k); }
                L(i, j) = (this->f_matrix_table[i][j] - sum)/L(j, j);
            }
        }

    }
    return L;

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

template<typename T>
std::vector<T> Matrix<T>::solve_by_cholesky_decomposition(const Matrix<T>& coeficient_matrix, std::vector<T>& b_vectors) {
    Matrix<T> L = coeficient_matrix.cholesky_decomposition();
    Matrix<T> L_T = coeficient_matrix.cholesky_decomposition().transpose();
    vector<T> y(b_vectors.size(), 0);
    for (int i=0; i < y.size(); ++i) {
        T sum = 0;
        for (int j=0; j<i; i++) {
            sum += y[j]*coeficient_matrix[i][j];
        }
        y[i] = (b_vectors[i] - sum)/coeficient_matrix[i][i];
    }
    // not finished, ... optional

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


#endif //MATRIX_MATRIX_H
int main() {
    // Cholesky decomposition test
    typedef int T;
    {
        T matrix[] = {  4,  12, -16,
                        12,  37, -43,
                        -16, -43,  98};
        Matrix<T> test_matrix(3, 3, matrix);
        std::cout << "Cholesky decomposition test: " << std::endl;
        std::cout << test_matrix << std::endl;
        std::cout << test_matrix.cholesky_decomposition() * test_matrix.cholesky_decomposition().transpose() << std::endl;

    }

    return 0;
}
