cmake_minimum_required(VERSION 3.17)
project(Matrix)

set(CMAKE_CXX_STANDARD 14)

add_executable(Matrix Matrix.h test_class_matrix.cpp)