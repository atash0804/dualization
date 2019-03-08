cmake_minimum_required(VERSION 3.13)
project(exp1)
project(exp2)

set(CMAKE_CXX_STANDARD 14)

add_executable(exp1 dualization.h experiment_ord_double.cpp)
add_executable(exp2 dualization.h experiment_exist_all.cpp)