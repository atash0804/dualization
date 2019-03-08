cmake_minimum_required(VERSION 3.13)
project(dualization)

set(CMAKE_CXX_STANDARD 14)

add_executable(dualization dualization.cpp experiment_ord_double.cpp dualization.h)