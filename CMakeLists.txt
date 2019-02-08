cmake_minimum_required(VERSION 3.13)
project(Numerical_Analysis)

set(CMAKE_CXX_STANDARD 17)

add_executable(Numerical_Analysis project/src/main.cpp project/src/OrdinaryDE.cpp project/include/OrdinaryDE.h)