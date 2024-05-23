set(CMAKE_MAKE_PROGRAM make)
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDART_REQUIRED true)

cmake_minimum_required(VERSION 3.15)
project(Genetic-algorithm CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_executable(TSP src/TSP.cpp src/TSP.hpp)
