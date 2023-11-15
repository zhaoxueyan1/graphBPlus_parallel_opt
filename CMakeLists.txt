project(graphBplus_OpenMP)

cmake_minimum_required(VERSION 2.8)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11 -fopenmp")

add_executable(graphBplus_OpenMP graphBplus_02.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -std=c++11")

add_executable(graphBplus_serial graphBplus_02.cpp)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

add_executable(graphBplus_serial_unfold graphBplus_02_serial_loop-unfold.cpp)