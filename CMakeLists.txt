cmake_minimum_required(VERSION 3.16)
project(dids)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -mavx2 -pthread -fopenmp")
set(CMAKE_BUILD_TYPE Release)

FIND_PACKAGE(OpenMP REQUIRED)
add_compile_options(-fopenmp)

add_executable(example example.cpp random_data.cpp)
target_link_libraries(example tcmalloc gsl)