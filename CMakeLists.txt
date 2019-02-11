cmake_minimum_required(VERSION 3.10)
project(SPEC)

enable_language(CXX)
find_package(MPI REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_PREFIX_PATH  "~/libraries/")

set(CMAKE_SOURCE_DIR src)
set(CMAKE_INCLUDE_PATH include)

add_executable(SPEC ${CMAKE_SOURCE_DIR}/main.cpp ${CMAKE_SOURCE_DIR}/zupply.cpp ${CMAKE_SOURCE_DIR}/zupply.hpp include/zoltan_fn.hpp include/cell.hpp include/communication.hpp include/utils.hpp)