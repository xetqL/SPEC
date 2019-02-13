cmake_minimum_required(VERSION 3.10)
project(SPEC)

enable_language(CXX)
find_package(MPI REQUIRED)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_PREFIX_PATH  "~/libraries/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

set(HighFive_DIR "~/libraries/HighFive")

find_library(ZOLTAN
        NAMES zoltan libzoltan
        HINTS /usr/local ~/local/lib ${CMAKE_PREFIX_PATH}/Zoltan/lib
        PATHS ${ZOLTAN_LIBRARY_PATH})

find_path(ZOLTAN_INCLUDES
        NAMES zoltan zoltan.h
        HINTS ${CMAKE_PREFIX_PATH}/Zoltan/include)

find_package(HighFive 2.0 REQUIRED)

find_library(CNPY
        NAMES libcnpy.a libcnpy.so
        HINTS ${CMAKE_PREFIX_PATH}/cnpy/build/)

find_path(CNPY_INCLUDES
        NAMES cnpy.h
        HINTS ${CMAKE_PREFIX_PATH}/cnpy/build/)

set(CMAKE_SOURCE_DIR src)
set(CMAKE_INCLUDE_PATH include)

add_executable(SPEC ${CMAKE_SOURCE_DIR}/main.cpp ${CMAKE_SOURCE_DIR}/zupply.cpp ${CMAKE_SOURCE_DIR}/zupply.hpp include/zoltan_fn.hpp include/cell.hpp include/communication.hpp include/utils.hpp include/io.hpp)

target_include_directories(SPEC INTERFACE ${CMAKE_SOURCE_DIR}/include)
target_include_directories(SPEC PRIVATE ${ZOLTAN_INCLUDES})
target_include_directories(SPEC PUBLIC $<TARGET_PROPERTY:HighFive,INTERFACE_INCLUDE_DIRECTORIES>)
target_include_directories(SPEC PRIVATE ${CNPY_INCLUDES})

target_link_libraries(SPEC $<TARGET_PROPERTY:HighFive,INTERFACE_LINK_LIBRARIES>)
target_link_libraries(SPEC ${ZOLTAN})
target_link_libraries(SPEC ${CNPY})
