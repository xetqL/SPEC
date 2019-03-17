cmake_minimum_required(VERSION 3.10)
project(SPEC C CXX)

#enable_language(C CXX)
find_package(MPI REQUIRED)
set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++14 -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_PREFIX_PATH  "~/libraries/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg -O3")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg -O3")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg -O3")
set(HighFive_DIR "~/libraries/HighFive")
set(BINARY_NAME "SPEC" CACHE STRING "SPEC")

find_package (Threads)

if(DEFINED WITH_ZOLTAN)
    find_library(ZOLTAN
        NAMES zoltan libzoltan
        HINTS /usr/local ~/local/lib ${CMAKE_PREFIX_PATH}/Zoltan/lib ~/Zoltan_v3.83/build/lib
        PATHS ${ZOLTAN_LIBRARY_PATH})
    find_path(ZOLTAN_INCLUDES
        NAMES zoltan zoltan.h
        HINTS ${CMAKE_PREFIX_PATH}/Zoltan/include ~/Zoltan_v3.83/build/include)
endif()
#find_package(HighFive 2.0 REQUIRED)
if(DEFINED PRODUCE_OUTPUTS)
    find_library(CNPY
        NAMES libcnpy.a libcnpy.so
        HINTS ${CMAKE_PREFIX_PATH}/cnpy/build/)

    find_path(CNPY_INCLUDES
        NAMES cnpy.h
        HINTS ${CMAKE_PREFIX_PATH}/cnpy/)
    link_directories(${CMAKE_PREFIX_PATH}/cnpy/build/)

    message(${CNPY})
    message(${CNPY_INCLUDES})
    
endif()

set(CMAKE_SOURCE_DIR src)
set(CMAKE_INCLUDE_PATH include)

add_executable(SPEC ${CMAKE_SOURCE_DIR}/main.cpp ${CMAKE_SOURCE_DIR}/zupply.cpp ${CMAKE_SOURCE_DIR}/zupply.hpp include/cell.hpp include/communication.hpp include/utils.hpp include/io.hpp include/window.hpp include/CPULoadDatabase.hpp src/main.hpp include/time.hpp include/band_partitioner.hpp)

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${BINARY_NAME})

target_include_directories(SPEC INTERFACE ${CMAKE_SOURCE_DIR}/include)
if(DEFINED WITH_ZOLTAN)
    target_include_directories(SPEC PRIVATE ${ZOLTAN_INCLUDES})
endif()
#target_include_directories(SPEC PUBLIC $<TARGET_PROPERTY:HighFive,INTERFACE_INCLUDE_DIRECTORIES>)

#target_link_libraries(SPEC $<TARGET_PROPERTY:HighFive,INTERFACE_LINK_LIBRARIES>)
if(DEFINED WITH_ZOLTAN)
    target_link_libraries(SPEC ${ZOLTAN})
endif()

target_link_libraries(SPEC z)
target_link_libraries(SPEC ${CMAKE_THREAD_LIBS_INIT})

set(ignoreMe "${PRODUCE_OUTPUTS}")
if(DEFINED PRODUCE_OUTPUTS)
    message("Linking with io library")
    target_include_directories(SPEC PRIVATE ${CNPY_INCLUDES})
    target_link_libraries(SPEC ${CNPY})
    target_compile_definitions(SPEC PRIVATE PRODUCE_OUTPUTS=1)
else()
    message("Will not produce outputs")
endif()

if(DEFINED DISABLE_LB)
    target_compile_definitions(SPEC PRIVATE LB_ENABLED=1)
endif()

if(DEFINED LB_METHOD)
    message("Compiled with: Load Balancing Flags = ${LB_METHOD}")
    target_compile_definitions(SPEC PRIVATE LB_METHOD=${LB_METHOD})
else()
    message("Compiled with: No Load Balancing Flags!")
endif()
if(DEFINED WITH_ZOLTAN)
    message("Compiled with Zoltan")
    target_compile_definitions(SPEC PRIVATE WITH_ZOLTAN=1)
else()
    message("Compiled without Zoltan")
endif()