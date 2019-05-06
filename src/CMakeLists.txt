cmake_minimum_required(VERSION 3.10)
project(SPEC C CXX)

#enable_language(C CXX)
find_package(MPI REQUIRED)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_PREFIX_PATH  "~/libraries/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS}")
set(HighFive_DIR "~/libraries/HighFive")
set(BINARY_NAME "SPEC" CACHE STRING "SPEC")
set(LOAD_BALANCING_CALL "AUTONOMIC" CACHE STRING "CACA")
message(${LOAD_BALANCING_CALL})
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

add_executable(SPEC main.cpp zupply.cpp SimulatedLBM.cpp LoadBalancer.cpp Cell.cpp Utils.cpp ../include/ZoltanLoadBalancer.hpp ../include/SimulationParams.hpp ../include/Simflow.hpp Simflow.cpp CLIParser.cpp ../include/CLIParser.hpp GossipPropagator.cpp ../include/GossipPropagator.hpp)

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME ${BINARY_NAME})

if(DEFINED WITH_ZOLTAN)
    message("Compiled with Zoltan")
    target_include_directories(SPEC PRIVATE ${ZOLTAN_INCLUDES})
    target_link_libraries(SPEC ${ZOLTAN})
    target_compile_definitions(SPEC PRIVATE WITH_ZOLTAN=1)
else()
    message("Compiled without Zoltan")
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

if(LOAD_BALANCING_CALL STREQUAL "NONE")
    message("Load balancing call strategy = No load balancing")
elseif(LOAD_BALANCING_CALL STREQUAL "CYCLIC")
    message("Load balancing call strategy = Cyclic")
    message("Cycle size = 100")
    target_compile_definitions(SPEC PRIVATE LB_CYCLE_SIZE=100)
    target_compile_definitions(SPEC PRIVATE CYCLIC_LOAD_BALANCING=1)
else()
    message("Load balancing call strategy = Autonomic")
    target_compile_definitions(SPEC PRIVATE AUTONOMIC_LOAD_BALANCING=1)
endif()