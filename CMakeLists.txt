cmake_minimum_required(VERSION 3.9)
project(sparse_solvers)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(WITH_MEMORY_PRINT    "Enable printing C++ objects' lifetime"         OFF)
option(WITH_NUMERICAL_PRINT "Enable printing numerical related issues"      OFF)
option(WITH_DEBUG           "Enable debugging"                              OFF)
option(WITH_TESTS           "Enable testing"                                ON)

if(WITH_MEMORY_PRINT)
    add_definitions(-DMEMORY_PRINT)
endif()
if(WITH_NUMERICAL_PRINT)
    add_definitions(-DNUMERICAL_PRINT)
endif()

add_compile_options(-Wno-writable-strings)

if(WITH_DEBUG)
    add_compile_options(-O0)
    add_definitions(-DDEBUG)
    add_definitions(${CMAKE_CXX_FLAGS} "-g")
    add_definitions(${CMAKE_CXX_FLAGS} "-Wall")
else()
    add_compile_options(-O3)
    add_definitions(-DNDEBUG)
    add_definitions(${CMAKE_CXX_FLAGS} "-march=native")
endif()

# OpenMP
SET(OpenMP_C_FLAGS          "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
SET(OpenMP_C_LIB_NAMES      "omp")
SET(OpenMP_CXX_FLAGS        "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include")
SET(OpenMP_CXX_LIB_NAMES    "omp")
SET(OpenMP_omp_LIBRARY      "/usr/local/opt/libomp/lib/libomp.dylib")

find_package(OpenMP REQUIRED)

# submodule dependencies
add_subdirectory(Catch2)
add_subdirectory(eigen)

# minitrace
add_definitions(-DMTR_ENABLED)

file(GLOB SRCFILES 
    src/*.h
    src/*.cpp)

add_executable(main main.cpp ${SRCFILES})
target_link_libraries(main Eigen3::Eigen OpenMP::OpenMP_CXX)

if (WITH_TESTS)
    add_executable(tests tests.cpp  ${SRCFILES})
    target_link_libraries(tests Catch2::Catch2 Eigen3::Eigen OpenMP::OpenMP_CXX)
endif()
