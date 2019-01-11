cmake_minimum_required(VERSION 3.1)
project(sparse_solvers)

set (CMAKE_CXX_STANDARD 14)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(WITH_MEMORY_PRINT "Enable printing memory management of C++ objects" OFF)
option(WITH_DEBUG "Enable debugging" ON)
option(WITH_TESTS "Enable testing" ON)

if(WITH_MEMORY_PRINT)
    add_definitions(-DMEMORY_PRINT)
endif()


if(WITH_DEBUG)
    add_compile_options(-O0)
    add_definitions(-DDEBUG)
    add_definitions(${CMAKE_CXX_FLAGS} "-g")
    add_definitions(${CMAKE_CXX_FLAGS} "-Wall")
else()
    add_compile_options(-O3)
    add_definitions(-DNDEBUG)
endif()

add_compile_options(-Wno-writable-strings)
# Catch2 
add_subdirectory(Catch2)
# minitrace
add_definitions(-DMTR_ENABLED)
# Eigen
add_subdirectory(eigen)

file(GLOB SRCFILES 
    src/*.h
    src/*.cpp)

add_executable(main main.cpp ${SRCFILES})
target_link_libraries(main Eigen3::Eigen)

if (WITH_TESTS)
    add_executable(tests tests.cpp  ${SRCFILES})
    target_link_libraries(tests Catch2::Catch2 Eigen3::Eigen)
endif()
