cmake_minimum_required(VERSION 3.1)

# Include library headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eigen/)

# Find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Example executables
add_executable(example_batch example_batch.cpp)
add_executable(example_fitting example_fitting.cpp)
add_executable(example_comparison example_comparison.cpp)
add_executable(example_export example_export.cpp)

set_property(TARGET example_batch PROPERTY CXX_STANDARD 11)
set_property(TARGET example_fitting PROPERTY CXX_STANDARD 11)
set_property(TARGET example_comparison PROPERTY CXX_STANDARD 11)
set_property(TARGET example_export PROPERTY CXX_STANDARD 11)
