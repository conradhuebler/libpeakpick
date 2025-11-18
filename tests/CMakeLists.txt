# Test CMakeLists.txt for libpeakpick
cmake_minimum_required(VERSION 3.1)

# Enable testing
enable_testing()

# Include the library headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../eigen/)

# Find OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Set C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define test executables
add_executable(test_spectrum test_spectrum.cpp)
add_executable(test_mathhelper test_mathhelper.cpp)
add_executable(test_analyse test_analyse.cpp)
add_executable(test_baseline test_baseline.cpp)
add_executable(test_savitzky test_savitzky.cpp)
add_executable(test_peakpick test_peakpick.cpp)
add_executable(test_integration test_integration.cpp)
add_executable(test_advanced test_advanced.cpp)

# Set C++11 for all test targets
set_property(TARGET test_spectrum PROPERTY CXX_STANDARD 11)
set_property(TARGET test_mathhelper PROPERTY CXX_STANDARD 11)
set_property(TARGET test_analyse PROPERTY CXX_STANDARD 11)
set_property(TARGET test_baseline PROPERTY CXX_STANDARD 11)
set_property(TARGET test_savitzky PROPERTY CXX_STANDARD 11)
set_property(TARGET test_peakpick PROPERTY CXX_STANDARD 11)
set_property(TARGET test_integration PROPERTY CXX_STANDARD 11)
set_property(TARGET test_advanced PROPERTY CXX_STANDARD 11)

# Register tests with CTest
add_test(NAME SpectrumTests COMMAND test_spectrum)
add_test(NAME MathHelperTests COMMAND test_mathhelper)
add_test(NAME AnalyseTests COMMAND test_analyse)
add_test(NAME BaselineTests COMMAND test_baseline)
add_test(NAME SavitzkyTests COMMAND test_savitzky)
add_test(NAME PeakPickTests COMMAND test_peakpick)
add_test(NAME IntegrationTests COMMAND test_integration)
add_test(NAME AdvancedTests COMMAND test_advanced)
