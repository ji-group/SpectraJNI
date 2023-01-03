cmake_minimum_required(VERSION 3.5)
project(SpectraJNI)

set(CMAKE_CXX_STANDARD 17)

add_library(SpectraJNI SHARED library.cpp)

target_link_libraries(SpectraJNI Eigen3::Eigen)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (Spectra REQUIRED NO_MODULE)

if (NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "No Eigen3 Found")
endif (NOT TARGET Eigen3::Eigen)

add_subdirectory(test)