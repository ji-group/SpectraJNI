cmake_minimum_required(VERSION 3.5)
project(SpectraJNI)

set(CMAKE_CXX_STANDARD 17)

add_library(spectra-jni SHARED
        dr_evomodel_substmodel_spectra_SpectraJNIWrapper.cpp
        dr_evomodel_substmodel_spectra_SpectraJNIWrapper.h)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
find_package (Spectra REQUIRED NO_MODULE)
target_link_libraries(spectra-jni Eigen3::Eigen)
target_link_libraries(spectra-jni Spectra::Spectra)
if (NOT TARGET Eigen3::Eigen)
    message(FATAL_ERROR "No Eigen3 Found")
endif (NOT TARGET Eigen3::Eigen)
if (NOT TARGET Spectra::Spectra)
    message(FATAL_ERROR "No Spectra Found")
endif (NOT TARGET Spectra::Spectra)

find_package(JNI)
if(JNI_FOUND)
    message (STATUS "JAVA_HOME=$ENV{JAVA_HOME}")
    message (STATUS "JNI_INCLUDE_DIRS=${JNI_INCLUDE_DIRS}")
    message (STATUS "JNI_LIBRARIES=${JNI_LIBRARIES}")

    include_directories(
            $ENV{JAVA_HOME}/include
            ${JNI_INCLUDE_DIRS}
    )
else(JNI_FOUND)
    if(BUILD_JNI)
        message (FATAL_ERROR "No JNI includes and libraries found")
    endif(BUILD_JNI)
endif(JNI_FOUND)

if(APPLE)
    set_target_properties(spectra-jni PROPERTIES PREFIX "lib" SUFFIX ".jnilib")
elseif(WIN32)
    set_target_properties(spectra-jni PROPERTIES PREFIX "" SUFFIX "64.dll")
elseif(UNIX)
    set_target_properties(spectra-jni PROPERTIES PREFIX "lib" SUFFIX ".so")
endif(APPLE)

add_subdirectory(test)

set(BEAGLE_INSTALL_DIR lib)
set(BEAGLE_HEADER_DIR include/libhmsbeagle-1/libhmsbeagle)
install(TARGETS spectra-jni
        DESTINATION ${BEAGLE_INSTALL_DIR}
        )