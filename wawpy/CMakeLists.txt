cmake_minimum_required(VERSION 2.8.9)

# Set the project name and version
cmake_policy(SET CMP0048 NEW)
project(waveform_pylib VERSION 2.1)


# Set the programming languages
# C++
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall")
# Fortran
enable_language(Fortran)
set(CMAKE_Fortran_FLAGS "-Wall")

# Check all needed libraries: TBD
find_library(HABEMUS_GSL gsl)
if(NOT HABEMUS_GSL)
  message(FATAL_ERROR "gsl library not found.")
endif()

find_library(HABEMUS_GSLCBLAS gslcblas)
if(NOT HABEMUS_GSLCBLAS)
  message(FATAL_ERROR "gslcblas library not found.")
endif()

find_library(HABEMUS_FFTW3 fftw3)
if(NOT HABEMUS_FFTW3)
  message(FATAL_ERROR "fftw3 library not found.")
endif()

find_package (Python3 COMPONENTS Development)
if (NOT Python3_FOUND)
  message(WARNING "PYTHON libraries not found. Plotting functions will not be available.")
endif()
#find_library(HABEMUS_MATPLOTLIB )

# SWIG TBD


set(CMAKE_BUILD_TYPE Debug)



# Define installation directories
set(CMAKE_INSTALL_LIBDIR     "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/lib/")
set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/include/")
set(CMAKE_INSTALL_DATADIR    "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/share/")
set(CMAKE_INSTALL_BINDIR     "${CMAKE_INSTALL_PREFIX}/${PROJECT_NAME}/bin/")
set(EXAMPLES_DIR             "${CMAKE_INSTALL_BINDIR}/examples/")
set(PROVA_DIR                "${CMAKE_INSTALL_BINDIR}/prova/")

set(LOCAL_DATA_PATH "${CMAKE_INSTALL_DATADIR}")

# The library
add_subdirectory(src)

# Testing 
add_subdirectory(test)

# The examples
# add_subdirectory(examples)

# BoostPython (to have the python interface)
add_subdirectory(boostpython)

# Proves per a comprovar la classes visitant (Process)
add_subdirectory(prova)

#target_include_directories(waveform_pylib PUBLIC "${PROJECT_BINARY_DIR}") 

# The documentation
add_subdirectory(docs)

