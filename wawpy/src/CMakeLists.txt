# Revisar aquesta part...


# Define the sources for the Library
file(GLOB SOURCES "*.cpp" "*.f90")
# Generate the shared library from the sources
add_library(waveform_pylib SHARED ${SOURCES})

configure_file(../waveform_pylib_config.h.in waveform_pylib_config.h)


# Include directories
target_include_directories(waveform_pylib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/inc PRIVATE ${Python3_INCLUDE_DIRS})

# Linking to other libraries
target_link_libraries(waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
# Installation 
#Set the location for library installation -- i.e., /usr/lib in this case
#install(TARGETS waveform_pylib DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME}")
#install(DIRECTORY "inc/" DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
#install(DIRECTORY "data_share/" DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}")
install(TARGETS waveform_pylib DESTINATION "${CMAKE_INSTALL_LIBDIR}")
install(DIRECTORY "inc/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
#install(DIRECTORY "data_share/" DESTINATION "${CMAKE_INSTALL_DATADIR}")
