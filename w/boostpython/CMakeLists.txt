FIND_PACKAGE(PythonInterp 3)
message(STATUS "PYTHON_VERSION_MAJOR: ${PYTHON_VERSION_MAJOR}")
message(STATUS "PYTHON_VERSION_MINOR: ${PYTHON_VERSION_MINOR}")

FIND_PACKAGE(PythonLibs 3 REQUIRED)
message(STATUS "PYTHON_LIBRARY: ${PYTHON_LIBRARY}")
message(STATUS "PYTHON_INCLUDE_DIR: ${PYTHON_INCLUDE_DIR}")

#set(PYTHON_VERSION_SUFFIX "3") #TODO This should be detected automatically
#set(NUMPY_VERSION_SUFFIX "") #TODO This should be detected automatically
#message(STATUS "PYTHON_VERSION_SUFFIX: ${PYTHON_VERSION_SUFFIX}")

if (PYTHONINTERP_FOUND)
  if (UNIX AND NOT APPLE)
    if (PYTHON_VERSION_MAJOR EQUAL 3)
      FIND_PACKAGE(Boost 1.65 COMPONENTS python${PYTHON_VERSION_MAJOR} numpy${PYTHON_VERSION_MAJOR} )
      
      #message(STATUS "NUMPY_VERSION_SUFFIX: ${NUMPY_VERSION_SUFFIX}")
      #message(STATUS "NUMPY_VERSION_MAJOR: ${NUMPY_VERSION_MAJOR}")
      #message(STATUS "NUMPY_VERSION_MINOR: ${NUMPY_VERSION_MINOR}")

    else()
      message(WARNING "Python 3 required: Interface to python will not be compiled")
    endif()
  else()
    message(WARNING "Interface to python not implemented yet for operating systems different from LINUX")
  endif()
endif()

message(STATUS "PYTHON_LIBRARIES = ${PYTHON_LIBRARIES}")
message(STATUS "PYTHON_EXECUTABLE = ${PYTHON_EXECUTABLE}")
message(STATUS "PYTHON_INCLUDE_DIRS = ${PYTHON_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}")


#ENABLE_TESTING()
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/../src/inc)
LINK_LIBRARIES(${Boost_LIBRARIES} ${PYTHON_LIBRARIES} waveform_pylib) # Deprecated but so convenient!


PYTHON_ADD_MODULE(referenceframe ReferenceFrame_boostpython.cpp ndarray_gsl.cpp)

PYTHON_ADD_MODULE(ctor ctor_boostpython.cpp)
PYTHON_ADD_MODULE(ctor2 ctor2_boostpython.cpp)

set(PYTHON_MODULES_DIR "${CMAKE_INSTALL_LIBDIR}wp${PYTHON_VERSION_SUFFIX}/")
message(STATUS "PYTHON_MODULES_DIR for installation: ${PYTHON_MODULES_DIR}")

install(TARGETS referenceframe DESTINATION "${PYTHON_MODULES_DIR}")
install(TARGETS ctor DESTINATION "${PYTHON_MODULES_DIR}")
install(TARGETS ctor2 DESTINATION "${PYTHON_MODULES_DIR}")

#FILE(COPY ReferenceFrame.py DESTINATION .)
# Revisar aquesta part...


# Define the sources for the wavpy interface
##file(GLOB SOURCES_WAVPY "*.cpp")
##message(STATUS "boostpython sources: " ${SOURCES_WAVPY} ) 

# Generate the shared library from the sources

##add_library(wavpy SHARED ${SOURCES_WAVPY} ${SOURCES} )


# Include directories
##target_include_directories(wavpy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../inc)
##message(STATUS "CMAKE_CURRENT_SOURCE_DIR=" ${CMAKE_CURRENT_SOURCE_DIR})

# Installation 
#Set the location for library installation -- i.e., /usr/lib in this case
#install(TARGETS waveform_pylib DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME}")
#install(DIRECTORY "inc/" DESTINATION "${CMAKE_INSTALL_PREFIX}/include/${PROJECT_NAME}")
#install(DIRECTORY "data_share/" DESTINATION "${CMAKE_INSTALL_PREFIX}/share/${PROJECT_NAME}")
###install(TARGETS waveform_pylib DESTINATION "${CMAKE_INSTALL_LIBDIR}")
###install(DIRECTORY "inc/" DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
#install(DIRECTORY "data_share/" DESTINATION "${CMAKE_INSTALL_DATADIR}")

