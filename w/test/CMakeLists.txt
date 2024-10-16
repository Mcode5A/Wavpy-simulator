enable_testing()

# Boost version
#set(Boost_ADDITIONAL_VERSIONS "1.58")

# find Boost, trigger error otherwise
find_package(Boost COMPONENTS unit_test_framework OPTIONAL_COMPONENTS)
#find_library(HABEMUS_BOOST_TEST boost_unit_test_framework)
#if (HABEMUS_BOOST_TEST)
#  message("Boost test available")
#else()
#  message( "Boost test not available")
#endif()


if (Boost_FOUND)
  message("    Boost Found. Test suite code will be compiled, and testing will be possible.")
  # create the executable
  add_executable(testexample testexample.cpp)
  add_executable(test_ReferenceFrame test_ReferenceFrame.cpp)
  add_executable(test_ComplexSignal test_ComplexSignal.cpp)
  add_executable(test_ComplexWaveform test_ComplexWaveform.cpp)
  add_executable(boost_test_functions boost_test_functions.cpp)
  
  # indicate the include paths
  target_include_directories(testexample           PRIVATE ${Boost_INCLUDE_DIRS} PUBLIC "${CMAKE_SOURCE_DIR}/src/inc")
  target_include_directories(test_ReferenceFrame   PRIVATE ${Boost_INCLUDE_DIRS} PUBLIC "${CMALE_SOURCE_DIR}/src/inc")
  target_include_directories(test_ComplexSignal    PRIVATE ${Boost_INCLUDE_DIRS} PUBLIC "${CMALE_SOURCE_DIR}/src/inc")
  target_include_directories(test_ComplexWaveform  PRIVATE ${Boost_INCLUDE_DIRS} PUBLIC "${CMALE_SOURCE_DIR}/src/inc")
  target_include_directories(boost_test_functions  PRIVATE ${Boost_INCLUDE_DIRS} PUBLIC "${CMALE_SOURCE_DIR}/src/inc")
  
  # indicates the shared library variant
  target_compile_definitions(testexample           PRIVATE "BOOST_TEST_DYN_LINK=1")
  target_compile_definitions(test_ReferenceFrame   PRIVATE "BOOST_TEST_DYN_LINK=1")
  target_compile_definitions(test_ComplexSignal    PRIVATE "BOOST_TEST_DYN_LINK=1")
  target_compile_definitions(test_ComplexWaveform  PRIVATE "BOOST_TEST_DYN_LINK=1")
  target_compile_definitions(boost_test_functions  PRIVATE "BOOST_TEST_DYN_LINK=1")
  
  # indicates the link paths
  target_link_libraries(testexample          LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Python3::Python)
  target_link_libraries(test_ReferenceFrame  LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Python3::Python)
  target_link_libraries(test_ComplexSignal   LINK_PUBLIC waveform_pylib gsl gslcblas stdc++ ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Python3::Python)
  target_link_libraries(test_ComplexWaveform LINK_PUBLIC waveform_pylib gsl gslcblas stdc++ ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Python3::Python)
  target_link_libraries(boost_test_functions LINK_PUBLIC waveform_pylib gsl gslcblas stdc++ ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} Python3::Python)
  
  #message("    Hola ${Boost_INCLUDE_DIRS} // ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY} ")
  # declares a test with our executable
  add_test(NAME TestExample      COMMAND testexample )
  add_test(NAME ReferenceFrame   COMMAND test_ReferenceFrame )
  add_test(NAME ComplexSignal    COMMAND test_ComplexSignal )
  add_test(NAME ComplexWaveform  COMMAND test_ComplexWaveform )
  add_test(NAME boost_test_functions  COMMAND boost_test_functions )
  
else()
  message("    Boost NOT Found. Testsuit code will not be compiled. Testing will not be possible")
endif()



