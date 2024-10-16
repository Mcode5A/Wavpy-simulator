# Generate the example programs
add_executable(prova   prova.cpp)
add_executable(prova_bibaspirfile prova_bibaspirfile.cpp)
add_executable(prova_ComplexSignal prova_ComplexSignal.cpp)
add_executable(prova_ComplexWaveform prova_ComplexWaveform.cpp)
add_executable(plotsomething plotsomething.cpp)
add_executable(check_Mallorca check_Mallorca.cpp)
add_executable(check_Mallorca_ddm check_Mallorca_ddm.cpp)
add_executable(check_Mallorca_CR check_Mallorca_CR.cpp)
add_executable(check_Mallorca_peak check_Mallorca_peak.cpp)
add_executable(check_Mallorca_waveforms check_Mallorca_waveforms.cpp)

# Include directories
target_include_directories(prova PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(prova_bibaspirfile PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(prova_ComplexSignal PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(prova_ComplexWaveform PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(plotsomething PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(check_Mallorca PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(check_Mallorca_ddm PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(check_Mallorca_CR PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(check_Mallorca_peak PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(check_Mallorca_waveforms PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")

# Link the tests to the waveform_pylib
target_link_libraries(prova                 LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(prova_bibaspirfile    LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(prova_ComplexSignal   LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(prova_ComplexWaveform LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(plotsomething         LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(check_Mallorca        LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(check_Mallorca_ddm    LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(check_Mallorca_CR     LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(check_Mallorca_peak   LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)
target_link_libraries(check_Mallorca_waveforms LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++ Python3::Python)

# Installation for the examples
install(TARGETS prova prova_bibaspirfile prova_ComplexSignal prova_ComplexWaveform plotsomething
                check_Mallorca check_Mallorca_ddm check_Mallorca_CR check_Mallorca_peak check_Mallorca_waveforms
        DESTINATION "${PROVA_DIR}")

