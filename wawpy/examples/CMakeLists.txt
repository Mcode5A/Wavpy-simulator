# Generate the example programs
add_executable(ComputeSPIR_realdataAutocorr   ComputeSPIR_realdataAutocorr.cc)
add_executable(GenerateGoldRTR_dataset        GenerateGoldRTR_dataset.cc)
add_executable(GenerateSPIR_dataset           GenerateSPIR_dataset.cc)
add_executable(GenerateSPIR_dataset_CR_GPS_L1 GenerateSPIR_dataset_CR_GPS_L1.cc)
#add_executable(ModelSPIR_altimetry            ModelSPIR_altimetry.cc)
#add_executable(ProcessDELTEBRE_phase          ProcessDELTEBRE_phase.cc)
#add_executable(ProcessLongIntegration         ProcessLongIntegration.cc)
#add_executable(ProcessSPIR_altimetry          ProcessSPIR_altimetry.cc)
#add_executable(ProcessTIGRISfreq_GOLD         ProcessTIGRISfreq_GOLD.cc)
#add_executable(ProcessTIGRIS_GOLD             ProcessTIGRIS_GOLD.cc)
#add_executable(ProcessTIGRIS_PIR              ProcessTIGRIS_PIR.cc)
add_executable(example_geometry_class         example_geometry_class.cc)
add_executable(example_surface_class          example_surface_class.cc)
add_executable(example_wavcomplex_class       example_wavcomplex_class.cc)
add_executable(example_waveform_class         example_waveform_class.cc)
#add_executable(Visibility_GNSS                Visibility_GNSS.cc)
add_executable(Visibility_track_GNSS          Visibility_track_GNSS.cc)

# Include directories
target_include_directories(example_wavcomplex_class PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
target_include_directories(example_waveform_class   PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/")
#target_include_directories(Visibility_GNSS       PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/inc/" "${CMAKE_CURRENT_SOURCE_DIR}/recipes_cpp")

# Link the tests to the waveform_pylib
target_link_libraries(ComputeSPIR_realdataAutocorr      LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(GenerateGoldRTR_dataset           LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(GenerateSPIR_dataset              LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(GenerateSPIR_dataset_CR_GPS_L1    LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ModelSPIR_altimetry               LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ProcessDELTEBRE_phase             LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ProcessLongIntegration            LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ProcessSPIR_altimetry             LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ProcessTIGRISfreq_GOLD            LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ProcessTIGRIS_GOLD                LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(ProcessTIGRIS_PIR                 LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(example_geometry_class            LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(example_surface_class             LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(example_wavcomplex_class          LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(example_waveform_class            LINK_PUBLIC waveform_pylib) # gsl gslcblas m fftw3 gfortran stdc++)
#target_link_libraries(Visibility_GNSS                   LINK_PUBLIC waveform_pylib gsl gslcblas m fftw3 gfortran stdc++)
target_link_libraries(Visibility_track_GNSS             LINK_PUBLIC waveform_pylib) #gsl gslcblas m fftw3 gfortran stdc++)


# Installation for the examples
install(TARGETS ComputeSPIR_realdataAutocorr
                GenerateGoldRTR_dataset 
                GenerateSPIR_dataset
                GenerateSPIR_dataset_CR_GPS_L1
#                ModelSPIR_altimetry
#                ProcessDELTEBRE_phase
#                ProcessLongIntegration
#                ProcessSPIR_altimetry
#                ProcessTIGRISfreq_GOLD
#                ProcessTIGRIS_GOLD
#                ProcessTIGRIS_PIR
                example_geometry_class
                example_surface_class
                example_wavcomplex_class
                example_waveform_class
#                Visibility_GNSS
                Visibility_track_GNSS
        DESTINATION "${EXAMPLES_DIR}")

