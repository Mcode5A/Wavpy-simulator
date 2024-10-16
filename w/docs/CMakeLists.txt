# indicate the documentation build as an option that is set to ON by default
option(BUILD_DOC "Build documentation" OFF)

# Check if doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
  # set input and output files
  set(DOXYGEN_IN  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
  set(DOXYGEN_OUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)

  # request to configure file
  configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
  message("Doxygen build started")

  # note the option ALL which allows to build the docs together with the applications
  if (BUILD_DOC)
    add_custom_target( doc_doxygen ALL
                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generate API documentation with Doxygen"
                       VERBATIM )
  else(BUILD_DOC)
    add_custom_target( doc_doxygen
                       COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                       COMMENT "Generate API documentation with Doxygen"
                       VERBATIM )
  endif (BUILD_DOC)
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation.")
  message("No documentation will be generated")
endif (DOXYGEN_FOUND)

