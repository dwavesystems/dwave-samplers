find_package(Boost 1.37 REQUIRED)
find_package(MATLAB COMPONENTS mex REQUIRED)

set(MEX_EXTRA_ARGS "" CACHE STRING "Extra arguments for mex command")

set(ORANG_MINSUM_MEX ${CMAKE_CURRENT_BINARY_DIR}/orang_minsum.${MATLAB_MEXEXT})
set(ORANG_MINCOUNT_MEX
  ${CMAKE_CURRENT_BINARY_DIR}/orang_mincount.${MATLAB_MEXEXT})
set(ORANG_SAMPLE_MEX ${CMAKE_CURRENT_BINARY_DIR}/orang_sample.${MATLAB_MEXEXT})
set(ORANG_GREEDYVARORDER_MEX
    ${CMAKE_CURRENT_BINARY_DIR}/orang_greedyvarorder.${MATLAB_MEXEXT})
set(ORANG_MEX_OBJ
    ${CMAKE_CURRENT_BINARY_DIR}/orang_mex${CMAKE_CXX_OUTPUT_EXTENSION})

add_custom_command(
  OUTPUT ${ORANG_MEX_OBJ}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/orang_mex.cpp
  IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/orang_mex.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${MATLAB_MEX_EXECUTABLE} -c -largeArrayDims
          ${CMAKE_CURRENT_SOURCE_DIR}/orang_mex.cpp
          -I${CMAKE_SOURCE_DIR}/src -I${Boost_INCLUDE_DIR}
          -DBOOST_ALL_NO_LIB
          ${MEX_EXTRA_ARGS} VERBATIM)

add_custom_command(
  OUTPUT ${ORANG_MINSUM_MEX}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/orang_minsum.cpp ${ORANG_MEX_OBJ}
  IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/orang_minsum.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${MATLAB_MEX_EXECUTABLE} -largeArrayDims
          ${CMAKE_CURRENT_SOURCE_DIR}/orang_minsum.cpp ${ORANG_MEX_OBJ}
          -DBOOST_ALL_NO_LIB
          -I${CMAKE_SOURCE_DIR}/src -I${Boost_INCLUDE_DIR} ${MEX_EXTRA_ARGS}
          VERBATIM)

add_custom_command(
  OUTPUT ${ORANG_MINCOUNT_MEX}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/orang_mincount.cpp ${ORANG_MEX_OBJ}
  IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/orang_mincount.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${MATLAB_MEX_EXECUTABLE} -largeArrayDims
          ${CMAKE_CURRENT_SOURCE_DIR}/orang_mincount.cpp ${ORANG_MEX_OBJ}
          -DBOOST_ALL_NO_LIB
          -I${CMAKE_SOURCE_DIR}/src -I${Boost_INCLUDE_DIR} ${MEX_EXTRA_ARGS}
          VERBATIM)

add_custom_command(OUTPUT ${ORANG_SAMPLE_MEX}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/orang_sample.cpp ${ORANG_MEX_OBJ}
  IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/orang_sample.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${MATLAB_MEX_EXECUTABLE} -largeArrayDims
          ${CMAKE_CURRENT_SOURCE_DIR}/orang_sample.cpp ${ORANG_MEX_OBJ}
          -DBOOST_ALL_NO_LIB
          -I${CMAKE_SOURCE_DIR}/src -I${Boost_INCLUDE_DIR} ${MEX_EXTRA_ARGS}
          VERBATIM)

add_custom_command(OUTPUT ${ORANG_GREEDYVARORDER_MEX}
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/orang_greedyvarorder.cpp ${ORANG_MEX_OBJ}
  IMPLICIT_DEPENDS CXX ${CMAKE_CURRENT_SOURCE_DIR}/orang_greedyvarorder.cpp
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMAND ${MATLAB_MEX_EXECUTABLE} -largeArrayDims
          ${CMAKE_CURRENT_SOURCE_DIR}/orang_greedyvarorder.cpp ${ORANG_MEX_OBJ}
          -DBOOST_ALL_NO_LIB
          -I${CMAKE_SOURCE_DIR}/src -I${Boost_INCLUDE_DIR} ${MEX_EXTRA_ARGS}
          VERBATIM)

add_custom_target(orang_minsum ALL DEPENDS ${ORANG_MINSUM_MEX})
add_custom_target(orang_mincount ALL DEPENDS ${ORANG_MINCOUNT_MEX})
add_custom_target(orang_sample ALL DEPENDS ${ORANG_SAMPLE_MEX})
add_custom_target(orang_greedyvarorder ALL DEPENDS ${ORANG_GREEDYVARORDER_MEX})

#add_dependencies(orang_greedyvarorder orang_sample)
#add_dependencies(orang_sample orang_minsum)

install(FILES ${ORANG_MINSUM_MEX}
              ${ORANG_MINCOUNT_MEX}
              ${ORANG_SAMPLE_MEX}
              ${ORANG_GREEDYVARORDER_MEX}
              orang_minsum.m
              orang_mincount.m
              orang_sample.m
              orang_greedyvarorder.m
              quboTables.m
              isingTables.m
              adjacencyTables.m
              chimeraVarOrder.m
              fattreeVarOrder.m
        DESTINATION matlab)

