set(SWIG_WRAPPER orang_wrap.cpp)
set(SWIG_WRAPPER_PY orang.py)
set(SWIG_WRAPPER_PYC orang.pyc)

find_package(Boost 1.37 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
add_definitions(-DBOOST_ALL_NO_LIB)

add_custom_command(OUTPUT ${SWIG_WRAPPER} ${SWIG_WRAPPER_PY}
  COMMAND ${SWIG_EXECUTABLE} -c++ -python -o ${SWIG_WRAPPER}
    ${CMAKE_CURRENT_SOURCE_DIR}/orang.i
  DEPENDS orang.i python-api.h VERBATIM
  COMMENT "Generating orang SWIG wrapper")

add_custom_command(OUTPUT ${SWIG_WRAPPER_PYC}
  COMMAND ${PYTHON_EXECUTABLE} -m py_compile ${SWIG_WRAPPER_PY}
  DEPENDS ${SWIG_WRAPPER_PY} VERBATIM
  COMMENT "Compiling orang Python bytecode")

find_path(NUMPY_INCLUDEDIR numpy/arrayobject.h)
include_directories(${NUMPY_INCLUDEDIR} ${Boost_INCLUDE_DIR})

add_library(_orang MODULE
  ${SWIG_WRAPPER}
  conversions.cpp
  solve.cpp
  sample.cpp)

set_target_properties(_orang PROPERTIES
  ${PYTHON_EXTENSION_PROPERTIES}
  PREFIX ""
  COMPILE_FLAGS "${PYTHON_CXXFLAGS}"
  LINK_FLAGS "${PYTHON_LDFLAGS}")
target_link_libraries(_orang ${PYTHON_LIBRARIES})

if(CMAKE_BUILD_TYPE MATCHES "Release|MinSizeRel")
  set(STRIP_ARGS)
  if(CMAKE_CXX_COMPILER_ID STREQUAL Clang)
    set(STRIP_ARGS "-x")
  endif()
  if(CMAKE_STRIP)
    add_custom_command(TARGET _orang POST_BUILD
      COMMAND ${CMAKE_STRIP} ${STRIP_ARGS} $<TARGET_FILE:_orang>
      COMMENT "Stripping orang module"
      VERBATIM)
  endif()
endif()

execute_process(COMMAND ${PYTHON_EXECUTABLE}
  ${CMAKE_CURRENT_SOURCE_DIR}/build_egg.py --egg-name
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  RESULT_VARIABLE EGGFILE_RESULT
  OUTPUT_VARIABLE EGGFILE
  ERROR_VARIABLE EGGFILE_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE)
if(EGGFILE_RESULT EQUAL 0)
  add_custom_command(OUTPUT ${EGGFILE}
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/build_egg.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS _orang ${SWIG_WRAPPER_PY} ${SWIG_WRAPPER_PYC}
    VERBATIM COMMENT "Building orang egg")
  add_custom_target(python-egg ALL DEPENDS ${EGGFILE})
else()
  message(SEND_ERROR "Error determining egg file name: ${EGGFILE_ERROR}")
endif()
