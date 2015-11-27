cmake_minimum_required(VERSION 2.6)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules)
project(orang CXX)

set(Boost_ADDITIONAL_VERSIONS "1.46.1")

add_subdirectory(src)

OPTION(ENABLE_TESTS "Build unit tests?" ON)
if(ENABLE_TESTS)
  add_subdirectory(test)
endif()

OPTION(ENABLE_MATLAB_INTERFACE "Build MATLAB interface?" OFF)
if(ENABLE_MATLAB_INTERFACE)
  add_subdirectory(matlab)
endif()

OPTION(ENABLE_PYTHON "Build Python interface?" OFF)
if(ENABLE_PYTHON)
  include(python)
  add_subdirectory(python)
endif()

OPTION(ENABLE_C "Build C interface?" OFF)
if(ENABLE_C)
  add_subdirectory(c)
endif()
