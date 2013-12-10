cmake_minimum_required(VERSION 2.6)
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

