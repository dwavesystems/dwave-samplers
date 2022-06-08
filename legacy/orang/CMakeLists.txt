cmake_minimum_required(VERSION 3.0)

project(orang CXX)


OPTION(ENABLE_TESTS "Build unit tests?" ON)
if(ENABLE_TESTS)
  add_subdirectory(tests/cpp)
endif()
