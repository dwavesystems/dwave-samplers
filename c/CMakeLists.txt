find_package(Boost 1.37 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)

add_library(orang SHARED interface.cc)

add_executable(orang_call example.cc)

target_link_libraries(orang_call orang)