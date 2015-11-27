find_package(Boost 1.37 REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)

add_library(orang SHARED interface.cc mempool.cc)

add_executable(run_orang example.cc)

target_link_libraries(run_orang orang)