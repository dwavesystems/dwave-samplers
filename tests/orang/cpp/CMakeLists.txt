find_package(Boost 1.42 COMPONENTS unit_test_framework REQUIRED)

add_executable(run_tests EXCLUDE_FROM_ALL 
               test-main.cpp
               test-table.cpp
               test-graph.cpp
               test-treedecomp.cpp
               test-task.cpp
               test-merger.cpp
               test-op-min.cpp
               test-op-count.cpp
               test-op-logsumprod.cpp
               test-buckettree.cpp
               test-varorder.cpp
               util.cpp)
include_directories(${CMAKE_SOURCE_DIR}/orang/src/include ${Boost_INCLUDE_DIRS})
add_definitions(-DBOOST_TEST_DYN_LINK -DBOOST_ALL_NO_LIB)
target_link_libraries(run_tests ${Boost_LIBRARIES})

option(TEST_XML_OUTPUT OFF)

if(TEST_XML_OUTPUT)
  set(TEST_ARGS --output_format=xml --report_level=no --log_level=all
                --log_sink=${CMAKE_CURRENT_BINARY_DIR}/test-results.xml)
else()
  set(TEST_ARGS)
endif()

add_custom_target(check COMMAND run_tests ${TEST_ARGS} DEPENDS run_tests)

