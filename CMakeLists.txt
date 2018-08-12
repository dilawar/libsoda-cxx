cmake_minimum_required(VERSION 2.8)
add_definitions( -std=c++11 )
add_library( lsoda ${CMAKE_CURRENT_SOURCE_DIR}/LSODA.cpp)

if( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )
    set_target_properties( ksolve lsoda PROPERTIES COMPILE_FLAGS "-Werror -Wall")
endif()

add_executable( test_lsoda ${CMAKE_CURRENT_SOURCE_DIR}/test_LSODA.cpp)
target_link_libraries( test_lsoda lsoda )

add_executable( benchmark_lsoda ${CMAKE_CURRENT_SOURCE_DIR}/benchmark_LSODA.cpp)
target_link_libraries( benchmark_lsoda lsoda )

set_target_properties(test_lsoda PROPERTIES COMPILE_FLAGS "-g -UNDEBUG" )
set_target_properties(benchmark_lsoda PROPERTIES COMPILE_FLAGS "-O3" ) 

enable_testing()
add_test(NAME test_lsoda COMMAND $<TARGET_FILE:test_lsoda>)
