cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Release)
add_definitions( -std=c++11 ) 

add_library( lsoda ${CMAKE_CURRENT_SOURCE_DIR}/LSODA.cpp)

if( "${CMAKE_BUILD_TYPE}" STREQUAL "Release" )
    message( STATUS "Building for release" )
    set_target_properties( lsoda PROPERTIES COMPILE_FLAGS "-Werror -Wall")
else()
    message( STATUS "Building for DEBUG mode" )
endif()

find_package(Threads REQUIRED)

add_executable( test_lsoda ${CMAKE_CURRENT_SOURCE_DIR}/test_LSODA.cpp)
target_link_libraries( test_lsoda lsoda )

add_executable( benchmark_lsoda ${CMAKE_CURRENT_SOURCE_DIR}/benchmark_LSODA.cpp)
target_link_libraries( benchmark_lsoda lsoda )

set_target_properties(test_lsoda PROPERTIES COMPILE_FLAGS "-g -UNDEBUG" )
set_target_properties(benchmark_lsoda 
    PROPERTIES COMPILE_FLAGS "-O3" LINK_FLAGS "-pthread"
    ) 
target_link_libraries(benchmark_lsoda ${CMAKE_THREAD_LIBS_INIT})

enable_testing()
add_test(NAME test_lsoda COMMAND $<TARGET_FILE:test_lsoda>)
add_test(NAME test_benchmark COMMAND $<TARGET_FILE:benchmark_lsoda>)

