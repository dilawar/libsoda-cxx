cmake_minimum_required(VERSION 3.12)

project(libsoda-cxx VERSION 0.1.0)

# compiler
set(CMAKE_CXX_STANDARD 11)

if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra)
endif()

# Configuration.
find_package(Threads REQUIRED)

# targets
add_library( lsoda ${CMAKE_CURRENT_SOURCE_DIR}/LSODA.cpp)

add_executable(test_lsoda ${CMAKE_CURRENT_SOURCE_DIR}/test_LSODA.cpp)
target_link_libraries(test_lsoda lsoda)

add_executable(benchmark_lsoda ${CMAKE_CURRENT_SOURCE_DIR}/benchmark_LSODA.cpp)
target_link_libraries(benchmark_lsoda lsoda)
target_link_libraries(benchmark_lsoda ${CMAKE_THREAD_LIBS_INIT})

# Tests and benchmark
enable_testing()
add_test(NAME test_lsoda COMMAND $<TARGET_FILE:test_lsoda>)
add_test(NAME test_benchmark COMMAND $<TARGET_FILE:benchmark_lsoda>)
