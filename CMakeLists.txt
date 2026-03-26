cmake_minimum_required(VERSION 3.12)
project(libsoda-cxx VERSION 0.1.2)

# compiler
set(CMAKE_CXX_STANDARD 17)


#
# Setup clang-tidy for this project
#
find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
if(CLANG_TIDY_EXE)
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}" "--checks=-*,cppcoreguidelines-*")
    set(CMAKE_CXX_CLANG_TIDY_EXPORT_FIXES_DIR "${CMAKE_BINARY_DIR}/__clang-tidy-fixes")
else()
    message(WARNING "clang-tidy not found, continuing without linting.")
endif()

if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra)
endif()

# Configuration.
add_subdirectory(src)

#
# Tests
#
find_package(Threads REQUIRED)

add_executable(test_lsoda ${CMAKE_SOURCE_DIR}/tests/test_LSODA.cpp)
target_link_libraries(test_lsoda lsoda)

add_executable(benchmark_lsoda ${CMAKE_SOURCE_DIR}/tests/benchmark_LSODA.cpp)
target_link_libraries(benchmark_lsoda lsoda)
target_link_libraries(benchmark_lsoda ${CMAKE_THREAD_LIBS_INIT})

# Tests and benchmark
enable_testing()
add_test(NAME test_lsoda
    COMMAND $<TARGET_FILE:test_lsoda>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})

add_test(NAME test_benchmark
    COMMAND $<TARGET_FILE:benchmark_lsoda>
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
