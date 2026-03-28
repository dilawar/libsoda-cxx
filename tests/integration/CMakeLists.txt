cmake_minimum_required(VERSION 3.21)
project(lsoda_integration_test LANGUAGES CXX)

find_package(lsoda REQUIRED CONFIG)

add_executable(integration_test main.cpp)
target_link_libraries(integration_test PRIVATE lsoda::lsoda)

enable_testing()
add_test(NAME integration_test COMMAND integration_test)
