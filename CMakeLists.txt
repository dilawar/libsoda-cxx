cmake_minimum_required(VERSION 3.21)
project(libsoda-cxx VERSION 0.1.2 LANGUAGES CXX)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

# Clang-tidy: only in standalone builds to avoid leaking into consumers
if(PROJECT_IS_TOP_LEVEL)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
    if(CLANG_TIDY_EXE)
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}" "--checks=-*,cppcoreguidelines-*")
        set(CMAKE_CXX_CLANG_TIDY_EXPORT_FIXES_DIR "${CMAKE_BINARY_DIR}/__clang-tidy-fixes")
    else()
        message(WARNING "clang-tidy not found, continuing without linting.")
    endif()
endif()

if(MSVC)
    add_compile_options(/W3)
else()
    add_compile_options(-Wall -Wextra)
endif()

add_subdirectory(src)

# --- Install / export ---

write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/lsodaConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/lsodaConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/lsodaConfig.cmake"
    INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lsoda"
)

install(EXPORT lsodaTargets
    FILE      lsodaTargets.cmake
    NAMESPACE lsoda::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lsoda"
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/lsodaConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/lsodaConfigVersion.cmake"
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/lsoda"
)

# --- Tests: only when building standalone ---
if(PROJECT_IS_TOP_LEVEL)
    find_package(Threads REQUIRED)

    add_executable(test_lsoda tests/test_LSODA.cpp)
    target_link_libraries(test_lsoda PRIVATE lsoda)

    add_executable(benchmark_lsoda tests/benchmark_LSODA.cpp)
    target_link_libraries(benchmark_lsoda PRIVATE lsoda ${CMAKE_THREAD_LIBS_INIT})

    enable_testing()
    add_test(NAME test_lsoda
        COMMAND $<TARGET_FILE:test_lsoda>
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
    add_test(NAME test_benchmark
        COMMAND $<TARGET_FILE:benchmark_lsoda>
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
endif()
