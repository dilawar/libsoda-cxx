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

    # Integration test: install into a temp prefix then build+run a standalone
    # consumer that uses find_package(lsoda) to verify the full install works.
    set(_integ_prefix "${CMAKE_CURRENT_BINARY_DIR}/integ_install")
    set(_integ_build  "${CMAKE_CURRENT_BINARY_DIR}/integ_build")

    add_test(NAME integration_install
        COMMAND ${CMAKE_COMMAND} --install ${CMAKE_BINARY_DIR}
                --prefix "${_integ_prefix}")
    add_test(NAME integration_configure
        COMMAND ${CMAKE_COMMAND}
                -S "${CMAKE_SOURCE_DIR}/tests/integration"
                -B "${_integ_build}"
                -DCMAKE_PREFIX_PATH=${_integ_prefix}
                -DCMAKE_BUILD_TYPE=Release)
    add_test(NAME integration_build
        COMMAND ${CMAKE_COMMAND} --build "${_integ_build}" --config Release)
    # Use ctest inside the consumer build directory so the executable path is
    # resolved correctly on all platforms (including MSVC multi-config generators).
    add_test(NAME integration_run
        COMMAND ${CMAKE_CTEST_COMMAND}
                --test-dir "${_integ_build}"
                --build-config Release
                --output-on-failure)

    set_tests_properties(integration_configure PROPERTIES DEPENDS integration_install)
    set_tests_properties(integration_build     PROPERTIES DEPENDS integration_configure)
    set_tests_properties(integration_run       PROPERTIES DEPENDS integration_build)
endif()
