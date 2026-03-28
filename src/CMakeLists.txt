add_library(lsoda LSODA.cpp)
add_library(lsoda::lsoda ALIAS lsoda)

# Mirror the installed layout (<prefix>/include/lsoda/LSODA.h) in the build
# tree so add_subdirectory/FetchContent consumers use the same #include spelling
# as find_package consumers: #include <lsoda/LSODA.h>
configure_file(LSODA.h  "${CMAKE_CURRENT_BINARY_DIR}/lsoda/LSODA.h"  COPYONLY)
configure_file(helper.h "${CMAKE_CURRENT_BINARY_DIR}/lsoda/helper.h" COPYONLY)

target_include_directories(lsoda PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_compile_features(lsoda PUBLIC cxx_std_17)

install(TARGETS lsoda
    EXPORT  lsodaTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(FILES LSODA.h helper.h
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/lsoda
)
