add_library(lsoda LSODA.cpp)
add_library(lsoda::lsoda ALIAS lsoda)

target_include_directories(lsoda PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
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
