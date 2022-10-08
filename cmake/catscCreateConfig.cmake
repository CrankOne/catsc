# - pkg-config config
configure_file( cmake/catsc.pc.in ${CMAKE_CURRENT_BINARY_DIR}/catsc.pc )
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/catsc.pc
         DESTINATION lib/pkgconfig )

# - CMake Config
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "catscConfigVersion.cmake"
    VERSION ${catsc_VERSION}
    COMPATIBILITY SameMajorVersion)

install( FILES "cmake/catscConfig.cmake"
               "${CMAKE_CURRENT_BINARY_DIR}/catscConfigVersion.cmake"
         DESTINATION lib/cmake/catsc
         )
