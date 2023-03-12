cmake_minimum_required(VERSION 3.1)

#
# Main library

project(catsc)
set( catsc_VERSION_MAJOR 1 )
set( catsc_VERSION_MINOR 0 )
set( catsc_VERSION ${catsc_VERSION_MAJOR}.${catsc_VERSION_MINOR} )

#option(BUILD_STATIC_LIBS "Build the static library" ON)
#option(BUILD_SHARED_LIBS "Build the shared library" ON)
#option(BUILD_TESTS "Build test programs" OFF)

set(CATS_COORDINATE_TYPE float CACHE STRING "C type used to describe spatial coordinates")
set(CATS_HIT_ID_TYPE "unsigned int" CACHE STRING "C type used to identify hit on layer")

configure_file(include/catsc/config.h.in ${CMAKE_BINARY_DIR}/include/catsc/config.h)

set( catsc_LIB catsc )
set( catsc_SOURCES src/cats.c
                   src/cats.cc )

add_library(${catsc_LIB} SHARED ${catsc_SOURCES})
target_include_directories(${catsc_LIB} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include> )
set_target_properties( ${catsc_LIB} PROPERTIES VERSION ${catsc_VERSION}
                       SOVERSION ${catsc_VERSION_MAJOR} )

#include(GenerateExportHeader)
#generate_export_header(catsc)

#
# Functional tests (require ROOT)
find_package(ROOT COMPONENTS Core)
if(ROOT_FOUND)
    include( ${ROOT_USE_FILE} )

    add_executable( catsc-test test/main-gpl.cc )
    target_link_libraries( catsc-test ${catsc_LIB} ${ROOT_LIBRARIES} )

    add_executable( catsc-test-mult test/main-gpl-mult.cc )
    target_link_libraries( catsc-test-mult ${catsc_LIB} ${ROOT_LIBRARIES} )

    add_executable( catsc-test-hard test/main-hard.cc )
    target_link_libraries( catsc-test-hard ${catsc_LIB} ${ROOT_LIBRARIES} )

    #separate_arguments(ROOT_CXX_FLAGS_ UNIX_COMMAND ${ROOT_CXX_FLAGS})
    #list( REMOVE_ITEM ROOT_CXX_FLAGS_ -std=c++11 -std=c++14 )
    #target_compile_options(catsc-test-mult PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${ROOT_CXX_FLAGS_}>)
    #list(APPEND CMAKE_BUILD_RPATH "${ROOT_LIBRARY_DIR}")
endif(ROOT_FOUND)

install( TARGETS ${catsc_LIB}
         EXPORT catscTargets
         #FILE_SET HEADERS
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib
         RUNTIME DESTINATION bin
         INCLUDES DESTINATION include
         )
install( EXPORT catscTargets
         FILE catscTargets.cmake
         NAMESPACE catsc::
         DESTINATION lib/cmake/catsc
         )
install( FILES include/catsc/cats.h
               include/catsc/cats.hh
               ${CMAKE_CURRENT_BINARY_DIR}/include/catsc/config.h
               #${CMAKE_CURRENT_BINARY_DIR}/catsc_export.h
         DESTINATION include/catsc
         COMPONENT Devel
         )

include(cmake/catscCreateConfig.cmake)

# Uninstall target (TODO: leaves empty dirs)
add_custom_target(uninstall
    "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_SOURCE_DIR}/cmake/uninstall.cmake"
)
