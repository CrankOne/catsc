cmake_minimum_required(VERSION 3.1)

#
# Main library

project(catsc)
set( CATSC_VERSION_MAJOR 1 )
set( CATSC_VERSION_MINOR 0 )
set( CATSC_VERSION ${CATSC_VERSION_MAJOR}.${CATSC_VERSION_MINOR} )

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
set_target_properties( ${catsc_LIB} PROPERTIES VERSION ${CATSC_VERSION}
                       SOVERSION ${CATSC_VERSION_MAJOR} )

install( TARGETS ${catsc_LIB} EXPORT ${catsc_LIB} LIBRARY
         PUBLIC_HEADER DESTINATION include/catsc )
install( FILES include/catsc/cats.h
               include/catsc/cats.hh
               ${CMAKE_CURRENT_BINARY_DIR}/include/catsc/config.h
         DESTINATION include/catsc )

#
# Functional tests (require ROOT)
find_package(ROOT COMPONENTS Core)
if(ROOT_FOUND)
    include( ${ROOT_USE_FILE} )

    add_executable( catsc-test main-gpl-test.cc )
    target_link_libraries( catsc-test ${catsc_LIB} ${ROOT_LIBRARIES} )

    add_executable( catsc-test-mult main-gpl-test-mult.cc )
    target_link_libraries( catsc-test-mult ${catsc_LIB} ${ROOT_LIBRARIES} )
endif(ROOT_FOUND)
