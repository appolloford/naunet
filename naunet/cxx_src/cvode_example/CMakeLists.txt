cmake_minimum_required(VERSION 3.10)

# set the project name
project(NAUNET VERSION 0.1.0)

include(CMakePrintHelpers)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

IF (CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DDEBUG)
ENDIF()


set(SUNDIALS_MIN_VERSION "5.6.1")
find_package(sundials ${SUNDIALS_MIN_VERSION} REQUIRED CONFIG)

# IF (CMAKE_BUILD_TYPE STREQUAL Debug)
#     cmake_print_properties(TARGETS SUNDIALS::cvode_static PROPERTIES
#                            LOCATION INTERFACE_INCLUDE_DIRECTORIES)
#     cmake_print_properties(TARGETS SUNDIALS::nvecserial_static PROPERTIES
#                            LOCATION INTERFACE_INCLUDE_DIRECTORIES)
# ENDIF()


add_subdirectory(src)
add_subdirectory(test)

