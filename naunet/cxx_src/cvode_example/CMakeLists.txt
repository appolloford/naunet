cmake_minimum_required(VERSION 3.10)

include(CMakePrintHelpers)
# set the project name
project(NAUNET VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

IF (CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DDEBUG)
ENDIF()

SET(SUNDIALS_MIN_VERSION "5.6.1")

find_package(sundials ${SUNDIALS_MIN_VERSION} REQUIRED CONFIG)

IF (CMAKE_BUILD_TYPE STREQUAL Debug)
    cmake_print_properties(TARGETS SUNDIALS::cvode_static PROPERTIES
                        LOCATION INTERFACE_INCLUDE_DIRECTORIES)
    cmake_print_properties(TARGETS SUNDIALS::nvecserial_static PROPERTIES
                        LOCATION INTERFACE_INCLUDE_DIRECTORIES)
ENDIF()

add_library(naunet_ode STATIC src/naunet_ode.cpp)
target_include_directories(naunet_ode PUBLIC include)
target_link_libraries(naunet_ode SUNDIALS::cvode_static SUNDIALS::nvecserial_static)

add_library(naunet STATIC src/naunet.cpp)
target_include_directories(naunet PUBLIC include)
target_link_libraries(naunet SUNDIALS::cvode_static SUNDIALS::nvecserial_static)

add_executable(naunet.x src/main.cpp)

target_link_libraries(naunet.x naunet naunet_ode SUNDIALS::cvode_static SUNDIALS::nvecserial_static)

