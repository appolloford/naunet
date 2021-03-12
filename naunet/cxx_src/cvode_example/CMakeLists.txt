cmake_minimum_required(VERSION 3.10)

# set the project name
project(NAUNET VERSION 0.1.0)

set (CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

IF (CMAKE_BUILD_TYPE STREQUAL Debug)
    ADD_DEFINITIONS(-DDEBUG)
ENDIF()

find_library(CVODE sundials_cvode)
if (NOT CVODE)
    message(FATAL_ERROR "cvode not found")
endif()

find_library(NVECTOR_S sundials_nvecserial)
if (NOT NVECTOR_S)
    message(FATAL_ERROR "sundials_nvecserial not found")
endif()

add_library(fex STATIC src/fex.cpp)
target_include_directories(fex PUBLIC include)
target_include_directories(fex PUBLIC /usr/local/include)

add_library(jtv STATIC src/jtv.cpp)
target_include_directories(jtv PUBLIC include)
target_include_directories(jtv PUBLIC /usr/local/include)

add_library(naunet STATIC src/naunet.cpp)
target_include_directories(naunet PUBLIC include)
target_include_directories(naunet PUBLIC /usr/local/include)

add_executable(naunet.x src/main.cpp)

target_include_directories(naunet.x PRIVATE include)
target_include_directories(naunet.x PRIVATE /usr/local/include)
#link_directories(/usr/local/lib)
#target_link_libraries(naunet sundials_cvode sundials_nvecserial)

target_link_libraries(naunet.x fex jtv naunet ${CVODE} ${NVECTOR_S})

