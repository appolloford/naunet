cmake_minimum_required(VERSION 3.10)

# set the project name
project(NAUNET VERSION 0.1.0)

find_library(CVODE sundials_cvode)
if (NOT CVODE)
    message(FATAL_ERROR "cvode not found")
endif()

find_library(NVECTOR_S sundials_nvecserial)
if (NOT NVECTOR_S)
    message(FATAL_ERROR "sundials_nvecserial not found")
endif()



add_executable(naunet main.cpp)

include_directories(/usr/local/include)
#link_directories(/usr/local/lib)
#target_link_libraries(naunet sundials_cvode sundials_nvecserial)

target_link_libraries(naunet ${CVODE} ${NVECTOR_S})

