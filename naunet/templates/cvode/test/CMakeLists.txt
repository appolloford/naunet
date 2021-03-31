
list(APPEND NAUNETLIB
    naunet_ode
    naunet
)

add_executable(naunet.x main.cpp)

target_include_directories(naunet.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(naunet.x ${NAUNETLIB})
file(COPY timeres.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

