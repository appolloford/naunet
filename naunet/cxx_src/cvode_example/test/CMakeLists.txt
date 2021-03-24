
list(APPEND NAUNETLIB
    naunet_ode
    naunet
)

add_executable(naunet.x main.cpp)

target_link_libraries(naunet.x ${NAUNETLIB})
file(COPY timeres.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

