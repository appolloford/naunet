
add_executable(naunet.x main.cpp)

# target_link_libraries(naunet.x naunet naunet_ode SUNDIALS::cvode_static SUNDIALS::nvecserial_static)
target_link_libraries(naunet.x naunet naunet_ode)
file(COPY timeres.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

