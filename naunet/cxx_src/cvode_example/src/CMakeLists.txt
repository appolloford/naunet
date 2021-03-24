

add_library(naunet_ode STATIC naunet_ode.cpp)
target_include_directories(naunet_ode PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(naunet_ode SUNDIALS::cvode_static SUNDIALS::nvecserial_static)

add_library(naunet STATIC naunet.cpp)
target_include_directories(naunet PUBLIC ${CMAKE_SOURCE_DIR}/include)
# for klu.h
target_include_directories(naunet PUBLIC /usr/local/include)
target_link_libraries(naunet SUNDIALS::cvode_static SUNDIALS::nvecserial_static)
target_link_libraries(naunet SUNDIALS::sunlinsolklu_static SUNDIALS::sunmatrixsparse_static)
target_link_libraries(naunet m)


