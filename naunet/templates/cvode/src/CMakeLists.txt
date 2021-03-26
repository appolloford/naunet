
list(APPEND SUNDIALSLIB
    SUNDIALS::cvode
    SUNDIALS::nvecserial
    SUNDIALS::sunlinsolklu
    SUNDIALS::sunmatrixsparse
)

add_library(naunet_ode STATIC naunet_ode.cpp)
target_include_directories(naunet_ode PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(naunet_ode ${SUNDIALSLIB})

add_library(naunet STATIC naunet.cpp)
target_include_directories(naunet PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(naunet ${SUNDIALSLIB})


