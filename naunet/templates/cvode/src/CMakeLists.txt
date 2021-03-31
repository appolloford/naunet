
list(APPEND SUNDIALSLIB
    SUNDIALS::cvode
    SUNDIALS::nvecserial
    SUNDIALS::sunlinsolklu
    SUNDIALS::sunmatrixsparse
    SUNDIALS::nveccuda
    SUNDIALS::sunmatrixcusparse
    SUNDIALS::sunlinsolcusolversp
)

add_library(naunet_ode STATIC naunet_ode.cu)
set_target_properties(naunet_ode PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(naunet_ode PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(naunet_ode PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(naunet_ode ${SUNDIALSLIB})

add_library(naunet STATIC naunet.cpp)
set_target_properties(naunet PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_include_directories(naunet PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_include_directories(naunet PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(naunet ${SUNDIALSLIB})


