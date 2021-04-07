list(APPEND NAUNETLIB
    naunet
    naunet_ode
)

add_executable(naunet.x main.cpp)

get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
if("CUDA" IN_LIST languages)
    target_include_directories(naunet.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

if(USE_BOOST)
    target_include_directories(naunet.x PRIVATE ${BOOST_ROOT}/include)
endif()

target_link_libraries(naunet.x ${NAUNETLIB})
file(COPY timeres.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

