list(APPEND NAUNETLIB
    naunet
    naunet_ode
)

# check enabled languages
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

add_executable(test_singlegrid.x singlegrid.cpp)

if("CUDA" IN_LIST languages)
    target_include_directories(test_singlegrid.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# for odeint test
if(USE_BOOST)
    target_include_directories(test_singlegrid.x PRIVATE ${BOOST_ROOT}/include)
endif()

target_link_libraries(test_singlegrid.x ${NAUNETLIB})

# if CUDA is on, add parallel test
if("CUDA" IN_LIST languages)
    add_executable(test_parallel.x parallel.cpp)
    target_include_directories(test_parallel.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(test_parallel.x ${NAUNETLIB})
endif()

file(COPY timeres.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

