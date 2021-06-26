if(MAKE_SHARED)
    set(NAUNETLIB ${NAUNETLIB_SHARED})
elseif(MAKE_STATIC)
    set(NAUNETLIB ${NAUNETLIB_STATIC})
endif()

# check enabled languages
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

add_executable(test_singlegrid.x singlegrid.cpp)

add_executable(test_serialdata.x serialdata.cpp)

if("CUDA" IN_LIST languages)
    target_include_directories(test_singlegrid.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_include_directories(test_serialdata.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# for odeint test
if(USE_BOOST)
    target_include_directories(test_singlegrid.x PRIVATE ${Boost_INCLUDE_DIRS})
    target_include_directories(test_serialdata.x PRIVATE ${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(test_singlegrid.x ${NAUNETLIB})
target_link_libraries(test_serialdata.x ${NAUNETLIB})

if("CUDA" IN_LIST languages)
    # if CUDA is on, parallel test should be included
    add_executable(test_multiplegrid.x multiplegrid.cpp)
    target_include_directories(test_multiplegrid.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(test_multiplegrid.x ${NAUNETLIB})

    add_executable(test_paralleldata.x paralleldata.cpp)
    target_include_directories(test_paralleldata.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(test_paralleldata.x ${NAUNETLIB})

    # test if gpu device is found when CUDA is enabled
    if(${NVIDIA_GPU_COUNT})
        add_test(
            NAME test_singlegrid
            COMMAND $<TARGET_FILE:test_singlegrid.x>
        )
        add_test(
            NAME test_serialdata
            COMMAND $<TARGET_FILE:test_serialdata.x>
        )
        add_test(
            NAME test_multiplegrid
            COMMAND $<TARGET_FILE:test_multiplegrid.x>
        )
        add_test(
            NAME test_paralleldata
            COMMAND $<TARGET_FILE:test_paralleldata.x>
        )
    endif()
else()
    # if no CUDA, only test the simple case
    add_test(
        NAME test_singlegrid
        COMMAND $<TARGET_FILE:test_singlegrid.x>
    )
    add_test(
        NAME test_serialdata
        COMMAND $<TARGET_FILE:test_serialdata.x>
    )
endif()

if(MAKE_PYTHON)
    configure_file(test.py.in test.py)

    if( ${NVIDIA_GPU_COUNT} OR NOT ("CUDA" IN_LIST languages) )
        add_test(
            NAME
                test_pymodule
            COMMAND
                ${CMAKE_COMMAND} -E env 
                    CHEMISTRY_MODULE_PATH=$<TARGET_FILE_DIR:${PYMODNAME}>
                    NAUNET_CONFIG_PATH=${CMAKE_SOURCE_DIR}
                ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/test.py
        )
    endif()
endif()

file(COPY timeres.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY grids.dat DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

