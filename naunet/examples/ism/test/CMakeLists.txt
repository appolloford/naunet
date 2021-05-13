if(MAKE_SHARED)
    set(NAUNETLIB ${NAUNETLIB_SHARED})
elseif(MAKE_STATIC)
    set(NAUNETLIB ${NAUNETLIB_STATIC})
endif()

# check enabled languages
get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)

add_executable(test_singlegrid.x singlegrid.cpp)

if("CUDA" IN_LIST languages)
    target_include_directories(test_singlegrid.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
endif()

# for odeint test
if(USE_BOOST)
    target_include_directories(test_singlegrid.x PRIVATE ${Boost_INCLUDE_DIRS})
endif()

target_link_libraries(test_singlegrid.x ${NAUNETLIB})

if("CUDA" IN_LIST languages)
    # if CUDA is on, parallel test should be included
    add_executable(test_parallel.x parallel.cpp)
    target_include_directories(test_parallel.x PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
    target_link_libraries(test_parallel.x ${NAUNETLIB})

    # test if gpu device is found when CUDA is enabled
    if(${NVIDIA_GPU_COUNT})
        add_test(
            NAME test_singlegrid
            COMMAND $<TARGET_FILE:test_singlegrid.x>
        )
        add_test(
            NAME test_parallel
            COMMAND $<TARGET_FILE:test_parallel.x>
        )
    endif()
else()
    # if no CUDA, only test the simple case
    add_test(
        NAME test_singlegrid
        COMMAND $<TARGET_FILE:test_singlegrid.x>
    )
endif()

