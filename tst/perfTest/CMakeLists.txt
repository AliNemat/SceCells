ENABLE_TESTING()

set(PROJECT_TEST_NAME PerfTest)
MESSAGE("Now building performance tests")
include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${CPU_SOURCE_DIR} ${GPU_SOURCE_DIR})

file(GLOB TEST_SRC_FILES ${CMAKE_PERFTEST_DIR}/*.cpp ${CMAKE_PERFTEST_DIR}/*.cu)
# DEBUGGING: 
#file(GLOB TEST_SRC_FILES ${CMAKE_UNITTEST_DIR}/*.cpp )
#MESSAGE(${TEST_SRC_FILES})
cuda_add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} gtest)
add_dependencies(${PROJECT_TEST_NAME} gtest_main)
add_dependencies(${PROJECT_TEST_NAME} SceGPU)

if(NOT WIN32)
        #target_link_libraries(${PROJECT_TEST_NAME} pthread)
        target_link_libraries(${PROJECT_TEST_NAME} SceGPU)
        target_link_libraries(${PROJECT_TEST_NAME} SceInit)
        target_link_libraries(${PROJECT_TEST_NAME} SceUtils)
        target_link_libraries(${PROJECT_TEST_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )
endif()

#add_test(performanceTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_TEST_NAME})
