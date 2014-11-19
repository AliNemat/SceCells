ENABLE_TESTING()

set(PROJECT_TEST_NAME ${PROJECT_NAME_STR}_unit_test)
#MESSAGE(${PROJECT_TEST_NAME})
MESSAGE("Now building unit tests")
include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${CPU_SOURCE_DIR} ${GPU_SOURCE_DIR})

file(GLOB TEST_SRC_FILES ${CMAKE_UNITTEST_DIR}/*.cpp ${CMAKE_UNITTEST_DIR}/*.cu)
# DEBUGGING: 
#file(GLOB TEST_SRC_FILES ${CMAKE_UNITTEST_DIR}/*.cpp )
#MESSAGE(${TEST_SRC_FILES})
cuda_add_executable(${PROJECT_TEST_NAME} ${TEST_SRC_FILES})
add_dependencies(${PROJECT_TEST_NAME} googletest)
add_dependencies(${PROJECT_TEST_NAME} SceGPU)

if(NOT WIN32)
        target_link_libraries(${PROJECT_TEST_NAME} pthread)
        target_link_libraries(${PROJECT_TEST_NAME} SceGPU)
        target_link_libraries(${PROJECT_TEST_NAME} SceInit)
        target_link_libraries(${PROJECT_TEST_NAME} SceUtils)
        target_link_libraries(${PROJECT_TEST_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )
endif()

add_test(unitTest ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/${PROJECT_TEST_NAME})
