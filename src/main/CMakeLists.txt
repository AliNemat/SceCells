set(BEAK_PROJECT_NAME runBeakSimulation)
MESSAGE("Now building " ${BEAK_PROJECT_NAME})
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} ${GPU_SOURCE_DIR} ${CPU_SOURCE_DIR})

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
        include(${CGAL_USE_FILE})
endif (CGAL_FOUND)

file(GLOB SRC_FILES ./beakMain.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${BEAK_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${BEAK_PROJECT_NAME} googletest)
add_dependencies(${BEAK_PROJECT_NAME} SceGPU)
#add_dependencies(${BEAK_PROJECT_NAME} SceUtilsGPU)

if(NOT WIN32)
        target_link_libraries(${BEAK_PROJECT_NAME} pthread)
        target_link_libraries(${BEAK_PROJECT_NAME} SceGPU)
        target_link_libraries(${BEAK_PROJECT_NAME} SceInit)
        target_link_libraries(${BEAK_PROJECT_NAME} SceUtils)
        #target_link_libraries(${BEAK_PROJECT_NAME} SceUtilsGPU)
        target_link_libraries(${BEAK_PROJECT_NAME}
                ${GTEST_LIBS_DIR}/libgtest.a
                ${GTEST_LIBS_DIR}/libgtest_main.a
        )
endif()

set(DISK_PROJECT_NAME runDiskSimulation)
MESSAGE("Now building " ${DISK_PROJECT_NAME})
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})

file(GLOB SRC_FILES ./diskMain.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${DISK_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${DISK_PROJECT_NAME} googletest)
add_dependencies(${DISK_PROJECT_NAME} SceGPU)


if(NOT WIN32)
        target_link_libraries(${DISK_PROJECT_NAME} pthread)
        target_link_libraries(${DISK_PROJECT_NAME} SceGPU)
        target_link_libraries(${DISK_PROJECT_NAME} SceInit)
        target_link_libraries(${DISK_PROJECT_NAME} SceUtils)
        #target_link_libraries(${DISK_PROJECT_NAME} SceUtilsGPU)

        target_link_libraries(${DISK_PROJECT_NAME}
                ${GTEST_LIBS_DIR}/libgtest.a
                ${GTEST_LIBS_DIR}/libgtest_main.a
        )
        target_link_libraries(${DISK_PROJECT_NAME}
                        ${CGAL_LIBRARY}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

set(MODEL_TEST_NAME runModelTest)
MESSAGE("Now building " ${MODEL_TEST_NAME})
include_directories(${GTEST_INCLUDE_DIRS} ${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})

file(GLOB SRC_FILES ./modelTest.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${MODEL_TEST_NAME} ${SRC_FILES})
add_dependencies(${MODEL_TEST_NAME} googletest)
add_dependencies(${MODEL_TEST_NAME} SceGPU)

if(NOT WIN32)
        target_link_libraries(${MODEL_TEST_NAME} pthread)
        target_link_libraries(${MODEL_TEST_NAME} SceGPU)
        target_link_libraries(${MODEL_TEST_NAME} SceInit)
        target_link_libraries(${MODEL_TEST_NAME} SceUtils)
        #target_link_libraries(${MODEL_TEST_NAME} SceUtilsGPU)

        target_link_libraries(${MODEL_TEST_NAME}
                ${GTEST_LIBS_DIR}/libgtest.a
                ${GTEST_LIBS_DIR}/libgtest_main.a
        )
        target_link_libraries(${MODEL_TEST_NAME}
                        ${CGAL_LIBRARY}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
endif()


