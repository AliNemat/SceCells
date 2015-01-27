set(BEAK_PROJECT_NAME runBeakSimulation)
MESSAGE("Now building " ${BEAK_PROJECT_NAME})
include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${GPU_SOURCE_DIR} ${CPU_SOURCE_DIR})

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
        include(${CGAL_USE_FILE})
endif (CGAL_FOUND)

file(GLOB SRC_FILES ./beakMain.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${BEAK_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${BEAK_PROJECT_NAME} gtest)
add_dependencies(${BEAK_PROJECT_NAME} gtest_main)
add_dependencies(${BEAK_PROJECT_NAME} SceGPU)
#add_dependencies(${BEAK_PROJECT_NAME} SceUtilsGPU)

if(NOT WIN32)
        #target_link_libraries(${BEAK_PROJECT_NAME} pthread)
        target_link_libraries(${BEAK_PROJECT_NAME} SceGPU)
        target_link_libraries(${BEAK_PROJECT_NAME} SceInit)
        target_link_libraries(${BEAK_PROJECT_NAME} SceUtils)
        #target_link_libraries(${BEAK_PROJECT_NAME} SceUtilsGPU)
        target_link_libraries(${BEAK_PROJECT_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )
endif()

set(DISK_PROJECT_NAME runDiscSimulation)
MESSAGE("Now building " ${DISK_PROJECT_NAME})
include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})

file(GLOB SRC_FILES ./discMain.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${DISK_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${DISK_PROJECT_NAME} gtest)
add_dependencies(${DISK_PROJECT_NAME} gtest_main)
add_dependencies(${DISK_PROJECT_NAME} SceGPU)

if(NOT WIN32)
        #target_link_libraries(${DISK_PROJECT_NAME} pthread)
        target_link_libraries(${DISK_PROJECT_NAME} SceGPU)
        target_link_libraries(${DISK_PROJECT_NAME} SceInit)
        target_link_libraries(${DISK_PROJECT_NAME} SceUtils)
        #target_link_libraries(${DISK_PROJECT_NAME} SceUtilsGPU)

        target_link_libraries(${DISK_PROJECT_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )
        target_link_libraries(${DISK_PROJECT_NAME}
                        ${CGAL_LIBRARY}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

set(MODEL_TEST_NAME runModelTest)
MESSAGE("Now building " ${MODEL_TEST_NAME})
include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})

file(GLOB SRC_FILES ./modelTest.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${MODEL_TEST_NAME} ${SRC_FILES})
add_dependencies(${MODEL_TEST_NAME} gtest)
add_dependencies(${MODEL_TEST_NAME} gtest_main)
add_dependencies(${MODEL_TEST_NAME} SceGPU)

if(NOT WIN32)
        #target_link_libraries(${MODEL_TEST_NAME} pthread)
        target_link_libraries(${MODEL_TEST_NAME} SceGPU)
        target_link_libraries(${MODEL_TEST_NAME} SceInit)
        target_link_libraries(${MODEL_TEST_NAME} SceUtils)
        #target_link_libraries(${MODEL_TEST_NAME} SceUtilsGPU)

        target_link_libraries(${MODEL_TEST_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                 ${GOOGLETEST_DIR}/libgtest_main.a
        )
        target_link_libraries(${MODEL_TEST_NAME}
                        ${CGAL_LIBRARY}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

set(LASER_PROJECT_NAME LaserAblationTest)
MESSAGE("Now building " ${LASER_PROJECT_NAME})
include_directories(${GOOGLETEST_DIR}/include/ ${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})

file(GLOB SRC_FILES ./laserAblationTest.cpp)
#MESSAGE("Compiling these files: "${SRC_FILES})
cuda_add_executable(${LASER_PROJECT_NAME} ${SRC_FILES})
add_dependencies(${LASER_PROJECT_NAME} gtest)
add_dependencies(${LASER_PROJECT_NAME} gtest_main)
add_dependencies(${LASER_PROJECT_NAME} SceGPU)

if(NOT WIN32)
        target_link_libraries(${LASER_PROJECT_NAME} SceGPU)
        target_link_libraries(${LASER_PROJECT_NAME} SceInit)
        target_link_libraries(${LASER_PROJECT_NAME} SceUtils)

        target_link_libraries(${LASER_PROJECT_NAME}
                ${GOOGLETEST_DIR}/libgtest.a
                ${GOOGLETEST_DIR}/libgtest_main.a
        )
        target_link_libraries(${LASER_PROJECT_NAME}
                        ${CGAL_LIBRARY}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
endif()

