include_directories(${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})
add_library( SceCPU
        CellInitHelper.h
        CellInitHelper.cpp
)

add_subdirectory("./Util")
