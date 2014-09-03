#MESSAGE(${GPU_SOURCE_DIR})
#include_directories(${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})

include_directories(${COMMON_INCLUDES})

CUDA_ADD_LIBRARY(SceGPU
	GrowthDistriMap.h
	GrowthDistriMap.cu
	SceNodes.h
	SceNodes.cu
	SceCells.h
	SceCells.cu
	SimulationDomainGPU.h
	SimulationDomainGPU.cu
        ${COMMON_INCLUDES}/commonData.h
	#CellInitHelper.cpp
        #CellInitHelper.h
)
#target_link_libraries(SceGPU ${CPU_SOURCE_DIR}/libSceCPU.a)
#target_link_libraries(SceGPU ${CPU_SOURCE_DIR}/Utils/libSceUtils.a)
