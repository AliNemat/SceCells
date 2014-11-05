#MESSAGE(${GPU_SOURCE_DIR})
#include_directories(${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES} )

include_directories(${COMMON_INCLUDES} ${GTEST_INCLUDE_DIRS})

CUDA_ADD_LIBRARY(SceGPU
    Cartilage.h
    Cartilage.cu
	GrowthDistriMap.h
	GrowthDistriMap.cu
	SceNodes.h
	SceNodes.cu
	SceCells.h
	SceCells.cu
	SimulationDomainGPU.h
	SimulationDomainGPU.cu
    ${COMMON_INCLUDES}/commonData.h
    ${COMMON_INCLUDES}/commonGPUData.h
    #${COMMON_INCLUDES}/commonGPUData.cu
)
#target_link_libraries(SceGPU SceUtilsGPU)
#target_link_libraries(SceGPU ${COMMON_INCLUDES}/libSceUtilsGPU.a)
#target_link_libraries(SceGPU ${CPU_SOURCE_DIR}/libSceCPU.a)
#target_link_libraries(SceGPU ${CPU_SOURCE_DIR}/Utils/libSceUtils.a)
