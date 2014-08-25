include_directories(${COMMON_INCLUDES} ${COMMON_CUDA_INCLUDES})
add_library( SceCPU
        CellInitHelper.h
        CellInitHelper.cpp
)
#target_link_libraries(SceCPU ./Util/libSceUtils.a)
#target_link_libraries(SceCPU SceUtils)
#target_link_libraries(SceCPU  
#                        ${CGAL_LIBRARIES}
#                        ${CGAL_3RD_PARTY_LIBRARIES} 
#)

add_subdirectory("./Util")
