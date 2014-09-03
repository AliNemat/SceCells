find_package(CGAL REQUIRED)
if(CGAL_FOUND)
        include(${CGAL_USE_FILE})
endif (CGAL_FOUND)
#MESSAGE(${CGAL_INCLUDE_DIRS})
#include_directories(${CGAL_INCLUDE_DIRS})

include_directories(${COMMON_INCLUDES})

add_library(SceInit
        MeshGen.cpp
        MeshGen.h
        UnstructMesh2D.cpp
        UnstructMesh2D.h
        CellInitHelper.h
        CellInitHelper.cpp
        ${COMMON_INCLUDES}/Point2D.cpp
        ${COMMON_INCLUDES}/Point2D.h
        ${COMMON_INCLUDES}/ConfigParser.h
        ${COMMON_INCLUDES}/ConfigParser.cpp
        ${COMMON_INCLUDES}/GeoVector.h
        ${COMMON_INCLUDES}/GeoVector.cpp
)

target_link_libraries(SceInit 
                        ${CGAL_LIBRARIES}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
