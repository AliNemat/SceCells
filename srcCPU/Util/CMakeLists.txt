find_package(CGAL REQUIRED)
if(CGAL_FOUND)
        include_directories(${CGAL_USE_FILE})
endif (CGAL_FOUND)

add_library(SceUtils
        MeshGen.cpp
        MeshGen.h
        Point2D.cpp
        Point2D.h
        UnstructMesh2D.cpp
        UnstructMesh2D.h
        ConfigParser.h
        ConfigParser.cpp
        GeoVector.h
        GeoVector.cpp
)

target_link_libraries(SceUtils 
                        ${CGAL_LIBRARIES}
                        ${CGAL_3RD_PARTY_LIBRARIES} )
