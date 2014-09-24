/*
 * MeshGen.h
 *
 *  Created on: Aug 20, 2014
 *      Author: wsun2
 */

#ifndef MESHGEN_H_
#define MESHGEN_H_

#include <iostream>

#include "UnstructMesh2D.h"
#include "GeoVector.h"
#include <assert.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Delaunay_mesh_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> Tds;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds> CDT;
typedef CGAL::Delaunay_mesh_size_criteria_2<CDT> Criteria;
typedef CDT::Vertex_handle Vertex_handle;
typedef CDT::Point Point;

namespace GEOMETRY {

class MeshInput {
public:
	std::vector<std::vector<CVector> > bdryPts;
	std::vector<CVector> seeds;

	double criteria_aspect_bound;
	double criteria_size_bound;
};

class MeshInputReader {
	static std::vector<CVector> readPointVec(fstream &fs);
	static CVector readPoint(fstream &fs);
public:
	static GEOMETRY::MeshInput readFile(std::string &fileName);
};

class MeshGen {
public:
	static std::list<Point> default_list_of_seeds;
	static Criteria default_criteria;
	MeshGen();
	std::vector<Point2D> readBdryPointsFromFile(std::string &fileName);
	std::vector<Point2D> createBdryPointsOnCircle(double r, int n);
	UnstructMesh2D generateMesh2D(std::vector<Point2D> &boundaryPoints,
			std::list<Point> list_of_seeds = default_list_of_seeds,
			Criteria criteria = default_criteria);
	UnstructMesh2D generateMesh2DFromFile(std::string &fileName);
	virtual ~MeshGen();
};

} /* namespace GEOMETRY */

#endif /* MESHGEN_H_ */
