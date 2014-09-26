/*
 * MeshGen.h
 *
 *  Created on: Aug 20, 2014
 *      Author: wsun2
 */

#ifndef MESHGEN_H_
#define MESHGEN_H_

#include <iostream>
#include <assert.h>
#include <algorithm>

#include "ConfigParser.h"
#include "UnstructMesh2D.h"
#include "GeoVector.h"

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

struct DistToPtComp {
	CVector _sourcePt;
	DistToPtComp(CVector pt) :
			_sourcePt(pt) {
	}
	bool operator()(const GEOMETRY::Point2D& pt1,
			const GEOMETRY::Point2D& pt2) {
		CVector vector1(pt1.getX() - _sourcePt.GetX(),
				pt1.getY() - _sourcePt.GetY(), 0);
		CVector vector2(pt2.getX() - _sourcePt.GetX(),
				pt2.getY() - _sourcePt.GetY(), 0);
		return (vector1.getModul() < vector2.getModul());
	}
};

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
	double delta1, delta2;
	UnstructMesh2D generateMeshGivenInput(MeshInput input);
	std::vector<GEOMETRY::Point2D> orderPointsOnLine(UnstructMesh2D mesh,
			CVector pt1, CVector pt2);
	std::vector<GEOMETRY::Point2D> getPointsOnLine(
			std::vector<GEOMETRY::Point2D> bdryPts, CVector pt1, CVector pt2);
	bool isBetweenPoints(GEOMETRY::Point2D point, CVector pt1, CVector pt2);
public:
	static std::list<Point> default_list_of_seeds;
	static Criteria default_criteria;
	MeshGen();
	std::vector<Point2D> createBdryPointsOnCircle(double r, int n);
	UnstructMesh2D generateMesh2D(std::vector<Point2D> &boundaryPoints,
			std::list<Point> list_of_seeds = default_list_of_seeds,
			Criteria criteria = default_criteria);

	std::vector<GEOMETRY::Point2D> obtainOrderedBdryPoints(UnstructMesh2D &mesh,
			MeshInput &input);

	UnstructMesh2D generateMesh2DFromFile(std::string &fileName);
	/**
	 * Generating finer mesh than the meshinputfile.
	 * Used to generate fine boundary nodes with everything else remain the same.
	 */
	UnstructMesh2D generateMesh2DFromFile(std::string &fileName, double ratio);
	virtual ~MeshGen();
};

} /* namespace GEOMETRY */

#endif /* MESHGEN_H_ */
