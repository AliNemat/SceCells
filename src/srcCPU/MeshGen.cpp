/*
 * MeshGen.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: wsun2
 */

#include "MeshGen.h"
#include <tr1/unordered_map>
#include <tr1/unordered_set>
#include <tr1/functional>

struct HandleHashFunc {
	std::size_t operator()(const CDT::Vertex_handle& vHandle) const {
		double xCoord = vHandle->point().x();
		double yCoord = vHandle->point().y();
		std::tr1::hash<double> hashFunc;
		return hashFunc(xCoord) ^ hashFunc(yCoord);
	}
};

struct IntPairHashFunc {
	std::size_t operator()(const std::pair<int, int>& myPair) const {
		return (myPair.first + myPair.second + (myPair.first ^ myPair.second));
	}
};

typedef std::tr1::unordered_map<CDT::Vertex_handle, int, HandleHashFunc> VertexHashMap;
typedef std::tr1::unordered_set<std::pair<int, int>, IntPairHashFunc> EdgeHashSet;

namespace GEOMETRY {

MeshGen::MeshGen() {

}

UnstructMesh2D MeshGen::generateMesh2D(std::vector<Point2D> &boundaryPoints) {
	UnstructMesh2D result;
	CDT cdt;

	int pointCount = boundaryPoints.size();
	std::vector<Vertex_handle> vertexHandles(pointCount);
	for (int i = 0; i < pointCount; i++) {
		vertexHandles[i] = cdt.insert(
				Point(boundaryPoints[i].getX(), boundaryPoints[i].getY()));
	}

	for (int i = 0; i < pointCount; i++) {
		cdt.insert_constraint(vertexHandles[i],
				vertexHandles[(i + 1) % pointCount]);
	}

	std::cout << "Number of vertices: " << cdt.number_of_vertices()
			<< std::endl;
	std::cout << "Meshing the triangulation ..." << std::endl;

	std::list<Point> list_of_seeds;
	list_of_seeds.push_back(Point(9999999, 0));
	CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(),
			list_of_seeds.end(), Criteria(0.125, 0.5));
	std::cout << "Number of vertices: " << cdt.number_of_vertices()
			<< std::endl;
	VertexHashMap vertexHandleToIndex;
	EdgeHashSet edgeSet;
	int vertexIndex = 0;
	for (CDT::Finite_vertices_iterator vit = cdt.finite_vertices_begin();
			vit != cdt.finite_vertices_end(); ++vit) {
		CDT::Vertex_handle vertexHandle = vit->handle();
		VertexHashMap::iterator it = vertexHandleToIndex.find(vertexHandle);

		if (it == vertexHandleToIndex.end()) {
			vertexHandleToIndex.insert(
					std::pair<CDT::Vertex_handle, int>(vertexHandle,
							vertexIndex));
			vertexIndex++;
		} else {
			std::string errMessage(
					"Unexpected duplicate while iterate vertices");
			throw GeoException(errMessage);
		}
		Point2D pt(vit->point().x(), vit->point().y());
		result.insertVertex(pt);
	}

	for (CDT::Finite_faces_iterator fit = cdt.finite_faces_begin();
			fit != cdt.finite_faces_end(); ++fit) {
		CDT::Vertex_handle v1 = fit->vertex(0);
		CDT::Vertex_handle v2 = fit->vertex(1);
		CDT::Vertex_handle v3 = fit->vertex(2);
		VertexHashMap::iterator it1 = vertexHandleToIndex.find(v1);
		VertexHashMap::iterator it2 = vertexHandleToIndex.find(v2);
		VertexHashMap::iterator it3 = vertexHandleToIndex.find(v3);

		if (it1 == vertexHandleToIndex.end() || it2 == vertexHandleToIndex.end()
				|| it3 == vertexHandleToIndex.end()) {
			std::string errMessage(
					"at least one triangle vertex was not found in vertex pool");
			throw GeoException(errMessage);
		} else {
			if (fit->is_in_domain()) {
				std::vector<int> triIndicies;
				triIndicies.push_back(it1->second);
				triIndicies.push_back(it2->second);
				triIndicies.push_back(it3->second);
				result.insertTriangle(triIndicies);

				// sorting in order to make sure there is no duplicates
				// when counting edges
				std::sort(triIndicies.begin(), triIndicies.end());
				std::pair<int, int> edge1 = std::pair<int, int>(triIndicies[0],
						triIndicies[1]);
				std::pair<int, int> edge2 = std::pair<int, int>(triIndicies[0],
						triIndicies[2]);
				std::pair<int, int> edge3 = std::pair<int, int>(triIndicies[1],
						triIndicies[2]);
				if (edgeSet.find(edge1) == edgeSet.end()) {
					edgeSet.insert(edge1);
					result.insertEdge(edge1);
				}
				if (edgeSet.find(edge2) == edgeSet.end()) {
					edgeSet.insert(edge2);
					result.insertEdge(edge2);
				}
				if (edgeSet.find(edge3) == edgeSet.end()) {
					edgeSet.insert(edge3);
					result.insertEdge(edge3);
				}
			}
		}
	}
	return result;
}

std::vector<GEOMETRY::Point2D> MeshGen::readBdryPointsFromFile(
		std::string& fileName) {
	std::vector<GEOMETRY::Point2D> result;
	double x, y, z;
	std::ifstream infile(fileName.c_str());
	if (infile.is_open()) {
		while (!infile.eof()) {
			infile >> x >> y >> z;
			GEOMETRY::Point2D pt(x, y);
			result.push_back(pt);
		}
	} else {
		throw GEOMETRY::GeoException("Boundary points input file not found!");
	}
	return result;
}

std::vector<Point2D> MeshGen::createBdryPointsOnCircle(double r, int n) {
	std::vector<GEOMETRY::Point2D> result;
	const double PI = acos(-1.0);
	double unitAngle = 2 * PI / n;
	for (int i = 0; i < n; i++) {
		GEOMETRY::Point2D pt(r * sin(i * unitAngle), r * cos(i * unitAngle));
		result.push_back(pt);
	}
	return result;
}

MeshGen::~MeshGen() {
// TODO Auto-generated destructor stub
}

} /* namespace GEOMETRY */

