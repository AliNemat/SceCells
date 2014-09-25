/*
 * Mesh.h
 *
 *  Created on: Aug 20, 2014
 *      Author: wsun2
 */

#ifndef MESH_H_
#define MESH_H_

#include <vector>
#include <string>
#include <fstream>

#include "Point2D.h"
#include "commonData.h"

namespace GEOMETRY {

struct GeoException: public std::exception {
	std::string errorMessage;
	GeoException(std::string errMsg) :
			errorMessage(errMsg) {
	}
	~GeoException() throw () {
	}
	const char* what() const throw () {
		return errorMessage.c_str();
	}
};

/**
 * Help class for further processing generated mesh data.
 *
 */
class UnstructMesh2D {

	std::vector<std::vector<int> > triangles;
	std::vector<std::pair<int, int> > edges;
	std::vector<GEOMETRY::Point2D> points;

public:
	void insertVertex(const GEOMETRY::Point2D& point2D);
	void insertTriangle(const std::vector<int>& triangle);
	void insertEdge(const std::pair<int, int>& edge);
	void setPointAsBdry(int index);

	UnstructMesh2D();
	void outputVtkFile(std::string outputFileName);
	std::vector<GEOMETRY::Point2D> outputTriangleCenters();
	std::vector<GEOMETRY::Point2D> outputTriangleVerticies();
	std::vector<GEOMETRY::Point2D> getAllInsidePoints();
	std::vector<GEOMETRY::Point2D> getAllBdryPoints();
	virtual ~UnstructMesh2D();
};

} /* namespace GEOMETRY */

#endif /* MESH_H_ */
