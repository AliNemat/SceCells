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
#include <assert.h>

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

	/**
	 * By default, boundary points are not generated in order.
	 * Therefore an external method was used to order boundary points.
	 */
	std::vector<GEOMETRY::Point2D> orderedBdryPts;

	std::vector<GEOMETRY::Point2D> finalBdryPts;
	std::vector<GEOMETRY::Point2D> finalProfilePts;

	/**
	 * profile nodes can be fetched from part of the bdry node array.
	 */
	uint beginIndexProfile, endIndexProfile;
	uint findIndexGivenPos(CVector pos);
public:
	void insertVertex(const GEOMETRY::Point2D& point2D);
	void insertTriangle(const std::vector<int>& triangle);
	void insertEdge(const std::pair<int, int>& edge);
	void setPointAsBdry(int index);

	UnstructMesh2D();
	virtual ~UnstructMesh2D();
	void outputVtkFile(std::string outputFileName);
	std::vector<GEOMETRY::Point2D> outputTriangleCenters();
	std::vector<GEOMETRY::Point2D> outputTriangleVerticies();
	std::vector<GEOMETRY::Point2D> getAllInsidePoints();
	std::vector<GEOMETRY::Point2D> getAllBdryPoints();

	void generateFinalBdryAndProfilePoints(CVector posBegin, CVector posEnd);

	const std::vector<GEOMETRY::Point2D>& getOrderedBdryPts() const {
		return orderedBdryPts;
	}

	void setOrderedBdryPts(
			const std::vector<GEOMETRY::Point2D>& orderedBdryPts) {
		this->orderedBdryPts = orderedBdryPts;
	}

	const std::vector<GEOMETRY::Point2D>& getFinalBdryPts() const {
		return finalBdryPts;
	}

	const std::vector<GEOMETRY::Point2D>& getFinalProfilePts() const {
		return finalProfilePts;
	}
};

} /* namespace GEOMETRY */

#endif /* MESH_H_ */
