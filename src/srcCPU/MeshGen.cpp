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
#include <list>
#include "commonData.h"

// C++ sucks; declaring these static class members here is so counter-intuitive.
//std::list<Point> initList = std::list<Point>();
//initList.push_back(Point(9999999, 0));
std::list<Point> GEOMETRY::MeshGen::default_list_of_seeds = std::list<Point>();
//GEOMETRY::MeshGen::default_list_of_seeds.
//.push_back(Point(9999999, 0));
Criteria GEOMETRY::MeshGen::default_criteria = Criteria(0.125, 0.5);

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
	if (default_list_of_seeds.size() == 0) {
		default_list_of_seeds.push_back(Point(9999999, 0));
	}

	delta1 = globalConfigVars.getConfigValue("MeshGen_Delta1").toDouble();
	delta2 = globalConfigVars.getConfigValue("MeshGen_Delta2").toDouble();

	std::string bdryInputFileName = globalConfigVars.getConfigValue(
			"Bdry_InputFileName").toString();

	meshInput = GEOMETRY::MeshInputReader::readFile(bdryInputFileName);
	//default_list_of_seeds.push_back(Point(9999999, 0));
	//default_criteria = Criteria(0.125, 0.5);
}

UnstructMesh2D MeshGen::generateMesh2D(std::vector<Point2D> &boundaryPoints,
		std::list<Point> list_of_seeds, Criteria criteria) {
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

	//std::list<Point> list_of_seeds;
	//list_of_seeds.push_back(Point(9999999, 0));
	CGAL::refine_Delaunay_mesh_2(cdt, list_of_seeds.begin(),
			list_of_seeds.end(), criteria);
	//CGAL::refine_Delaunay_mesh_2(cdt, criteria);
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
	for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
			eit != cdt.finite_edges_end(); ++eit) {

		const CDT::Face_handle& fh = eit->first;

		int ctr = 0;
		if (fh->is_in_domain()) {
			ctr++;
		}
		if (fh->neighbor(eit->second)->is_in_domain()) {
			ctr++;
		}

		// this means boundary edges
		if (ctr == 1) {
			int i = eit->second;
			CDT::Vertex_handle vs = fh->vertex(fh->cw(i));
			CDT::Vertex_handle vt = fh->vertex(fh->ccw(i));
			VertexHashMap::iterator it = vertexHandleToIndex.find(vs);
			result.setPointAsBdry(it->second);
			it = vertexHandleToIndex.find(vt);
			result.setPointAsBdry(it->second);
		}
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

UnstructMesh2D MeshGen::generateMeshGivenInput(MeshInput input) {
	Criteria criteria(input.criteria_aspect_bound, input.criteria_size_bound);

	UnstructMesh2D result;
	CDT cdt;

	uint ptCount = 0;
	for (uint i = 0; i < input.bdryPts.size(); i++) {
		ptCount += input.bdryPts[i].size();
	}

	std::cout.flush();
	std::vector<Vertex_handle> vertexHandles(ptCount);
	uint vertexIndex = 0;
	uint beginIndexOfLevel = 0;

	for (uint i = 0; i < input.bdryPts.size(); i++) {

		for (uint j = 0; j < input.bdryPts[i].size(); j++) {
			vertexHandles[vertexIndex] = cdt.insert(
					Point(input.bdryPts[i][j].GetX(),
							input.bdryPts[i][j].GetY()));
			vertexIndex++;
		}

		for (uint j = 0; j < input.bdryPts[i].size(); j++) {
			if (j != input.bdryPts[i].size() - 1) {
				cdt.insert_constraint(vertexHandles[beginIndexOfLevel + j],
						vertexHandles[beginIndexOfLevel + j + 1]);
			} else {
				cdt.insert_constraint(vertexHandles[beginIndexOfLevel + j],
						vertexHandles[beginIndexOfLevel]);
			}
		}
		beginIndexOfLevel += input.bdryPts[i].size();
	}

	std::list<Point> seedList;
	for (uint i = 0; i < input.seeds.size(); i++) {
		seedList.push_back(Point(input.seeds[i].GetX(), input.seeds[i].GetY()));
	}

	CGAL::refine_Delaunay_mesh_2(cdt, seedList.begin(), seedList.end(),
			criteria);

	VertexHashMap vertexHandleToIndex;
	EdgeHashSet edgeSet;
	vertexIndex = 0;
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
	for (CDT::Finite_edges_iterator eit = cdt.finite_edges_begin();
			eit != cdt.finite_edges_end(); ++eit) {

		const CDT::Face_handle& fh = eit->first;

		int ctr = 0;
		if (fh->is_in_domain()) {
			ctr++;
		}
		if (fh->neighbor(eit->second)->is_in_domain()) {
			ctr++;
		}

		// this means boundary edges
		if (ctr == 1) {
			int i = eit->second;
			CDT::Vertex_handle vs = fh->vertex(fh->cw(i));
			CDT::Vertex_handle vt = fh->vertex(fh->ccw(i));
			VertexHashMap::iterator it = vertexHandleToIndex.find(vs);
			result.setPointAsBdry(it->second);
			it = vertexHandleToIndex.find(vt);
			result.setPointAsBdry(it->second);
		}
	}

	std::vector<GEOMETRY::Point2D> orderedBdryPts = obtainOrderedBdryPoints(
			result, input);
	result.setOrderedBdryPts(orderedBdryPts);

	return result;
}

UnstructMesh2D MeshGen::generateMesh2DFromFile(std::string& fileName) {
	//MeshInput input = GEOMETRY::MeshInputReader::readFile(fileName);
	return generateMeshGivenInput(meshInput);
}

UnstructMesh2D MeshGen::generateMesh2DFromFile(std::string &fileName,
		double ratio) {
	MeshInput input = meshInput;
	input.criteria_size_bound = input.criteria_size_bound / ratio;
	return generateMeshGivenInput(input);
}

std::vector<GEOMETRY::Point2D> MeshGen::obtainOrderedBdryPoints(
		UnstructMesh2D& mesh, MeshInput& input) {
	std::vector<GEOMETRY::Point2D> result;

	for (uint i = 0; i < input.bdryPts.size(); i++) {
		uint cirSize = input.bdryPts[i].size();
		for (uint j = 0; j < cirSize; j++) {
			CVector pt1 = input.bdryPts[i][j];
			// % for the last point, the neighbor point should be the first point.
			CVector pt2 = input.bdryPts[i][(j + 1) % cirSize];
			// orderPointsOnLine method should not only find points on the line but also
			std::vector<GEOMETRY::Point2D> orderedPtOnLine = orderPointsOnLine(
					mesh, pt1, pt2);
			// not inserting the entire vector to avoid duplication.
			result.insert(result.end(), orderedPtOnLine.begin(),
					orderedPtOnLine.end() - 1);
		}
	}
	return result;
}

CartilageRawData MeshGen::obtainCartilageData(UnstructMesh2D& mesh,
		MeshInput& input) {
	CartilageRawData result;
	/**
	 * this method is valid only for two-loop boundary.
	 */
	assert(input.bdryPts.size() == 2);
	/**
	 * this method is valid for four-point second boundary only.
	 */
	assert(input.bdryPts[1].size() == 4);
	/**
	 * again, this method is not robust enough. it assumes several
	 */
	assert(input.bdryPts[1][1].x < input.bdryPts[1][0].x);
	assert(input.bdryPts[1][0].x < input.bdryPts[1][2].x);
	assert(input.bdryPts[1][2].x < input.bdryPts[1][3].x);

	CVector pt1_1 = input.bdryPts[1][3];
	CVector pt1_2 = input.bdryPts[1][0];
	std::vector<CVector> orderedPtOnLine1 = orderPointsOnLine_vec3(mesh, pt1_1,
			pt1_2);
	CVector pt2_1 = input.bdryPts[1][0];
	CVector pt2_2 = input.bdryPts[1][1];
	std::vector<CVector> orderedPtOnLine2 = orderPointsOnLine_vec3(mesh, pt2_1,
			pt2_2);
	CVector pt3_1 = input.bdryPts[1][1];
	CVector pt3_2 = input.bdryPts[1][2];
	std::vector<CVector> orderedPtOnLine3 = orderPointsOnLine_vec3(mesh, pt3_1,
			pt3_2);

	// the +1 and -1 index are intentional.
	result.nonTipVerticies.insert(result.nonTipVerticies.end(),
			orderedPtOnLine1.begin() + 1, orderedPtOnLine1.end() - 1);
	// the -1 index is intentional.
	result.nonTipVerticies.insert(result.nonTipVerticies.end(),
			orderedPtOnLine2.begin(), orderedPtOnLine2.end() - 1);
	// the -1 index is intentional.
	result.nonTipVerticies.insert(result.nonTipVerticies.end(),
			orderedPtOnLine3.begin(), orderedPtOnLine3.end() - 1);

	CVector pt4_1 = input.bdryPts[1][2];
	CVector pt4_2 = input.bdryPts[1][3];
	std::vector<CVector> orderedPtOnLine4 = orderPointsOnLine_vec3(mesh, pt4_1,
			pt4_2);

	result.tipVerticies.insert(result.tipVerticies.end(),
			orderedPtOnLine4.begin(), orderedPtOnLine4.end());

	result.growNode1Index_on_tip = findClosestArrIndexGivenPos(
			result.tipVerticies, pt4_1);
	assert(result.growNode1Index_on_tip == 0);
	result.growNode2Index_on_tip = findClosestArrIndexGivenPos(
			result.tipVerticies, pt4_2);

	result.growNodeBehind1Index = findClosestArrIndexGivenPos(
			result.nonTipVerticies, pt4_1);
	result.growNodeBehind2Index = findClosestArrIndexGivenPos(
			result.nonTipVerticies, pt4_2);

	// pt2_2 is the left node of side 2
	result.pivotNode1Index = findClosestArrIndexGivenPos(result.nonTipVerticies,
			pt2_2);
	result.pivotNode2Index = findClosestArrIndexGivenPos(result.nonTipVerticies,
			pt2_1);

	return result;
}

std::vector<GEOMETRY::Point2D> MeshGen::orderPointsOnLine(UnstructMesh2D &mesh,
		CVector pt1, CVector pt2) {
	//std::vector<GEOMETRY::Point2D> result;
	std::vector<GEOMETRY::Point2D> bdryPts = mesh.getAllBdryPoints();
	std::vector<GEOMETRY::Point2D> pointsOnLine = getPointsOnLine(bdryPts, pt1,
			pt2);
	//CVector pos(pt1.GetX(), pt1.GetY(), 0);
	sort(pointsOnLine.begin(), pointsOnLine.end(), DistToPtComp(pt1));
	return pointsOnLine;
}

std::vector<CVector> MeshGen::orderPointsOnLine_vec3(UnstructMesh2D &mesh,
		CVector pt1, CVector pt2) {
	std::vector<CVector> result;
	std::vector<GEOMETRY::Point2D> vec2 = orderPointsOnLine(mesh, pt1, pt2);
	for (uint i = 0; i < vec2.size(); i++) {
		result.push_back(CVector(vec2[i].getX(), vec2[i].getY(), 0));
	}
	return result;
}

std::vector<GEOMETRY::Point2D> MeshGen::getPointsOnLine(
		std::vector<GEOMETRY::Point2D> bdryPts, CVector pt1, CVector pt2) {
	std::vector<GEOMETRY::Point2D> result;
	for (uint i = 0; i < bdryPts.size(); i++) {
		if (isBetweenPoints(bdryPts[i], pt1, pt2)) {
			result.push_back(bdryPts[i]);
		}
	}
	//cout << "inside getPointsOnline, result size = " << result.size() << endl;
	return result;
}

bool MeshGen::isBetweenPoints(GEOMETRY::Point2D point, CVector pt1,
		CVector pt2) {
	bool isParallel, isBounded;
	CVector pt(point.getX(), point.getY(), 0);
	CVector vec1 = pt - pt1;
	CVector vec2 = pt1 - pt2;
	double crossProduct = vec1.GetX() * vec2.GetY() - vec1.GetY() * vec2.GetX();
	if (fabs(crossProduct) < delta1) {
		isParallel = true;
	} else {
		isParallel = false;
	}
	double minX = std::min(pt1.GetX(), pt2.GetX()) - delta2;
	double maxX = std::max(pt1.GetX(), pt2.GetX()) + delta2;
	double minY = std::min(pt1.GetY(), pt2.GetY()) - delta2;
	double maxY = std::max(pt1.GetY(), pt2.GetY()) + delta2;
	if (pt.GetX() >= minX && pt.GetX() <= maxX && pt.GetY() >= minY
			&& pt.GetY() <= maxY) {
		isBounded = true;
	} else {
		isBounded = false;
	}
	if (isParallel && isBounded) {
		return true;
	} else {
		return false;
	}
}

MeshInput MeshGen::obtainMeshInput() {
	return meshInput;
}

MeshGen::~MeshGen() {
// TODO Auto-generated destructor stub
}

} /* namespace GEOMETRY */

std::vector<CVector> GEOMETRY::MeshInputReader::readPointVec(fstream& fs) {
	std::vector<CVector> result;
	char specialChar;
	fs >> specialChar;
	assert(specialChar == '{');
	fs >> specialChar;
	assert(specialChar == '$');
	int numPoints;
	fs >> numPoints;
	assert(numPoints > 0);
	for (int i = 0; i < numPoints; i++) {
		CVector point = readPoint(fs);
		result.push_back(point);
	}
	fs >> specialChar;
	assert(specialChar == '}');
	return result;
}

CVector GEOMETRY::MeshInputReader::readPoint(fstream& fs) {
	char specialChar;
	fs >> specialChar;
	assert(specialChar == '(');
	double x, y, z;
	fs >> x >> y >> z;
	fs >> specialChar;
	assert(specialChar == ')');
	CVector result(x, y, z);
	return result;
}

GEOMETRY::MeshInput GEOMETRY::MeshInputReader::readFile(std::string& fileName) {
	GEOMETRY::MeshInput meshInput;
	fstream fs(fileName.c_str());
	char specialChar;
	fs >> specialChar;
	assert(specialChar == '#');

	int count;
	fs >> count;

	assert(count > 0);
	for (int i = 0; i < count; i++) {
		std::vector<CVector> vec = readPointVec(fs);
		meshInput.bdryPts.push_back(vec);
	}

	fs >> specialChar;
	assert(specialChar == '#');
	fs >> count;
	assert(count == 1);
	std::vector<CVector> vec = readPointVec(fs);
	meshInput.seeds = vec;

	fs >> specialChar;
	assert(specialChar == '#');
	fs >> count;
	assert(count == 1);
	std::vector<CVector> vec2 = readPointVec(fs);
	meshInput.internalBdryPts = vec2;

	fs >> specialChar;
	assert(specialChar == '#');
	fs >> count;
	assert(count == 1);
	std::vector<CVector> endPts = readPointVec(fs);
	assert(endPts.size() == 2);
	meshInput.profileBeginPos = endPts[0];
	meshInput.profileEndPos = endPts[1];

	fs >> specialChar;
	assert(specialChar == '[');
	double aspectRatio, size;
	fs >> aspectRatio >> size;
	fs >> specialChar;
	assert(specialChar == ']');

	meshInput.criteria_aspect_bound = aspectRatio;
	meshInput.criteria_size_bound = size;

	return meshInput;
}

GEOMETRY::UnstructMesh2D GEOMETRY::MeshGen::generateMesh2DWithProfile(
		std::string& fileName, double ratio, bool isInnerBdryIncluded) {
	MeshInput input = meshInput;
	if (!isInnerBdryIncluded) {
		input.bdryPts.erase(input.bdryPts.begin() + 1);
		input.seeds.erase(input.seeds.begin());
	}
	input.criteria_size_bound = input.criteria_size_bound / ratio;
	GEOMETRY::UnstructMesh2D mesh = generateMeshGivenInput(input);
	mesh.generateFinalBdryAndProfilePoints(input.profileBeginPos,
			input.profileEndPos);
	return mesh;
}

