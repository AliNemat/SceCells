//============================================================================
// Name        : MeshGen.cpp
// Author      : Wenzhao Sun
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include "MeshGen.h"

using namespace std;

int main() {
	GEOMETRY::MeshGen meshGen;
	//std::string meshBdryPointsFile("boundaryPoints.dat");
	//std::vector<GEOMETRY::Point2D> points = meshGen.readBdryPointsFromFile(
	//		meshBdryPointsFile);

	std::vector<GEOMETRY::Point2D> points = meshGen.createBdryPointsOnCircle(3,
			10);

	//points.push_back(GEOMETRY::Point2D(3, 3));
	//points.push_back(GEOMETRY::Point2D(-3, 3));
	//points.push_back(GEOMETRY::Point2D(-3, -3));
	//points.push_back(GEOMETRY::Point2D(3, -3));
	GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2D(points);
	std::string vtkFileName = "testVtk";
	mesh.outputVtkFile(vtkFileName);
	return 0;
}
