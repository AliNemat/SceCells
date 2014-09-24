//============================================================================
// Name        : MeshGen.cpp
// Author      : Wenzhao Sun
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include "MeshGen.h"
#include "commonData.h"
#include "CellInitHelper.h"
#include <vector>
#include "SimulationDomainGPU.h"

using namespace std;

GlobalConfigVars globalConfigVars;

void generateStringInputs(std::string &loadMeshInput,
		std::string &animationInput, std::string &animationFolder,
		std::vector<std::string> &boundaryMeshFileNames) {
	std::string meshLocation =
			globalConfigVars.getConfigValue("MeshLocation").toString();
	std::string meshName =
			globalConfigVars.getConfigValue("MeshName").toString();
	std::string meshExtention =
			globalConfigVars.getConfigValue("MeshExtention").toString();
	loadMeshInput = meshLocation + meshName + meshExtention;

	animationFolder =
			globalConfigVars.getConfigValue("AnimationFolder").toString();
	animationInput = animationFolder
			+ globalConfigVars.getConfigValue("AnimationName").toString();

	std::string boundaryMeshLocation = globalConfigVars.getConfigValue(
			"BoundaryMeshLocation").toString();
	std::string boundaryMeshName = globalConfigVars.getConfigValue(
			"BoundaryMeshName").toString();
	std::string boundaryMeshExtention = globalConfigVars.getConfigValue(
			"BoundaryMeshExtention").toString();
	std::string boundaryMeshInput = boundaryMeshLocation + boundaryMeshName
			+ boundaryMeshExtention;
	boundaryMeshFileNames.push_back(boundaryMeshInput);
}

int main() {
	srand(time(NULL));
	ConfigParser parser;
	std::string configFileName = "./resources/modelTest.cfg";
	globalConfigVars = parser.parseConfigFile(configFileName);

	/*
	 double SimulationTotalTime = globalConfigVars.getConfigValue(
	 "SimulationTotalTime").toDouble();
	 double SimulationTimeStep = globalConfigVars.getConfigValue(
	 "SimulationTimeStep").toDouble();
	 int TotalNumOfOutputFrames = globalConfigVars.getConfigValue(
	 "TotalNumOfOutputFrames").toInt();

	 std::string loadMeshInput;
	 std::string animationInput;
	 std::vector<std::string> boundaryMeshFileNames;
	 std::string animationFolder;
	 generateStringInputs(loadMeshInput, animationInput, animationFolder,
	 boundaryMeshFileNames);

	 const double simulationTime = SimulationTotalTime;
	 const double dt = SimulationTimeStep;
	 const int numOfTimeSteps = simulationTime / dt;
	 const int totalNumOfOutputFrame = TotalNumOfOutputFrames;
	 const int outputAnimationAuxVarible = numOfTimeSteps
	 / totalNumOfOutputFrame;

	 AnimationCriteria aniCri;
	 aniCri.defaultEffectiveDistance = globalConfigVars.getConfigValue(
	 "IntraLinkDisplayRange").toDouble();
	 aniCri.isStressMap = true;

	 CellInitHelper initHelper;

	 SimulationDomainGPU simuDomain;
	 SimulationInitData initData = initHelper.generateDiskInput(loadMeshInput);
	 simuDomain.initialize_V2(initData);

	 simuDomain.checkIfAllDataFieldsValid();

	 */
	GEOMETRY::MeshGen meshGen;
	//std::vector<GEOMETRY::Point2D> points = meshGen.createBdryPointsOnCircle(7,
	//		8);
	//Criteria criteria(0.125, 2.0);
	//GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2D(points,
	//		GEOMETRY::MeshGen::default_list_of_seeds, criteria);
	std::string testString = "./resources/BdryData_unit_test.txt";
	GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2DFromFile(testString);
	std::vector<GEOMETRY::Point2D> centerPoints = mesh.outputTriangleCenters();
	mesh.outputVtkFile("modelTest.vtk");

	return 0;
}
