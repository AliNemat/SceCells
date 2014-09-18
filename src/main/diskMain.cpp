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
	std::string configFileName = "./resources/disk.cfg";
	globalConfigVars = parser.parseConfigFile(configFileName);

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

	CellInitHelper initHelper;

	SimulationDomainGPU simuDomain;
	SimulationInitData initData = initHelper.generateDiskInput(loadMeshInput);
	simuDomain.initialize_V2(initData);

	simuDomain.checkIfAllDataFieldsValid();

	for (int i = 0; i <= numOfTimeSteps; i++) {
		cout << "step number = " << i << endl;
		if (i % outputAnimationAuxVarible == 0) {
			//simuDomain.outputVtkFilesWithColor_v2(animationInput, i);
			//simuDomain.outputVtkFilesWithColor_v2_stress(animationInput, i);
			simuDomain.outputVtkFilesWithColor_v3(animationInput, i);
			cout << "finished output Animation" << endl;
		}
		simuDomain.runAllLogic(dt);
	}
	return 0;
}
