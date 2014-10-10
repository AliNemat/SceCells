#include <iostream>
#include "SimulationDomainGPU.h"
#include "CellInitHelper.h"
#include <stdlib.h>
using namespace std;

const int myDeviceID = 0;
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

/**
 * Generate cell init
 */
void generateCellInitNodeInfo(std::string meshInput, vector<CVector> &initPos) {
	fstream fs;
	uint NNum, LNum, ENum;
	uint i;
	double tmp;
	std::cout << "mesh input is :" << meshInput << std::endl;
	fs.open(meshInput.c_str(), ios::in);
	if (!fs.is_open()) {
		std::string errString =
				"Unable to load mesh in string input mode, meshname: "
						+ meshInput
						+ " ,possible reason is the file is not located in the project folder \n";
		std::cout << errString;
	}
	fs >> NNum >> LNum >> ENum;
	cout << "NNum = " << NNum << std::endl;

	initPos.resize(NNum);

	for (i = 0; i < NNum; i++) {
		fs >> initPos[i].x >> initPos[i].y >> tmp;
		initPos[i].x = initPos[i].x / 20.0;
		initPos[i].y = initPos[i].y / 10.0;
	}
	fs.close();
}

/**
 * Depreciated -- Previous test input setup.
 * This method is highly not robust and is for test purpose only.
 * This method assumes that the size of the computation domain
 * is 50 X 50 and left bottom point is (0,0). it initializes cells of a vertical line.
 * starting from position (5,5) with step of 3.0
 */
void generateCellInitInfo(std::string meshInput, uint numberOfInitCells,
		std::vector<double> &initCellNodePosX,
		std::vector<double> &initCellNodePosY, std::vector<double> &centerPosX,
		std::vector<double> &centerPosY) {
	centerPosX.resize(numberOfInitCells);
	centerPosY.resize(numberOfInitCells);

	fstream fs;
	uint NNum, LNum, ENum;
	uint i;
	double tmp;
	std::cout << "mesh input is :" << meshInput << std::endl;
	fs.open(meshInput.c_str(), ios::in);
	if (!fs.is_open()) {
		std::string errString =
				"Unable to load mesh in string input mode, meshname: "
						+ meshInput
						+ " ,possible reason is the file is not located in the project folder \n";
		std::cout << errString;
	}
	fs >> NNum >> LNum >> ENum;
	cout << "NNum = " << NNum << std::endl;

	initCellNodePosX.resize(NNum);
	initCellNodePosY.resize(NNum);

	for (i = 0; i < NNum; i++) {
		fs >> initCellNodePosX[i] >> initCellNodePosY[i] >> tmp;
		initCellNodePosX[i] = initCellNodePosX[i] / 20.0;
		initCellNodePosY[i] = initCellNodePosY[i] / 10.0;
	}
	fs.close();
	double startCenterPosX = 10.0;
	double startCenterPosY = 5.0;
	double currentCenterPosX = startCenterPosX;
	double currentCenterPosY = startCenterPosY;
	double stepIncrease = 0.99;
	for (i = 0; i < numberOfInitCells; i++) {
		centerPosX[i] = currentCenterPosX;
		centerPosY[i] = currentCenterPosY;
		currentCenterPosY = currentCenterPosY + stepIncrease;
	}
}

int main() {
	srand(time(NULL));
	ConfigParser parser;
	cudaSetDevice(myDeviceID);
	std::string configFileName = "./resources/beak.cfg";
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
	 SimulationInitData initData = initHelper.generateInput(loadMeshInput);
	 simuDomain.initialize(initData);

	 simuDomain.checkIfAllDataFieldsValid();

	 for (int i = 0; i <= numOfTimeSteps; i++) {
	 cout << "step number = " << i << endl;
	 if (i % outputAnimationAuxVarible == 0) {
	 //simuDomain.outputVtkFilesWithColor_v2(animationInput, i);
	 simuDomain.outputVtkFilesWithColor(animationInput, i, aniCri);
	 cout << "finished output Animation" << endl;
	 }
	 simuDomain.runAllLogic(dt);
	 }

	 */

	std::string loadMeshInput;
	std::string animationInput;
	std::vector<std::string> boundaryMeshFileNames;
	std::string animationFolder;
	generateStringInputs(loadMeshInput, animationInput, animationFolder,
			boundaryMeshFileNames);

	double SimulationTotalTime = globalConfigVars.getConfigValue(
			"SimulationTotalTime").toDouble();
	double SimulationTimeStep = globalConfigVars.getConfigValue(
			"SimulationTimeStep").toDouble();
	int TotalNumOfOutputFrames = globalConfigVars.getConfigValue(
			"TotalNumOfOutputFrames").toInt();

	const double simulationTime = SimulationTotalTime;
	const double dt = SimulationTimeStep;
	const int numOfTimeSteps = simulationTime / dt;
	const int totalNumOfOutputFrame = TotalNumOfOutputFrames;
	const int outputAnimationAuxVarible = numOfTimeSteps
			/ totalNumOfOutputFrame;

	CellInitHelper initHelper;

	RawDataInput rawInput = initHelper.generateRawInput_stab(loadMeshInput);

	SimulationInitData simuData = initHelper.initInputsV2(rawInput);

	SimulationDomainGPU simuDomain;

	std::vector<CVector> stabilizedCenters;
	// TODO: These initialization statments are removed for debugging purpose.
	//stabilizedCenters = simuDomain.stablizeCellCenters(
	//		simuData);

	RawDataInput rawInput2 = initHelper.generateRawInput_V2(stabilizedCenters);

	SimulationInitData_V2 simuData2 = initHelper.initInputsV3(rawInput2);

	simuDomain.initialize_v2(simuData2);

	AnimationCriteria aniCri;
	aniCri.defaultEffectiveDistance = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();
	aniCri.isStressMap = false;

	uint aniFrame = 0;
	for (int i = 0; i <= numOfTimeSteps; i++) {
		cout << "step number = " << i << endl;
		if (i % outputAnimationAuxVarible == 0) {
			simuDomain.outputVtkFilesWithColor(animationInput, aniFrame,
					aniCri);
			aniFrame++;
		}
		simuDomain.runAllLogic(SimulationTimeStep);
	}
}

//vector<CVector> initNodePosFromMesh;
//generateCellInitNodeInfo(loadMeshInput, initNodePosFromMesh);

//double Cell_Center_Interval = globalConfigVars.getConfigValue(
//		"Cell_Center_Interval").toDouble();
//double bdryNodeInterval = globalConfigVars.getConfigValue(
//		"Bdry_Node_Interval").toDouble();
//double Cell_deformationRatio = globalConfigVars.getConfigValue(
//		"Cell_deformationRatio").toDouble();

//vector<CellPlacementInfo> cellInfoArray = initHelper.obtainCellInfoArray(
//		1.0, 0.5/sqrt(55));
//vector<CellPlacementInfo> cellInfoArray =
//		initHelper.obtainPreciseCellInfoArray(Cell_Center_Interval,
//				Cell_deformationRatio);

//vector<CVector> bdryNodes, FNMCellCenters, MXCellCenters;
//initHelper.generateThreeInputCellInfoArrays(bdryNodes, FNMCellCenters,
//MXCellCenters, Cell_Center_Interval, bdryNodeInterval
//);

//std::vector<CellType> cellTypes;
//std::vector<uint> numOfInitNodesOfCells;
//std::vector<double> initBdryCellNodePosX;
//std::vector<double> initBdryCellNodePosY;
//std::vector<double> initFNMCellNodePosX;
//std::vector<double> initFNMCellNodePosY;
//std::vector<double> initMXCellNodePosX;
//std::vector<double> initMXCellNodePosY;

//initHelper.initInputsFromCellInfoArray(cellTypes, numOfInitNodesOfCells,
//		initBdryCellNodePosX, initBdryCellNodePosY, initFNMCellNodePosX,
//		initFNMCellNodePosY, initMXCellNodePosX, initMXCellNodePosY,
//		bdryNodes, FNMCellCenters, MXCellCenters, initNodePosFromMesh);

// This is the previous version of cell initialization. depreciated.
//simuDomain.initializeCells(initCellNodePosX, initCellNodePosY, centerPosX,
//		centerPosY, numberOfCellSpaceReserveForBdry);
//simuDomain.initializeCellTypes(initTypes);

//simuDomain.initialCellsOfThreeTypes(cellTypes, numOfInitNodesOfCells,
//		initBdryCellNodePosX, initBdryCellNodePosY, initFNMCellNodePosX,
//		initFNMCellNodePosY, initMXCellNodePosX, initMXCellNodePosY);

//cout << "number of cell spaces reserved:" << cellTypes.size() << endl;

//int jj;
//cin >> jj;

//double lengthDiff = simuDomain.cells.lengthDifference[20];
//double expectedLen = simuDomain.cells.expectedLength[20];
//double growthProgress = simuDomain.cells.growthProgress[20];
//double isGoingToDivide = simuDomain.cells.isDivided[20];
//cout << "length difference: " << lengthDiff << endl;
//cout << "expected length: " << expectedLen << endl;
//cout << "growth progress: " << growthProgress << endl;
//cout << "schedule to divide? " << isGoingToDivide << endl;

//double cellGrowProgress = simuDomain.cells.growthProgress[0];
//bool cellIsCheduledToGrow = simuDomain.cells.isScheduledToGrow[0];
//double lastCheckPoint = simuDomain.cells.lastCheckPoint[0];
//if (cellIsCheduledToGrow) {
//	cout << "growth progress = " << cellGrowProgress << endl;
//	cout << "Is scheduled to grow?" << cellIsCheduledToGrow << endl;
//	cout << "last check point is " << lastCheckPoint << endl;
//}

