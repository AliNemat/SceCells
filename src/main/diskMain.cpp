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

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort =
		true) {
	if (code != cudaSuccess) {
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file,
				line);
		if (abort)
			exit(code);
	}
}

int main() {
	srand(time(NULL));
	ConfigParser parser;
	std::string configFileName = "./resources/disk.cfg";
	globalConfigVars = parser.parseConfigFile(configFileName);

	// set GPU device.
	int myDeviceID = globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
	gpuErrchk(cudaSetDevice(myDeviceID));

	std::string animationInput = globalConfigVars.getConfigValue(
			"AnimationFolder").toString()
			+ globalConfigVars.getConfigValue("AnimationName").toString();
	std::string dataOutput =
			globalConfigVars.getConfigValue("DataOutputFolder").toString()
					+ globalConfigVars.getConfigValue("PolygonStatFileName").toString();
	std::string imgOutput =
			globalConfigVars.getConfigValue("DataOutputFolder").toString()
					+ globalConfigVars.getConfigValue("ImgFileNameBase").toString();

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

	AnimationCriteria aniCri;
	aniCri.defaultEffectiveDistance = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();
	aniCri.isStressMap = true;

	std::string dataFolder =
			globalConfigVars.getConfigValue("DataFolder").toString();
	std::string dataName = dataFolder
			+ globalConfigVars.getConfigValue("DataName").toString();
	PixelizePara pixelPara;
	pixelPara.initFromConfigFile();

	CellInitHelper initHelper;

	RawDataInput rawInput = initHelper.generateRawInput_stab();

	SimulationInitData simuData = initHelper.initInputsV2(rawInput);

	SimulationDomainGPU simuDomain;

	std::vector<CVector> stabilizedCenters;

	stabilizedCenters = simuDomain.stablizeCellCenters(simuData);
	std::cout << "begin generating raw input data" << std::endl;
	std::cout.flush();
	RawDataInput rawInput2 = initHelper.generateRawInput_V2(stabilizedCenters);
	std::cout << "finished generating raw input data" << std::endl;
	std::cout.flush();
	SimulationInitData_V2 simuData2 = initHelper.initInputsV3(rawInput2);

	std::cout << "finished generating simulation init data V2" << std::endl;
	std::cout.flush();
	simuDomain.initialize_v2(simuData2);

	uint aniFrame = 0;
	for (int i = 0; i <= numOfTimeSteps; i++) {
		cout << "step number = " << i << endl;
		if (i % outputAnimationAuxVarible == 0) {
			cout << "started to output Animation" << endl;
			simuDomain.outputVtkFilesWithColor(animationInput, aniFrame,
					aniCri);
			cout << "finished output Animation" << endl;
			cout << "started writing label matrix" << endl;
			vector<vector<int> > labelMatrix = simuDomain.outputLabelMatrix(
					dataName, aniFrame, pixelPara);
			simuDomain.analyzeLabelMatrix(labelMatrix, aniFrame, imgOutput,
					dataOutput);
			cout << "finished writing label matrix" << endl;
			aniFrame++;
		}
		simuDomain.runAllLogic(dt);
	}

	return 0;
}
