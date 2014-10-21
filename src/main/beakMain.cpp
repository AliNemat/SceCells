#include <iostream>
#include "SimulationDomainGPU.h"
#include "CellInitHelper.h"
#include <stdlib.h>
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
	// initialize config file.
	srand(time(NULL));
	ConfigParser parser;
	std::string configFileName = "./resources/beak.cfg";
	globalConfigVars = parser.parseConfigFile(configFileName);

	// set GPU device.
	int myDeviceID = globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
	gpuErrchk(cudaSetDevice(myDeviceID));

	std::string animationInput = globalConfigVars.getConfigValue(
			"AnimationFolder").toString()
			+ globalConfigVars.getConfigValue("AnimationName").toString();

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

	RawDataInput rawInput = initHelper.generateRawInput_stab();

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

