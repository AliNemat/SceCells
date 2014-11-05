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

	SimulationGlobalParameter mainPara;
	mainPara.initFromConfig();

	CellInitHelper initHelper;
	SimulationDomainGPU simuDomain;

	SimulationInitData_V2 initData = initHelper.initStabInput();
	std::vector<CVector> stabilizedCenters = simuDomain.stablizeCellCenters(
			initData);

	SimulationInitData_V2 simuData = initHelper.initSimuInput(
			stabilizedCenters);
	simuDomain.initialize_v2(simuData);

	uint aniFrame = 0;
	for (int i = 0; i <= mainPara.totalTimeSteps; i++) {
		cout << "step number = " << i << endl;
		if (i % mainPara.aniAuxVar == 0) {
			simuDomain.outputVtkFilesWithColor(mainPara.animationNameBase,
					aniFrame, mainPara.aniCri);
			aniFrame++;
		}
		simuDomain.runAllLogic(mainPara.dt);
	}
}

