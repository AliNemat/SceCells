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

	SimulationGlobalParameter mainPara;
	mainPara.initFromConfig();

	PixelizePara pixelPara;
	pixelPara.initFromConfigFile();

	CellInitHelper initHelper;
	SimulationDomainGPU simuDomain;

	SimulationInitData_V2 initData = initHelper.initStabInput();
	std::vector<CVector> stabilizedCenters = simuDomain.stablizeCellCenters(
			initData);

	SimulationInitData_V2 simuData2 = initHelper.initSimuInput(
			stabilizedCenters);
	simuDomain.initialize_v2(simuData2);

	uint aniFrame = 0;
	std::remove(mainPara.dataOutput.c_str());
	for (int i = 0; i <= mainPara.totalTimeSteps; i++) {
		cout << "step number = " << i << endl;
		if (i % mainPara.aniAuxVar == 0) {
			simuDomain.outputVtkFilesWithColor(mainPara.animationNameBase,
					aniFrame, mainPara.aniCri);
			cout << "finished output Animation" << endl;
			vector<vector<int> > labelMatrix = simuDomain.outputLabelMatrix(
					mainPara.dataName, aniFrame, pixelPara);
			cout << "finished writing label matrix" << endl;
			simuDomain.analyzeLabelMatrix(labelMatrix, aniFrame,
					mainPara.imgOutput, mainPara.dataOutput);
			cout << "finished output matrix analysis" << endl;
			aniFrame++;
		}
		simuDomain.runAllLogic(mainPara.dt);
	}

	return 0;
}
