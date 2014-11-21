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
ConfigVarsCollection configCollection;

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

void initModelTestConfigCollection() {
	// read configuration collection.
	ConfigParser parser;
	std::string configCollectionFile = "./resources/modelTestCollection.cfg";
	configCollection = parser.parseConfigCollection(configCollectionFile);

	std::cout << "Now printing all parameter sets for model testing"
			<< std::endl;
	std::cout << "Number of parameter set = "
			<< configCollection.configVarSets.size() << std::endl;

	for (uint i = 0; i < configCollection.configVarSets.size(); i++) {
		configCollection.configVarSets[i].printAll();
		std::cout << std::endl;
	}

	std::cout << "Finished printing parameter sets. Press ENter to continue."
			<< std::endl;
	getchar();

}

void readModelTestConfig(uint i) {
	ConfigParser parser;
	std::string configFileNameDefault = "./resources/modelTest.cfg";
	globalConfigVars = parser.parseConfigFile(configFileNameDefault);

	globalConfigVars.updateFromConfig(configCollection.configVarSets[i]);

	// set GPU device.
	int myDeviceID = globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
	gpuErrchk(cudaSetDevice(myDeviceID));

	//globalConfigVars.printAll();
	//getchar();
}

int main(int argc, char* argv[]) {
	// initialize random seed.
	srand(time(NULL));

	initModelTestConfigCollection();

	for (uint ii = 0; ii < configCollection.configVarSets.size(); ii++) {
		readModelTestConfig(ii);

		// initialize simulation control related parameters from config file.
		SimulationGlobalParameter mainPara;
		mainPara.initFromConfig();

		// initialize post-processing related parameters from config file.
		PixelizePara pixelPara;
		pixelPara.initFromConfigFile();

		// initialize simulation initialization helper.
		CellInitHelper initHelper;
		// initialize simulation domain.
		SimulationDomainGPU simuDomain;

		// Generate initial inputs for simulation domain.
		SimulationInitData_V2 simuData = initHelper.initSingleCellTest();
		// initialize domain based on initial inputs.
		simuDomain.initialize_v2(simuData);

		// delete old data file.
		std::remove(mainPara.dataOutput.c_str());

		// preparation.
		uint aniFrame = 0;
		// main simulation steps.
		for (int i = 0; i <= mainPara.totalTimeSteps; i++) {
			cout << "step number = " << i << endl;
			if (i % mainPara.aniAuxVar == 0) {
				simuDomain.outputVtkFilesWithColor(mainPara.animationNameBase,
						aniFrame, mainPara.aniCri);
				cout << "finished output Animation" << endl;
				//vector<vector<int> > labelMatrix = simuDomain.outputLabelMatrix(
				//		mainPara.dataName, aniFrame, pixelPara);
				//cout << "finished writing label matrix" << endl;
				//simuDomain.analyzeLabelMatrix(labelMatrix, aniFrame,
				//		mainPara.imgOutput, mainPara.dataOutput);
				//cout << "finished output matrix analysis" << endl;
				aniFrame++;
			}
			// for each step, run all logics of the domain.
			simuDomain.runAllLogic(mainPara.dt);
		}
	}

	return 0;
}
