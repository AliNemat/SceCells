//============================================================================
// Name        : Main.cpp
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
//test here
void initializeSlurmConfig(int argc, char* argv[]) {
	// read configuration.
	ConfigParser parser;
	std::string configFileNameDefault = "./resources/disc_M.cfg";
	globalConfigVars = parser.parseConfigFile(configFileNameDefault);
	std::string configFileNameBaseL = "./resources/disc_";
	std::string configFileNameBaseR = ".cfg";

	// Unknown number of input arguments.
	if (argc != 1 && argc != 3) {
		std::cout << "ERROR: Incorrect input argument count.\n"
				<< "Expect either no command line argument or three arguments"
				<< std::endl;
		exit(0);
	}
	// one input argument. It has to be "-slurm".
	else if (argc == 3) {
		if (strcmp(argv[1], "-slurm") != 0) {
			std::cout
					<< "ERROR: one argument received from commandline but it's not recognized.\n"
					<< "Currently, the argument value must be -slurm"
					<< std::endl;
			exit(0);
		} else {
			std::string configFileNameM(argv[2]);
			std::string configFileNameCombined = configFileNameBaseL
					+ configFileNameM + configFileNameBaseR;
			parser.updateConfigFile(globalConfigVars, configFileNameCombined);
		}
	}
	// no input argument. Take default.
	else if (argc == 1) {

		// set GPU device.
		int myDeviceID =
				globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
		gpuErrchk(cudaSetDevice(myDeviceID));
	}
}

void updateDivThres(double& curDivThred, uint& i, double& dt,
		double& decayCoeff, double& divThreshold) {
	double curTime = i * dt;
	double decay = exp(-curTime * decayCoeff);
	//curDivThred = 1.0 - (1.0 - divThreshold) * decay;
	curDivThred = divThreshold ;
}

int main(int argc, char* argv[]) {
	// initialize random seed.
	srand(time(NULL));

	// Slurm is computer-cluster management system.
	initializeSlurmConfig(argc, argv);

	// initialize simulation control related parameters from config file.
	SimulationGlobalParameter mainPara;
	mainPara.initFromConfig();

	// initialize simulation initialization helper.
	CellInitHelper initHelper;

	// initialize simulation domain.
	SimulationDomainGPU simuDomain;

	SimulationInitData_V2_M initData = initHelper.initInput_M();
	simuDomain.initialize_v2_M(initData);

	std::string polyStatFileNameBase = globalConfigVars.getConfigValue(
			"PolygonStatFileName").toString();
	std::string uniqueSymbol =
			globalConfigVars.getConfigValue("UniqueSymbol").toString();
	std::string polyStatFileName = polyStatFileNameBase + uniqueSymbol + ".txt";

	std::remove(polyStatFileName.c_str());

	std::string detailStatFileNameBase = globalConfigVars.getConfigValue(
			"DetailStatFileNameBase").toString() + uniqueSymbol;
	double divThreshold =
			globalConfigVars.getConfigValue("DivThreshold").toDouble();
	double decayCoeff =
			globalConfigVars.getConfigValue("ProlifDecayCoeff").toDouble();
	double curDivThred;

	int maxStepTraceBack =
			globalConfigVars.getConfigValue("MaxStepTraceBack").toInt();

	// preparation.
	uint aniFrame = 0;
	// main simulation steps.
	for (uint i = 0; i <= (uint) (mainPara.totalTimeSteps); i++) {
		if (i % mainPara.aniAuxVar == 0) {
			std::cout << "substep 1 " << std::endl;
			CellsStatsData polyData = simuDomain.outputPolyCountData();

			std::cout << "substep 2 " << std::endl;
			//////// update division threshold //////
			updateDivThres(curDivThred, i, mainPara.dt, decayCoeff,
					divThreshold);

			std::cout << "substep 3 " << std::endl;
			// prints brief polygon counting statistics to file
			polyData.printPolyCountToFile(polyStatFileName, curDivThred);
			// prints detailed individual cell statistics to file
			polyData.printDetailStatsToFile(detailStatFileNameBase, aniFrame);
			// prints the animation frames to file. They can be open by Paraview

			//std::cout << "substep 4 " << std::endl;
			//if (i != 0) {
				//simuDomain.processT1Info(maxStepTraceBack, polyData);
			//}

			std::cout << "substep 5 " << std::endl;
			//simuDomain.outputVtkFilesWithCri_M(mainPara.animationNameBase,
			//		aniFrame, mainPara.aniCri);
			//simuDomain.outputVtkColorByCell_T1(mainPara.animationNameBase,
			//		aniFrame, mainPara.aniCri);
			simuDomain.outputVtkColorByCell_polySide(mainPara.animationNameBase,
					aniFrame, mainPara.aniCri);
			// std::cout << "in ani step " << aniFrame << std::endl;
			std::cout << "substep 6 " << std::endl;
			aniFrame++;
		}
		simuDomain.runAllLogic_M(mainPara.dt);
	}

	return 0;
}
