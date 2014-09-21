#ifndef SimulationDomainGPU_H_
#define SimulationDomainGPU_H_

#include "SceNodes.h"
#include "SceCells.h"
//#include "CellInitHelper.h"
#include "commonData.h"

#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

/**
 * This class is responsible for domain-wise highest level logic, e.g. output animation.
 */
class SimulationDomainGPU {

	SceNodes nodes;
	SceCells_M cells_m;

	GrowthDistriMap growthMap; // first map
	GrowthDistriMap growthMap2; // second map

	SceMemPara memPara;
	SceDomainPara domainPara;
	SceChemPara chemPara;

	void initialCellsOfFiveTypes(std::vector<SceNodeType> &cellTypes,
			std::vector<uint> &numOfInitActiveNodesOfCells,
			std::vector<double> &initBdryCellNodePosX,
			std::vector<double> &initBdryCellNodePosY,
			std::vector<double> &initProfileNodePosX,
			std::vector<double> &initProfileNodePosY,
			std::vector<double> &initECMNodePosX,
			std::vector<double> &initECMNodePosY,
			std::vector<double> &initFNMCellNodePosX,
			std::vector<double> &initFNMCellNodePosY,
			std::vector<double> &initMXCellNodePosX,
			std::vector<double> &initMXCellNodePosY);

	void readMemPara();
	void readDomainPara();
	void readChemPara();
	void readAllParameters();

	void initializeGrowthMap();
public:
	/**
	 * Default constructor.
	 * Reads all the configuration values from file.
	 */
	SimulationDomainGPU();

	/**
	 * Domain initialization.
	 * Assigns values to the data fields in simulation domain.
	 * @param initData initial data set for simulation domain
	 */
	void initialize_V2(SimulationInitData &initData);

	/**
	 * Checks if all data fields are valid.
	 * This methods only loosely checks the data validity.
	 */
	void checkIfAllDataFieldsValid();

	/**
	 * Run one step of simulation in the domain.
	 * Contains cell level logics and node level logics.
	 * @param dt timestep
	 */
	void runAllLogic(double dt);

	/**
	 * Method that animates the domain to VTK format.
	 * @param scriptNameBase name of the vtk animation series.
	 * @param rank frame sequence in the vtk animation series.
	 * @param aniCri criteria for outputing animation.
	 */
	void outputVtkFilesWithColor_v3(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);
};

#endif
