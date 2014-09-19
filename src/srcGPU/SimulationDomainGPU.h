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
	uint maxCellInDomain;
	uint maxNodePerCell;
	uint maxECMInDomain;
	uint maxNodePerECM;
	//uint initECMCount;
	double FinalToInitProfileNodeCountRatio;
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
public:
	SceNodes nodes;
	//SceCells cells;
	SceCells_M cells_m;

	GrowthDistriMap growthMap; // first map

	GrowthDistriMap growthMap2; // second map

	// boundary nodes share same attribute with cell nodes
	// but nodes can't move.
	// uint cellSpaceForBdry;

	double intraLinkDisplayRange;

	double minX;
	double maxX;
	double minY;
	double maxY;
	double gridSpacing;

	uint growthGridXDim;
	uint growthGridYDim;
	double growthGridSpacing;
	double growthGridLowerLeftPtX;
	double growthGridLowerLeftPtY;

	// first morphogen distribution
	double growthMorCenterXCoord;
	double growthMorCenterYCoord;
	double growthMorHighConcen;
	double growthMorLowConcen;
	double growthMorDiffSlope;

	// second morphogen distribution
	double growthMorCenterXCoordMX;
	double growthMorCenterYCoordMX;
	double growthMorHighConcenMX;
	double growthMorLowConcenMX;
	double growthMorDiffSlopeMX;

	double compuDist(double x1, double y1, double z1, double x2, double y2,
			double z2) {
		return sqrt(
				(x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2)
						+ (z1 - z2) * (z1 - z2));
	}

	SimulationDomainGPU();

	void initialize_V2(SimulationInitData &initData);

	void runAllLogic(double dt);
	//void outputVtkFilesWithColor(std::string scriptNameBase, int rank);
	//void outputVtkFilesWithColor_v2(std::string scriptNameBase, int rank);
	//void outputVtkFilesWithColor_v3(std::string scriptNameBase, int rank);
	void outputVtkFilesWithColor_v3(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);
	//void outputVtkFilesWithColor_v2_stress(std::string scriptNameBase,
	//int rank);
	void checkIfAllDataFieldsValid();
};

#endif
