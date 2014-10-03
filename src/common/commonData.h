#include <vector>
#include "GeoVector.h"
#include <fstream>
#include <exception>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <map>
#include <sstream>
#include <iomanip>
#include <fstream>

#ifndef COMMONDATA_H_
#define COMMONDATA_H_

typedef unsigned int uint;
typedef std::map<uint, uint> IndexMap;

enum SceExceptionType {
	BaseException,
	InputInitException,
	ConfigFileNotFound,
	ConfigValueException,
	OutputAnalysisDataException
};

std::string toString(SceExceptionType type);
double compuDistHost(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2);

enum SceNodeType {
	Boundary, Profile, ECM, FNM, MX, Cart, Base
};
std::string toString(SceNodeType type);

double nodeTypeToScale(SceNodeType type);

class SceException: public std::exception {
private:
	std::string _message;
	SceExceptionType _exceptionType;
public:
	SceException(const std::string& message) :
			_message(message), _exceptionType(BaseException) {
	}
	SceException(const std::string& message, const SceExceptionType type) :
			_message(message), _exceptionType(type) {
	}
	~SceException() throw () {
	}
	virtual const char* what() const throw () {
		std::string fullErrorMsg = std::string(
				_message + ", Exception type: " + toString(_exceptionType));
		return fullErrorMsg.c_str();
	}
};

/**
 * This data structure contains mechanical parameters of the model.
 */
struct SceMechPara {
	double sceInterParaCPU[5];
	double sceIntraParaCPU[4];
	double sceInterDiffParaCPU[5];
	double sceProfileParaCPU[7];
	double sceECMParaCPU[5];
	double sceDiffParaCPU[5];
};

/**
 * This data structure contains chemical parameters of the model.
 */
struct SceChemPara {
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
};

/**
 * This data structure contains biology parameters of the model.
 */
struct SceBioPara {
	double cellInitLength;
	double cellFinalLength;
	double elongationCoefficient;
	double chemoCoefficient;
};

/**
 * This data structure contains miscellaneous parameters of the model.
 */
struct SceMiscPara {
	double growThreshold;
	double isDivideCriticalRatio;
	double addNodeDistance;
	double minDistanceToOtherNode;
};

/**
 * This data structure contains parameters about the memory layout of the simulation domain.
 */
struct SceMemPara {
	uint maxCellInDomain;
	uint maxNodePerCell;
	uint maxECMInDomain;
	uint maxNodePerECM;
	double FinalToInitProfileNodeCountRatio;
	double FinalToInitCartNodeCountRatio;
};

/**
 * This data structure contains parameters about the memory layout of the cells
 */
struct CellsMemPara {
	uint maxTotalCellNodeCount; // max possible node count for cell representation
	uint currentActiveCellCount; // the number of active cells would keep changing
	uint currentActiveECMCount;     // the number of active ECM might change.
	uint currectActiveProfileNode; // the number of epithilum nodes might change.
};

/**
 * This data structure contains parameters about the setup of the simulation domain.
 */
struct SceDomainPara {
	double minX;
	double maxX;
	double minY;
	double maxY;
	double minZ;
	double maxZ;
	double gridSpacing;
	uint numOfBucketsInXDim;
	uint numOfBucketsInYDim;
	uint totalBucketCount;
};

/**
 * This data structure contains parameters about the memory allocation in SceNodes.
 */
struct NodeAllocPara {
	// @maxNodeOfOneCell represents maximum number of nodes per cell
	uint maxNodeOfOneCell;
	// @maxCellCount represents maximum number of cells in the system
	uint maxCellCount;
	// @maxTotalNodeCount represents maximum total number of nodes of all cells
	// maxTotalCellNodeCount = maxNodeOfOneCell * maxCellCount;
	uint maxTotalCellNodeCount;
	// @currentActiveCellCount represents number of cells that are currently active.
	uint currentActiveCellCount;

	// @maxNodePerECM represents maximum number of nodes per ECM
	uint maxNodePerECM;
	// @maxECMCount represents maximum number of ECM
	uint maxECMCount;
	// @maxTotalECMNodeCount represents maximum total number of node of ECM
	uint maxTotalECMNodeCount;
	// @currentActiveECM represents number of ECM that are currently active.
	uint currentActiveECM;

	// epithilum might grow or might not. Set maxProfileNodeCount as the maximum possible node count
	uint maxProfileNodeCount;
	// no matter whether epithilum grow or not we need to track the cucrent count.
	uint currentActiveProfileNodeCount;

	// epithilum might grow or might not. Set maxProfileNodeCount as the maximum possible node count
	uint maxCartNodeCount;
	// no matter whether epithilum grow or not we need to track the cucrent count.
	uint currentActiveCartNodeCount;

	uint BdryNodeCount;

	uint startPosProfile;
	uint startPosCart;
	uint startPosECM;
	uint startPosCells;
};

/**
 * Generated Raw data that needs reformatting in order to be used for domain initialization.
 */
struct RawDataInput {
	std::vector<CVector> bdryNodes;
	std::vector<CVector> profileNodes;
	std::vector<CVector> FNMCellCenters;
	std::vector<CVector> MXCellCenters;
	std::vector<CVector> ECMCenters;
	std::vector<double> ECMAngles;
	std::vector<CVector> initCellNodePoss;
	std::vector<CVector> initECMNodePoss;
};

/**
 * a data structure that was specifically designed for Beak project.
 */
struct SimulationInitData {
	std::vector<SceNodeType> cellTypes;
	std::vector<uint> numOfInitActiveNodesOfCells;
	std::vector<double> initBdryCellNodePosX;
	std::vector<double> initBdryCellNodePosY;
	std::vector<double> initProfileNodePosX;
	std::vector<double> initProfileNodePosY;
	std::vector<double> initECMNodePosX;
	std::vector<double> initECMNodePosY;
	std::vector<double> initFNMCellNodePosX;
	std::vector<double> initFNMCellNodePosY;
	std::vector<double> initMXCellNodePosX;
	std::vector<double> initMXCellNodePosY;
};

/**
 * a data structure that was specifically designed for Beak project.
 */
struct SimulationInitData_V2 {
	std::vector<SceNodeType> cellTypes;
	std::vector<uint> numOfInitActiveNodesOfCells;
	std::vector<CVector> initBdryNodeVec;
	std::vector<CVector> initProfileNodeVec;
	std::vector<CVector> initCartNodeVec;
	std::vector<CVector> initECMNodeVec;
	std::vector<CVector> initFNMNodeVec;
	std::vector<CVector> initMXNodeVec;
};

std::vector<double> getArrayXComp(std::vector<CVector> &nodePosVec);
std::vector<double> getArrayYComp(std::vector<CVector> &nodePosVec);
std::vector<double> getArrayZComp(std::vector<CVector> &nodePosVec);

/**
 * This class is not used for now but might be useful in the future
 */
struct SceInputPoint {
	static const std::string delimiter;
	uint cellRank;
	SceNodeType cellType;
	double xCoord;
	double yCoord;
	double zCoord;
};

/**
 * This class is not used for nwo but might be useful in the future
 */
struct inputInitialData {
	SceMemPara memParas;
	SceMechPara mechParas;
	std::vector<SceInputPoint> inputPoints;
	void addNewPoints(std::vector<SceInputPoint> &newPoints);
};

/**
 * Data structure that controls the animation criteria.
 */
struct AnimationCriteria {
	// If this varible is set to be true, output stress map;
	// otherwise, output normal animation.
	bool isStressMap;
	// We will only animate links that are close enough.
	double defaultEffectiveDistance;
	// determines if a potential pair is qualified for animation.
	bool isPairQualify(uint seq1, uint seq2, double x1, double y1, double z1,
			SceNodeType t1, uint r1, double x2, double y2, double z2,
			SceNodeType t2, uint r2);

};

/**
 * Necessary information to animate a point in VTK.
 */
struct PointAniData {
	// position of the node.
	CVector pos;
	// In VTK animation software, color scale represents the relative distance to red and blue;
	// bigger value means close to red. smaller value means close to blue.
	double colorScale;
};

/**
 * Necessary information to animate a link in VTK.
 */
struct LinkAniData {
	uint node1Index, node2Index;
};

/**
 * contains all information need to be
 */
struct VtkAnimationData {
	std::vector<PointAniData> pointsAniData;
	std::vector<LinkAniData> linksAniData;
	void outputVtkAni(std::string scriptNameBase, int rank);
};

#endif /* COMMONDATA_H_ */
