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
	FileIOException,
	OutputAnalysisDataException,
	MemoryInvalidAccess,
	InvalidInput
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

enum SimulationType {
	Beak, Disc
};

SimulationType parseTypeFromConfig(int configValue);

/**
 * This data structure contains mechanical parameters of the model.
 */
struct SceMechPara {
	double sceInterParaCPU[5];
	double sceIntraParaCPU[4];
	double sceCartParaCPU[5];
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
	bool isStab;
	SimulationType simuType;
	uint maxCellInDomain;
	uint maxNodePerCell;
	uint maxECMInDomain;
	uint maxNodePerECM;
	double FinalToInitProfileNodeCountRatio;
	//double FinalToInitCartNodeCountRatio;
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
 * contains raw information in order to initialize cartilage.
 */
struct CartilageRawData {
	uint pivotNode1Index;
	uint pivotNode2Index;

	uint growNodeBehind1Index;
	uint growNodeBehind2Index;

	// these two indices are in tip.
	uint growNode1Index_on_tip;
	uint growNode2Index_on_tip;

	std::vector<CVector> nonTipVerticies;
	std::vector<CVector> tipVerticies;
};

/**
 * parameters for cartilage.
 */
class CartPara {
public:
	/**
	 * pivot node index should not change after initialization
	 */
	uint pivotNode1Index;
	uint pivotNode2Index;

	/**
	 * grow node 1 and 2 index should be 0 and 1 respectively,
	 * and they should not change.
	 */
	uint growNode1Index;
	uint growNode2Index;

	/**
	 * these two indicies should be always changing.
	 */
	uint growNodeBehind1Index;
	uint growNodeBehind2Index;

	/**
	 * memory allocation related parameter.
	 * value should be 2.
	 */
	uint tipNodeStartPos;

	/**
	 * this value changes with the cartilage grows.
	 */
	uint tipNodeIndexEnd;

	/**
	 * memory allocation related parameter.
	 */
	uint nonTipNodeStartPos;

	/**
	 * this value changes with the cartilage grows.
	 */
	uint nodeIndexEnd;

	/**
	 * means maximum number of nodes in Cartilage.
	 */
	uint nodeIndexTotal;

	CVector fixedPt;
	CVector growthDir;

	CVector node1GrowthDir;
	CVector node2GrowthDir;

	CVector growNode1;
	CVector growNode2;

	double growthSpeedNode1;
	double growthSpeedNode2;

	double currentLength;
	/**
	 * moment of intertia.
	 */
	double moInertia;
	//uint activeCartilageNodeCount;

	double totalTorque;
	double torqueFromEpi;
	double angularSpeed;

	double growthThreshold;
};

/**
 * Generated Raw data that needs reformatting in order to be used for domain initialization.
 */
struct RawDataInput {
	bool isStab;
	CartilageRawData cartilageData;
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
	bool isStab;
	/**
	 * This parameter is necessary for Cartilage, because this CartPara
	 * cannot be generated from config file directly.
	 */
	CartPara cartPara;
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

template<class T>
void printMatrixToFile(vector<vector<T> >& matrix, std::string &fileName) {
	ofstream ofs(fileName.c_str());
	if (ofs.fail()) {
		throw SceException("unable to open file for writing", FileIOException);
	}
	for (uint i = 0; i < matrix.size(); i++) {
		for (uint j = 0; j < matrix[i].size(); j++) {
			ofs << matrix[i][j] << " ";
		}
		ofs << std::endl;
	}
	ofs.close();
}

uint findClosestArrIndexGivenPos(std::vector<CVector> &vecArr, CVector &pos);
#endif /* COMMONDATA_H_ */
