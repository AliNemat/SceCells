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
#include <set>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <algorithm>

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
	InvalidInput,
	AlgorithmBug
};

enum ECellType {notActive, pouch, peri, bc} ; 
std::string toString(SceExceptionType type);
double compuDistHost(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2);

enum SceNodeType {
	Boundary, Profile, ECM, FNM, MX, Cart, Base, CellIntnl, CellMembr
};
std::string toString(SceNodeType type);

double nodeTypeToScale(SceNodeType type);

bool valueToType(int value);

/**
 * There are only two possible states of switch:
 * On and Off
 */
enum SwitchState {
	ON, OFF
};

class SceException: public std::exception {
private:
	std::string _message;
	SceExceptionType _exceptionType;
	std::string _msg_combine;
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
		std::string _msg_combine = _message + ", "
				+ std::string(", Exception type: " + toString(_exceptionType));
		return _msg_combine.c_str();
	}
};

enum SimulationType {
	Beak, Disc, SingleCellTest, Disc_M
};

SimulationType parseTypeFromConfig(int configValue);

enum AniType {
	CellType, ForceAbsVal, Force, Tension, T1Tran, PolySide
};

class BondInfo {
public:
	CVector pos1, pos2;
	uint cellRank1, cellRank2;
	double val1, val2;
};

AniType parseAniTpFromConfig(int configValue);

struct ControlSwitchs {
	SwitchState stab;
	SwitchState outputLabelMatrix;
	SwitchState outputBmpImg;
	SwitchState outputVtkFile;
	SwitchState outputStat;
};

struct ControlPara {
	SimulationType simuType;
	ControlSwitchs controlSwitchs;
};

/**
 * This data structure contains mechanical parameters of the model.
 */
struct SceMechPara {

	/*
	 * These parameters are useful in all versions.
	 */
	double sceIntraParaCPU[5];
	double sceIntraParaDivCPU[5];

	/*
	 * This parameter is useful in original sce model.
	 */
	double sceInterParaCPU[5];

	/*
	 * There parameters are useful modified disc model
	 */
	double sceBdryBdryParaCPU[5];
	double sceBdryInterParaCPU[5];

	/*
	 * These parameters are for beak project only
	 */
	double sceInterDiffParaCPU[5];
	double sceCartParaCPU[5];
	double sceProfileParaCPU[7];
	double sceECMParaCPU[5];
	double sceDiffParaCPU[5];
};

/**
 * This data structure contains mechanical parameters of the model.
 */
struct SceMechPara_M {
	double sceInterBParaCPU_M[5];
        int    sceInterBParaCPU_Jones_On_M  ; //Ali
	double sceInterBParaCPU_Jones_M[3]; //Ali

	double sceIntnlBParaCPU_M[5];
	double sceIntraParaCPU_M[5];
	double sceIntraParaDivCPU_M[5];
	double growthPrgrCriValCPU_M;
	double maxAdhBondLenCPU_M;
	double minAdhBondLenCPU_M;
	double bondStiffCPU_M;
	double bondStiffCPU_Mitotic;//Ali June 16
	double bondAdhCriLenCPU_M;
};

/**
 * This data structure contains chemical parameters of the model.
 * used only in beak project.
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
 * These parameters controls cell polarity.
 */
struct SceBioPara {
	double cellInitLength;
	double cellFinalLength;
	double elongationCoefficient;
	/*
	 * chemical related parameter is
	 */
	double chemoCoefficient;
};

/**
 * This data structure contains miscellaneous parameters of the model.
 */
struct SceMiscPara {
	double growThreshold;
	double isDivideCriticalRatio;
	double addNodeDistance;
	/*
	 * prevents numerical instability.
	 */
	double minDistanceToOtherNode;

	/**
	 * proliferation slows down with time.
	 * this parameter controls exponential decay.
	 */
	double prolifDecayCoeff;
};

/**
 * This data structure contains parameters about the memory layout of the simulation domain.
 */
struct SceMemPara {
	bool isStab;
	SimulationType simuType;
	uint maxCellInDomain;
	uint maxNodePerCell;

	// these parameters are useful in Disc_M
	uint maxMembrNodePerCell;
	uint maxIntnlNodePerCell;
	uint maxAllNodePerCell;

	// these parameters are useful in Beak project
	uint maxECMInDomain;
	uint maxNodePerECM;
	double FinalToInitProfileNodeCountRatio;
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
 * This data structure contains parameters about the memory layout of the cells
 */
struct CellsMemPara_M {
	uint maxEpiNodeCount; // max possible node count for epithilium nodes
	uint maxInternalNodeCount; //max possible node count for internal nodes
	uint maxTotalCellNodeCount; // max possible node count for epithilium and internal node combined
	uint currentActiveCellCount; // the number of active cells would keep changing
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
	uint XBucketSize;
	uint YBucketSize;
	uint ZBucketSize;
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
 * This data structure contains parameters about the memory allocation in SceNodes_M.
 */
struct NodeAllocPara_M {
	uint bdryNodeCount;
	// @maxNodeOfOneCell represents maximum number of nodes per cell
	uint maxAllNodePerCell;
	uint maxMembrNodePerCell;
	uint maxIntnlNodePerCell;
	// @maxCellCount represents maximum number of cells in the system
	uint maxCellCount;
	// @maxTotalNodeCount represents maximum total number of nodes of all cells
	// maxTotalCellNodeCount = maxNodeOfOneCell * maxCellCount;
	uint maxTotalNodeCount;
	// @currentActiveCellCount represents number of cells that are currently active.
	uint currentActiveCellCount;
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
	SimulationType simuType;
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
 * Generated Raw data that needs reformatting in order to be used for domain initialization.
 */
struct RawDataInput_M {
	bool isStab;
	SimulationType simuType;
	std::vector<CVector> bdryNodes;
	std::vector<CVector> initCellCenters;
	std::vector<double> cellGrowProgVec;
	std::vector<ECellType> cellsTypeCPU; //Ali 
	std::vector<std::vector<CVector> > initIntnlNodePoss;
	std::vector<std::vector<CVector> > initMembrNodePoss;
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
 *
 */
struct SimulationInitData_V2 {
	bool isStab;
	SimulationType simuType;
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

/**
 * a data structure that was specifically designed for Disc project (modified).
 */
struct SimulationInitData_V2_M {
	bool isStab;
	SimulationType simuType;
	std::vector<uint> initActiveMembrNodeCounts;
	std::vector<uint> initActiveIntnlNodeCounts;
	std::vector<double> initGrowProgVec;
	std::vector<SceNodeType> nodeTypes;
	std::vector<CVector> initNodeVec;
	std::vector<bool> initIsActive;
	std::vector<ECellType> eCellTypeV1;  //Ali 
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
 * This class is not used for now but might be useful in the future
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
	AniType animationType;
	// We will only animate links that are close enough.
	double pairDisplayDist;

	double threshold;
	double arrowLength;
	// determines if a potential pair is qualified for animation.
	bool isPairQualify(uint seq1, uint seq2, double x1, double y1, double z1,
			SceNodeType t1, uint r1, double x2, double y2, double z2,
			SceNodeType t2, uint r2);

	bool isPairQualify_M(double x1, double y1, double x2, double y2);

};

/**
 * Necessary information to animate a point in VTK.
 */
struct PointAniData {
	// position of the node.
	CVector pos;
	// In VTK animation software, color scale represents the relative distance to red and blue;
	// bigger value means close to red. smaller value means close to blue.
	CVector dir;
	CVector F_MI_M; //AliE
	double F_MI_M_MagN_Int; //AliE
	CVector extForce;//AAMIRI
	double colorScale;
	double colorScale2;//AAMIRI //curvature
	double colorScale3;//Ali  //membrane tension
	double colorScale4;//Ali //actin 
	int rankScale;//AAMIRI
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
	bool isArrowIncluded;
	std::vector<BondInfo> bondsInfo;
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
	std::cout << "Now printing matrix to file " << fileName.c_str()
			<< std::endl;
	for (uint i = 0; i < matrix.size(); i++) {
		for (uint j = 0; j < matrix[i].size(); j++) {
			ofs << matrix[i][j] << " ";
		}
		ofs << std::endl;
	}
	ofs.close();
}

uint findClosestArrIndexGivenPos(std::vector<CVector> &vecArr, CVector &pos);

struct AblaInfo {
	uint cellNum;
	std::vector<uint> nodeNums;
};

class AblationEvent {
public:
	uint timeStep;
	std::vector<AblaInfo> ablationCells;
	void printInfo();
};

AblationEvent readAblationEvent(std::string inputName);

struct AniRawData {
	std::vector<CVector> aniNodePosArr;
	std::vector<CVector> aniNodeF_MI_M;//AAMIRI // AliE
	std::vector<double> aniNodeF_MI_M_MagN_Int; //AliE
	std::vector<CVector> aniNodeExtForceArr;//AAMIRI
	std::vector<double> aniNodeVal;
	std::vector<double> aniNodeCurvature;//AAMIRI
	std::vector<double> aniNodeMembTension ;//Ali 
	std::vector<double> aniNodeActinLevel ;//Ali 
	std::vector<int> aniNodeRank;//AAMIRI
	std::vector<LinkAniData> memLinks;
	std::vector<LinkAniData> internalLinks;
	std::vector<BondInfo> bondsArr;
};

struct VecVal {
	CVector vec;
	double val;
	bool operator <(const VecVal& other) const {
		return (val < other.val);
	}
};

std::vector<CVector> obtainPtsBetween(CVector& pt1, CVector& pt2,
		double& spacing, uint maxNewMembrNodeCount);

struct CellStats {
	uint cellRank;
	double cellGrowthProgress;
	bool isBdryCell;
	uint numNeighbors;
        double cellCenterX ;  //Ali	
        double membrGrowthProgress;
        double cellPerim;//AAMIRI
	double cellArea;
        int cellNeighborStrength[10]; //Ali
	std::set<int> neighborVec; 
	std::vector<int> neighborVecV; //Ali
	uint currentActiveIntnlNodes;
	uint currentActiveMembrNodes;
	CVector cellCenter;
	void printToFile(ofstream& ofs);
};

struct CountEntry {
	uint numOfNeighbor;
	uint count;
	bool operator <(const CountEntry& other) const {
		return (numOfNeighbor < other.numOfNeighbor);
	}
};

class CellsStatsData {


public:
        //Ali
        double Cells_Extrem_Loc[4] ;
        double F_Ext_Out ; //Ali  
        //Ali 
        double MaxDistanceX ; //Ali 
	std::vector<CellStats> cellsStats;
	void printPolyCountToFile(std::string fileName, double divThreshold);
	void printDetailStatsToFile(std::string fileNameBase, int timestep);
	vector<double> outputPolySides();
        void printStressStrain(std::string FileName1,double curTime,double Init_Displace );   //Ali
        void printStressStrain_Ini(std::string FileName1); // Ali
};

void insertCount(uint numNeighbor, std::map<uint, uint>& count);
void printCountsToFile(std::string fileName, std::map<uint, uint>& countNormal,
		std::map<uint, uint>& countDiv, std::map<uint, uint>& countBdry);

std::vector<CountEntry> processCountMap(std::map<uint, uint>& countMap);

void printEntriesToFile(ofstream& fs, std::vector<CountEntry>& countEntries);
//Ali
//Ali 
#endif /* COMMONDATA_H_ */
