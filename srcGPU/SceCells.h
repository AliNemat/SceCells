#ifndef SCECELLS_H_
#define SCECELLS_H_

#include "SceNodes.h"
#include "GrowthDistriMap.h"
#include <time.h>
#include <thrust/tabulate.h>

typedef thrust::tuple<double, double, CellType> CVec2Type;
typedef thrust::tuple<bool, CellType> boolType;
typedef thrust::tuple<double, double, bool, CellType, uint> Vel2DActiveTypeRank;

/**
 * Functor for divide operation.
 * @param dividend divisor for divide operator.
 * @param input1: number to be divided
 * @return output: result from division
 */
struct DivideFunctor: public thrust::unary_function<uint, uint> {
	uint dividend;__host__ __device__
	DivideFunctor(uint dividendInput) :
			dividend(dividendInput) {
	}
	__host__ __device__
	uint operator()(const uint &num) {
		return num / dividend;
	}
};

/**
 * Functor for modulo operation.
 * @param dividend divisor for modulo operator.
 * @param input1: number to be moduled
 * @return output: result from modulo
 */
struct ModuloFunctor: public thrust::unary_function<uint, uint> {
	uint dividend;__host__ __device__
	ModuloFunctor(uint dividendInput) :
			dividend(dividendInput) {
	}
	__host__ __device__
	uint operator()(const uint &num) {
		return num % dividend;
	}
};

/**
 * Functor predicate see if a boolean varible is true(seems unnecessary but still required).
 */
struct isTrue {
	__host__ __device__
	bool operator()(bool b) {
		if (b == true) {
			return true;
		} else {
			return false;
		}
	}
};

struct isActiveNoneBdry {
	__host__ __device__
	bool operator()(boolType b) {
		bool isActive = thrust::get<0>(b);
		CellType type = thrust::get<1>(b);
		if (isActive && type != Boundary) {
			return true;
		} else {
			return false;
		}
	}
};

/**
 * Functor for add two three dimensional vectors.
 * @param input1 first three dimensional vector to add
 * @param input2 second three dimensional vector to add
 * @return output result of addition
 */
struct CVec3Add: public thrust::binary_function<CVec3, CVec3, CVec3> {
	__host__ __device__
	CVec3 operator()(const CVec3 &vec1, const CVec3 &vec2) {
		return thrust::make_tuple(thrust::get<0>(vec1) + thrust::get<0>(vec2),
				thrust::get<1>(vec1) + thrust::get<1>(vec2),
				thrust::get<2>(vec1) + thrust::get<2>(vec2));
	}
};

/**
 * Divide three inputs by one same number.
 * @param input1 first number to be divide \n
 *        input2 second number to be divide \n
 *        input3 third number to be divide \n
 * @return output1 first division result \n
 *        output2 second division result \n
 *        output3 third division result \n
 */
struct CVec3Divide: public thrust::binary_function<CVec3, double, CVec3> {
	__host__ __device__
	CVec3 operator()(const CVec3 &vec1, const double &divisor) {
		return thrust::make_tuple(thrust::get<0>(vec1) / divisor,
				thrust::get<1>(vec1) / divisor, thrust::get<2>(vec1) / divisor);
	}
};

/**
 * Obtain growth speed and direction given node position.
 * @param _gridDimensionX number of grid points in x direction
 * @param _gridDimensionY number of grid points in y direction
 * @param _gridSpacing spacing of the chemical signal mesh.
 * @param _gridMagValue begin address of growth speed vector
 * @param _gridDirXCompValue begin address of growth direction x component vector
 * @param _gridDirYCompValue begin address of growth direction y component vector
 *
 * @param input1 x coordinate of node position
 * @param input2 y coordinate of node position
 *
 * @return output1 growth speed \n
 *         output2 x component of growth direction \n
 *         output3 y component of growth direction \n
 *
 */
struct LoadGridDataToNode: public thrust::unary_function<CVec2, CVec3> {
	uint _gridDimensionX;
	uint _gridDimensionY;
	double _gridSpacing;
	double* _gridMagValue;
	double* _gridDirXCompValue;
	double* _gridDirYCompValue;__host__ __device__
	LoadGridDataToNode(uint gridDimensionX, uint gridDimensionY,
			double gridSpacing, double* gridMagValue, double* gridDirXCompValue,
			double* gridDirYCompValue) :
			_gridDimensionX(gridDimensionX), _gridDimensionY(gridDimensionY), _gridSpacing(
					gridSpacing), _gridMagValue(gridMagValue), _gridDirXCompValue(
					gridDirXCompValue), _gridDirYCompValue(gridDirYCompValue) {
	}
	__host__ __device__
	CVec3 operator()(const CVec2 &d2) const {
		double xCoord = thrust::get<0>(d2);
		double yCoord = thrust::get<1>(d2);
		uint gridLoc = (uint) (xCoord / _gridSpacing)
				+ (uint) (yCoord / _gridSpacing) * _gridDimensionX;
		double magRes = _gridMagValue[gridLoc];
		double xDirRes = _gridDirXCompValue[gridLoc];
		double yDirRes = _gridDirYCompValue[gridLoc];
		return thrust::make_tuple(magRes, xDirRes, yDirRes);
	}
};

/**
 * Obtain growth speed and direction given node position and chemical field.
 * @param _gridDimensionX number of grid points in x direction
 * @param _gridDimensionY number of grid points in y direction
 * @param _gridSpacing spacing of the chemical signal mesh.
 * @param _gridMagValue begin address of growth speed vector
 * @param _gridDirXCompValue begin address of growth direction x component vector
 * @param _gridDirYCompValue begin address of growth direction y component vector
 *
 * @param input1 x coordinate of node position
 * @param input2 y coordinate of node position
 * @param input3 type of the cell. Could be boundary, FNM or MX
 *
 * @return output1 growth speed \n
 *         output2 x component of growth direction \n
 *         output3 y component of growth direction \n
 *
 */
struct LoadChemDataToNode: public thrust::unary_function<CVec2Type, CVec3> {
	uint _gridDimensionX;
	uint _gridDimensionY;
	double _gridSpacing;
	double* _gridMagValue;
	double* _gridDirXCompValue;
	double* _gridDirYCompValue;
	uint _gridDimensionX2;
	uint _gridDimensionY2;
	double _gridSpacing2;
	double* _gridMagValue2;
	double* _gridDirXCompValue2;
	double* _gridDirYCompValue2;

	__host__ __device__
	LoadChemDataToNode(uint gridDimensionX, uint gridDimensionY,
			double gridSpacing, double* gridMagValue, double* gridDirXCompValue,
			double* gridDirYCompValue, uint gridDimensionX2,
			uint gridDimensionY2, double gridSpacing2, double* gridMagValue2,
			double* gridDirXCompValue2, double* gridDirYCompValue2) :
			_gridDimensionX(gridDimensionX), _gridDimensionY(gridDimensionY), _gridSpacing(
					gridSpacing), _gridMagValue(gridMagValue), _gridDirXCompValue(
					gridDirXCompValue), _gridDirYCompValue(gridDirYCompValue), _gridDimensionX2(
					gridDimensionX2), _gridDimensionY2(gridDimensionY2), _gridSpacing2(
					gridSpacing2), _gridMagValue2(gridMagValue2), _gridDirXCompValue2(
					gridDirXCompValue2), _gridDirYCompValue2(gridDirYCompValue2) {
	}
	__host__ __device__
	CVec3 operator()(const CVec2Type &d2) const {
		double xCoord = thrust::get<0>(d2);
		double yCoord = thrust::get<1>(d2);
		CellType type = thrust::get<2>(d2);
		uint gridLoc = (uint) (xCoord / _gridSpacing)
				+ (uint) (yCoord / _gridSpacing) * _gridDimensionX;
		if (type == FNM) {
			double magRes = _gridMagValue[gridLoc];
			double xDirRes = _gridDirXCompValue[gridLoc];
			double yDirRes = _gridDirYCompValue[gridLoc];
			return thrust::make_tuple(magRes, xDirRes, yDirRes);
		} else if (type == MX) {
			double magRes = _gridMagValue2[gridLoc];
			double xDirRes = _gridDirXCompValue2[gridLoc];
			double yDirRes = _gridDirYCompValue2[gridLoc];
			return thrust::make_tuple(magRes, xDirRes, yDirRes);
		} else {
			return thrust::make_tuple(0.0, 1.0, 0.0);
		}

	}
};

/**
 * One dimensional version of a*X plus Y.
 * @param input1 X
 * @param input2 Y
 *
 * @return output1 a*X+Y
 */
struct SaxpyFunctor: public thrust::binary_function<double, double, double> {
	double _dt;__host__ __device__
	SaxpyFunctor(double dt) :
			_dt(dt) {
	}
	__host__ __device__
	double operator()(const double &x, const double &y) {
		return _dt * x + y;
	}
};

/**
 * One dimensional version of a*X plus Y, return one if result is larger than one.
 * @param input1 X
 * @param input2 Y
 *
 * @return output1 a*X+Y
 */
struct SaxpyFunctorWithMaxOfOne: public thrust::binary_function<double, double,
		double> {
	double _dt;__host__ __device__
	SaxpyFunctorWithMaxOfOne(double dt) :
			_dt(dt) {
	}
	__host__ __device__
	double operator()(const double &x, const double &y) {
		double result = _dt * x + y;
		if (result > 1.0) {
			return 1.0;
		} else {
			return result;
		}
	}
};

/**
 * Two dimensional version of a*X plus Y.
 * @param input1 x and y components of X
 * @param input2 x and y components of Y
 *
 * @return output1 x and y compoents of result
 */
struct SaxpyFunctorDim2: public thrust::binary_function<CVec2, CVec2, CVec2> {
	double _dt;__host__ __device__
	SaxpyFunctorDim2(double dt) :
			_dt(dt) {
	}
	__host__ __device__
	CVec2 operator()(const CVec2 &vec1, const CVec2 &vec2) {
		double xRes = thrust::get<0>(vec1) * _dt + thrust::get<0>(vec2);
		double yRes = thrust::get<1>(vec1) * _dt + thrust::get<1>(vec2);
		return thrust::make_tuple(xRes, yRes);
	}
};

/**
 * Point condition operater, decide if cell is ready to add a new point.
 * @param _threshold threshold value for difference of current progress and last checkpoint.
 * if difference is bigger than threshold then the cell is ready for adding a new node.
 * @param input1 growth progress \n
 * @param input2 last check point \n
 * @return output1 is the cell going to add one more node?
 * @return output2 updated check point value (change or unchanged)
 */
struct PtCondiOp: public thrust::unary_function<CVec2, bool> {
	double _threshold;__host__ __device__
	PtCondiOp(double threshold) :
			_threshold(threshold) {
	}
	__host__ __device__
	bool operator()(const CVec2 &d2) const {
		double progress = thrust::get<0>(d2);
		double lastCheckPoint = thrust::get<1>(d2);
		bool resBool = false;
		if (progress == 1.0 && lastCheckPoint < 1.0) {
			resBool = true;
		}
		if (progress - lastCheckPoint >= _threshold) {
			resBool = true;
		}
		return resBool;
	}
};

/**
 * return zero given a celltype
 */
struct GetZero: public thrust::unary_function<CellType, double> {
	__host__ __device__
	double operator()(const CellType &type) {
		return 0.0;
	}
};
/**
 * determines if cell type is boundary.
 */
struct IsBoundary: public thrust::unary_function<CellType, bool> {
	__host__ __device__
	uint operator()(const CellType &type) {
		if (type == Boundary) {
			return true;
		} else {
			return false;
		}
	}
};

/**
 * Unary opterator for adding new node in the cell.
 * BoolUIDDUI consists of the following:
 *
 * (1) - (Bool,bool) is this cell scheduled to grow?
 *
 * (2) - (UI,unsigned integer) how many active nodes are there in this cell?
 * we need this input to decide where should we place the coordinate information of the new node
 *
 * (3) - (D, double) x coordinate of the cell center
 *
 * (4) - (D, double) y coordinate of the cell center
 * we need these two cell center coordinates because cuda only has a pseduo-random number generator,
 * so we need to obtain a good seed to generate a random number. Here we choose center position of the cell.
 *
 * (5) - (UI, unsigned integer) rank of the cell
 *
 * BoolUI consists of the following:
 *
 * (1) - (Bool,bool) if operation succeed, this will return 0. otherwise, return 1
 *
 * (2) - (UI,unsigned integer) how many active nodes are there in this cell? if operation succeed,
 * this will input active node count + 1. otherwise, return input active node count
 *
 * @param _maxNodeOfOneCell Maximum node count of a cell.
 * @param _addNodeDistance  While adding a node, we need to set a fixed distance as radius of the circle
 *        that we would like to add a point.
 * @param _minDistanceToOtherNode Minimum distance of the newly added point to any other node.
 *        If the distance of the newly added node is greater than this min distance,
 *        the add operation will fail and the method will change nothing.
 * @param _nodeIsActiveAddress pointer to the begining of vector nodeIsActive of SceNodes
 * @param _nodeXPosAddress pointer to the begining of vector nodeLocX of SceNodes
 * @param _nodeYPosAddress pointer to the begining of vector nodeLocY of SceNodes
 */
struct AddPtOp: thrust::unary_function<BoolUIDDUID, BoolUID> {
	uint _maxNodeOfOneCell;
	double _addNodeDistance;
	double _minDistanceToOtherNode;
	bool* _nodeIsActiveAddress;
	double* _nodeXPosAddress;
	double* _nodeYPosAddress;
	double _growThreshold;

	unsigned int m_seed;

	__host__ __device__
	AddPtOp(uint maxNodeOfOneCell, double addNodeDistance,
			double minDistanceToOtherNode, bool* nodeIsActiveAddress,
			double* nodeXPosAddress, double* nodeYPosAddress, uint seed,
			double growThreshold) :
			_maxNodeOfOneCell(maxNodeOfOneCell), _addNodeDistance(
					addNodeDistance), _minDistanceToOtherNode(
					minDistanceToOtherNode), _nodeIsActiveAddress(
					nodeIsActiveAddress), _nodeXPosAddress(nodeXPosAddress), _nodeYPosAddress(
					nodeYPosAddress), _growThreshold(growThreshold), m_seed(
					seed) {
	}
	__host__ __device__
	BoolUID operator()(const BoolUIDDUID &biddi) {
		const double pI = acos(-1.0);
		bool isScheduledToGrow = thrust::get<0>(biddi);
		uint activeNodeCountOfThisCell = thrust::get<1>(biddi);
		double lastCheckPoint = thrust::get<5>(biddi);
		if (isScheduledToGrow == false) {
			return thrust::make_tuple(isScheduledToGrow,
					activeNodeCountOfThisCell, lastCheckPoint);
		}
		double cellCenterXCoord = thrust::get<2>(biddi);
		double cellCenterYCoord = thrust::get<3>(biddi);
		uint cellRank = thrust::get<4>(biddi);

		bool isSuccess = true;

		thrust::default_random_engine rng(m_seed);

		// discard n numbers to avoid correlation
		rng.discard(cellRank);

		thrust::uniform_real_distribution<double> u0Pi(0, 2.0 * pI);
		double randomAngle = u0Pi(rng);
		double xOffset = _addNodeDistance * cos(randomAngle);
		double yOffset = _addNodeDistance * sin(randomAngle);
		double xCoordNewPt = cellCenterXCoord + xOffset;
		double yCoordNewPt = cellCenterYCoord + yOffset;
		uint cellNodeStartPos = cellRank * _maxNodeOfOneCell;
		uint cellNodeEndPos = cellNodeStartPos + activeNodeCountOfThisCell;
		for (uint i = cellNodeStartPos; i < cellNodeEndPos; i++) {
			double distance = sqrt(
					(xCoordNewPt - _nodeXPosAddress[i])
							* (xCoordNewPt - _nodeXPosAddress[i])
							+ (yCoordNewPt - _nodeYPosAddress[i])
									* (yCoordNewPt - _nodeYPosAddress[i]));
			if (distance < _minDistanceToOtherNode) {
				isSuccess = false;
				break;
			}
		}

		if (isSuccess) {
			_nodeXPosAddress[cellNodeEndPos] = xCoordNewPt;
			_nodeYPosAddress[cellNodeEndPos] = yCoordNewPt;
			isScheduledToGrow = false;
			activeNodeCountOfThisCell = activeNodeCountOfThisCell + 1;
			lastCheckPoint = lastCheckPoint + _growThreshold;
			if (lastCheckPoint > 1.0) {
				lastCheckPoint = 1.0;
			}
		}
		return thrust::make_tuple(isScheduledToGrow, activeNodeCountOfThisCell,
				lastCheckPoint);
	}

};

/**
 * Compute the target length of a cell given growth progress.
 * @param _cellInitLength initial length of a cell. (when growth progress = 0)
 * @param _cellFinalLength final length of a cell. (when growth progress =1)
 * @param input1 progress cell growth progress.
 * @return cell expected length
 */
struct CompuTarLen: thrust::unary_function<double, double> {
	double _cellInitLength, _cellFinalLength;__host__ __device__
	CompuTarLen(double initLen, double finalLen) :
			_cellInitLength(initLen), _cellFinalLength(finalLen) {
	}
	__host__ __device__
	double operator()(const double &progress) {
		return _cellInitLength + progress * (_cellFinalLength - _cellInitLength);
	}
};

/**
 * Compute the distance of a node to its corresponding center, return 0 if node is inactive.
 * @param input1 x component of center position of cell center
 * @param input2 y component of center position of cell center
 * @param input3 x component of cell growth direction
 * @param input4 y component of cell growth direction
 * @param input5 x component of node location
 * @param input6 y component of node location
 * @param input7 flag for node activeness
 */
struct CompuDist: thrust::unary_function<CVec6Bool, double> {
	__host__ __device__
	double operator()(const CVec6Bool &vec6b) {
		double centerXPos = thrust::get<0>(vec6b);
		double centerYPos = thrust::get<1>(vec6b);
		double growthXDir = thrust::get<2>(vec6b);
		double growthYDir = thrust::get<3>(vec6b);
		double nodeXPos = thrust::get<4>(vec6b);
		double nodeYPos = thrust::get<5>(vec6b);
		bool nodeIsActive = thrust::get<6>(vec6b);
		if (nodeIsActive == false) {
			// this is not true. but those nodes that are inactive will be omitted.
			// I choose 0 because 0 will not be either maximum or minimum
			return 0;
		} else {
			double dirModule = sqrt(
					growthXDir * growthXDir + growthYDir * growthYDir);
			return ((nodeXPos - centerXPos) * (growthXDir)
					+ (nodeYPos - centerYPos) * growthYDir) / dirModule;
		}
	}
};

/**
 * Compute difference of cell expected length and current length.
 * @param input1 expected length of the cell
 * @param input2 minimum distance of nodes of the cell to its corresponding center along growth direction
 * @param input3 maximum distance of nodes of the cell to its corresponding center along growth direction
 * @return difference of expected and current length.
 */
struct CompuDiff: thrust::unary_function<CVec3, double> {
	__host__ __device__
	double operator()(const CVec3 &vec3) {
		double expectedLen = thrust::get<0>(vec3);
		// minimum distance of node to its corresponding center along growth direction
		double minDistance = thrust::get<1>(vec3);
		double maxDistance = thrust::get<2>(vec3);
		return (expectedLen - (maxDistance - minDistance));
	}
};

/**
 * Apply stretch force to all cell nodes.
 * @param _elongationCoefficient elongationForce = _elongationCoefficient*distInElongationDirection
 * 			* elongateDirection;
 * @param input1 distToCenterAlongGrowDir distance of a node to the corresponding cell center along growth direction
 * @param input2 lengthDifference length difference of the expected length of a cell and currentl length of the same cell.
 * @param input3 x component of growth direction.
 * @param input4 y component of growth direction.
 * @param input5 x direction of original velocity
 * @param input6 y direction of original velocity
 */
struct ApplyStretchForce: thrust::unary_function<CVec6, CVec2> {
	double _elongationCoefficient;__host__ __device__
	ApplyStretchForce(double elongationCoefficient) :
			_elongationCoefficient(elongationCoefficient) {
	}
	__host__ __device__
	CVec2 operator()(const CVec6 &vec6) {
		double distToCenterAlongGrowDir = thrust::get<0>(vec6);
		// minimum distance of node to its corresponding center along growth direction
		double lengthDifference = thrust::get<1>(vec6);
		double growthXDir = thrust::get<2>(vec6);
		double growthYDir = thrust::get<3>(vec6);
		double originalVelX = thrust::get<4>(vec6);
		double originalVelY = thrust::get<5>(vec6);
		double xRes = lengthDifference * _elongationCoefficient
				* distToCenterAlongGrowDir * growthXDir;
		xRes = xRes + originalVelX;
		double yRes = lengthDifference * _elongationCoefficient
				* distToCenterAlongGrowDir * growthYDir;
		yRes = yRes + originalVelY;
		return thrust::make_tuple(xRes, yRes);
	}
};

struct ApplyChemoVel: thrust::unary_function<CVec5, CVec2> {
	double _chemoCoefficient;__host__ __device__
	ApplyChemoVel(double chemoCoefficient) :
			_chemoCoefficient(chemoCoefficient) {
	}
	__host__ __device__
	CVec2 operator()(const CVec5 &vec5) {
		double growthSpeed = thrust::get<0>(vec5);
		double growthXDir = thrust::get<1>(vec5);
		double growthYDir = thrust::get<2>(vec5);
		double originalVelX = thrust::get<3>(vec5);
		double originalVelY = thrust::get<4>(vec5);
		double xRes = growthSpeed * _chemoCoefficient * growthXDir;
		xRes = xRes + originalVelX;
		double yRes = growthSpeed * _chemoCoefficient * growthYDir;
		yRes = yRes + originalVelY;
		return thrust::make_tuple(xRes, yRes);
	}
};

/**
 * compute the left shifted global position of a node.
 * @param _shiftLeftOffset number of spaces the node should left shift \n
 * @param input original global position of a node \n
 * @return output shifted global position of a node.\n
 */
struct LeftShiftFunctor: thrust::unary_function<uint, uint> {
	uint _shiftLeftOffset;__host__ __device__
	LeftShiftFunctor(uint maxNodeOfOneCell) :
			_shiftLeftOffset(maxNodeOfOneCell / 2) {
	}
	__host__ __device__
	uint operator()(const uint &position) {
		uint result;
		if (position < _shiftLeftOffset) {
			// could be 0, because these region will actually never be used
			result = 0;
		} else {
			result = position - _shiftLeftOffset;
		}
		return result;
	}
};

/**
 * decide if a node, given by its global rank, is on the right side of a cell.
 * @param _maxNodeCountPerCell maximum number of nodes per cell \n
 * @param _halfMaxNode half of maximum number of nodes per cell \n
 * @param nodeGlobalRank global rank of a node \n
 * @return IsRightSide : true if is on the left side.\n
 */
struct IsRightSide: thrust::unary_function<uint, bool> {
	uint _maxNodeCountPerCell;
	uint _halfMaxNode;__host__ __device__
	IsRightSide(uint maxNodeOfOneCell) :
			_maxNodeCountPerCell(maxNodeOfOneCell), _halfMaxNode(
					maxNodeOfOneCell / 2) {
	}
	__host__ __device__
	bool operator()(const uint &position) {
		if (position % _maxNodeCountPerCell < _halfMaxNode) {
			return false;
		} else {
			return true;
		}
	}
};

/**
 * decide if a node, given by its global rank, is on the left side of a cell.
 * @param _maxNodeCountPerCell maximum number of nodes per cell \n
 * @param _halfMaxNode half of maximum number of nodes per cell \n
 * @param nodeGlobalRank global rank of a node \n
 * @return IsLeftSide : true if is on the left side.\n
 */
struct IsLeftSide: thrust::unary_function<uint, bool> {
	uint _maxNodeCountPerCell;
	uint _halfMaxNode;__host__ __device__
	IsLeftSide(uint maxNodeOfOneCell) :
			_maxNodeCountPerCell(maxNodeOfOneCell), _halfMaxNode(
					maxNodeOfOneCell / 2) {
	}
	__host__ __device__
	bool operator()(const uint &position) {
		if (position % _maxNodeCountPerCell < _halfMaxNode) {
			return true;
		} else {
			return false;
		}
	}
};

/**
 * Given rank of a node inside a cell and rank of the cell, get the global rank of the node.
 * @param _maxNodeCountPerCell maximum number of nodes of a cell
 * @param vec first input: rank of a node inside a cell. \n
 * second input: rank of the cell \n
 * @return nodePosition global rank of a node
 */
struct CompuPos: thrust::unary_function<Tuint2, uint> {
	uint _maxNodeCountPerCell;__host__ __device__
	CompuPos(uint maxNodeOfOneCell) :
			_maxNodeCountPerCell(maxNodeOfOneCell) {
	}
	__host__ __device__
	uint operator()(const Tuint2 &vec) {
		uint rankInCell = thrust::get<0>(vec) % _maxNodeCountPerCell;
		uint cellRank = thrust::get<1>(vec);
		return (cellRank * _maxNodeCountPerCell + rankInCell);
	}
};

/**
 * struct for decide if a cell is ready to divide.
 * @param _isDivideCriticalRatio If the length difference to expected length
 *     is less than this critical ratio and growth progress is equal or bigger than 1.0
 *     it means the cell is ready to divide.
 * @param vec first input : length difference of current length and expected length \n
 * second input: expected length \n
 * thrid input: growth progress. should be 0.0 to 1.0. \n
 * @return isGoingToDivide result that indicates whether a cell is ready to divide.
 */
struct CompuIsDivide: thrust::unary_function<CVec3Int, bool> {
	uint _isDivideCriticalRatio;
	uint _maxNodePerCell;__host__ __device__
	CompuIsDivide(double isDivideCriticalRatio, uint maxNodePerCell) :
			_isDivideCriticalRatio(isDivideCriticalRatio), _maxNodePerCell(
					maxNodePerCell) {
	}
	__host__ __device__
	uint operator()(const CVec3Int &vec) {
		double lengthDifference = thrust::get<0>(vec);
		double expectedLength = thrust::get<1>(vec);
		double currentLength = expectedLength - lengthDifference;
		double growthProgress = thrust::get<2>(vec);
		int nodeCount = thrust::get<3>(vec);
		if (currentLength / expectedLength > _isDivideCriticalRatio
				&& growthProgress >= 1.0 && nodeCount == _maxNodePerCell) {
			return true;
		} else {
			return false;
		}
	}
};

/**
 * Functor for modify veolcities of all nodes given node type and isActive
 * @param input1: take velocity , type and isActive info of node
 * @return output: modified velocity.
 */
struct VelocityModifier: public thrust::unary_function<Vel2DActiveTypeRank,
		CVec2> {
	uint beginPosOfProfileNodes;
	uint currentActiveProfileNodes;__host__ __device__
	VelocityModifier(uint beginPos, uint currentProfileNodeCount) {
		beginPosOfProfileNodes = beginPos;
		currentActiveProfileNodes = currentProfileNodeCount;
	}
	__host__ __device__
	CVec2 operator()(const Vel2DActiveTypeRank &nodeInfo) {
		double velX = thrust::get<0>(nodeInfo);
		double velY = thrust::get<1>(nodeInfo);
		bool isActive = thrust::get<2>(nodeInfo);
		CellType type = thrust::get<3>(nodeInfo);
		uint nodeRank = thrust::get<4>(nodeInfo);
		// boundary nodes should not move. Also, inactive nodes should not move.
		if (type == Boundary || isActive == false) {
			return thrust::make_tuple(0.0, 0.0);
		}
		// The end nodes of the profile should be fixed.
		if (type == Profile
				&& (nodeRank == beginPosOfProfileNodes
						|| nodeRank
								== (beginPosOfProfileNodes
										+ currentActiveProfileNodes - 1))) {
			return thrust::make_tuple(0.0, 0.0);
		}
		return thrust::make_tuple(velX, velY);
	}
};

/**
 * Modified implementation of subcellular element.
 * Important component to process cell growth and division.
 * Boundary, epithilum layer and ECM are also treated as cells but computed seperately.
 * @param beginPosOfBdry   represents begining position of boundary.
 * @param maxNodeOfBdry    represents maximum number of nodes of boundary.
 * @param beginPosOfEpi    represents begining position of epithilum layer.
 * @param maxNodeOfEpi     represents maximum number of nodes of epithilum layer.
 * @param maxNodeOfECM     represents maximum number of nodes per ECM
 * @param beginPosOfECM    represents begining position of ECM.
 * @param maxECMCount      represents maximum number of ECM.
 * @param maxNodeOfOneCell represents maximum number of nodes per cell
 * @param beginPosOfCells  represents begining position of cells.
 * @param maxCellCount     represents maximum number of cells.
 * @param isDivideCrticalRatio If the current cell length divide
 * CellFinalLength is larger than this ratio and the cell growth progress
 *  is complete then we set cell ready to divide
 */
class SceCells_M {
public:

	uint beginPosOfBdry;     // represents begining position of boundary.
	uint maxNodeOfBdry;      // represents maximum number of nodes of boundary.

	uint beginPosOfEpi; // represents begining position of epithilum layer (in cell perspective)
	uint beginPosOfEpiNode; // represents begining position of epithilum layer (in node perspective)
	uint maxNodeOfEpi; // represents maximum number of nodes of epithilum layer.

	uint maxNodeOfECM;       // represents maximum number of nodes per ECM
	uint beginPosOfECM; // represents begining position of ECM (in cell perspective)
	uint beginPosOfECMNode; // represents begining position of ECM (in node perspective)
	uint maxECMCount;        // represents maximum number of ECM.

	uint beginPosOfCells; // represents begining position of cells (in cell perspective)
	uint beginPosOfCellsNode; // represents begining position of cells (in node perspective)
	uint maxNodeOfOneCell;   // represents maximum number of nodes per cell
	uint maxCellCount;       // represents maximum number of cells in the system

	uint maxCellCountAll;    // represents maximum number of all cells
							 // (including pesudo-cells) in the system

	uint maxTotalNodeCountCellOnly; // max possible node count for cell representation
	uint currentActiveCellCount; // the number of active cells would keep changing
	uint currentActiveECMCount;     // the number of active ECM might change.
	uint currectActiveEpiNode;    // the number of epithilum nodes might change.

	// if growthProgress[i] - lastCheckPoint[i] > growThreshold then isScheduledToGrow[i] = true;
	double growThreshold;

	double isDivideCriticalRatio;

	double addNodeDistance;
	double minDistanceToOtherNode;
	double cellInitLength;
	double cellFinalLength;
	double elongationCoefficient;
	double chemoCoefficient;

	SceNodes* nodes;

	/**
	 * @param growthProgress is a vector of size maxCellCount.
	 * In each cell, \n
	 * progress == 0 means recently divided
	 * progress == 1 means ready to divide
	 */
	thrust::device_vector<double> growthProgress;
	/**
	 *
	 */
	thrust::device_vector<double> expectedLength;
	thrust::device_vector<double> currentLength;
	thrust::device_vector<double> lengthDifference;
	// array to hold temp value computed in growth phase.
	// obtained by smallest value of vector product of (a cell node to the cell center)
	// and (growth direction). This will be a negative value
	thrust::device_vector<double> smallestDistance;
	// biggest value of vector product of (a cell node to the cell center)
	// and (growth direction). This will be a positive value
	thrust::device_vector<double> biggestDistance;
	thrust::device_vector<uint> activeNodeCountOfThisCell;
	thrust::device_vector<double> lastCheckPoint;
	thrust::device_vector<bool> isDivided;
	// This cell type array should be initialized together with the host class.
	thrust::device_vector<CellType> cellTypes;
	thrust::device_vector<bool> isScheduledToGrow;
	thrust::device_vector<double> centerCoordX;
	thrust::device_vector<double> centerCoordY;
	thrust::device_vector<double> centerCoordZ;

	thrust::device_vector<double> growthSpeed;
	thrust::device_vector<double> growthXDir;
	thrust::device_vector<double> growthYDir;

	// some memory for holding intermediate values instead of dynamically allocating.
	thrust::device_vector<uint> cellRanksTmpStorage;

	// these tmp coordinates will be temporary storage for division info
	// their size will be the same with maximum node count.
	//thrust::device_vector<double> xCoordTmp;
	//thrust::device_vector<double> yCoordTmp;
	//thrust::device_vector<double> zCoordTmp;

	// some memory for holding intermediate values instead of dynamically allocating.
	thrust::device_vector<uint> cellRanks;
	thrust::device_vector<double> activeXPoss;
	thrust::device_vector<double> activeYPoss;
	thrust::device_vector<double> activeZPoss;

	// temp value for holding a node direction to its corresponding center
	thrust::device_vector<double> distToCenterAlongGrowDir;

	// ************************* these parameters are used for cell growth **************************
	double* growthFactorMagAddress;
	double* growthFactorDirXAddress;
	double* growthFactorDirYAddress;

	// obtain pointer address for second region
	double* growthFactorMagAddress2;
	double* growthFactorDirXAddress2;
	double* growthFactorDirYAddress2;

	uint totalNodeCountForActiveCells;

	bool* nodeIsActiveAddress;
	double* nodeXPosAddress;
	double* nodeYPosAddress;

	double dt;
	// ************************* these parameters are used for cell growth **************************

	// ************************ these parameters are used for cell division *************************
	// sum all bool values which indicate whether the cell is going to divide.
	// toBeDivideCount is the total number of cells going to divide.
	uint toBeDivideCount;
	uint nodeStorageCount;
	thrust::device_vector<bool> tmpIsActiveHold1;
	thrust::device_vector<double> tmpDistToCenter1;
	thrust::device_vector<uint> tmpCellRankHold1;
	thrust::device_vector<double> tmpXValueHold1;
	thrust::device_vector<double> tmpYValueHold1;
	thrust::device_vector<double> tmpZValueHold1;

	thrust::device_vector<bool> tmpIsActiveHold2;
	thrust::device_vector<double> tmpDistToCenter2;
	thrust::device_vector<double> tmpXValueHold2;
	thrust::device_vector<double> tmpYValueHold2;
	thrust::device_vector<double> tmpZValueHold2;

	thrust::device_vector<CellType> tmpCellTypes;
	// ************************ these parameters are used for cell division *************************

	// in this class, there are a lot of arrays that store information for each cell
	// this counting iterator is used for all these arrays indicating the begining.
	thrust::counting_iterator<uint> countingBegin;
	thrust::constant_iterator<uint> initCellCount;
	thrust::constant_iterator<double> initGrowthProgress;

	SceCells_M(SceNodes* nodesInput);
	SceCells_M() {
	}

	void distributeIsActiveInfo();

	void distributeBdryIsActiveInfo();
	void distributeProfileIsActiveInfo();
	void distributeECMIsActiveInfo();
	void distributeCellIsActiveInfo();

	void distributeIsCellRank();
	void runAllCellLevelLogics(double dt, GrowthDistriMap &region1,
			GrowthDistriMap &region2);
	void allComponentsMove();
	void grow2DTwoRegions(double dt, GrowthDistriMap &region1,
			GrowthDistriMap &region2);

	void computeCenterPos();
	void processDivisionInfoAndAddNewCells();
	void divide2DSimplified();

	/**
	 * we need to copy the growth information from grid for chemical to cell nodes.
	 */
	void copyGrowInfoFromGridToCells(GrowthDistriMap &region1,
			GrowthDistriMap &region2);

	/**
	 * Use the growth magnitude and dt to update growthProgress.
	 */
	void updateGrowthProgress();

	/**
	 * Decide if the cells are going to add a node or not.
	 * Use lastCheckPoint and growthProgress to decide whether add point or not
	 */
	void decideIsScheduleToGrow();

	/**
	 * Calculate target length of cell given the cell growth progress.
	 * length is along the growth direction.
	 */
	void computeCellTargetLength();

	/**
	 * Compute distance of each node to its corresponding cell center.
	 * The distantce could be either positive or negative, depending on the pre-defined
	 * growth direction.
	 */
	void computeDistToCellCenter();

	/**
	 * For nodes of each cell, find the maximum and minimum distance to the center.
	 * We will then calculate the current length of a cell along its growth direction
	 * using max and min distance to the center.
	 */
	void findMinAndMaxDistToCenter();

	/**
	 * Compute the difference for cells between their expected length and current length.
	 */
	void computeLenDiffExpCur();

	/**
	 * Use the difference that just computed and growthXDir&growthYDir
	 * to apply stretching force (velocity) on nodes of all cells
	 */
	void stretchCellGivenLenDiff();

	/**
	 * This is just an attempt. Cells move according to chemicals.
	 */
	void cellChemotaxis();

	/**
	 * Adjust the velocities of nodes.
	 * For example, velocity of boundary nodes must be zero.
	 */
	void adjustNodeVel();

	/**
	 * Move nodes according to the velocity we just adjusted.
	 */
	void moveNodes();

	/**
	 * Add a point to a cell if it is scheduled to grow.
	 * This step does not guarantee success ; If adding new point failed, it will not change
	 * isScheduleToGrow and activeNodeCount;
	 */
	void addPointIfScheduledToGrow();

	thrust::device_vector<CellType> getCellTypes() const {
		return cellTypes;
	}

	/**
	 * 2D version of cell division.
	 * Division process is done by creating two temporary vectors to hold the node information
	 * that are going to divide.
	 *
	 * step 1: based on lengthDifference, expectedLength and growthProgress,
	 *     this process determines whether a certain cell is ready to divide and then assign
	 *     a boolean value to isDivided.
	 *
	 * step 2. copy those cells that will divide in to the temp vectors created
	 *
	 * step 3. For each cell in the temp vectors, we sort its nodes by its distance to the
	 * corresponding cell center.
	 * This step is not very effcient when the number of cells going to divide is big.
	 * but this is unlikely to happen because cells will divide according to external chemical signaling
	 * and each will have different divide progress.
	 *
	 * step 4. copy the right part of each cell of the sorted array (temp1) to left part of each cell of
	 * another array
	 *
	 * step 5. transform isActive vector of both temp1 and temp2, making only left part of each cell active.
	 *
	 * step 6. insert temp2 to the end of the cell array
	 *
	 * step 7. copy temp1 to the previous position of the cell array.
	 *
	 * step 8. add activeCellCount of the system.
	 *
	 * step 9. mark isDivide of all cells to false.
	 */

	void decideIfGoingToDivide();

	void copyCellsPreDivision();

	void sortNodesAccordingToDist();

	void copyLeftAndRightToSeperateArrays();

	void transformIsActiveArrayOfBothArrays();

	void addSecondArrayToCellArray();

	void copyFirstArrayToPreviousPos();

	void updateActiveCellCount();

	void markIsDivideFalse();

	/**
	 * This is different from the default set function.
	 */
	void setCellTypes(thrust::device_vector<CellType> cellTypesInput) {
		thrust::copy(cellTypesInput.begin(), cellTypesInput.end(),
				cellTypes.begin());
	}

	uint getMaxTotalNodeCountCellOnly() const {
		return maxTotalNodeCountCellOnly;
	}

	void setMaxTotalNodeCountCellOnly(uint maxTotalNodeCountCellOnly) {
		this->maxTotalNodeCountCellOnly = maxTotalNodeCountCellOnly;
	}
};

#endif /* SCECELLS_H_ */
