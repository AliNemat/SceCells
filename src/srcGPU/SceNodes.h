#ifndef SCENODES_H_
#define SCENODES_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/copy.h>
#include <thrust/sort.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/tabulate.h>
#include <thrust/binary_search.h>
#include <thrust/tuple.h>
#include <thrust/random.h>
#include <thrust/equal.h>
#include <thrust/inner_product.h>
#include <thrust/transform_reduce.h>

#include <thrust/gather.h>
#include <thrust/scan.h>
#include <thrust/fill.h>
#include <thrust/count.h>
#include <thrust/unique.h>
#include <thrust/extrema.h>

#include <thrust/for_each.h>
#include <cuda_runtime.h>
#include "ConfigParser.h"
#include <assert.h>
#include "commonData.h"
#include "ResAnalysisHelper.h"

#include <ctime>

// include google test files in order to test private functions.
#include "gtest/gtest_prod.h"
#include "SceNodes.h"
//#include "SceCells.h"

// I wish I could include some c++ 11 data structure here but it seems
// Thrust is not compatible with c++ 11.
// #include <unordered_map>
// TODO: maybe it's because of the -ansi compiler option instead of Thrust?

/**
 * @mainpage Subcellular Element Method for development biology
 *
 * @author Wenzhao Sun wsun2@nd.edu
 *
 *  This simulation package is highly dependent on Thrust library, which is an official
 *  fast - GPU development toolkit released by NVIDIA research.
 *  The code requires CMake to build. When making project, html style documentation will be
 *  automatically generated in "html" folder in project root directory.
 *
 *  1) The goal of the project is to simulate a developmental biology phenomenon
 *  that the shape of the beak is changed from normal triangular shape to
 *  a hooked shape after being treated by VPA, which will support our hypothesis
 *  that the final shape of the chicken beak is determined by the mechanical interaction between cells.
 *
 *  2) There are five essential parts in our simulation model -
 *  FNM bud, MX bud, ECM, epithelium, and external chemical, which is FGF signal.
 *
 *  3) The approximate time scale of our simulation is several days. (day 8 - day 11)?
 */

typedef thrust::tuple<double, double> CVec2;
typedef thrust::tuple<double, double, bool> CVec2Bool;
typedef thrust::tuple<double, double, int> CVec2Int; //Ali 
typedef thrust::tuple<double, uint> DUi;
typedef thrust::tuple<double, uint, double, double> DUiDD;
typedef thrust::tuple<double, uint, double, double,double> DUiDDD;  //Ali 
typedef thrust::tuple<bool, double> BoolD;
typedef thrust::tuple<double, double, double, uint> DDDUi;//AAMIRI
typedef thrust::tuple<bool, int> BoolInt;
typedef thrust::tuple<uint, bool> UiB;
typedef thrust::tuple<bool, uint, double> BoolUID;
typedef thrust::tuple<bool, int, uint, double, double> BoolIUiDD;
typedef thrust::tuple<bool, uint, double, double> BoolUiDD; //Ali
typedef thrust::tuple<bool, uint, double, double, uint, double> BoolUIDDUID;
typedef thrust::tuple<bool, uint, double, double, uint, uint, bool, double> BoolUIDDUIUIBoolD;//AAMIRI
typedef thrust::tuple<uint, uint, bool, double> UiUiBD;//AAMIRI
typedef thrust::tuple<double, double, double> CVec3;
typedef thrust::tuple<double, double, double, uint> CVec3Int;
typedef thrust::tuple<double, double, double, bool> CVec3Bool;
typedef thrust::tuple<double, double, double, bool, uint> CVec3BoolInt;
typedef thrust::tuple<ECellType,double, double, double, bool, uint> TypeCVec3BoolInt; //Ali
typedef thrust::tuple<double, double, double, double> CVec4;
typedef thrust::tuple<double, double, double, double, bool> CVec4Bool;
typedef thrust::tuple<double, double, double, double, double> CVec5;
typedef thrust::tuple<double, double, double, double, double, double> CVec6;
typedef thrust::tuple<double, double, double, double, double, double, double,
		double, double, double> CVec10;
typedef thrust::tuple<double, double, double, double, double, double, bool> CVec6Bool;
typedef thrust::tuple<double, double, double, double, double, double, uint> CVec6UI;
typedef thrust::tuple<int, int> Int2;
typedef thrust::tuple<uint, uint> Tuint2;
typedef thrust::tuple<uint, uint, uint> Tuint3;
typedef thrust::tuple<uint, uint, uint, double, double, double> Tuuuddd;
typedef thrust::tuple<uint, uint, uint, int, int> TuuuII; // Ali
typedef thrust::tuple<uint, uint, uint, double, double> Tuuudd;
typedef thrust::tuple<uint, uint, uint, double, double,int ,int > TuuuddII;


struct SubApicalInfoEachCell{

	int nodeIdFront[10];
	int nodeIdBehind[10] ; 

	SubApicalInfoEachCell() {
		for (int i=0 ; i<10 ; i++) {
			nodeIdFront[i]= 0 ; 
			nodeIdBehind[i]= 0 ; 
		}
	}
}; 


// special datatype required for Thrust minmax element function.
typedef thrust::pair<thrust::device_vector<int>::iterator,
		thrust::device_vector<int>::iterator> MinMaxRes;

__device__
bool bothInternal(uint nodeGlobalRank1, uint nodeGlobalRank2);

__device__
bool bothMembr(uint nodeGlobalRank1, uint nodeGlobalRank2);

__device__
bool bothMembrDiffCell(uint nodeGlobalRank1, uint nodeGlobalRank2);

__device__
bool sameCellMemIntnl(uint nodeGlobalRank1, uint nodeGlobalRank2);

//Ali
__device__
bool Is_Lennard_Jones();
//Ali
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

struct NanCount: public thrust::binary_function<double, double, CVec3> {
	__device__
	int operator()(const double& num) {
		if (isnan(num)) {
			return 1;
		} else {
			return 0;
		}
	}
};

/**
 * Functor predicate see if a boolean varible is true(seems unnecessary but still required).
 */
struct ActiveAndIntnl {
	__device__
	bool operator()(const thrust::tuple<bool, SceNodeType> &bt) {
		bool isActive = thrust::get<0>(bt);
		SceNodeType type = thrust::get<1>(bt);
		if (isActive == true && type == CellIntnl) {
			return true;
		} else {
			return false;
		}
	}
};

//Ali 
struct ActiveAndApical {
	__device__
	bool operator()(const thrust::tuple<bool, MembraneType1> &bm) {
		bool isActive = thrust::get<0>(bm);
		MembraneType1  type = thrust::get<1>(bm);
		if (isActive == true && type == apical1) {
			return true;
		} else {
			return false;
		}
	}
};


/**
 * Functor predicate see if a boolean varible is true(seems unnecessary but still required).
 */
struct ActiveAndMembr {
	__device__
	bool operator()(const thrust::tuple<bool, SceNodeType> &bt) {
		bool isActive = thrust::get<0>(bt);
		SceNodeType type = thrust::get<1>(bt);
		if (isActive == true && type == CellMembr) {
			return true;
		} else {
			return false;
		}
	}
};

/**
 * functor to add two three dimensional vectors.
 */
struct AddFunctor: public thrust::binary_function<CVec3, CVec3, CVec3> {
	double _dt;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__ AddFunctor(double dt) :
			_dt(dt) {
	}
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__ CVec3 operator()(const CVec3 &vel, const CVec3 &loc) {
		double xMoveDist = thrust::get<0>(vel) * _dt;
		double yMoveDist = thrust::get<1>(vel) * _dt;
		double zMoveDist = thrust::get<2>(vel) * _dt;
		double xPosAfterMove = xMoveDist + thrust::get<0>(loc);
		double yPosAfterMove = yMoveDist + thrust::get<1>(loc);
		double zPosAfterMove = zMoveDist + thrust::get<2>(loc);
		return thrust::make_tuple(xPosAfterMove, yPosAfterMove, zPosAfterMove);
	}
};

struct AdjustAdh: public thrust::unary_function<Int2, int> {
	int* _adhArrAddr;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	AdjustAdh(int* adhArrAddr) :
			_adhArrAddr(adhArrAddr) {
	}
	__host__ __device__
	int operator()(const Int2& myInt2) const {
		int nodeIndx = thrust::get<0>(myInt2);
		int adhIndx = thrust::get<1>(myInt2);
		if (adhIndx != -1) {
			if (_adhArrAddr[adhIndx] != nodeIndx) {
				return -1;
			} else {
				return adhIndx;
			}
		} else {
			return -1;
		}
	}
};

/**
 * random number generation engine.
 */
struct Prg {
	double a, b;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__ Prg(double _a = 0.f, double _b = 1.f) :
			a(_a), b(_b) {
	}
	__host__ __device__
	double operator()(const unsigned int n) const {
		thrust::default_random_engine rng;
		thrust::uniform_real_distribution<double> dist(a, b);
		rng.discard(n);
		return dist(rng);
	}
};

/**
 * map a node coordinate to its bucket index.
 *
 */
struct pointToBucketIndex2D: public thrust::unary_function<CVec3BoolInt, Tuint2> {
	double _minX;
	double _maxX;
	double _minY;
	double _maxY;

	double _bucketSize;
	unsigned int width;

	__host__ __device__ pointToBucketIndex2D(double minX, double maxX,
			double minY, double maxY, double bucketSize) :
			_minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY), _bucketSize(
					bucketSize), width((maxX - minX) / bucketSize + 1) {
	}

	__host__ __device__ Tuint2 operator()(const CVec3BoolInt& v) const {
		// find the raster indices of p's bucket
		if (thrust::get<3>(v) == true) {
			unsigned int x = static_cast<unsigned int>((thrust::get<0>(v)
					- _minX) / _bucketSize);
			unsigned int y = static_cast<unsigned int>((thrust::get<1>(v)
					- _minY) / _bucketSize);
			// return the bucket's linear index and node's global index
			return thrust::make_tuple(y * width + x, thrust::get<4>(v));
		} else {
			// return UINT_MAX to indicate the node is inactive and its value should not be used
			return thrust::make_tuple(UINT_MAX, UINT_MAX);
		}
	}
};

struct BucketIndexer3D: public thrust::unary_function<CVec3BoolInt, Tuint2> {
	double _minX;
	double _maxX;
	double _minY;
	double _maxY;
	double _minZ;
	double _maxZ;

	double _unitLen;
	uint _XSize;
	uint _YSize;

	__host__ __device__ BucketIndexer3D(double minX, double maxX, double minY,
			double maxY, double minZ, double maxZ, double unitLen) :
			_minX(minX), _maxX(maxX), _minY(minY), _maxY(maxY), _minZ(minZ), _maxZ(
					maxZ), _unitLen(unitLen), _XSize(
					(maxX - minX) / _unitLen + 1), _YSize(
					(maxY - minY) / _unitLen + 1) {
	}

	__host__ __device__ Tuint2 operator()(const CVec3BoolInt& v) const {
		// find the raster indices of p's bucket
		if (thrust::get<3>(v) == true) {
			uint x = static_cast<uint>((thrust::get<0>(v) - _minX) / _unitLen);
			uint y = static_cast<uint>((thrust::get<1>(v) - _minY) / _unitLen);
			uint z = static_cast<uint>((thrust::get<2>(v) - _minZ) / _unitLen);
			// return the bucket's linear index and node's global index
			return thrust::make_tuple(z * _XSize * _YSize + y * _XSize + x,
					thrust::get<4>(v));
		} else {
			// return UINT_MAX to indicate the node is inactive and its value should not be used
			return thrust::make_tuple(UINT_MAX, UINT_MAX);
		}
	}
};

/**
 * Functor to compute neighbor buckets(center bucket included) of a node.
 * @param input1 bucket index of node
 * @param input2 pick from the sequence, which is also global rank of the node
 *
 * @return output1 bucket indices of node ( all neighbors and the center bucket) of node
 * @return output2 global rank of node
 */
struct NeighborFunctor2D: public thrust::unary_function<Tuint2, Tuint2> {
	uint _numOfBucketsInXDim;
	uint _numOfBucketsInYDim;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__ NeighborFunctor2D(uint numOfBucketsInXDim,
			uint numOfBucketsInYDim) :
			_numOfBucketsInXDim(numOfBucketsInXDim), _numOfBucketsInYDim(
					numOfBucketsInYDim) {
	}
	__host__ __device__ Tuint2 operator()(const Tuint2 &v) {
		uint relativeRank = thrust::get<1>(v) % 9;
		uint xPos = thrust::get<0>(v) % _numOfBucketsInXDim;
		uint yPos = thrust::get<0>(v) / _numOfBucketsInXDim;
		switch (relativeRank) {
		case 0:
			return thrust::make_tuple(thrust::get<0>(v), thrust::get<1>(v));
		case 1:
			if (xPos > 0 && yPos > 0) {
				uint topLeft = (xPos - 1) + (yPos - 1) * _numOfBucketsInXDim;
				return thrust::make_tuple(topLeft, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 2:
			if (yPos > 0) {
				uint top = xPos + (yPos - 1) * _numOfBucketsInXDim;
				return thrust::make_tuple(top, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 3:
			if (yPos > 0 && xPos < _numOfBucketsInXDim - 1) {
				uint topRight = xPos + 1 + (yPos - 1) * _numOfBucketsInXDim;
				return thrust::make_tuple(topRight, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 4:
			if (xPos < _numOfBucketsInXDim - 1) {
				uint right = xPos + 1 + yPos * _numOfBucketsInXDim;
				return thrust::make_tuple(right, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 5:
			if (xPos < _numOfBucketsInXDim - 1
					&& yPos < _numOfBucketsInYDim - 1) {
				uint bottomRight = xPos + 1 + (yPos + 1) * _numOfBucketsInXDim;
				return thrust::make_tuple(bottomRight, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 6:
			if (yPos < _numOfBucketsInYDim - 1) {
				uint bottom = xPos + (yPos + 1) * _numOfBucketsInXDim;
				return thrust::make_tuple(bottom, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 7:
			if (xPos > 0 && yPos < _numOfBucketsInYDim - 1) {
				uint bottomLeft = (xPos - 1) + (yPos + 1) * _numOfBucketsInXDim;
				return thrust::make_tuple(bottomLeft, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		case 8:
			if (xPos > 0) {
				uint left = xPos - 1 + yPos * _numOfBucketsInXDim;
				return thrust::make_tuple(left, thrust::get<1>(v));
			} else {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		default:
			return thrust::make_tuple(UINT_MAX, UINT_MAX);
		}
	}
};

struct NgbrFunc3D: public thrust::unary_function<Tuint2, Tuint2> {
	uint _XSize;
	uint _YSize;
	uint _ZSize;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	NgbrFunc3D(uint XSize, uint YSize, uint ZSize) :
			_XSize(XSize), _YSize(YSize), _ZSize(ZSize) {
	}
	__host__ __device__ Tuint2 operator()(const Tuint2 &v) {
		uint relativeRank = thrust::get<1>(v) % 27;
		uint XYSize = _XSize * _YSize;
		uint zPos = thrust::get<0>(v) / XYSize;
		uint xyPos = thrust::get<0>(v) % XYSize;
		uint xPos = xyPos % _XSize;
		uint yPos = xyPos / _XSize;

		int rankZ = relativeRank / 9;
		int rankTmp = relativeRank % 9;
		int rankY = rankTmp / 3;
		int rankX = rankTmp % 3;

		rankX--;
		rankY--;
		rankZ--;

		if (rankX == -1) {
			// try to find left neighbor, but already in leftmost position
			if (xPos == 0) {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		} else if (rankX == 1) {
			// try to find right neighbor, but already in rightmost position
			if (xPos == _XSize - 1) {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		}
		if (rankY == -1) {
			if (yPos == 0) {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		} else if (rankY == 1) {
			if (yPos == _YSize - 1) {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		}
		if (rankZ == -1) {
			if (zPos == 0) {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		} else if (rankZ == 1) {
			if (zPos == _ZSize - 1) {
				return thrust::make_tuple(UINT_MAX, thrust::get<1>(v));
			}
		}

		uint ngbrPos = (xPos + rankX) + (yPos + rankY) * _XSize
				+ (zPos + rankZ) * XYSize;
		return thrust::make_tuple(ngbrPos, thrust::get<1>(v));
	}
};

__device__
double computeDist(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2);

__device__
double computeDist2D(double &xPos, double &yPos, double &xPos2, double &yPos2);

__device__
void calAndAddInter_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes);
__device__
void calAndAddInter_M2(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes);

__device__
void calAndAddIntraB_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& xRes, double& yRes);

__device__
void calAndAddIntraDiv_M(double& xPos, double& yPos, double& xPos2,
		double& yPos2, double& growPro, double& xRes, double& yRes);

__device__
void calculateAndAddInterForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes);

__device__
void calculateAndAddIntraForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes);

__device__
void calAndAddIntraForceDiv(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &growPro,
		double &xRes, double &yRes, double &zRes);

__device__
void calculateAndAddInterForceDiffType(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes);
__device__ bool bothNodesCellNode(uint nodeGlobalRank1, uint nodeGlobalRank2,
		uint cellNodesThreshold);

__device__ bool isSameCell(uint nodeGlobalRank1, uint nodeGlobalRank2);

__device__ bool ofSameType(uint cellType1, uint cellType2);

__device__
void handleSceForceNodesBasic(uint &nodeRank1, uint &nodeRank2, double &xPos,
		double &yPos, double &zPos, double &xPos2, double &yPos2, double &zPos2,
		double &xRes, double &yRes, double &zRes, double* _nodeLocXAddress,
		double* _nodeLocYAddress, double* _nodeLocZAddress);

__device__
void handleSceForceNodesDisc(uint& nodeRank1, uint& nodeRank2, double& xPos,
		double& yPos, double& zPos, double& xPos2, double& yPos2, double& zPos2,
		double& xRes, double& yRes, double& zRes, double& interForceX,
		double& interForceY, double& interForceZ, double* _nodeLocXAddress,
		double* _nodeLocYAddress, double* _nodeLocZAddress,
		double* _nodeGrowProAddr);

__device__
void handleSceForceNodesDisc_M(uint& nodeRank1, uint& nodeRank2, double& xPos,
		double& yPos, double& xPos2, double& yPos2, double& xRes, double& yRes,
		double* _nodeLocXAddress, double* _nodeLocYAddress,
		double* _nodeGrowProAddr);

__device__
void handleAdhesionForce_M(int& adhereIndex, double& xPos, double& yPos,
		double& curAdherePosX, double& curAdherePosY, double& xRes,
		double& yRes, double& alpha);

//Ali for adhesion reaction force
__device__
void handleAdhesionForce_M2(double& xPos, double& yPos,
		double& curAdherePosX, double& curAdherePosY, double& xRes,
		double& yRes, double& alpha);

__device__
double getMitoticAdhCoef(double& growProg, double& growProgNeigh);

__device__
void attemptToAdhere(bool& isSuccess, uint& index, double& dist,
		uint& nodeRank2, double& xPos1, double& yPos1, double& xPos2,
		double& yPos2);

__device__
void calculateForceBetweenLinkNodes(double &xLoc, double &yLoc, double &zLoc,
		double &xLocLeft, double &yLocLeft, double &zLocLeft, double &xLocRight,
		double &yLocRight, double &zLocRight, double &xVel, double &yVel,
		double &zVel);

/**
 * A compute functor designed for computing SceForce for Disc project.
 */
struct AddSceForceDisc: public thrust::unary_function<Tuuuddd, CVec6> {
	uint* _extendedValuesAddress;
	double* _nodeLocXAddress;
	double* _nodeLocYAddress;
	double* _nodeLocZAddress;
	double* _nodeGroProAddr;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	AddSceForceDisc(uint* valueAddress, double* nodeLocXAddress,
			double* nodeLocYAddress, double* nodeLocZAddress,
			double* nodeGrowProAddr) :
			_extendedValuesAddress(valueAddress), _nodeLocXAddress(
					nodeLocXAddress), _nodeLocYAddress(nodeLocYAddress), _nodeLocZAddress(
					nodeLocZAddress), _nodeGroProAddr(nodeGrowProAddr) {
	}
	__device__
	CVec6 operator()(const Tuuuddd &u3d3) const {
		double xRes = 0.0;
		double yRes = 0.0;
		double zRes = 0.0;

		double interForceX = 0.0;
		double interForceY = 0.0;
		double interForceZ = 0.0;

		uint begin = thrust::get<0>(u3d3);
		uint end = thrust::get<1>(u3d3);
		uint myValue = thrust::get<2>(u3d3);
		double xPos = thrust::get<3>(u3d3);
		double yPos = thrust::get<4>(u3d3);
		double zPos = thrust::get<5>(u3d3);

		for (uint i = begin; i < end; i++) {
			uint nodeRankOfOtherNode = _extendedValuesAddress[i];
			if (nodeRankOfOtherNode == myValue) {
				continue;
			}
			handleSceForceNodesDisc(myValue, nodeRankOfOtherNode, xPos, yPos,
					zPos, _nodeLocXAddress[nodeRankOfOtherNode],
					_nodeLocYAddress[nodeRankOfOtherNode],
					_nodeLocZAddress[nodeRankOfOtherNode], xRes, yRes, zRes,
					interForceX, interForceY, interForceZ, _nodeLocXAddress,
					_nodeLocYAddress, _nodeLocZAddress, _nodeGroProAddr);
		}
		return thrust::make_tuple(xRes, yRes, zRes, interForceX, interForceY,
				interForceZ);
	}
};

/**
 * A compute functor designed for computing SceForce for Disc project.
 */
struct AddForceDisc_M: public thrust::unary_function<Tuuudd, CVec2> {
	uint* _extendedValuesAddress;
	double* _nodeLocXAddress;
	double* _nodeLocYAddress;
	int* _nodeAdhereIndex;
	int* _membrIntnlIndex;
	double* _nodeGroProAddr;
	bool _adhNotSet ; 
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	AddForceDisc_M(uint* valueAddress, double* nodeLocXAddress,
			double* nodeLocYAddress, int* nodeAdhereIndex, int* membrIntnlIndex,
			double* nodeGrowProAddr,bool adhNotSet) :
			_extendedValuesAddress(valueAddress), _nodeLocXAddress(
					nodeLocXAddress), _nodeLocYAddress(nodeLocYAddress), _nodeAdhereIndex(
					nodeAdhereIndex), _membrIntnlIndex(membrIntnlIndex), _nodeGroProAddr(
					nodeGrowProAddr),_adhNotSet(adhNotSet) {
	}
	__device__
	CVec2 operator()(const Tuuudd &u3d2) const {
		double xRes = 0.0;
		double yRes = 0.0;

		uint begin = thrust::get<0>(u3d2);
		uint end = thrust::get<1>(u3d2);
		uint myValue = thrust::get<2>(u3d2);
		double xPos = thrust::get<3>(u3d2);
		double yPos = thrust::get<4>(u3d2);

		bool isSuccess = false;
		uint index;
		double dist;
                bool  Lennard_Jones =Is_Lennard_Jones() ;
//		if (_adhNotSet){
	//	_nodeAdhereIndex[myValue] = -1 ;  Ali commented to deactive this part of the code
//		}
		for (uint i = begin; i < end; i++) {
			uint nodeRankOther = _extendedValuesAddress[i];
			if (nodeRankOther == myValue) {
				continue;
			}
			if (bothMembrDiffCell(myValue, nodeRankOther)) {
                            if (Lennard_Jones) {
				calAndAddInter_M2(xPos, yPos, _nodeLocXAddress[nodeRankOther],
						_nodeLocYAddress[nodeRankOther], xRes, yRes);
                                               }
                            else               {    
				calAndAddInter_M(xPos, yPos, _nodeLocXAddress[nodeRankOther],
						_nodeLocYAddress[nodeRankOther], xRes, yRes);
                                               }
		//	if(_adhNotSet){
				//if (_nodeAdhereIndex[myValue] == -1) {
	//				attemptToAdhere(isSuccess, index, dist, nodeRankOther, xPos,
	//						yPos, _nodeLocXAddress[nodeRankOther],
	//						_nodeLocYAddress[nodeRankOther]);
			//	}
//Ali

		//	}
			}
		}
	//	if (_adhNotSet) {
	//		if (isSuccess) {
		//			_nodeAdhereIndex[myValue] = index;  //Ali commented to deactive this part of the code
		//		_nodeAdhereIndex[index] = myValue; Ali added and then commentted out
	//		}
	//	}
		return thrust::make_tuple(xRes, yRes);
	}
};

/**
 * a complicated data structure for adding subcellular element force to cell nodes.
 * This data structure is designed in an unconventional way because of performance considerations.
 */
struct AddSceForceBasic: public thrust::unary_function<Tuuuddd, CVec3> {
	uint* _extendedValuesAddress;
	double* _nodeLocXAddress;
	double* _nodeLocYAddress;
	double* _nodeLocZAddress;
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	AddSceForceBasic(uint* valueAddress, double* nodeLocXAddress,
			double* nodeLocYAddress, double* nodeLocZAddress) :
			_extendedValuesAddress(valueAddress), _nodeLocXAddress(
					nodeLocXAddress), _nodeLocYAddress(nodeLocYAddress), _nodeLocZAddress(
					nodeLocZAddress) {
	}
	__device__ CVec3 operator()(const Tuuuddd &u3d3) const {
		double xRes = 0.0;
		double yRes = 0.0;
		double zRes = 0.0;

		uint begin = thrust::get<0>(u3d3);
		uint end = thrust::get<1>(u3d3);
		uint myValue = thrust::get<2>(u3d3);
		double xPos = thrust::get<3>(u3d3);
		double yPos = thrust::get<4>(u3d3);
		double zPos = thrust::get<5>(u3d3);

		for (uint i = begin; i < end; i++) {
			uint nodeRankOfOtherNode = _extendedValuesAddress[i];
			if (nodeRankOfOtherNode == myValue) {
				continue;
			}

			handleSceForceNodesBasic(myValue, nodeRankOfOtherNode, xPos, yPos,
					zPos, _nodeLocXAddress[nodeRankOfOtherNode],
					_nodeLocYAddress[nodeRankOfOtherNode],
					_nodeLocZAddress[nodeRankOfOtherNode], xRes, yRes, zRes,
					_nodeLocXAddress, _nodeLocYAddress, _nodeLocZAddress);

		}

		return thrust::make_tuple(xRes, yRes, zRes);
	}
};

struct ApplyAdh: public thrust::unary_function<BoolIUiDD, CVec2> {
	double* _nodeLocXArrAddr;
	double* _nodeLocYArrAddr;
	double* _nodeGrowProAddr;
	int* _nodeAdhAddr ; 
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	ApplyAdh(double* nodeLocXArrAddr, double* nodeLocYArrAddr, double* nodeGrowProAddr, int* nodeAdhAddr  ) :
			_nodeLocXArrAddr(nodeLocXArrAddr), _nodeLocYArrAddr(nodeLocYArrAddr), _nodeGrowProAddr(nodeGrowProAddr), _nodeAdhAddr(nodeAdhAddr) {
	}
	__device__
	CVec2 operator()(const BoolIUiDD& adhInput) const {
		bool isActive = thrust::get<0>(adhInput);
		int adhIndx = thrust::get<1>(adhInput);
		uint nodeIndx = thrust::get<2>(adhInput);
		double oriVelX = thrust::get<3>(adhInput);
		double oriVelY = thrust::get<4>(adhInput);
		double growProg = _nodeGrowProAddr[nodeIndx];
		double growProgNeigh = _nodeGrowProAddr[adhIndx];
		//bool adhSkipped = false;	
		double alpha = getMitoticAdhCoef(growProg, growProgNeigh);//to adjust the mitotic values of stiffness
		/*int maxNodePerCell=680  ; 
		int cellRank=nodeIndx/maxNodePerCell ;
		int nodeRank=nodeIndx%maxNodePerCell ;
		int activeMembCount=600 ; 

		int indexLeft=nodeRank-1;
		if (indexLeft==-1){
			indexLeft=activeMembCount-1 ;
		}
		indexLeft=indexLeft+cellRank*maxNodePerCell ;
		int indexRight=nodeRank+1 ;
		if (indexRight==activeMembCount){
			indexRight=0 ; 
		}
		indexRight=indexRight+cellRank*maxNodePerCell ;
		*/
		if (adhIndx == -1 || !isActive) {
			return thrust::make_tuple(oriVelX, oriVelY);
		} 
		else {	//else if  (_nodeAdhAddr[indexLeft]==-1 || _nodeAdhAddr[indexRight]==-1){ 
				double locX = _nodeLocXArrAddr[nodeIndx];
				double locY = _nodeLocYArrAddr[nodeIndx];
				double adhLocX = _nodeLocXArrAddr[adhIndx];
				double adhLocY = _nodeLocYArrAddr[adhIndx];
				handleAdhesionForce_M(adhIndx, locX, locY, adhLocX, adhLocY,
					oriVelX, oriVelY, alpha);
				return thrust::make_tuple(oriVelX, oriVelY);
		}
	//	else { 
	//		return thrust::make_tuple(oriVelX, oriVelY);
	//	}
	}
};

/**
 * calculate force in epithilum.
 */
/*
struct ApplyAdhReaction: public thrust::unary_function<BoolUiDD, CVec2> {
	double* _nodeLocXArrAddr;
	double* _nodeLocYArrAddr;
	double* _nodeGrowProAddr;
	int* _nodeAdhIndx ; 
	uint _maxTotalNode ; 
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	ApplyAdhReaction(double* nodeLocXArrAddr, double* nodeLocYArrAddr, double* nodeGrowProAddr, int*nodeAdhIndx,uint maxTotalNode) :
			_nodeLocXArrAddr(nodeLocXArrAddr), _nodeLocYArrAddr(nodeLocYArrAddr), _nodeGrowProAddr(nodeGrowProAddr),_nodeAdhIndx(nodeAdhIndx),_maxTotalNode(maxTotalNode) {
	}
	__device__
	CVec2 operator()(const BoolUiDD& adhInput2) const {
		bool isActive = thrust::get<0>(adhInput2);
		uint nodeIndx = thrust::get<1>(adhInput2);
		double oriVelX = thrust::get<2>(adhInput2);
		double oriVelY = thrust::get<3>(adhInput2);
		//bool adhSkipped = false;	
		double sumOriVelX= 0 ;
		double sumOriVelY=0 ; 
		double oriVelXTmp= 0 ;
		double oriVelYTmp= 0 ; 

		double growProg = _nodeGrowProAddr[nodeIndx];
		double locX = _nodeLocXArrAddr[nodeIndx];
		double locY = _nodeLocYArrAddr[nodeIndx];

		if ( isActive==false) {
			return thrust::make_tuple(oriVelX, oriVelY);
		}
		else {
			for (int i=0 ; i<_maxTotalNode; i++) {
				if (_nodeAdhIndx[i] ==nodeIndx) {
					double growProgNeigh = _nodeGrowProAddr[i];
					double alpha = getMitoticAdhCoef(growProg, growProgNeigh);//to adjust the mitotic values of stiffness
					double adhLocX = _nodeLocXArrAddr[i];
					double adhLocY = _nodeLocYArrAddr[i];
					handleAdhesionForce_M2(locX, locY, adhLocX, adhLocY,oriVelXTmp, oriVelYTmp, alpha);
					sumOriVelX=sumOriVelX+oriVelXTmp ; 
					sumOriVelY=sumOriVelY+oriVelYTmp ;
					if (_nodeAdhIndx[nodeIndx]==-1){
						_nodeAdhIndx[nodeIndx]=i ; 
					}
					 
				}
			}
			return thrust::make_tuple(sumOriVelX+oriVelX, sumOriVelY+oriVelY);

		}
	}
};

*/

struct AddLinkForces: public thrust::unary_function<uint, CVec3> {
	double* _nodeLocXLinkBeginAddress;
	double* _nodeLocYLinkBeginAddress;
	double* _nodeLocZLinkBeginAddress;
	double* _nodeVelXLinkBeginAddress;
	double* _nodeVelYLinkBeginAddress;
	double* _nodeVelZLinkBeginAddress;
	uint _maxNumberOfNodes;
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__ AddLinkForces(double* nodeLocXLinkBeginAddress,
			double* nodeLocYLinkBeginAddress, double* nodeLocZLinkBeginAddress,
			double* nodeVelXLinkBeginAddress, double* nodeVelYLinkBeginAddress,
			double* nodeVelZLinkBeginAddress, uint maxNumberOfNodes) :
			_nodeLocXLinkBeginAddress(nodeLocXLinkBeginAddress), _nodeLocYLinkBeginAddress(
					nodeLocYLinkBeginAddress), _nodeLocZLinkBeginAddress(
					nodeLocZLinkBeginAddress), _nodeVelXLinkBeginAddress(
					nodeVelXLinkBeginAddress), _nodeVelYLinkBeginAddress(
					nodeVelYLinkBeginAddress), _nodeVelZLinkBeginAddress(
					nodeVelZLinkBeginAddress), _maxNumberOfNodes(
					maxNumberOfNodes) {
	}
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__device__ CVec3 operator()(const uint &pos) const {
		// nodes on two ends should be fixed.
		// therefore, velocity is set to 0.
		if (pos == 0 || pos == _maxNumberOfNodes - 1) {
			return thrust::make_tuple(0.0, 0.0, 0.0);
		}

		uint left = pos - 1;
		uint right = pos + 1;

		double xLoc = _nodeLocXLinkBeginAddress[pos];
		double yLoc = _nodeLocYLinkBeginAddress[pos];
		double zLoc = _nodeLocZLinkBeginAddress[pos];

		double xLocLeft = _nodeLocXLinkBeginAddress[left];
		double yLocLeft = _nodeLocYLinkBeginAddress[left];
		double zLocLeft = _nodeLocZLinkBeginAddress[left];

		double xLocRight = _nodeLocXLinkBeginAddress[right];
		double yLocRight = _nodeLocYLinkBeginAddress[right];
		double zLocRight = _nodeLocZLinkBeginAddress[right];

		double xVel = _nodeVelXLinkBeginAddress[pos];
		double yVel = _nodeVelYLinkBeginAddress[pos];
		double zVel = _nodeVelZLinkBeginAddress[pos];

		calculateForceBetweenLinkNodes(xLoc, yLoc, zLoc, xLocLeft, yLocLeft,
				zLocLeft, xLocRight, yLocRight, zLocRight, xVel, yVel, zVel);

		return thrust::make_tuple(xVel, yVel, zVel);
	}
};

/**
 * Store information for each node, including position and velocity.
 */
struct NodeInfoVecs {
public:
// this vector is used to indicate whether a node is active or not.
// E.g, max number of nodes of a cell is 100 maybe the first 75 nodes are active.
// The value of this vector will be changed by external process.
	thrust::device_vector<bool> nodeIsActive;
	thrust::host_vector<bool> nodeIsActiveHost; //Ali
// X locations of nodes
	thrust::device_vector<double> nodeLocX;
	thrust::host_vector<double> nodeLocXHost; // Ali
// Y locations of nodes
	thrust::device_vector<double> nodeLocY;
	thrust::host_vector<double> nodeLocYHost; // Ali
// Z locations of nodes
	thrust::device_vector<double> nodeLocZ;
// X velocities of nodes
	thrust::device_vector<double> nodeVelX;
// Y velocities of nodes
	thrust::device_vector<double> nodeVelY;
// Z velocities of nodes
	thrust::device_vector<double> nodeVelZ;
// For from internal ndoes to each membrane node in x direction
	thrust::device_vector<double> nodeF_MI_M_x; //Ali
// For from internal ndoes to each membrane node in y  direction
	thrust::device_vector<double> nodeF_MI_M_y; //Ali 

// For from internal ndoes to each membrane node in Tangential  direction
	thrust::device_vector<double> nodeF_MI_M_T; //Ali 

// For from internal ndoes to each membrane node in Normal direction
	thrust::device_vector<double> nodeF_MI_M_N; //Ali 
// Tangent to the nodes
	thrust::device_vector<double> nodeVelTangent;//AAMIRI
// Normal to the nodes
	thrust::device_vector<double> nodeVelNormal;//AAMIRI
//Sitffness of the cells representing actin level
	thrust::device_vector<double> nodeActinLevel;//Ali 
// Curvature at the node
	thrust::device_vector<double> nodeCurvature;//AAMIRI

//External forces on nodes in x dir
	thrust::device_vector<double> nodeExtForceX;//AAMIRI

//External forces on nodes in y dir
	thrust::device_vector<double> nodeExtForceY;//AAMIRI

//External forces on nodes in y dir
	thrust::device_vector<double> nodeExtForceTangent;//AAMIRI

//External forces on nodes in y dir
	thrust::device_vector<double> nodeExtForceNormal;//AAMIRI
// is subApical node , for adhesion purpose
	thrust::device_vector<bool> isSubApicalJunction;//Ali 
	thrust::host_vector<bool> isSubApicalJunctionHost;//Ali 

// represents nodes's stress level.
	thrust::device_vector<double> nodeMaxForce;

// growth progress of the cell that the node belongs to.
	thrust::device_vector<double> nodeGrowPro;

	thrust::device_vector<double> nodeInterForceX;

	thrust::device_vector<double> nodeInterForceY;

	thrust::device_vector<double> nodeInterForceZ;

	thrust::host_vector<double> nodeAdhMinDist; // Ali

	thrust::device_vector<MembraneType1> memNodeType1; // Ali
	thrust::host_vector<MembraneType1> memNodeType1Host; // Ali

// in order to represent differential adhesion, we also need an vector
// for each cell node to identify the cell type.
	thrust::device_vector<SceNodeType> nodeCellType;
// for each node, we need to identify which cell it belongs to.
	thrust::device_vector<uint> nodeCellRank;

// only for modified version
	thrust::device_vector<int> nodeAdhereIndex;
	thrust::host_vector<int> nodeAdhIndxHostCopy;
	thrust::host_vector<int> nodeAdhereIndexHost;    //Ali
	thrust::device_vector<int> membrIntnlIndex;

	thrust::device_vector<double> membrTensionMag;
	thrust::device_vector<double> membrTenMagRi;
	thrust::device_vector<double> membrDistToRi;//AAMIRI
	thrust::device_vector<double> membrLinkRiMidX;
	thrust::device_vector<double> membrLinkRiMidY;

	thrust::device_vector<double> membrBendLeftX;
	thrust::device_vector<double> membrBendLeftY;
	thrust::device_vector<double> membrBendRightX;
	thrust::device_vector<double> membrBendRightY;
	//thrust::device_vector<bool> nodeIsBasalMem;//Ali
	//thrust::device_vector<bool> nodeIsLateralMem;//Ali
	thrust::device_vector<int> nodeIsApicalMem;//Ali it only gets 0 and 1
	//thrust::device_vector<bool> nodeIsLateralMemHost;//Ali
	thrust::device_vector<int>  nodeCellRankFront;//Ali it is cell size
	thrust::device_vector<int>  nodeCellRankBehind;//Ali it is cell size
	thrust::device_vector<int>  nodeCellRankFrontOld;//Ali it is cell size
	thrust::device_vector<int>  nodeCellRankBehindOld;//Ali it is cell size
	thrust::host_vector<int>  nodeCellRankFrontHost;//Ali it is cell size
	thrust::host_vector<int>  nodeCellRankBehindHost;//Ali it is cell size
};

/**
 * Store temporary values while computing for Sub-cellular element method.
 */
struct NodeAuxVecs {
// bucket key means which bucket ID does a certain point fit into
	thrust::device_vector<uint> bucketKeys;
// bucket value means global rank of a certain point
	thrust::device_vector<uint> bucketValues;
// bucket key expanded means what are the bucket IDs are the neighbors of a certain point
	thrust::device_vector<uint> bucketKeysExpanded;
// bucket value expanded means each point ( represented by its global rank) will have multiple copies
	thrust::device_vector<uint> bucketValuesIncludingNeighbor;
// begin position of a keys in bucketKeysExpanded and bucketValuesIncludingNeighbor
	thrust::device_vector<uint> keyBegin;
// end position of a keys in bucketKeysExpanded and bucketValuesIncludingNeighbor
	thrust::device_vector<uint> keyEnd;
};

/**
 * Data structure for calculating nearby node interaction.
 * To maximize GPU performance, I choose to implement Structure of Arrays(SOA)
 * instead of Array Of Structure, which is commonly used in serial algorithms.
 * This data structure is a little bit counter-intuitive, but would expect to
 * give 2X - 3X performance gain.
 * To gain a better performance on GPU, we pre-allocate memory for everything that
 * will potentially
 * This data structure need to be initialized with following parameters:
 * 1) maximum number of cells to be simulated in the simulation domain
 * 2) maximum number of nodes inside a cell
 * 3) maximum number of uintra-cellular links per node
 * 4) maximum number of uinter-cellular links per node
 */
class SceNodes {
//	SceCells* cells ;
	bool adhNotSet ; 
	SceDomainPara domainPara;
	SceMechPara mechPara;
	NodeAllocPara allocPara;
	//NodeInfoVecs infoVecs;
	NodeAuxVecs auxVecs;
	ControlPara controlPara;

	NodeAllocPara_M allocPara_M;
	SceMechPara_M mechPara_M;  
	/**
	 * reads domain related parameters.
	 */
	void readDomainPara();

	/**
	 * reads mechanics related parameters.
	 */
	void readMechPara();

	void readParas_M();

	void initNodeAllocPara(uint totalBdryNodeCount, uint maxProfileNodeCount,
			uint maxCartNodeCount, uint maxTotalECMCount, uint maxNodeInECM,
			uint maxTotalCellCount, uint maxNodeInCell);

	void initNodeAllocPara_M(uint totalBdryNodeCount, uint maxTotalCellCount,
			uint maxEpiNodePerCell, uint maxInternalNodePerCell);

	/**
	 * This function copies parameters to GPU constant memory.
	 */
	void copyParaToGPUConstMem();

	void copyParaToGPUConstMem_M();

	void allocSpaceForNodes(uint maxTotalNodeCount, uint maxNumCells, uint currentActiveCellCount); //Ali updated
	/**
	 * this method maps the points to their bucket ID.
	 * writes data to thrust::device_vector<uint> bucketValues and
	 * thrust::device_vector<uint> bucketValues;
	 * each data in bucketValues will represent
	 */
	void buildBuckets2D();

	void buildBuckets2D_M();
	/**
	 * this method extends the previously obtained vector of keys and values to a map
	 * that each point will map to all bucket IDs that are near the specific point.
	 */
	void extendBuckets2D();

	void extendBuckets2D_M();

	/**
	 * this method prepares data for apply Sce forces.
	 */
	void findBucketBounds();

	void findBucketBounds_M();

	/**
	 * Basic Sce method for performance testing purpose.
	 */
	void applySceForcesBasic();

	void applySceForcesBasic_M();

	/**
	 * Sce method specifically for wing disc.
	 */
	void applySceForcesDisc();

	void applySceForcesDisc_M();

	/**
	 * This method outputs a vector of possible neighbor pairs.
	 * Reason why this method exist is that outputting animation frame
	 * is really slow using previous version of animation function.
	 * This new method significantly improves the speed of producing each
	 * animation frame.
	 */
	std::vector<std::pair<uint, uint> > obtainPossibleNeighborPairs();

	void initControlPara(bool isStab);

	void debugNAN();

	void keepAdhIndxCopyInHost_M();
	void processMembrAdh_M();
	void removeInvalidPairs_M();
	void applyMembrAdh_M();

	void copyExtForces_M();//AAMIRI

	uint endIndx_M;
	uint endIndxExt_M;
	uint endIndxExtProc_M;

// friend unit test so these it can test private functions
	FRIEND_TEST(SceNodeTest, FindingPossiblePairsTest);
// friend unit test so these it can test private functions
	FRIEND_TEST(SceNodeTest, outputAnimationLinks);
public:
    
	bool adhUpdate ; //Ali 
	bool isInitPhase ; //Ali 

	NodeInfoVecs infoVecs; //Ali 
	/**
	 * Default constructor -- explicit usage is discouraged.
	 */
	SceNodes();

	/**
	 */
	//SceNodes(uint maxTotalCellCount, uint maxAllNodePerCell);
	SceNodes(uint maxTotalCellCount, uint maxAllNodePerCell, uint currentActiveCellCount);

	/**
	 * recommended constructor for beak project.
	 */
	SceNodes(uint totalBdryNodeCount, uint maxProfileNodeCount,
			uint maxCartNodeCount, uint maxTotalECMCount, uint maxNodeInECM,
			uint maxTotalCellCount, uint maxNodeInCell, bool isStab);

	/**
	 * Override dimension data introduced by config files.
	 */
	void initDimension(double minX, double maxX, double minY, double maxY,
			double bucketSize);

	/**
	 * initialize data fields.
	 */
	void initValues(std::vector<CVector> &initBdryCellNodePos,
			std::vector<CVector> &initProfileNodePos,
			std::vector<CVector> &initCartNodePos,
			std::vector<CVector> &initECMNodePos,
			std::vector<CVector> &initFNMCellNodePos,
			std::vector<CVector> &initMXCellNodePos);

	/**
	 * initialize data fields.
	 */
	void initValues_M(std::vector<bool>& initIsActive,
			std::vector<CVector> &initCellNodePos,
			std::vector<SceNodeType>& nodeTypes);

	/**
	 * this method contains all preparation work for SCE force calculation.
	 */
	void prepareSceForceComputation();

	/**
	 * this method contains all preparation work for SCE force calculation.
	 */
	void prepareSceForceComputation_M();

	void buildBuckets3D();
	void extendBuckets3D();
	void findBucketBounds3D();
	/**
	 * this method contains all preparation work for SCE force calculation
	 * in 3D.
	 */
	void prepareSceForceComputation3D();

	/**
	 * wrap apply forces methods together.
	 * This is for performance testing only.
	 */
	void sceForcesPerfTesting();

	void sceForcesPerfTesting_M();

	/**
	 * wrap apply forces methods together.
	 * This is for Wing Disc development only.
	 */
	void sceForcesDisc();

	void sceForcesDisc_M();

	/**
	 * add maxNodeOfOneCell
	 */
	void addNewlyDividedCells(thrust::device_vector<double> &nodeLocXNewCell,
			thrust::device_vector<double> &nodeLocYNewCell,
			thrust::device_vector<double> &nodeLocZNewCell,
			thrust::device_vector<bool> &nodeIsActiveNewCell,
			thrust::device_vector<SceNodeType> &nodeCellTypeNewCell);

	/**
	 * This method outputs a data structure for animation.
	 */
	VtkAnimationData obtainAnimationData(AnimationCriteria aniCri);

	VtkAnimationData obtainAnimationData_M(AnimationCriteria aniCri);

	/**
	 * method that outputs label matrix.
	 */
	std::vector<std::vector<int> > obtainLabelMatrix(PixelizePara &pixelPara);

	void removeNodes(int cellRank, vector<uint> &removeSeq);

	const SceDomainPara&getDomainPara() const;
	void setDomainPara(const SceDomainPara& domainPara);

	const NodeAllocPara& getAllocPara() const;
	void setAllocPara(const NodeAllocPara& allocPara);

	const NodeAuxVecs& getAuxVecs() const;
	void setAuxVecs(const NodeAuxVecs& auxVecs);

	/**
	 * This getter does not contain "const", meaning the value can be changed inside this getter function.
	 * This setting might not fit code design principles but might be necessary for performance considerations.
	 */
	NodeInfoVecs& getInfoVecs();
	void setInfoVecs(const NodeInfoVecs& infoVecs);

	const SceMechPara& getMechPara() const {
		return mechPara;
	}

	ControlPara getControlPara() const {
		return controlPara;
	}

	void setControlPara(ControlPara controlPara) {
		this->controlPara = controlPara;
	}

	double getMaxEffectiveRange();

	const NodeAllocPara_M& getAllocParaM() const {
		return allocPara_M;
	}

	void setAllocParaM(const NodeAllocPara_M& allocParaM) {
		allocPara_M = allocParaM;
	}

	std::vector<std::pair<uint, uint> > obtainPossibleNeighborPairs_M();

	const SceMechPara_M& getMechParaM() const {
		return mechPara_M;
	}

	void setMechParaM(const SceMechPara_M& mechParaM) {
		mechPara_M = mechParaM;
	}

	void setActiveCellCount(uint activeCellCount) {
		allocPara_M.currentActiveCellCount = activeCellCount;
	}
};

#endif /* SCENODES_H_ */
