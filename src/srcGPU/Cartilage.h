/*
 * Cartilage.h
 *
 *  Created on: Sep 29, 2014
 *      Author: wenzhao
 */

#include "GeoVector.h"
#include "SceNodes.h"

#ifndef CARTILAGE_H_
#define CARTILAGE_H_

__device__
void crossProduct(double &ax, double &ay, double &az, double &bx, double &by,
		double &bz, double &cx, double &cy, double &cz);

__device__
double crossProduct2D(double& ax, double& ay, double& bx, double& by);

__device__
double dotProduct(double &ax, double &ay, double &az, double &bx, double &by,
		double &bz);

/**
 * vector based method to determine if a node is close enough to cartilage.
 */
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
__device__
bool closeToCart(double &xPos, double &yPos);

/**
 * if a node is adhered to cartilage. we need to adjust the velocity.
 */
// comment prevents bad formatting issues of __host__ and __device__ in Nsight
__device__
void modifyVelNoSlip(double &xPos, double &yPos, double &xVel, double &yVel);

struct TorqueCompute: public thrust::binary_function<CVec2Bool, CVec2, double> {
	double _fixedPtX, _fixedPtY, _growthDirX, _growthDirY, _normalX, _normalY;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	TorqueCompute(double fixedPtX, double fixePtY, double growthDirX,
			double growthDirY) :
			_fixedPtX(fixedPtX), _fixedPtY(fixePtY) {
		double length = sqrt(growthDirX * growthDirX + growthDirY * growthDirY);
		_growthDirX = growthDirX / length;
		_growthDirY = growthDirY / length;

		// normal direction is counter-clockwise rotation for 90 degrees.
		_normalX = -_growthDirY;
		_normalY = _growthDirX;
	}
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	double operator()(const CVec2Bool &loc, const CVec2 &force) {
		double posX = thrust::get<0>(loc);
		double posY = thrust::get<1>(loc);
		bool isActive = thrust::get<2>(loc);
		if (!isActive) {
			return 0;
		}
		double dirX = posX - _fixedPtX;
		double dirY = posY - _fixedPtY;
		double length = dirX * _growthDirX + dirY * _growthDirY;

		double forceDirX = thrust::get<0>(force);
		double forceDirY = thrust::get<1>(force);
		double perpForce = forceDirX * _normalX + forceDirY * _normalY;
		double torque = length * perpForce;

		return torque;
	}
};

struct Rotation2D: public thrust::unary_function<CVec2Bool, CVec2> {
	double _fixedPtX, _fixedPtY, _angle;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	Rotation2D(double fixedPtX, double fixedPtY, double angle) :
			_fixedPtX(fixedPtX), _fixedPtY(fixedPtY), _angle(angle) {
	}
	__host__ __device__
	CVec2 operator()(const CVec2Bool &loc) {
		double posX = thrust::get<0>(loc);
		double posY = thrust::get<1>(loc);
		bool isActive = thrust::get<2>(loc);
		if (!isActive) {
			thrust::make_tuple(posX, posY);
		}

		double dirToFixX = posX - _fixedPtX;
		double dirToFixY = posY - _fixedPtY;
		double dirToFixXNew = dirToFixX * cos(_angle) - dirToFixY * sin(_angle);
		double dirToFixYNew = dirToFixX * sin(_angle) + dirToFixY * cos(_angle);
		double xPosNew = _fixedPtX + dirToFixXNew;
		double yPosNew = _fixedPtY + dirToFixYNew;
		return thrust::make_tuple(xPosNew, yPosNew);
	}
};

struct NoSlipHandler: public thrust::unary_function<CVec4Bool, CVec2Bool> {
	__host__ __device__
	NoSlipHandler() {

	}
	__device__
	CVec2Bool operator()(const CVec4Bool &locVel) {
		double xPos = thrust::get<0>(locVel);
		double yPos = thrust::get<1>(locVel);
		double xVel = thrust::get<2>(locVel);
		double yVel = thrust::get<3>(locVel);
		bool isActive = thrust::get<4>(locVel);

		if (!isActive) {
			return thrust::make_tuple(0.0, 0.0, false);
		}
		bool isCloseToCart = false;
		if (closeToCart(xPos, yPos)) {
			modifyVelNoSlip(xPos, yPos, xVel, yVel);
			isCloseToCart = true;
		}
		return thrust::make_tuple(xVel, yVel, isCloseToCart);
	}
};

struct MinCount: public thrust::binary_function<BoolInt, BoolInt, int> {
	__host__ __device__
	MinCount() {

	}
	__device__
	int operator()(const BoolInt &i1, const BoolInt &i2) {
		bool bool1 = thrust::get<0>(i1);
		int index1 = thrust::get<1>(i1);
		bool bool2 = thrust::get<0>(i2);
		int index2 = thrust::get<1>(i2);
		if (bool1 & bool2) {
			return min(index1, index2);
		} else if ((!bool1) && (!bool2)) {
			return INT_MAX;
		} else if (bool1) {
			return index1;
		} else {
			return index2;
		}
	}
};

struct MaxCount: public thrust::binary_function<BoolInt, BoolInt, int> {
	__host__ __device__
	MaxCount() {

	}
	__device__
	int operator()(const BoolInt &i1, const BoolInt &i2) {
		bool bool1 = thrust::get<0>(i1);
		int index1 = thrust::get<1>(i1);
		bool bool2 = thrust::get<0>(i2);
		int index2 = thrust::get<1>(i2);
		if (bool1 & bool2) {
			return max(index1, index2);
		} else if ((!bool1) && (!bool2)) {
			return INT_MIN;
		} else if (bool1) {
			return index1;
		} else {
			return index2;
		}
	}
};

__device__
double calDist(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2);

struct ComputeLinkLen: public thrust::unary_function<uint, double> {
	double* _nodeLocXLinkBeginAddress;
	double* _nodeLocYLinkBeginAddress;
	double* _nodeLocZLinkBeginAddress;

	__host__ __device__
	ComputeLinkLen(double* nodeLocXLinkBeginAddress,
			double* nodeLocYLinkBeginAddress, double* nodeLocZLinkBeginAddress) :
			_nodeLocXLinkBeginAddress(nodeLocXLinkBeginAddress), _nodeLocYLinkBeginAddress(
					nodeLocYLinkBeginAddress), _nodeLocZLinkBeginAddress(
					nodeLocZLinkBeginAddress) {
	}

	__device__
	double operator()(const uint &index) {
		double xLoc = _nodeLocXLinkBeginAddress[index];
		double yLoc = _nodeLocYLinkBeginAddress[index];
		double zLoc = _nodeLocZLinkBeginAddress[index];

		double xLocN = _nodeLocXLinkBeginAddress[index + 1];
		double yLocN = _nodeLocYLinkBeginAddress[index + 1];
		double zLocN = _nodeLocZLinkBeginAddress[index + 1];

		return calDist(xLoc, yLoc, zLoc, xLocN, yLocN, zLocN);
	}
};

class Cartilage {
	bool isInitialized;
	bool isParaInitialized;
	SceNodes* nodes;
	CartPara cartPara;

	/**
	 * the first two will be mainly for the two ends.
	 */
	//std::vector<uint> tipNodeIndicies;
	void calculateGrowthDir();
	void addPoint1(CVector &nodeBehindPos);
	void addPoint2(CVector &nodeBehindPos);
	void grow1(double dt);
	void grow2(double dt);
	void runGrowthLogics1(double dt);
	void runGrowthLogics2(double dt);

	/**
	 * TODO
	 * really bad implementation for now.
	 */
	void updateTipNodes(double dt);
	void runGrowthLogics(double dt);
	void calculateTotalTorque();
	void move(double dt);

	/**
	 * distributes is active information to nodes.
	 */
	void initIsActive();

	void readValuesFromConfig();

	/**
	 * Per conversation with Dr.Newman, we learned that cartilage would not slip over the epithelium,
	 * This function was created to handle this special condition.
	 * Basic idea is to limit the velocity of the epithelium nodes that are contacting cartilage.
	 * velocity that is along the cartilage growth direction will not be changed.
	 * velocity that is perpendicular to cartilage growth direction will be changed to
	 * accomendate cartilage angular speed.
	 */
	void handleCartNoSlippage();

public:
	Cartilage();
	void initializeMem(SceNodes* nodeInput);

	void runAllLogics(double dt);

	const CartPara& getCartPara() const {
		return cartPara;
	}

	void setCartPara(const CartPara& cartPara) {
		this->cartPara = cartPara;
		isParaInitialized = 1;
	}
};

#endif /* CARTILAGE_H_ */
