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

struct TorqueCompute: public thrust::binary_function<CVec2Bool, CVec2, double> {
	double _fixedPtX, _fixedPtY, _growthDirX, _growthDirY, _normalX, _normalY;
	// comment prevents bad formatting issues of __host__ and __device__ in Nsight
	__host__ __device__
	TorqueCompute(double fixedPtX, double fixePtY, double growthDirX,
			double growthDirY) :
			_fixedPtX(fixedPtX), _fixedPtY(fixePtY), _growthDirX(growthDirX), _growthDirY(
					growthDirY) {
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
		double xPosNew = posX * cos(_angle) - posY * sin(_angle);
		double yPosNew = posX * sin(_angle) + posY * cos(_angle);
		return thrust::make_tuple(xPosNew, yPosNew);
	}
};

class Cartilage {
	SceNodes* nodes;
	CartPara cartPara;

	/**
	 * the first two will be mainly for the two ends.
	 */
	std::vector<uint> tipNodeIndicies;

	void calculateGrowthDir();
	void addPoint1(CVector &nodeBehindPos);
	void addPoint2(CVector &nodeBehindPos);
	void grow1(double dt);
	void grow2(double dt);
	void runGrowthLogics1(double dt);
	void runGrowthLogics2(double dt);
	void updateTipNodes();
	void runGrowthLogics(double dt);
	void calculateTotalTorque();
	void move(double dt);

	/**
	 * Initialize the values using processed information.
	 * Need to change isActive and
	 */
	void initializeNodes(std::vector<CVector> &initFixedNodes,
			std::vector<CVector> &initTipNodes, CVector &initGrowingNode1,
			CVector &initGrowingNode2, CVector& initPivotNode1,
			CVector& initPivotNode2);

public:
	Cartilage();
	void initializeMem(SceNodes* nodeInput);

	/**
	 * Initialize the values.
	 * @param initialNodes contains positions of all initial nodes
	 * @param tipNode1 is the location of node1 on the tip of growth
	 * @param tipNode2 is the location of node2 on the tip of growth
	 * @param pivotNode1 is used for determine the fixed point
	 * @param pivotNode2 is used for determine the fixed point
	 */
	void initializeVal(std::vector<CVector> &initialNodes, CVector &tipNode1,
			CVector &tipNode2, CVector &pivotNode1, CVector pivotNode2);
	void runAllLogics(double dt);
};

#endif /* CARTILAGE_H_ */
