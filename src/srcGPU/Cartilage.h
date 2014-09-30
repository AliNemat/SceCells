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

/**
 * parameters for cartilage.
 */
class CartPara {
public:
	CVector fixedPt;
	CVector GrowthDir;

	CVector growNode1;
	CVector growNode2;

	uint growNode1Index;
	uint growNode2Index;

	double growthSpeedNode1;
	double growthSpeedNode2;
	/**
	 * moment of intertia.
	 */
	double moInertia;
	//uint activeCartilageNodeCount;

	double totalGrowthSinceLastAdd;
	double growthThreshold;
};

class Cartilage {
	SceNodes* nodes;
	CartPara cartPara;

	std::vector<uint> tipIndicies;

	void calculateGrowthDir();
	void addPoint();
	void grow(double dt);
	void calculateTotalTorque();
	void move(double dt);
public:
	Cartilage();
	void initialize(SceNodes* nodeInput);
	void runAllLogics(double dt);
};

#endif /* CARTILAGE_H_ */
