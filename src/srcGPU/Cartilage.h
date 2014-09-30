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

	void calculateGrowthDir();
	void addPoint();
	void grow();
	void calculateTotalTorque();
	void move(double dt);
public:
	void initialize(SceNodes* nodeInput);
	void runAllLogics(double dt);
};

#endif /* CARTILAGE_H_ */
