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
	CVector growthDir;

	CVector growthDirNode1;
	CVector growthDirNode2;

	CVector growNode1;
	CVector growNode2;

	uint pivotNode1Index;
	uint pivotNode2Index;

	uint growNode1Index;
	uint growNode2Index;

	uint growNodeBehind1Index;
	uint growNodeBehind2Index;

	/**
	 * memory allocation related parameter.
	 * number of spaces allocated for tip nodes.
	 * The are designed to fit the first part of Cart nodes.
	 */
	uint tipNodeIndexTotal;

	/**
	 * memory allocation related parameter.
	 * should be the same with the corresponding parameter in SceNodes.
	 */
	uint totalAvailableIndicies;

	/**
	 * this value changes with the cartilage grows.
	 */
	uint nodeIndexEnd;

	double growthSpeedNode1;
	double growthSpeedNode2;

	/**
	 * moment of intertia.
	 */
	double moInertia;
	//uint activeCartilageNodeCount;

	double totalTorque;
	double angularSpeed;

	double growthThreshold;
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
