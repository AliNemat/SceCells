#include "Cartilage.h"

void Cartilage::calculateGrowthDir() {
	// first need to obtain the location of first node
	// finding index of node 1 in node location vector.
	uint globalIndex1 = cartPara.growNode1Index
			+ nodes->getAllocPara().startPosCart;
	// update node location 1 given index 1
	cartPara.growNode1 = CVector(nodes->getInfoVecs().nodeLocX[globalIndex1],
			nodes->getInfoVecs().nodeLocY[globalIndex1],
			nodes->getInfoVecs().nodeLocZ[globalIndex1]);

	// finding index of node 2 in node location vector.
	uint globalIndex2 = cartPara.growNode2Index
			+ nodes->getAllocPara().startPosCart;
	// update node location 2 given index 2
	cartPara.growNode2 = CVector(nodes->getInfoVecs().nodeLocX[globalIndex2],
			nodes->getInfoVecs().nodeLocY[globalIndex2],
			nodes->getInfoVecs().nodeLocZ[globalIndex2]);

	// compute position of midpoint.
	CVector midPoint = cartPara.growNode1 + cartPara.growNode2;
	midPoint = midPoint / 2;

	// compute current growth direction.
	CVector direction = midPoint - cartPara.fixedPt;
	cartPara.GrowthDir = direction.getUnitVector();
}

void Cartilage::addPoint() {
	//
}

void Cartilage::grow(double dt) {
	//

}

void Cartilage::calculateTotalTorque() {

}

void Cartilage::move(double dt) {
}

Cartilage::Cartilage() {
	nodes = NULL;
}

void Cartilage::initialize(SceNodes* nodeInput) {

}

void Cartilage::runAllLogics(double dt) {

}
