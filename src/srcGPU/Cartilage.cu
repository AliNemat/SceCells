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
	cartPara.growthDir = direction.getUnitVector();

	cartPara.growthDirNode1 = cartPara.growthDir;
	cartPara.growthDirNode2 = cartPara.growthDir;
}

void Cartilage::calculateTotalTorque() {
	thrust::plus<double> plusOp;
	uint indexBegin = nodes->getAllocPara().startPosCart;
	uint indexEnd = indexBegin + cartPara.nodeIndexEnd;
	cartPara.totalTorque = thrust::inner_product(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin()))
					+ indexBegin,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin()))
					+ indexEnd,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())), 0.0,
			plusOp,
			TorqueCompute(cartPara.fixedPt.x, cartPara.fixedPt.y,
					cartPara.growthDir.x, cartPara.growthDir.y));
}

void Cartilage::move(double dt) {
	uint indexBegin = nodes->getAllocPara().startPosCart;
	uint indexEnd = indexBegin + cartPara.nodeIndexEnd;
	cartPara.angularSpeed = cartPara.totalTorque / cartPara.moInertia;
	// angle is counter-clock wise.
	double angle = cartPara.angularSpeed * dt;
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin()))
					+ indexBegin,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin()))
					+ indexEnd,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin()))
					+ indexBegin,
			Rotation2D(cartPara.fixedPt.x, cartPara.fixedPt.y, angle));

}

Cartilage::Cartilage() {
	nodes = NULL;
}

void Cartilage::initializeMem(SceNodes* nodeInput) {
	nodes = nodeInput;
}

void Cartilage::addPoint1(CVector &nodeBehindPos) {
	CVector node1GrowthDir = cartPara.growNode1 - nodeBehindPos;
// get unit vector
	node1GrowthDir = node1GrowthDir.getModul();
	CVector incrementVec = node1GrowthDir * cartPara.growthThreshold;
	CVector posNew = nodeBehindPos + incrementVec;
	uint indexNew = cartPara.nodeIndexEnd + nodes->getAllocPara().startPosCart;
	nodes->getInfoVecs().nodeLocX[indexNew] = posNew.GetX();
	nodes->getInfoVecs().nodeLocY[indexNew] = posNew.GetY();
	nodes->getInfoVecs().nodeLocZ[indexNew] = posNew.GetZ();
	cartPara.growNodeBehind1Index = cartPara.nodeIndexEnd;
	cartPara.nodeIndexEnd++;
}

void Cartilage::addPoint2(CVector &nodeBehindPos) {
	CVector node2GrowthDir = cartPara.growNode2 - nodeBehindPos;
// get unit vector
	node2GrowthDir = node2GrowthDir.getModul();
	CVector incrementVec = node2GrowthDir * cartPara.growthThreshold;
	CVector posNew = nodeBehindPos + incrementVec;
	uint indexNew = cartPara.nodeIndexEnd + nodes->getAllocPara().startPosCart;
	nodes->getInfoVecs().nodeLocX[indexNew] = posNew.GetX();
	nodes->getInfoVecs().nodeLocY[indexNew] = posNew.GetY();
	nodes->getInfoVecs().nodeLocZ[indexNew] = posNew.GetZ();
	cartPara.growNodeBehind2Index = cartPara.nodeIndexEnd;
	cartPara.nodeIndexEnd++;
}

void Cartilage::grow1(double dt) {
	CVector movementVec = cartPara.growthDirNode1 * cartPara.growthSpeedNode1;
	uint node1GlobalIndex = cartPara.growNode1Index
			+ nodes->getAllocPara().startPosCart;
	nodes->getInfoVecs().nodeLocX[node1GlobalIndex] += movementVec.GetX();
	nodes->getInfoVecs().nodeLocY[node1GlobalIndex] += movementVec.GetY();
	nodes->getInfoVecs().nodeLocZ[node1GlobalIndex] += movementVec.GetZ();
}

void Cartilage::grow2(double dt) {
	CVector movementVec = cartPara.growthDirNode2 * cartPara.growthSpeedNode2;
	uint node2GlobalIndex = cartPara.growNode2Index
			+ nodes->getAllocPara().startPosCart;
	nodes->getInfoVecs().nodeLocX[node2GlobalIndex] += movementVec.GetX();
	nodes->getInfoVecs().nodeLocY[node2GlobalIndex] += movementVec.GetY();
	nodes->getInfoVecs().nodeLocZ[node2GlobalIndex] += movementVec.GetZ();
}

void Cartilage::runGrowthLogics1(double dt) {
	uint behind1Index = cartPara.growNodeBehind1Index
			+ nodes->getAllocPara().startPosCart;
	CVector node1BehindPos = CVector(
			nodes->getInfoVecs().nodeLocX[behind1Index],
			nodes->getInfoVecs().nodeLocY[behind1Index],
			nodes->getInfoVecs().nodeLocZ[behind1Index]);

	CVector node1GrowthDir = cartPara.growNode1 - node1BehindPos;
	double distFromBehind1 = node1GrowthDir.getModul();
	if (distFromBehind1 < cartPara.growthThreshold) {
		grow1(dt);
	} else {
		addPoint1(node1BehindPos);
		grow1(dt);
	}
}

void Cartilage::runGrowthLogics2(double dt) {
	uint behind2Index = cartPara.growNodeBehind2Index
			+ nodes->getAllocPara().startPosCart;
	CVector node2BehindPos = CVector(
			nodes->getInfoVecs().nodeLocX[behind2Index],
			nodes->getInfoVecs().nodeLocY[behind2Index],
			nodes->getInfoVecs().nodeLocZ[behind2Index]);

	CVector node2GrowthDir = cartPara.growNode2 - node2BehindPos;
	double distFromBehind2 = node2GrowthDir.getModul();
	if (distFromBehind2 < cartPara.growthThreshold) {
		grow2(dt);
	} else {
		addPoint2(node2BehindPos);
		grow2(dt);
	}
}

void Cartilage::runGrowthLogics(double dt) {
	runGrowthLogics1(dt);
	runGrowthLogics2(dt);
	updateTipNodes();
}

void Cartilage::initializeNodes(std::vector<CVector>& initFixedNodes,
		std::vector<CVector>& initTipNodes, CVector& initGrowingNode1,
		CVector& initGrowingNode2, CVector& initPivotNode1,
		CVector& initPivotNode2) {
}

void Cartilage::initializeVal(std::vector<CVector>& initialNodes,
		CVector& tipNode1, CVector& tipNode2, CVector &pivotNode1,
		CVector pivotNode2) {
	uint pivotIndex1 = cartPara.pivotNode1Index
			+ nodes->getAllocPara().startPosCart;
	CVector pivotNode1Pos = CVector(nodes->getInfoVecs().nodeLocX[pivotIndex1],
			nodes->getInfoVecs().nodeLocY[pivotIndex1],
			nodes->getInfoVecs().nodeLocZ[pivotIndex1]);

	uint pivotIndex2 = cartPara.pivotNode2Index
			+ nodes->getAllocPara().startPosCart;
	CVector pivotNode2Pos = CVector(nodes->getInfoVecs().nodeLocX[pivotIndex2],
			nodes->getInfoVecs().nodeLocY[pivotIndex2],
			nodes->getInfoVecs().nodeLocZ[pivotIndex2]);

	cartPara.fixedPt = pivotNode1Pos + pivotNode2Pos;
	cartPara.fixedPt = cartPara.fixedPt / 2;
}

void Cartilage::updateTipNodes() {

}

void Cartilage::runAllLogics(double dt) {
	calculateTotalTorque();
	move(dt);
	runGrowthLogics(dt);
}
