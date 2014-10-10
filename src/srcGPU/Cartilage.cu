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

	cartPara.node1GrowthDir = cartPara.growthDir;
	cartPara.node2GrowthDir = cartPara.growthDir;
}

void Cartilage::calculateTotalTorque() {
	thrust::plus<double> plusOp;
	uint indexBegin = nodes->getAllocPara().startPosCart;
	uint indexEnd = indexBegin + cartPara.nodeIndexEnd;
	std::cout << "growth dir = ";
	cartPara.growthDir.Print();
	std::cout << " fixed point = ";
	cartPara.fixedPt.Print();

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
							nodes->getInfoVecs().nodeVelY.begin()))
					+ indexBegin, 0.0, plusOp,
			TorqueCompute(cartPara.fixedPt.x, cartPara.fixedPt.y,
					cartPara.growthDir.x, cartPara.growthDir.y));

	std::cout << " total torque = " << cartPara.totalTorque << std::endl;
}

void Cartilage::move(double dt) {
	uint indexBegin = nodes->getAllocPara().startPosCart;
	uint indexEnd = indexBegin + cartPara.nodeIndexEnd;
	std::cout << " in move, torque  = " << cartPara.totalTorque << std::endl;
	std::cout << " in move, inertia  = " << cartPara.moInertia << std::endl;
	cartPara.angularSpeed = cartPara.totalTorque / cartPara.moInertia;
	std::cout << " angular speed  = " << cartPara.angularSpeed << std::endl;
	// angle is counter-clock wise.
	double angle = cartPara.angularSpeed * dt;
	std::cout << "angle = " << angle << std::endl;
	std::cout << "fixed point = ";
	cartPara.fixedPt.Print();
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
	isInitialized = false;
	isParaInitialized = false;
}

void Cartilage::initializeMem(SceNodes* nodeInput) {
	if (!isParaInitialized) {
		std::string errorMsg =
				"Error: Initialize memory is not allowed before initialize CartPara";
		throw SceException(errorMsg, InputInitException);
	}
	readValuesFromConfig();

	nodes = nodeInput;
	initIsActive();
	isInitialized = true;
	uint pivotIndex1 = cartPara.pivotNode1Index
			+ nodes->getAllocPara().startPosCart;
	CVector pivotNode1Pos = CVector(nodes->getInfoVecs().nodeLocX[pivotIndex1],
			nodes->getInfoVecs().nodeLocY[pivotIndex1],
			nodes->getInfoVecs().nodeLocZ[pivotIndex1]);

	std::cout << "pivot point 1 = ";
	pivotNode1Pos.Print();

	uint pivotIndex2 = cartPara.pivotNode2Index
			+ nodes->getAllocPara().startPosCart;
	CVector pivotNode2Pos = CVector(nodes->getInfoVecs().nodeLocX[pivotIndex2],
			nodes->getInfoVecs().nodeLocY[pivotIndex2],
			nodes->getInfoVecs().nodeLocZ[pivotIndex2]);

	std::cout << "pivot point 2 = ";
	pivotNode2Pos.Print();

	cartPara.fixedPt = pivotNode1Pos + pivotNode2Pos;
	cartPara.fixedPt = cartPara.fixedPt / 2;

	calculateGrowthDir();
}

void Cartilage::addPoint1(CVector &nodeBehindPos) {
	CVector node1GrowthDir = cartPara.growNode1 - nodeBehindPos;
	// get unit vector
	node1GrowthDir = node1GrowthDir.getUnitVector();
	CVector incrementVec = node1GrowthDir * cartPara.growthThreshold;
	CVector posNew = nodeBehindPos + incrementVec;
	uint indexNew = cartPara.nodeIndexEnd + nodes->getAllocPara().startPosCart;

	std::cout << "In run addPoint1, new point index = " << indexNew
			<< std::endl;
	std::cout << "grow node1 pos = ";
	cartPara.growNode1.Print();
	std::cout << "nodeBehindPos position = ";
	nodeBehindPos.Print();
	std::cout << "node1GrowthDir= ";
	node1GrowthDir.Print();
	std::cout << "increment vector = ";
	incrementVec.Print();

	nodes->getInfoVecs().nodeLocX[indexNew] = posNew.GetX();
	nodes->getInfoVecs().nodeLocY[indexNew] = posNew.GetY();
	nodes->getInfoVecs().nodeLocZ[indexNew] = posNew.GetZ();
	nodes->getInfoVecs().nodeIsActive[indexNew] = true;
	cartPara.growNodeBehind1Index = cartPara.nodeIndexEnd;
	cartPara.nodeIndexEnd++;
}

void Cartilage::addPoint2(CVector &nodeBehindPos) {
	CVector node2GrowthDir = cartPara.growNode2 - nodeBehindPos;
// get unit vector
	node2GrowthDir = node2GrowthDir.getUnitVector();
	CVector incrementVec = node2GrowthDir * cartPara.growthThreshold;
	CVector posNew = nodeBehindPos + incrementVec;
	uint indexNew = cartPara.nodeIndexEnd + nodes->getAllocPara().startPosCart;

	std::cout << "grow node1 pos = ";
	cartPara.growNode1.Print();
	std::cout << "nodeBehindPos position = ";
	nodeBehindPos.Print();
	std::cout << "node2GrowthDir= ";
	node2GrowthDir.Print();
	std::cout << "In run addPoint2, new point index = " << indexNew
			<< std::endl;
	std::cout << "increment vector = ";
	incrementVec.Print();

	nodes->getInfoVecs().nodeLocX[indexNew] = posNew.GetX();
	nodes->getInfoVecs().nodeLocY[indexNew] = posNew.GetY();
	nodes->getInfoVecs().nodeLocZ[indexNew] = posNew.GetZ();
	nodes->getInfoVecs().nodeIsActive[indexNew] = true;
	cartPara.growNodeBehind2Index = cartPara.nodeIndexEnd;
	cartPara.nodeIndexEnd++;
}

void Cartilage::grow1(double dt) {
	CVector movementVec = cartPara.node1GrowthDir * cartPara.growthSpeedNode1;
	uint node1GlobalIndex = cartPara.growNode1Index
			+ nodes->getAllocPara().startPosCart;

	std::cout << "In run grow1, node 1 index = " << node1GlobalIndex
			<< std::endl;
	std::cout << "movement vector = ";
	movementVec.Print();

	nodes->getInfoVecs().nodeLocX[node1GlobalIndex] += movementVec.GetX();
	nodes->getInfoVecs().nodeLocY[node1GlobalIndex] += movementVec.GetY();
	nodes->getInfoVecs().nodeLocZ[node1GlobalIndex] += movementVec.GetZ();
}

void Cartilage::grow2(double dt) {
	CVector movementVec = cartPara.node2GrowthDir * cartPara.growthSpeedNode2;
	uint node2GlobalIndex = cartPara.growNode2Index
			+ nodes->getAllocPara().startPosCart;

	std::cout << "In run grow2, node 2 index = " << node2GlobalIndex
			<< std::endl;
	std::cout << "movement vector = ";
	movementVec.Print();

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

	std::cout << "In run growth logics 1, behind 1 index = " << behind1Index
			<< std::endl;
	std::cout << "In run growth logics 1, behind 1 position = ";
	node1BehindPos.Print();

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

	std::cout << "In run growth logics 2 , behind 2 index = " << behind2Index
			<< std::endl;
	std::cout << "In run growth logics 2 , behind 2 position = ";
	node2BehindPos.Print();

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
	std::cout << "before run growth logic 1 " << std::endl;
	runGrowthLogics1(dt);
	std::cout << "before run growth logic 2 " << std::endl;
	runGrowthLogics2(dt);
	std::cout << "before updating tip nodes " << std::endl;
	updateTipNodes(dt);

	//int jj;
	//std::cin >> jj;
}

/*
 void Cartilage::initializeVal(std::vector<CVector>& initialNodes,
 CVector& tipNode1, CVector& tipNode2, CVector &pivotNode1,
 CVector pivotNode2) {

 }
 */

void Cartilage::updateTipNodes(double dt) {
	CVector movementVec = cartPara.node2GrowthDir * cartPara.growthSpeedNode2;
	for (uint i = cartPara.tipNodeStartPos; i < cartPara.tipNodeIndexEnd; i++) {
		uint globalIndex = i + nodes->getAllocPara().startPosCart;
		nodes->getInfoVecs().nodeLocX[globalIndex] += movementVec.GetX();
		nodes->getInfoVecs().nodeLocY[globalIndex] += movementVec.GetY();
		nodes->getInfoVecs().nodeLocZ[globalIndex] += movementVec.GetZ();
	}
}

void Cartilage::initIsActive() {
	for (uint i = 0; i < cartPara.nodeIndexEnd; i++) {
		uint gloablIndex = i + nodes->getAllocPara().startPosCart;
		bool isActive = false;
		if (i < cartPara.tipNodeIndexEnd
				|| (i >= cartPara.nonTipNodeStartPos
						&& i < cartPara.nodeIndexEnd)) {
			isActive = true;
		}
		nodes->getInfoVecs().nodeIsActive[gloablIndex] = isActive;
	}
}

void Cartilage::readValuesFromConfig() {
	cartPara.growthSpeedNode1 = globalConfigVars.getConfigValue(
			"GrowthSpeedNode1").toDouble();
	cartPara.growthSpeedNode2 = globalConfigVars.getConfigValue(
			"GrowthSpeedNode2").toDouble();
	cartPara.growthThreshold = globalConfigVars.getConfigValue(
			"GrowthThreshold").toDouble();
	cartPara.moInertia =
			globalConfigVars.getConfigValue("InitMemonetInertia").toDouble();
}

//void Cartilage::initializeNodes(CartilageRawData& rawData) {
//}

void Cartilage::runAllLogics(double dt) {
	if (isInitialized) {
		calculateGrowthDir();
		calculateTotalTorque();
		move(dt);
		runGrowthLogics(dt);
	}
}
