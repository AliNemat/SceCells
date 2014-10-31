#include "Cartilage.h"

__device__ __constant__ double cartFixedPtConst[3];
__device__ __constant__ double cartGrowDirConst[3];
__device__ __constant__ double cartAngSpeedConst;
__device__ __constant__ double cartCurLenConst;
__device__ __constant__ double effectiveRangeConst;

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
	cartPara.currentLength = direction.getModul();
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

	cartPara.totalTorque += cartPara.torqueFromEpi;
	cartPara.angularSpeed = cartPara.totalTorque / cartPara.moInertia;
	std::cout << " total torque = " << cartPara.totalTorque << std::endl;
}

void Cartilage::move(double dt) {
	uint indexBegin = nodes->getAllocPara().startPosCart;
	uint indexEnd = indexBegin + cartPara.nodeIndexEnd;
	std::cout << " in move, torque  = " << cartPara.totalTorque << std::endl;
	std::cout << " in move, inertia  = " << cartPara.moInertia << std::endl;

	std::cout << " angular speed  = " << cartPara.angularSpeed << std::endl;
	// angle is counter-clock wise.
	double angle = cartPara.angularSpeed * dt;
	std::cout << "angle = " << angle << std::endl;
	std::cout << "fixed point = ";
	//cartPara.fixedPt.Print();
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

	cartPara.torqueFromEpi = 0;

	cudaMemcpyToSymbol(effectiveRangeConst,
			&(nodes->getMechPara().sceCartParaCPU[4]), sizeof(double));

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

void Cartilage::handleCartNoSlippage() {
	CVector cartFixePt = cartPara.fixedPt;
	CVector cartGrowthDir = cartPara.growthDir;

	cartGrowthDir = cartGrowthDir.getUnitVector();
	double effectiveRange = nodes->getMechPara().sceCartParaCPU[4];

	double fixedPt[3];
	fixedPt[0] = cartFixePt.GetX();
	fixedPt[1] = cartFixePt.GetY();
	fixedPt[2] = cartFixePt.GetZ();
	cudaMemcpyToSymbol(cartFixedPtConst, fixedPt, 3 * sizeof(double));

	double growthDir[3];
	growthDir[0] = cartGrowthDir.GetX();
	growthDir[1] = cartGrowthDir.GetY();
	growthDir[2] = cartGrowthDir.GetZ();
	cudaMemcpyToSymbol(cartGrowDirConst, growthDir, 3 * sizeof(double));

	double angSpeed = cartPara.angularSpeed;
	cudaMemcpyToSymbol(cartAngSpeedConst, &angSpeed, sizeof(double));

	double curLen = cartPara.currentLength;
	cudaMemcpyToSymbol(cartCurLenConst, &curLen, sizeof(double));

	uint profileStartPos = nodes->getAllocPara().startPosProfile;
	uint profileEndPos = profileStartPos
			+ nodes->getAllocPara().currentActiveProfileNodeCount;

	uint activeEpiNodeCount =
			nodes->getAllocPara().currentActiveProfileNodeCount;
	thrust::device_vector<bool> isMergedToCart(activeEpiNodeCount);

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin()))
					+ profileStartPos,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin()))
					+ profileEndPos,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin()
									+ profileStartPos,
							nodes->getInfoVecs().nodeVelY.begin()
									+ profileStartPos, isMergedToCart.begin())),
			NoSlipHandler());

	thrust::counting_iterator<int> countingBegin(0);

	//thrust::copy_if(countingBegin);

	//FIXME: We probably need to do this on GPU. This problem seems to be easy but is actually non-

	int contactingNodeCount = thrust::reduce(isMergedToCart.begin(),
			isMergedToCart.end(), 0, thrust::plus<int>());

	thrust::device_vector<int> epiIndicies(contactingNodeCount);
	thrust::copy_if(countingBegin, countingBegin + activeEpiNodeCount,
			isMergedToCart.begin(), epiIndicies, isTrue());

	/*
	 int minIndex = thrust::reduce(
	 thrust::make_zip_iterator(
	 thrust::make_tuple(isMergedToCart.begin(), countingBegin)),
	 thrust::make_zip_iterator(
	 thrust::make_tuple(isMergedToCart.begin(), countingBegin))
	 + nodes->getAllocPara().currentActiveProfileNodeCount,
	 INT_MAX, MinCount());

	 */
	/*
	 int maxIndex = thrust::reduce(
	 thrust::make_zip_iterator(
	 thrust::make_tuple(isMergedToCart.begin(), countingBegin)),
	 thrust::make_zip_iterator(
	 thrust::make_tuple(isMergedToCart.begin(), countingBegin))
	 + nodes->getAllocPara().currentActiveProfileNodeCount, -1,
	 MaxCount());
	 */
}

void Cartilage::runAllLogics(double dt) {
	if (isInitialized) {
		calculateGrowthDir();
		calculateTotalTorque();
		handleCartNoSlippage();
		move(dt);
		runGrowthLogics(dt);
	}
}

__device__
double crossProduct2D(double& ax, double& ay, double& bx, double& by) {
	return (ax * by - ay * bx);
}

__device__
void crossProduct(double& ax, double& ay, double& az, double& bx, double& by,
		double& bz, double& cx, double& cy, double& cz) {
	cx = ay * bz - az * by;
	cy = az * bx - ax * bz;
	cz = ax * by - ay * bx;
}

__device__
double dotProduct(double& ax, double& ay, double& az, double& bx, double& by,
		double& bz) {
	return (ax * bx + ay * by + az * bz);
}

__device__
bool closeToCart(double& xPos, double& yPos) {

	double ax = xPos - cartFixedPtConst[0];
	double ay = yPos - cartFixedPtConst[1];
	double az = 0;
	double bx = cartGrowDirConst[0];
	double by = cartGrowDirConst[1];
	double bz = 0;
	double cx, cy, cz;
	crossProduct(ax, ay, az, bx, by, bz, cx, cy, cz);
	double distToFixedPt = dotProduct(ax, ay, az, bx, by, bz);
// fabs(cz) is actually the area. because growthDir is a unit vector,
// fabs(cz) is also the distance of node to a vector.
	if (fabs(cz) < effectiveRangeConst
			&& fabs(distToFixedPt - cartCurLenConst) < effectiveRangeConst) {
		// that means a node is very close to the cartilage bar.
		return true;
	} else {
		return false;
	}
}

__device__
void modifyVelNoSlip(double& xPos, double& yPos, double& xVel, double& yVel) {

// if an epithelium node is adhere to the cartilage, it will move together with cartilage.
// so as a first step we clear the previous result.
	double dummyz1 = 0, dummyz2 = 0;
	double product = dotProduct(xVel, yVel, dummyz1, cartGrowDirConst[0],
			cartGrowDirConst[1], dummyz2);
	xVel = cartGrowDirConst[0] * product;
	yVel = cartGrowDirConst[1] * product;

// we also want the attached epithelium nodes move along the cartilage.
	double dirToFixX = xPos - cartFixedPtConst[0];
	double dirToFixY = yPos - cartFixedPtConst[1];

	double angleSpeedX = -dirToFixY * cartAngSpeedConst;
	double angleSpeedY = dirToFixX * cartAngSpeedConst;
	xVel = xVel + angleSpeedX;
	yVel = yVel + angleSpeedY;

}
