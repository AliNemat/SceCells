#include "SceCells.h"
#include <cmath>

double epsilon = 1.0e-12;

__constant__ double membrEquLen;
__constant__ double membrStiff;
__constant__ double pI;
__constant__ double minLength;
__constant__ double minDivisor;
__constant__ uint maxAllNodePerCell;
__constant__ uint maxMembrPerCell;
__constant__ uint maxIntnlPerCell;
__constant__ double bendCoeff;

__constant__ double sceIB_M[5];
__constant__ double sceIBDiv_M[5];
__constant__ double sceII_M[5];
__constant__ double sceIIDiv_M[5];
__constant__ double grthPrgrCriEnd_M;
__constant__ double F_Ext_Incline_M2 ;  //Ali

__device__
double calMembrForce(double& length) {
//Ali	if (length < membrEquLen) {
//		return 0;
//	} else {
		return (length - membrEquLen) * membrStiff;
//	}
}
//Ali
__device__
double calExtForce(double& curTime) {
		return curTime * F_Ext_Incline_M2;
}
//Ali
__device__
double obtainRandAngle(uint& cellRank, uint& seed) {
	thrust::default_random_engine rng(seed);
	// discard n numbers to avoid correlation
	rng.discard(cellRank);
	thrust::uniform_real_distribution<double> u0Pi(0, 2.0 * pI);
	double randomAngle = u0Pi(rng);
	return randomAngle;
}

__device__ uint obtainNewIntnlNodeIndex(uint& cellRank, uint& curActiveCount) {
	return (cellRank * maxAllNodePerCell + maxMembrPerCell + curActiveCount);
}

__device__
bool isAllIntnlFilled(uint& currentIntnlCount) {
	if (currentIntnlCount < maxIntnlPerCell) {
		return false;
	} else {
		return true;
	}
}

__device__
bool longEnough(double& length) {
	if (length > minLength) {
		return true;
	} else {
		return false;
	}
}

__device__
double compDist2D(double &xPos, double &yPos, double &xPos2, double &yPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2));
}

void SceCells::distributeBdryIsActiveInfo() {
	thrust::fill(nodes->getInfoVecs().nodeIsActive.begin(),
			nodes->getInfoVecs().nodeIsActive.begin()
					+ allocPara.startPosProfile, true);
}

void SceCells::distributeProfileIsActiveInfo() {
	thrust::fill(
			nodes->getInfoVecs().nodeIsActive.begin()
					+ allocPara.startPosProfile,
			nodes->getInfoVecs().nodeIsActive.begin()
					+ allocPara.startPosProfile
					+ nodes->getAllocPara().currentActiveProfileNodeCount,
			true);
}

void SceCells::distributeECMIsActiveInfo() {
	uint totalNodeCountForActiveECM = allocPara.currentActiveECM
			* allocPara.maxNodePerECM;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveECM);
	thrust::fill(
			nodes->getInfoVecs().nodeIsActive.begin() + allocPara.startPosECM,
			nodes->getInfoVecs().nodeIsActive.begin()
					+ totalNodeCountForActiveECM + allocPara.startPosECM, true);
}

void SceCells::distributeCellIsActiveInfo() {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	thrust::transform(
			thrust::make_transform_iterator(countingBegin,
					ModuloFunctor(allocPara.maxNodeOfOneCell)),
			thrust::make_transform_iterator(countingEnd,
					ModuloFunctor(allocPara.maxNodeOfOneCell)),
			thrust::make_permutation_iterator(
					cellInfoVecs.activeNodeCountOfThisCell.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			nodes->getInfoVecs().nodeIsActive.begin() + allocPara.startPosCells,
			thrust::less<uint>());
}

void SceCells::distributeCellGrowthProgress() {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	thrust::copy(
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingEnd,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			nodes->getInfoVecs().nodeGrowPro.begin() + allocPara.startPosCells);
}

void MembrPara::initFromConfig() {
	membrEquLenCPU = globalConfigVars.getConfigValue("MembrEquLen").toDouble();
	membrStiffCPU = globalConfigVars.getConfigValue("MembrStiff").toDouble();
	membrGrowCoeff_Ori =
			globalConfigVars.getConfigValue("MembrGrowCoeff").toDouble();
	membrGrowLimit_Ori =
			globalConfigVars.getConfigValue("MembrGrowLimit").toDouble();
	membrGrowCoeff = membrGrowCoeff_Ori;
	membrGrowLimit = membrGrowLimit_Ori;
        //Ali
        F_Ext_Incline = 
			globalConfigVars.getConfigValue("FExtIncline").toDouble();
        //Ali
	membrBendCoeff =
			globalConfigVars.getConfigValue("MembrBenCoeff").toDouble();
	adjustLimit =
			globalConfigVars.getConfigValue("MembrAdjustLimit").toDouble();
	adjustCoeff =
			globalConfigVars.getConfigValue("MembrAdjustCoeff").toDouble();

	growthConst_N =
			globalConfigVars.getConfigValue("MembrGrowthConst").toDouble();
	initMembrCt_N =
			globalConfigVars.getConfigValue("InitMembrNodeCount").toInt();
	initIntnlCt_N =
			globalConfigVars.getConfigValue("InitCellNodeCount").toInt();
}

SceCells::SceCells() {
	curTime = 0;
}

void SceCells::growAtRandom(double d_t) {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;

	// randomly select growth direction and speed.
	randomizeGrowth();

	//std::cout << "after copy grow info" << std::endl;
	updateGrowthProgress();
	//std::cout << "after update growth progress" << std::endl;
	decideIsScheduleToGrow();
	//std::cout << "after decode os schedule to grow" << std::endl;
	computeCellTargetLength();
	//std::cout << "after compute cell target length" << std::endl;
	computeDistToCellCenter();
	//std::cout << "after compute dist to center" << std::endl;
	findMinAndMaxDistToCenter();
	//std::cout << "after find min and max dist" << std::endl;
	computeLenDiffExpCur();
	//std::cout << "after compute diff " << std::endl;
	stretchCellGivenLenDiff();
	//std::cout << "after apply stretch force" << std::endl;
	cellChemotaxis();
	//std::cout << "after apply cell chemotaxis" << std::endl;
	addPointIfScheduledToGrow();
	//std::cout << "after adding node" << std::endl;
}

/**
 * Use the growth magnitude and dt to update growthProgress.
 */
void SceCells::updateGrowthProgress() {
	thrust::transform(cellInfoVecs.growthSpeed.begin(),
			cellInfoVecs.growthSpeed.begin() + allocPara.currentActiveCellCount,
			cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin(), SaxpyFunctorWithMaxOfOne(dt));
}

/**
 * Decide if the cells are going to add a node or not.
 * Use lastCheckPoint and growthProgress to decide whether add point or not
 */
void SceCells::decideIsScheduleToGrow() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.isScheduledToGrow.begin(),
			PtCondiOp(miscPara.growThreshold));
}

/**
 * Calculate target length of cell given the cell growth progress.
 * length is along the growth direction.
 */
void SceCells::computeCellTargetLength() {
	thrust::transform(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.expectedLength.begin(),
			CompuTarLen(bioPara.cellInitLength, bioPara.cellFinalLength));
}

/**
 * Compute distance of each node to its corresponding cell center.
 * The distantce could be either positive or negative, depending on the pre-defined
 * growth direction.
 */
void SceCells::computeDistToCellCenter() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(), CompuDist());
}

/**
 * For nodes of each cell, find the maximum and minimum distance to the center.
 * We will then calculate the current length of a cell along its growth direction
 * using max and min distance to the center.
 */
void SceCells::findMinAndMaxDistToCenter() {
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.smallestDistance.begin(), thrust::equal_to<uint>(),
			thrust::minimum<double>());

	// for nodes of each cell, find the maximum distance from the node to the corresponding

	// cell center along the pre-defined growth direction.
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara.maxNodeOfOneCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.biggestDistance.begin(), thrust::equal_to<uint>(),
			thrust::maximum<double>());

}

/**
 * Compute the difference for cells between their expected length and current length.
 */
void SceCells::computeLenDiffExpCur() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin()))
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.lengthDifference.begin(), CompuDiff());
}

/**
 * Use the difference that just computed and growthXDir&growthYDir
 * to apply stretching force (velocity) on nodes of all cells
 */
void SceCells::stretchCellGivenLenDiff() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			ApplyStretchForce(bioPara.elongationCoefficient));
}
/**
 * This is just an attempt. Cells move according to chemicals.
 */
void SceCells::cellChemotaxis() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.growthSpeed.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.growthSpeed.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(
													allocPara.maxNodeOfOneCell))),
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara.startPosCells)),
			ApplyChemoVel(bioPara.chemoCoefficient));
}
/**
 * Adjust the velocities of nodes.
 * For example, velocity of boundary nodes must be zero.
 */
void SceCells::adjustNodeVel() {
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(
			totalNodeCountForActiveCells + allocPara.startPosCells);

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin(),
							countingIterBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin(),
							countingIterBegin)) + totalNodeCountForActiveCells
					+ allocPara.startPosCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			VelocityModifier(allocPara.startPosProfile,
					allocPara.currentActiveProfileNodeCount));
}
/**
 * Move nodes according to the velocity we just adjusted.
 */
void SceCells::moveNodes() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells + allocPara.startPosCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			SaxpyFunctorDim2(dt));
}
/**
 * Add a point to a cell if it is scheduled to grow.
 * This step does not guarantee success ; If adding new point failed, it will not change
 * isScheduleToGrow and activeNodeCount;
 */
void SceCells::addPointIfScheduledToGrow() {

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), countingBegin,
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), countingBegin,
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			AddPtOp(allocPara.maxNodeOfOneCell, miscPara.addNodeDistance,
					miscPara.minDistanceToOtherNode,
					growthAuxData.nodeIsActiveAddress,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress, time(NULL),
					miscPara.growThreshold));
}

SceCells::SceCells(SceNodes* nodesInput,
		std::vector<uint>& numOfInitActiveNodesOfCells,
		std::vector<SceNodeType>& cellTypes) :
		countingBegin(0), initIntnlNodeCount(
				nodesInput->getAllocPara().maxNodeOfOneCell / 2), initGrowthProgress(
				0.0) {
	curTime = 0.0;

	initialize(nodesInput);

	copyInitActiveNodeCount(numOfInitActiveNodesOfCells);

	thrust::device_vector<SceNodeType> cellTypesToPass = cellTypes;
	setCellTypes(cellTypesToPass);

	distributeIsActiveInfo();
}

SceCells::SceCells(SceNodes* nodesInput,
		std::vector<uint>& initActiveMembrNodeCounts,
		std::vector<uint>& initActiveIntnlNodeCounts,
		std::vector<double> &initGrowProgVec) {
	curTime = 0.0;
	tmpDebug = false;
	aniDebug = false;
	membrPara.initFromConfig();
	shrinkRatio = globalConfigVars.getConfigValue("ShrinkRatio").toDouble();
	centerShiftRatio =
			globalConfigVars.getConfigValue("CenterShiftRatio").toDouble();

	memNewSpacing = globalConfigVars.getConfigValue("MembrLenDiv").toDouble();

	initialize_M(nodesInput);
	copyToGPUConstMem();
	copyInitActiveNodeCount_M(initActiveMembrNodeCounts,
			initActiveIntnlNodeCounts, initGrowProgVec);
}

void SceCells::initCellInfoVecs() {
	cellInfoVecs.growthProgress.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.expectedLength.resize(allocPara.maxCellCount,
			bioPara.cellInitLength);
	cellInfoVecs.lengthDifference.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.smallestDistance.resize(allocPara.maxCellCount);
	cellInfoVecs.biggestDistance.resize(allocPara.maxCellCount);
	cellInfoVecs.activeNodeCountOfThisCell.resize(allocPara.maxCellCount);
	cellInfoVecs.lastCheckPoint.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.isDividing.resize(allocPara.maxCellCount);
	cellInfoVecs.cellTypes.resize(allocPara.maxCellCount, MX);
	cellInfoVecs.isScheduledToGrow.resize(allocPara.maxCellCount, false);
	cellInfoVecs.centerCoordX.resize(allocPara.maxCellCount);
	cellInfoVecs.centerCoordY.resize(allocPara.maxCellCount);
	cellInfoVecs.centerCoordZ.resize(allocPara.maxCellCount);
	cellInfoVecs.cellRanksTmpStorage.resize(allocPara.maxCellCount);
	cellInfoVecs.growthSpeed.resize(allocPara.maxCellCount, 0.0);
	cellInfoVecs.growthXDir.resize(allocPara.maxCellCount);
	cellInfoVecs.growthYDir.resize(allocPara.maxCellCount);
	cellInfoVecs.isRandGrowInited.resize(allocPara.maxCellCount, false);
}

void SceCells::initCellInfoVecs_M() {
	//std::cout << "max cell count = " << allocPara_m.maxCellCount << std::endl;
	cellInfoVecs.growthProgress.resize(allocPara_m.maxCellCount, 0.0);
        //Ali
	cellInfoVecs.Cell_Time.resize(allocPara_m.maxCellCount, 0.0);
        //Ali
	cellInfoVecs.expectedLength.resize(allocPara_m.maxCellCount,
			bioPara.cellInitLength);
	cellInfoVecs.lengthDifference.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.smallestDistance.resize(allocPara_m.maxCellCount);
	cellInfoVecs.biggestDistance.resize(allocPara_m.maxCellCount);
	cellInfoVecs.activeMembrNodeCounts.resize(allocPara_m.maxCellCount);
	cellInfoVecs.activeIntnlNodeCounts.resize(allocPara_m.maxCellCount);
	cellInfoVecs.lastCheckPoint.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.isDividing.resize(allocPara_m.maxCellCount);
	cellInfoVecs.isScheduledToGrow.resize(allocPara_m.maxCellCount, false);
	cellInfoVecs.centerCoordX.resize(allocPara_m.maxCellCount);
	cellInfoVecs.centerCoordY.resize(allocPara_m.maxCellCount);
	cellInfoVecs.centerCoordZ.resize(allocPara_m.maxCellCount);
	cellInfoVecs.cellRanksTmpStorage.resize(allocPara_m.maxCellCount);
	cellInfoVecs.growthSpeed.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.growthXDir.resize(allocPara_m.maxCellCount);
	cellInfoVecs.growthYDir.resize(allocPara_m.maxCellCount);
	cellInfoVecs.isRandGrowInited.resize(allocPara_m.maxCellCount, false);
	cellInfoVecs.isMembrAddingNode.resize(allocPara_m.maxCellCount, false);
	cellInfoVecs.maxTenIndxVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.maxTenRiVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.maxTenRiMidXVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.maxTenRiMidYVec.resize(allocPara_m.maxCellCount);
	cellInfoVecs.aveTension.resize(allocPara_m.maxCellCount);
	cellInfoVecs.membrGrowProgress.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.membrGrowSpeed.resize(allocPara_m.maxCellCount, 0.0);
	cellInfoVecs.cellAreaVec.resize(allocPara_m.maxCellCount, 0.0);
	//std::cout << "finished " << std::endl;
}

void SceCells::initCellNodeInfoVecs() {
	cellNodeInfoVecs.cellRanks.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.activeXPoss.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.activeYPoss.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.activeZPoss.resize(allocPara.maxTotalCellNodeCount);
	cellNodeInfoVecs.distToCenterAlongGrowDir.resize(
			allocPara.maxTotalCellNodeCount);
}

void SceCells::initCellNodeInfoVecs_M() {
	std::cout << "max total node count = " << allocPara_m.maxTotalNodeCount
			<< std::endl;
	cellNodeInfoVecs.cellRanks.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.activeXPoss.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.activeYPoss.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.activeZPoss.resize(allocPara_m.maxTotalNodeCount);
	cellNodeInfoVecs.distToCenterAlongGrowDir.resize(
			allocPara_m.maxTotalNodeCount);
}

void SceCells::initGrowthAuxData() {
	growthAuxData.nodeIsActiveAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[allocPara.startPosCells]));
	growthAuxData.nodeXPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[allocPara.startPosCells]));
	growthAuxData.nodeYPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[allocPara.startPosCells]));
	growthAuxData.randomGrowthSpeedMin = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMin").toDouble();
	growthAuxData.randomGrowthSpeedMax = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMax").toDouble();
	growthAuxData.randGenAuxPara = globalConfigVars.getConfigValue(
			"RandomGenerationAuxPara").toDouble();

	if (controlPara.simuType == SingleCellTest) {
		growthAuxData.fixedGrowthSpeed = globalConfigVars.getConfigValue(
				"FixedGrowthSpeed").toDouble();
	}
}

void SceCells::initGrowthAuxData_M() {
	growthAuxData.nodeIsActiveAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[allocPara_m.bdryNodeCount]));
	growthAuxData.nodeXPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[allocPara_m.bdryNodeCount]));
	growthAuxData.nodeYPosAddress = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[allocPara_m.bdryNodeCount]));
	growthAuxData.adhIndxAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeAdhereIndex[allocPara_m.bdryNodeCount]));
	growthAuxData.randomGrowthSpeedMin_Ori = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMin").toDouble();
	growthAuxData.randomGrowthSpeedMax_Ori = globalConfigVars.getConfigValue(
			"RandomGrowthSpeedMax").toDouble();
	growthAuxData.randomGrowthSpeedMin = growthAuxData.randomGrowthSpeedMin_Ori;
	growthAuxData.randomGrowthSpeedMax = growthAuxData.randomGrowthSpeedMax_Ori;
	growthAuxData.grthPrgrCriVal_M_Ori = globalConfigVars.getConfigValue(
			"GrowthPrgrCriVal").toDouble();
	growthAuxData.grthProgrEndCPU = globalConfigVars.getConfigValue(
			"GrowthPrgrValEnd").toDouble();
}

void SceCells::initialize(SceNodes* nodesInput) {
	nodes = nodesInput;
	controlPara = nodes->getControlPara();
	readMiscPara();
	readBioPara();

	allocPara = nodesInput->getAllocPara();

	// max internal node count must be even number.
	assert(allocPara_m.maxIntnlNodePerCell % 2 == 0);

	initCellInfoVecs();
	initCellNodeInfoVecs();
	initGrowthAuxData();

	distributeIsCellRank();
}

void SceCells::initialize_M(SceNodes* nodesInput) {
	std::cout << "Initializing cells ...... " << std::endl;
	//std::cout.flush();
	nodes = nodesInput;
	allocPara_m = nodesInput->getAllocParaM();
	// max internal node count must be even number.
	assert(allocPara_m.maxIntnlNodePerCell % 2 == 0);

	//std::cout << "break point 1 " << std::endl;
	//std::cout.flush();
	controlPara = nodes->getControlPara();
	//std::cout << "break point 2 " << std::endl;
	//std::cout.flush();
	readMiscPara_M();
	//std::cout << "break point 3 " << std::endl;
	//std::cout.flush();
	initCellInfoVecs_M();

	//std::cout << "break point 4 " << std::endl;
	//std::cout.flush();
	readBioPara();
	//std::cout << "break point 5 " << std::endl;
	//std::cout.flush();

	//std::cout << "break point 6 " << std::endl;
	//std::cout.flush();
	initCellNodeInfoVecs_M();
	//std::cout << "break point 7 " << std::endl;
	//std::cout.flush();
	initGrowthAuxData_M();
	//std::cout << "break point 8 " << std::endl;
	//std::cout.flush();

}

void SceCells::copyInitActiveNodeCount(
		std::vector<uint>& numOfInitActiveNodesOfCells) {
	thrust::copy(numOfInitActiveNodesOfCells.begin(),
			numOfInitActiveNodesOfCells.end(),
			cellInfoVecs.activeNodeCountOfThisCell.begin());
}

void SceCells::allComponentsMove() {
	adjustNodeVel();
	moveNodes();
}

/**
 * Mark cell node as either activdistributeIsActiveInfo()e or inactive.
 * left part of the node array will be active and right part will be inactive.
 * the threshold is defined by array activeNodeCountOfThisCell.
 * e.g. activeNodeCountOfThisCell = {2,3} and  maxNodeOfOneCell = 5
 */
void SceCells::distributeIsActiveInfo() {
//std::cout << "before distribute bdry isActive" << std::endl;
	distributeBdryIsActiveInfo();
//std::cout << "before distribute profile isActive" << std::endl;
	distributeProfileIsActiveInfo();
//std::cout << "before distribute ecm isActive" << std::endl;
	distributeECMIsActiveInfo();
//std::cout << "before distribute cells isActive" << std::endl;
	distributeCellIsActiveInfo();
}

void SceCells::distributeIsCellRank() {
	uint totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingCellEnd(
			totalNodeCountForActiveCells);
	std::cerr << "totalNodeCount for active cells "
			<< totalNodeCountForActiveCells << std::endl;

//thrust::counting_iterator<uint> countingECMEnd(countingECMEnd);

// only computes the cell ranks of cells. the rest remain unchanged.
	thrust::transform(countingBegin, countingCellEnd,
			nodes->getInfoVecs().nodeCellRank.begin() + allocPara.startPosCells,
			DivideFunctor(allocPara.maxNodeOfOneCell));
	std::cerr << "finished cellRank transformation" << std::endl;
}

/**
 * This method computes center of all cells.
 * more efficient then simply iterating the cell because of parallel reducing.
 */
void SceCells::computeCenterPos() {
	uint totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);
	uint totalNumberOfActiveNodes = thrust::reduce(
			cellInfoVecs.activeNodeCountOfThisCell.begin(),
			cellInfoVecs.activeNodeCountOfThisCell.begin()
					+ allocPara.currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().nodeIsActive.begin() + allocPara.startPosCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin(),
							cellNodeInfoVecs.activeZPoss.begin())), isTrue());

	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalNumberOfActiveNodes,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin(),
							cellNodeInfoVecs.activeZPoss.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin())),
			thrust::equal_to<uint>(), CVec3Add());
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin()))
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.activeNodeCountOfThisCell.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(),
							cellInfoVecs.centerCoordZ.begin())), CVec3Divide());
}

/**
 * 2D version of cell division.
 * Division process is done by creating two temporary vectors to hold the node information
 * that are going to divide.
 *
 * step 1: based on lengthDifference, expectedLength and growthProgress,
 *     this process determines whether a certain cell is ready to divide and then assign
 *     a boolean value to isDivided.
 *
 * step 2. copy those cells that will divide in to the temp vectors created
 *
 * step 3. For each cell in the temp vectors, we sort its nodes by its distance to the
 * corresponding cell center.
 * This step is not very effcient when the number of cells going to divide is big.
 * but this is unlikely to happen because cells will divide according to external chemical signaling
 * and each will have different divide progress.
 *
 * step 4. copy the right part of each cell of the sorted array (temp1) to left part of each cell of
 * another array
 *
 * step 5. transform isActive vector of both temp1 and temp2, making only left part of each cell active.
 *
 * step 6. insert temp2 to the end of the cell array
 *
 * step 7. copy temp1 to the previous position of the cell array.
 *
 * step 8. add activeCellCount of the system.
 *
 * step 9. mark isDivide of all cells to false.
 */

void SceCells::divide2DSimplified() {
	bool isDivisionPresent = decideIfGoingToDivide();
	if (!isDivisionPresent) {
		return;
	}
	copyCellsPreDivision();
	sortNodesAccordingToDist();
	copyLeftAndRightToSeperateArrays();
	transformIsActiveArrayOfBothArrays();
	addSecondArrayToCellArray();
	copyFirstArrayToPreviousPos();
	updateActiveCellCount();
	markIsDivideFalse();
}

bool SceCells::decideIfGoingToDivide() {
// step 1
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.lengthDifference.begin(),
							cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.lengthDifference.begin(),
							cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeNodeCountOfThisCell.begin()))
					+ allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isDividing.begin(),
							cellInfoVecs.growthProgress.begin())),
			CompuIsDivide(miscPara.isDivideCriticalRatio,
					allocPara.maxNodeOfOneCell));
// sum all bool values which indicate whether the cell is going to divide.
// toBeDivideCount is the total number of cells going to divide.
	divAuxData.toBeDivideCount = thrust::reduce(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin() + allocPara.currentActiveCellCount,
			(uint) (0));
	if (divAuxData.toBeDivideCount > 0) {
		return true;
	} else {
		return false;
	}
}

void SceCells::copyCellsPreDivision() {
// step 2 : copy all cell rank and distance to its corresponding center with divide flag = 1
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;

	divAuxData.nodeStorageCount = divAuxData.toBeDivideCount
			* allocPara.maxNodeOfOneCell;
	divAuxData.tmpIsActiveHold1 = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, true);
	divAuxData.tmpDistToCenter1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpCellRankHold1 = thrust::device_vector<uint>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpXValueHold1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYValueHold1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpZValueHold1 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

	divAuxData.tmpCellTypes = thrust::device_vector<SceNodeType>(
			divAuxData.nodeStorageCount);

	divAuxData.tmpIsActiveHold2 = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpDistToCenter2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpXValueHold2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYValueHold2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpZValueHold2 = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

// step 2 , continued
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeCellType.begin()
									+ allocPara.startPosCells)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(allocPara.maxNodeOfOneCell)),
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeCellType.begin()
									+ allocPara.startPosCells))
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(cellInfoVecs.isDividing.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara.maxNodeOfOneCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpCellRankHold1.begin(),
							divAuxData.tmpDistToCenter1.begin(),
							divAuxData.tmpXValueHold1.begin(),
							divAuxData.tmpYValueHold1.begin(),
							divAuxData.tmpZValueHold1.begin(),
							divAuxData.tmpCellTypes.begin())), isTrue());
}

/**
 * performance wise, this implementation is not the best because I can use only one sort_by_key
 * with speciialized comparision operator. However, This implementation is more robust and won't
 * compromise performance too much.
 */
void SceCells::sortNodesAccordingToDist() {
//step 3
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		thrust::sort_by_key(
				divAuxData.tmpDistToCenter1.begin()
						+ i * allocPara.maxNodeOfOneCell,
				divAuxData.tmpDistToCenter1.begin()
						+ (i + 1) * allocPara.maxNodeOfOneCell,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								divAuxData.tmpXValueHold1.begin()
										+ i * allocPara.maxNodeOfOneCell,
								divAuxData.tmpYValueHold1.begin()
										+ i * allocPara.maxNodeOfOneCell,
								divAuxData.tmpZValueHold1.begin()
										+ i * allocPara.maxNodeOfOneCell)));
	}
}

/**
 * scatter_if() is a thrust function.
 * inputIter1 first,
 * inputIter1 last,
 * inputIter2 map,
 * inputIter3 stencil
 * randomAccessIter output
 */
void SceCells::copyLeftAndRightToSeperateArrays() {
//step 4.
	thrust::scatter_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpXValueHold1.begin(),
							divAuxData.tmpYValueHold1.begin(),
							divAuxData.tmpZValueHold1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpXValueHold1.end(),
							divAuxData.tmpYValueHold1.end(),
							divAuxData.tmpZValueHold1.end())),
			make_transform_iterator(countingBegin,
					LeftShiftFunctor(allocPara.maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					IsRightSide(allocPara.maxNodeOfOneCell)),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpXValueHold2.begin(),
							divAuxData.tmpYValueHold2.begin(),
							divAuxData.tmpZValueHold2.begin())));
}

void SceCells::transformIsActiveArrayOfBothArrays() {
	thrust::transform(countingBegin,
			countingBegin + divAuxData.nodeStorageCount,
			divAuxData.tmpIsActiveHold1.begin(),
			IsLeftSide(allocPara.maxNodeOfOneCell));
	thrust::transform(countingBegin,
			countingBegin + divAuxData.nodeStorageCount,
			divAuxData.tmpIsActiveHold2.begin(),
			IsLeftSide(allocPara.maxNodeOfOneCell));
	if (divAuxData.toBeDivideCount != 0) {
		std::cout << "before insert, active cell count in nodes:"
				<< nodes->getAllocPara().currentActiveCellCount << std::endl;
	}
}

void SceCells::addSecondArrayToCellArray() {
/// step 6. call SceNodes function to add newly divided cells
	nodes->addNewlyDividedCells(divAuxData.tmpXValueHold2,
			divAuxData.tmpYValueHold2, divAuxData.tmpZValueHold2,
			divAuxData.tmpIsActiveHold2, divAuxData.tmpCellTypes);
}

void SceCells::copyFirstArrayToPreviousPos() {
	thrust::scatter(
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActiveHold1.begin(),
							divAuxData.tmpXValueHold1.begin(),
							divAuxData.tmpYValueHold1.begin(),
							divAuxData.tmpZValueHold1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActiveHold1.end(),
							divAuxData.tmpXValueHold1.end(),
							divAuxData.tmpYValueHold1.end(),
							divAuxData.tmpZValueHold1.end())),
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(countingBegin,
									divAuxData.tmpCellRankHold1.begin())),
					CompuPos(allocPara.maxNodeOfOneCell)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara.startPosCells,
							nodes->getInfoVecs().nodeLocZ.begin()
									+ allocPara.startPosCells)));

	/**
	 * after dividing, the cell should resume the initial
	 * (1) node count, which defaults to be half size of max node count
	 * (2) growth progress, which defaults to 0
	 * (3) last check point, which defaults to 0
	 */
	thrust::scatter_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(initIntnlNodeCount, initGrowthProgress,
							initGrowthProgress)),
			thrust::make_zip_iterator(
					thrust::make_tuple(initIntnlNodeCount, initGrowthProgress,
							initGrowthProgress))
					+ allocPara.currentActiveCellCount, countingBegin,
			cellInfoVecs.isDividing.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeNodeCountOfThisCell.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())), isTrue());

// TODO: combine this one with the previous scatter_if to improve efficiency.
	thrust::fill(
			cellInfoVecs.activeNodeCountOfThisCell.begin()
					+ allocPara.currentActiveCellCount,
			cellInfoVecs.activeNodeCountOfThisCell.begin()
					+ allocPara.currentActiveCellCount
					+ divAuxData.toBeDivideCount,
			allocPara.maxNodeOfOneCell / 2);

}

void SceCells::updateActiveCellCount() {
	allocPara.currentActiveCellCount = allocPara.currentActiveCellCount
			+ divAuxData.toBeDivideCount;
	NodeAllocPara para = nodes->getAllocPara();
	para.currentActiveCellCount = allocPara.currentActiveCellCount;
	nodes->setAllocPara(para);
}

void SceCells::markIsDivideFalse() {
	thrust::fill(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin() + allocPara.currentActiveCellCount,
			false);
}

void SceCells::readMiscPara() {
	miscPara.addNodeDistance = globalConfigVars.getConfigValue(
			"DistanceForAddingNode").toDouble();
	miscPara.minDistanceToOtherNode = globalConfigVars.getConfigValue(
			"MinDistanceToOtherNode").toDouble();
	miscPara.isDivideCriticalRatio = globalConfigVars.getConfigValue(
			"IsDivideCrticalRatio").toDouble();
// reason for adding a small term here is to avoid scenario when checkpoint might add many times
// up to 0.99999999 which is theoretically 1.0 but not in computer memory. If we don't include
// this small term we might risk adding one more node.
	int maxNodeOfOneCell =
			globalConfigVars.getConfigValue("MaxNodePerCell").toInt();
	miscPara.growThreshold = 1.0 / (maxNodeOfOneCell - maxNodeOfOneCell / 2)
			+ epsilon;
}

void SceCells::readMiscPara_M() {
	miscPara.addNodeDistance = globalConfigVars.getConfigValue(
			"DistanceForAddingNode").toDouble();
	miscPara.minDistanceToOtherNode = globalConfigVars.getConfigValue(
			"MinDistanceToOtherNode").toDouble();
	miscPara.isDivideCriticalRatio = globalConfigVars.getConfigValue(
			"IsDivideCrticalRatio").toDouble();
// reason for adding a small term here is to avoid scenario when checkpoint might add many times
// up to 0.99999999 which is theoretically 1.0 but not in computer memory. If we don't include
// this small term we might risk adding one more node.
	int maxIntnlNodePerCell = globalConfigVars.getConfigValue(
			"MaxIntnlNodeCountPerCell").toInt();
	miscPara.growThreshold = 1.0
			/ (maxIntnlNodePerCell - maxIntnlNodePerCell / 2) + epsilon;
	miscPara.prolifDecayCoeff = globalConfigVars.getConfigValue(
			"ProlifDecayCoeff").toDouble();
}

void SceCells::readBioPara() {
	if (controlPara.simuType != Disc_M) {
		bioPara.cellInitLength = globalConfigVars.getConfigValue(
				"CellInitLength").toDouble();
		std::cout << "break point 1 " << bioPara.cellInitLength << std::endl;
		std::cout.flush();
		bioPara.cellFinalLength = globalConfigVars.getConfigValue(
				"CellFinalLength").toDouble();
		std::cout << "break point 2 " << bioPara.cellFinalLength << std::endl;
		std::cout.flush();
		bioPara.elongationCoefficient = globalConfigVars.getConfigValue(
				"ElongateCoefficient").toDouble();

		std::cout << "break point 3 " << bioPara.elongationCoefficient
				<< std::endl;
		std::cout.flush();

	}

	if (controlPara.simuType == Beak) {
		std::cout << "break point 4 " << std::endl;
		std::cout.flush();
		bioPara.chemoCoefficient = globalConfigVars.getConfigValue(
				"ChemoCoefficient").toDouble();
	}
	//std::cin >> jj;
}

void SceCells::randomizeGrowth() {
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)) + allocPara.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin())),
			AssignRandIfNotInit(growthAuxData.randomGrowthSpeedMin,
					growthAuxData.randomGrowthSpeedMax,
					allocPara.currentActiveCellCount,
					growthAuxData.randGenAuxPara));
}

/**
 * To run all the cell level logics.
 * First step we got center positions of cells.
 * Grow.
 */
void SceCells::runAllCellLevelLogicsDisc(double dt) {
	this->dt = dt;
//std::cerr << "enter run all cell level logics" << std::endl;
	computeCenterPos();
//std::cerr << "after compute center position." << std::endl;

	if (nodes->getControlPara().controlSwitchs.stab == OFF) {
		growAtRandom(dt);
		//grow2DTwoRegions(dt, region1, region2);
		//std::cerr << "after grow cells" << std::endl;
		//distributeIsActiveInfo();
		//std::cerr << "after distribute is active info." << std::endl;
		divide2DSimplified();
		//std::cerr << "after divide 2D simplified." << std::endl;
		distributeIsActiveInfo();
		//std::cerr << "after distribute is active info." << std::endl;
		distributeCellGrowthProgress();
	}

	allComponentsMove();
	//std::cerr << "after all components move." << std::endl;
}

//Ali void SceCells::runAllCellLogicsDisc_M(double dt) {
void SceCells::runAllCellLogicsDisc_M(double dt, double Damp_Coef) {   //Ali
	std::cout << "     *** 1 ***" << endl;
	std::cout.flush();
	this->dt = dt;
        this->Damp_Coef=Damp_Coef ; //Ali 
	growthAuxData.prolifDecay = exp(-curTime * miscPara.prolifDecayCoeff);
	growthAuxData.randomGrowthSpeedMin = growthAuxData.prolifDecay
			* growthAuxData.randomGrowthSpeedMin_Ori;
	growthAuxData.randomGrowthSpeedMax = growthAuxData.prolifDecay
			* growthAuxData.randomGrowthSpeedMax_Ori;
	curTime = curTime + dt;

	std::cout << "     *** 2 ***" << endl;
	std::cout.flush();
	applySceCellDisc_M();
	std::cout << "     *** 3 ***" << endl;
	std::cout.flush();
//Ali        
	computeCenterPos_M();
	std::cout << "     *** 5 ***" << endl;
	std::cout.flush();
//Ali 
	applyMemForce_M();
	std::cout << "     *** 4 ***" << endl;
	std::cout.flush();
     //Ali cmment //
//	computeCenterPos_M();
	std::cout << "     *** 5 ***" << endl;
	std::cout.flush();
     //Ali cmment //
	growAtRandom_M(dt);
	std::cout << "     *** 6 ***" << endl;
	std::cout.flush();
	divide2D_M();
	std::cout << "     *** 7 ***" << endl;
	std::cout.flush();
	distributeCellGrowthProgress_M();
	std::cout << "     *** 8 ***" << endl;
	std::cout.flush();
	allComponentsMove_M();
	std::cout << "     *** 9 ***" << endl;
	std::cout.flush();
	handleMembrGrowth_M();
	std::cout << "     *** 10 ***" << endl;
	std::cout.flush();
}

void SceCells::runStretchTest(double dt) {
	this->dt = dt;
	computeCenterPos();
	growAlongX(false, dt);
	moveNodes();
}

void SceCells::growAlongX(bool isAddPt, double d_t) {
	totalNodeCountForActiveCells = allocPara.currentActiveCellCount
			* allocPara.maxNodeOfOneCell;

	setGrowthDirXAxis();

//std::cout << "after copy grow info" << std::endl;
	updateGrowthProgress();
//std::cout << "after update growth progress" << std::endl;
	decideIsScheduleToGrow();
//std::cout << "after decode os schedule to grow" << std::endl;
	computeCellTargetLength();
//std::cout << "after compute cell target length" << std::endl;
	computeDistToCellCenter();
//std::cout << "after compute dist to center" << std::endl;
	findMinAndMaxDistToCenter();
//std::cout << "after find min and max dist" << std::endl;
	computeLenDiffExpCur();
//std::cout << "after compute diff " << std::endl;
	stretchCellGivenLenDiff();

	if (isAddPt) {
		addPointIfScheduledToGrow();
	}
}

void SceCells::growWithStress(double d_t) {
}

std::vector<CVector> SceCells::getAllCellCenters() {
	thrust::host_vector<double> centerX = cellInfoVecs.centerCoordX;
	thrust::host_vector<double> centerY = cellInfoVecs.centerCoordY;
	thrust::host_vector<double> centerZ = cellInfoVecs.centerCoordZ;
	std::vector<CVector> result;
	for (uint i = 0; i < allocPara.currentActiveCellCount; i++) {
		CVector pos = CVector(centerX[i], centerY[i], centerZ[i]);
		result.push_back(pos);
	}
	return result;
}

void SceCells::setGrowthDirXAxis() {
	thrust::fill(cellInfoVecs.growthXDir.begin(),
			cellInfoVecs.growthXDir.begin() + allocPara.currentActiveCellCount,
			1.0);
	thrust::fill(cellInfoVecs.growthYDir.begin(),
			cellInfoVecs.growthYDir.begin() + allocPara.currentActiveCellCount,
			0.0);
	thrust::fill(cellInfoVecs.growthSpeed.begin(),
			cellInfoVecs.growthSpeed.begin() + allocPara.currentActiveCellCount,
			growthAuxData.fixedGrowthSpeed);
}

std::vector<double> SceCells::getGrowthProgressVec() {
	thrust::host_vector<double> growthProVec = cellInfoVecs.growthProgress;
	std::vector<double> result;
	for (uint i = 0; i < allocPara.currentActiveCellCount; i++) {
		result.push_back(growthProVec[i]);
	}
	return result;
}

void SceCells::copyCellsPreDivision_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;

	divAuxData.nodeStorageCount = divAuxData.toBeDivideCount
			* allocPara_m.maxAllNodePerCell;

	divAuxData.tmpIsActive_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, true);
	divAuxData.tmpNodePosX_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpNodePosY_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

	divAuxData.tmpCellRank_M = thrust::device_vector<uint>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpDivDirX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpDivDirY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpCenterPosX_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);
	divAuxData.tmpCenterPosY_M = thrust::device_vector<double>(
			divAuxData.toBeDivideCount, 0);

	divAuxData.tmpIsActive1_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpXPos1_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYPos1_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

	divAuxData.tmpIsActive2_M = thrust::device_vector<bool>(
			divAuxData.nodeStorageCount, false);
	divAuxData.tmpXPos2_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);
	divAuxData.tmpYPos2_M = thrust::device_vector<double>(
			divAuxData.nodeStorageCount, 0.0);

// step 2 , continued
	thrust::counting_iterator<uint> iStart(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(cellInfoVecs.isDividing.begin(),
					make_transform_iterator(iStart,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpIsActive_M.begin(),
							divAuxData.tmpNodePosX_M.begin(),
							divAuxData.tmpNodePosY_M.begin())), isTrue());
// step 3 , continued
	thrust::counting_iterator<uint> iBegin(0);
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin, cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isDividing.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(divAuxData.tmpCellRank_M.begin(),
							divAuxData.tmpDivDirX_M.begin(),
							divAuxData.tmpDivDirY_M.begin(),
							divAuxData.tmpCenterPosX_M.begin(),
							divAuxData.tmpCenterPosY_M.begin())), isTrue());
}

void SceCells::createTwoNewCellArr_M() {
	divAuxData.tmp1MemActiveCounts.clear();
	divAuxData.tmp1InternalActiveCounts.clear();
	divAuxData.tmp2MemActiveCounts.clear();
	divAuxData.tmp2InternalActiveCounts.clear();

	//divDebug();

	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		divAuxData.tmp1IntnlVec.clear();
		divAuxData.tmp2IntnlVec.clear();

		vector<CVector> membrNodes;
		vector<CVector> intnlNodes;
		obtainMembrAndIntnlNodes(i, membrNodes, intnlNodes);

		CVector oldCenter = obtainCenter(i);
		double lenAlongMajorAxis;
		CVector divDir = calDivDir_MajorAxis(oldCenter, membrNodes,
				lenAlongMajorAxis);

		std::vector<VecVal> tmp1Membr, tmp2Membr;
		CVector cell1Center, cell2Center;

		obtainTwoNewCenters(oldCenter, divDir, lenAlongMajorAxis, cell1Center,
				cell2Center);

		prepareTmpVec(i, divDir, oldCenter, tmp1Membr, tmp2Membr);

		processMemVec(tmp1Membr, tmp2Membr);

		shiftIntnlNodesByCellCenter(cell1Center, cell2Center);

		assembleVecForTwoCells(i);
	}
	//divDebug();
}

void SceCells::copyFirstCellArr_M() {
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		uint cellRank = divAuxData.tmpCellRank_M[i];
		uint nodeStartIndx = cellRank * maxAllNodePerCell
				+ allocPara_m.bdryNodeCount;
		uint tmpStartIndx = i * maxAllNodePerCell;
		uint tmpEndIndx = (i + 1) * maxAllNodePerCell;
		thrust::constant_iterator<int> noAdhesion(-1), noAdhesion2(-1);
		thrust::copy(
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos1_M.begin(),
								divAuxData.tmpYPos1_M.begin(),
								divAuxData.tmpIsActive1_M.begin(), noAdhesion,
								noAdhesion2)) + tmpStartIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos1_M.begin(),
								divAuxData.tmpYPos1_M.begin(),
								divAuxData.tmpIsActive1_M.begin(), noAdhesion,
								noAdhesion2)) + tmpEndIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								nodes->getInfoVecs().nodeLocX.begin(),
								nodes->getInfoVecs().nodeLocY.begin(),
								nodes->getInfoVecs().nodeIsActive.begin(),
								nodes->getInfoVecs().nodeAdhereIndex.begin(),
								nodes->getInfoVecs().membrIntnlIndex.begin()))
						+ nodeStartIndx);
		cellInfoVecs.activeIntnlNodeCounts[cellRank] =
				divAuxData.tmp1InternalActiveCounts[i];
		cellInfoVecs.activeMembrNodeCounts[cellRank] =
				divAuxData.tmp1MemActiveCounts[i];
		cellInfoVecs.growthProgress[cellRank] = 0;
		cellInfoVecs.membrGrowProgress[cellRank] = 0.0;
		cellInfoVecs.isRandGrowInited[cellRank] = false;
		cellInfoVecs.lastCheckPoint[cellRank] = 0;
	}
}

void SceCells::copySecondCellArr_M() {
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	for (uint i = 0; i < divAuxData.toBeDivideCount; i++) {
		uint cellRank = allocPara_m.currentActiveCellCount + i;
		uint nodeStartIndx = cellRank * maxAllNodePerCell
				+ allocPara_m.bdryNodeCount;
		uint tmpStartIndx = i * maxAllNodePerCell;
		uint tmpEndIndx = (i + 1) * maxAllNodePerCell;
		thrust::constant_iterator<int> noAdhesion(-1), noAdhesion2(-1);
		thrust::copy(
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos2_M.begin(),
								divAuxData.tmpYPos2_M.begin(),
								divAuxData.tmpIsActive2_M.begin(), noAdhesion,
								noAdhesion2)) + tmpStartIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(divAuxData.tmpXPos2_M.begin(),
								divAuxData.tmpYPos2_M.begin(),
								divAuxData.tmpIsActive2_M.begin(), noAdhesion,
								noAdhesion2)) + tmpEndIndx,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								nodes->getInfoVecs().nodeLocX.begin(),
								nodes->getInfoVecs().nodeLocY.begin(),
								nodes->getInfoVecs().nodeIsActive.begin(),
								nodes->getInfoVecs().nodeAdhereIndex.begin(),
								nodes->getInfoVecs().membrIntnlIndex.begin()))
						+ nodeStartIndx);
		cellInfoVecs.activeIntnlNodeCounts[cellRank] =
				divAuxData.tmp2InternalActiveCounts[i];
		cellInfoVecs.activeMembrNodeCounts[cellRank] =
				divAuxData.tmp2MemActiveCounts[i];
		cellInfoVecs.growthProgress[cellRank] = 0;
		cellInfoVecs.membrGrowProgress[cellRank] = 0;
		cellInfoVecs.isRandGrowInited[cellRank] = false;
		cellInfoVecs.lastCheckPoint[cellRank] = 0;
	}
}

void SceCells::updateActiveCellCount_M() {
	allocPara_m.currentActiveCellCount = allocPara_m.currentActiveCellCount
			+ divAuxData.toBeDivideCount;
	nodes->setActiveCellCount(allocPara_m.currentActiveCellCount);
}

void SceCells::markIsDivideFalse_M() {
	thrust::fill(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin()
					+ allocPara_m.currentActiveCellCount, false);
}

void SceCells::adjustNodeVel_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ allocPara_m.bdryNodeCount + totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			ForceZero());
}

void SceCells::moveNodes_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells + allocPara_m.bdryNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin())),
//Ali		SaxpyFunctorDim2(dt));
			SaxpyFunctorDim2_Damp(dt,Damp_Coef));   //Ali
}

void SceCells::applyMemForce_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0), iBegin1(0);
        //Ali
        thrust::fill(cellInfoVecs.Cell_Time.begin(),cellInfoVecs.Cell_Time.begin() +allocPara_m.currentActiveCellCount,curTime);
        
        
       //Ali 
        
        thrust::device_vector<double>::iterator  MinX_Itr=thrust::min_element(nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount,
                                              nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
        thrust::device_vector<double>::iterator  MaxX_Itr=thrust::max_element(nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount,
                                              nodes->getInfoVecs().nodeLocX.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
        thrust::device_vector<double>::iterator  MinY_Itr=thrust::min_element(nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount,
                                              nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
        thrust::device_vector<double>::iterator  MaxY_Itr=thrust::max_element(nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount,
                                              nodes->getInfoVecs().nodeLocY.begin()+ allocPara_m.bdryNodeCount+ totalNodeCountForActiveCells) ;
        MinX= *MinX_Itr ; 
        MaxX= *MaxX_Itr ; 
        MinY= *MinY_Itr ; 
        MaxY= *MaxY_Itr ;  
        cout<< "# of boundary nodes"<< allocPara_m.bdryNodeCount<<endl ;
        cout<< "# of total active nodes"<<totalNodeCountForActiveCells <<endl ;

        cout<<"The minimum location in X is="<<MinX<< endl;  
        cout<<"The maximum location in X is="<<MaxX<< endl;  
        cout<<"The minimum location in Y is="<<MinY<< endl;  
        cout<<"The maximum location in Y is="<<MaxY<< endl;  
        //Ali 
	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));
/**Ali Comment start

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().membrTensionMag.begin(),
							nodes->getInfoVecs().membrTenMagRi.begin(),
							nodes->getInfoVecs().membrLinkRiMidX.begin(),
							nodes->getInfoVecs().membrLinkRiMidY.begin(),
							nodes->getInfoVecs().membrBendLeftX.begin(),
							nodes->getInfoVecs().membrBendLeftY.begin(),
							nodes->getInfoVecs().membrBendRightX.begin(),
							nodes->getInfoVecs().membrBendRightY.begin()))
					+ allocPara_m.bdryNodeCount,
			AddMembrForce(allocPara_m.bdryNodeCount, maxAllNodePerCell,
					nodeLocXAddr, nodeLocYAddr, nodeIsActiveAddr));

**/
// Ali comment end
//Ali 

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                                                        thrust::make_permutation_iterator(
									cellInfoVecs.Cell_Time.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),

							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
                                                        thrust::make_permutation_iterator(
									cellInfoVecs.Cell_Time.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeVelY.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin(),
							nodes->getInfoVecs().membrTensionMag.begin(),
							nodes->getInfoVecs().membrTenMagRi.begin(),
							nodes->getInfoVecs().membrLinkRiMidX.begin(),
							nodes->getInfoVecs().membrLinkRiMidY.begin(),
							nodes->getInfoVecs().membrBendLeftX.begin(),
							nodes->getInfoVecs().membrBendLeftY.begin(),
							nodes->getInfoVecs().membrBendRightX.begin(),
							nodes->getInfoVecs().membrBendRightY.begin()))
					+ allocPara_m.bdryNodeCount,
			AddMembrForce(allocPara_m.bdryNodeCount, maxAllNodePerCell,
					nodeLocXAddr, nodeLocYAddr, nodeIsActiveAddr));


//Ali



	double* bendLeftXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendLeftX[0]));
	double* bendLeftYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendLeftY[0]));
	double* bendRightXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendRightX[0]));
	double* bendRightYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().membrBendRightY[0]));

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin1,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin1,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin1,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin1,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin1,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin1,
									ModuloFunctor(maxAllNodePerCell)),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			AddMembrBend(maxAllNodePerCell, nodeIsActiveAddr, bendLeftXAddr,
					bendLeftYAddr, bendRightXAddr, bendRightYAddr));
}

void SceCells::runAblationTest(AblationEvent& ablEvent) {
	for (uint i = 0; i < ablEvent.ablationCells.size(); i++) {
		int cellRank = ablEvent.ablationCells[i].cellNum;
		std::vector<uint> removeSeq = ablEvent.ablationCells[i].nodeNums;
		cellInfoVecs.activeNodeCountOfThisCell[cellRank] =
				cellInfoVecs.activeNodeCountOfThisCell[cellRank]
						- removeSeq.size();
		nodes->removeNodes(cellRank, removeSeq);
	}
}

void SceCells::computeCenterPos_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	//uint totalMembrActiveNodeCount = thrust::reduce(
	//		cellInfoVecs.activeMembrNodeCounts.begin(),
	//		cellInfoVecs.activeMembrNodeCounts.begin()
	//				+ allocPara_m.currentActiveCellCount);
	uint totalIntnlActiveNodeCount = thrust::reduce(
			cellInfoVecs.activeIntnlNodeCounts.begin(),
			cellInfoVecs.activeIntnlNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(iBegin,
									DivideFunctor(
											allocPara_m.maxAllNodePerCell)),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeCellType.begin()))
					+ allocPara_m.bdryNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.cellRanks.begin(),
							cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin())),
			ActiveAndIntnl());

	thrust::reduce_by_key(cellNodeInfoVecs.cellRanks.begin(),
			cellNodeInfoVecs.cellRanks.begin() + totalIntnlActiveNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellNodeInfoVecs.activeXPoss.begin(),
							cellNodeInfoVecs.activeYPoss.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin())),
			thrust::equal_to<uint>(), CVec2Add());
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.activeIntnlNodeCounts.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin())), CVec2Divide());
}

void SceCells::growAtRandom_M(double dt) {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;

	randomizeGrowth_M();

	updateGrowthProgress_M();

	decideIsScheduleToGrow_M();

	//computeCellTargetLength_M();

	//computeDistToCellCenter_M();

	//findMinAndMaxDistToCenter_M();

	//computeLenDiffExpCur_M();

	//stretchCellGivenLenDiff_M();

	addPointIfScheduledToGrow_M();

	adjustGrowthInfo_M();
}

void SceCells::divide2D_M() {
	bool isDivisionPresent = decideIfGoingToDivide_M();
	if (!isDivisionPresent) {
		return;
	}
	//aniDebug = true;
	copyCellsPreDivision_M();
	createTwoNewCellArr_M();
	copyFirstCellArr_M();
	copySecondCellArr_M();
	updateActiveCellCount_M();
	markIsDivideFalse_M();
	//divDebug();
}

void SceCells::distributeCellGrowthProgress_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);

	thrust::copy(
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			thrust::make_permutation_iterator(
					cellInfoVecs.growthProgress.begin(),
					make_transform_iterator(countingEnd,
							DivideFunctor(allocPara_m.maxAllNodePerCell))),
			nodes->getInfoVecs().nodeGrowPro.begin()
					+ allocPara_m.bdryNodeCount);
}

void SceCells::allComponentsMove_M() {
	moveNodes_M();
}

void SceCells::randomizeGrowth_M() {
	uint seed = time(NULL);
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin(),
							countingBegin))
					+ allocPara_m.currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthSpeed.begin(),
							cellInfoVecs.growthXDir.begin(),
							cellInfoVecs.growthYDir.begin(),
							cellInfoVecs.isRandGrowInited.begin())),
			RandomizeGrow_M(growthAuxData.randomGrowthSpeedMin,
					growthAuxData.randomGrowthSpeedMax, seed));
}

void SceCells::updateGrowthProgress_M() {
	thrust::transform(cellInfoVecs.growthSpeed.begin(),
			cellInfoVecs.growthSpeed.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin(), SaxpyFunctorWithMaxOfOne(dt));
}

void SceCells::decideIsScheduleToGrow_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isScheduledToGrow.begin(),
			PtCondiOp(miscPara.growThreshold));
}

void SceCells::computeCellTargetLength_M() {
	thrust::transform(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.expectedLength.begin(),
			CompuTarLen(bioPara.cellInitLength, bioPara.cellFinalLength));
}

void SceCells::computeDistToCellCenter_M() {
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(totalNodeCountForActiveCells);
	uint endIndx = allocPara_m.bdryNodeCount + totalNodeCountForActiveCells;
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							nodes->getInfoVecs().nodeLocX.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeLocY.begin()
									+ allocPara_m.bdryNodeCount,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ allocPara_m.bdryNodeCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(
									cellInfoVecs.centerCoordX.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.centerCoordY.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(
													allocPara_m.maxAllNodePerCell))),
							nodes->getInfoVecs().nodeLocX.begin() + endIndx,
							nodes->getInfoVecs().nodeLocY.begin() + endIndx,
							nodes->getInfoVecs().nodeIsActive.begin()
									+ endIndx)),
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(), CompuDist());
}

void SceCells::findMinAndMaxDistToCenter_M() {
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.smallestDistance.begin(), thrust::equal_to<uint>(),
			thrust::minimum<double>());

// for nodes of each cell, find the maximum distance from the node to the corresponding
// cell center along the pre-defined growth direction.

	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(allocPara_m.maxAllNodePerCell))
					+ totalNodeCountForActiveCells,
			cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.biggestDistance.begin(), thrust::equal_to<uint>(),
			thrust::maximum<double>());
}

void SceCells::computeLenDiffExpCur_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.expectedLength.begin(),
							cellInfoVecs.smallestDistance.begin(),
							cellInfoVecs.biggestDistance.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.lengthDifference.begin(), CompuDiff());
}

void SceCells::stretchCellGivenLenDiff_M() {
	uint count = allocPara_m.maxAllNodePerCell;
	uint bdry = allocPara_m.bdryNodeCount;
	uint actCount = totalNodeCountForActiveCells;
	uint all = bdry + actCount;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(actCount);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(count))),
							nodes->getInfoVecs().nodeVelX.begin() + bdry,
							nodes->getInfoVecs().nodeVelY.begin() + bdry,
							make_transform_iterator(iBegin,
									ModuloFunctor(count)))),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellNodeInfoVecs.distToCenterAlongGrowDir.begin()
									+ actCount,
							make_permutation_iterator(
									cellInfoVecs.lengthDifference.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthXDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(count))),
							make_permutation_iterator(
									cellInfoVecs.growthYDir.begin(),
									make_transform_iterator(iEnd,
											DivideFunctor(count))),
							nodes->getInfoVecs().nodeVelX.begin() + all,
							nodes->getInfoVecs().nodeVelY.begin() + all,
							make_transform_iterator(iEnd,
									ModuloFunctor(count)))),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().nodeVelX.begin() + bdry,
							nodes->getInfoVecs().nodeVelY.begin() + bdry)),
			ApplyStretchForce_M(bioPara.elongationCoefficient,
					allocPara_m.maxMembrNodePerCell));
}

void SceCells::addPointIfScheduledToGrow_M() {
	uint seed = time(NULL);
	uint activeCellCount = allocPara_m.currentActiveCellCount;
	thrust::counting_iterator<uint> iBegin(0);
	thrust::counting_iterator<uint> iEnd(activeCellCount);
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.centerCoordX.begin(),
							cellInfoVecs.centerCoordY.begin(), iBegin,
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.isScheduledToGrow.begin()
									+ activeCellCount,
							cellInfoVecs.activeIntnlNodeCounts.begin()
									+ activeCellCount,
							cellInfoVecs.centerCoordX.begin() + activeCellCount,
							cellInfoVecs.centerCoordY.begin() + activeCellCount,
							iEnd,
							cellInfoVecs.lastCheckPoint.begin()
									+ activeCellCount)),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.lastCheckPoint.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin())),
			AddPtOp_M(seed, miscPara.addNodeDistance, miscPara.growThreshold,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress,
					growthAuxData.nodeIsActiveAddress));
}

bool SceCells::decideIfGoingToDivide_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isDividing.begin(),
			CompuIsDivide_M(allocPara_m.maxIntnlNodePerCell));
	// sum all bool values which indicate whether the cell is going to divide.
	// toBeDivideCount is the total number of cells going to divide.
	divAuxData.toBeDivideCount = thrust::reduce(cellInfoVecs.isDividing.begin(),
			cellInfoVecs.isDividing.begin()
					+ allocPara_m.currentActiveCellCount, (uint) (0));
	if (divAuxData.toBeDivideCount > 0) {
		return true;
	} else {
		return false;
	}
}

AniRawData SceCells::obtainAniRawData(AnimationCriteria& aniCri) {
	uint activeCellCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	uint beginIndx = allocPara_m.bdryNodeCount;

	AniRawData rawAniData;
	//cout << "size of potential pairs = " << pairs.size() << endl;

// unordered_map is more efficient than map, but it is a c++ 11 feature
// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

	uint maxActiveNode = activeCellCount * maxNodePerCell;
	thrust::host_vector<double> hostTmpVectorLocX(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorLocY(maxActiveNode);
	thrust::host_vector<bool> hostIsActiveVec(maxActiveNode);
	thrust::host_vector<int> hostBondVec(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorTenMag(maxActiveNode);

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin()))
					+ maxActiveNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(hostTmpVectorLocX.begin(),
							hostTmpVectorLocY.begin(), hostIsActiveVec.begin(),
							hostBondVec.begin(), hostTmpVectorTenMag.begin())));

	thrust::host_vector<uint> curActiveMemNodeCounts =
			cellInfoVecs.activeMembrNodeCounts;

	CVector tmpPos;
	uint index1;
	int index2;
	std::vector<BondInfo> bondInfoVec;

	double node1X, node1Y;
	double node2X, node2Y;
	double aniVal;

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < maxMemNodePerCell; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (hostIsActiveVec[index1] == true) {
				index2 = hostBondVec[index1];
				if (index2 > index1 && index2 != -1) {
					BondInfo bond;
					bond.cellRank1 = i;
					bond.pos1 = CVector(hostTmpVectorLocX[index1],
							hostTmpVectorLocY[index1], 0);
					bond.cellRank2 = (index2 - beginIndx) / maxNodePerCell;
					bond.pos2 = CVector(hostTmpVectorLocX[index2],
							hostTmpVectorLocY[index2], 0);
					bondInfoVec.push_back(bond);
				}
			}
		}
	}

	rawAniData.bondsArr = bondInfoVec;

	uint curIndex = 0;

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < curActiveMemNodeCounts[i]; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (j == curActiveMemNodeCounts[i] - 1) {
				index2 = beginIndx + i * maxNodePerCell;
			} else {
				index2 = beginIndx + i * maxNodePerCell + j + 1;
			}

			if (hostIsActiveVec[index1] == true
					&& hostIsActiveVec[index2] == true) {
				node1X = hostTmpVectorLocX[index1];
				node1Y = hostTmpVectorLocY[index1];
				node2X = hostTmpVectorLocX[index2];
				node2Y = hostTmpVectorLocY[index2];
				IndexMap::iterator it = locIndexToAniIndexMap.find(index1);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index1, curIndex));
					curIndex++;
					tmpPos = CVector(node1X, node1Y, 0);
					//aniVal = hostTmpVectorNodeType[index1];
					aniVal = hostTmpVectorTenMag[index1];
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);
				}
				it = locIndexToAniIndexMap.find(index2);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index2, curIndex));
					curIndex++;
					tmpPos = CVector(node2X, node2Y, 0);
					//aniVal = hostTmpVectorNodeType[index2];
					aniVal = hostTmpVectorTenMag[index2];
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);
				}

				it = locIndexToAniIndexMap.find(index1);
				uint aniIndex1 = it->second;
				it = locIndexToAniIndexMap.find(index2);
				uint aniIndex2 = it->second;

				LinkAniData linkData;
				linkData.node1Index = aniIndex1;
				linkData.node2Index = aniIndex2;
				rawAniData.memLinks.push_back(linkData);
			}
		}
	}

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < allocPara_m.maxIntnlNodePerCell; j++) {
			for (uint k = j + 1; k < allocPara_m.maxIntnlNodePerCell; k++) {
				index1 = i * maxNodePerCell + maxMemNodePerCell + j;
				index2 = i * maxNodePerCell + maxMemNodePerCell + k;
				if (hostIsActiveVec[index1] && hostIsActiveVec[index2]) {
					node1X = hostTmpVectorLocX[index1];
					node1Y = hostTmpVectorLocY[index1];
					node2X = hostTmpVectorLocX[index2];
					node2Y = hostTmpVectorLocY[index2];
					if (aniCri.isPairQualify_M(node1X, node1Y, node2X,
							node2Y)) {
						IndexMap::iterator it = locIndexToAniIndexMap.find(
								index1);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index1, curIndex));
							curIndex++;
							tmpPos = CVector(node1X, node1Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = -1;
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}
						it = locIndexToAniIndexMap.find(index2);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index2, curIndex));
							curIndex++;
							tmpPos = CVector(node2X, node2Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = -1;
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}

						it = locIndexToAniIndexMap.find(index1);
						uint aniIndex1 = it->second;
						it = locIndexToAniIndexMap.find(index2);
						uint aniIndex2 = it->second;

						LinkAniData linkData;
						linkData.node1Index = aniIndex1;
						linkData.node2Index = aniIndex2;
						rawAniData.internalLinks.push_back(linkData);
					}
				}
			}
		}
	}
	return rawAniData;
}

AniRawData SceCells::obtainAniRawDataGivenCellColor(vector<double>& cellColors,
		AnimationCriteria& aniCri) {

	uint activeCellCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	uint beginIndx = allocPara_m.bdryNodeCount;

	assert(cellColors.size() >= activeCellCount);

	AniRawData rawAniData;
	//cout << "size of potential pairs = " << pairs.size() << endl;

	// unordered_map is more efficient than map, but it is a c++ 11 feature
	// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

	uint maxActiveNode = activeCellCount * maxNodePerCell;
	thrust::host_vector<double> hostTmpVectorLocX(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorLocY(maxActiveNode);
	thrust::host_vector<bool> hostIsActiveVec(maxActiveNode);
	thrust::host_vector<int> hostBondVec(maxActiveNode);
	thrust::host_vector<double> hostTmpVectorTenMag(maxActiveNode);

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeLocX.begin(),
							nodes->getInfoVecs().nodeLocY.begin(),
							nodes->getInfoVecs().nodeIsActive.begin(),
							nodes->getInfoVecs().nodeAdhereIndex.begin(),
							nodes->getInfoVecs().membrTensionMag.begin()))
					+ maxActiveNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(hostTmpVectorLocX.begin(),
							hostTmpVectorLocY.begin(), hostIsActiveVec.begin(),
							hostBondVec.begin(), hostTmpVectorTenMag.begin())));

	thrust::host_vector<uint> curActiveMemNodeCounts =
			cellInfoVecs.activeMembrNodeCounts;

	CVector tmpPos;
	uint index1;
	int index2;
	std::vector<BondInfo> bondInfoVec;

	double node1X, node1Y;
	double node2X, node2Y;
	double aniVal;

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < maxMemNodePerCell; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (hostIsActiveVec[index1] == true) {
				index2 = hostBondVec[index1];
				if (index2 > index1 && index2 != -1) {
					BondInfo bond;
					bond.cellRank1 = i;
					bond.pos1 = CVector(hostTmpVectorLocX[index1],
							hostTmpVectorLocY[index1], 0);
					bond.cellRank2 = (index2 - beginIndx) / maxNodePerCell;
					bond.pos2 = CVector(hostTmpVectorLocX[index2],
							hostTmpVectorLocY[index2], 0);
					bondInfoVec.push_back(bond);
				}
			}
		}
	}

	rawAniData.bondsArr = bondInfoVec;

	uint curIndex = 0;

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < curActiveMemNodeCounts[i]; j++) {
			index1 = beginIndx + i * maxNodePerCell + j;
			if (j == curActiveMemNodeCounts[i] - 1) {
				index2 = beginIndx + i * maxNodePerCell;
			} else {
				index2 = beginIndx + i * maxNodePerCell + j + 1;
			}

			if (hostIsActiveVec[index1] == true
					&& hostIsActiveVec[index2] == true) {
				node1X = hostTmpVectorLocX[index1];
				node1Y = hostTmpVectorLocY[index1];
				node2X = hostTmpVectorLocX[index2];
				node2Y = hostTmpVectorLocY[index2];
				IndexMap::iterator it = locIndexToAniIndexMap.find(index1);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index1, curIndex));
					curIndex++;
					tmpPos = CVector(node1X, node1Y, 0);
					//aniVal = hostTmpVectorNodeType[index1];
					aniVal = cellColors[i];
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);
				}
				it = locIndexToAniIndexMap.find(index2);
				if (it == locIndexToAniIndexMap.end()) {
					locIndexToAniIndexMap.insert(
							std::pair<uint, uint>(index2, curIndex));
					curIndex++;
					tmpPos = CVector(node2X, node2Y, 0);
					//aniVal = hostTmpVectorNodeType[index2];
					aniVal = cellColors[i];
					rawAniData.aniNodePosArr.push_back(tmpPos);
					rawAniData.aniNodeVal.push_back(aniVal);
				}

				it = locIndexToAniIndexMap.find(index1);
				uint aniIndex1 = it->second;
				it = locIndexToAniIndexMap.find(index2);
				uint aniIndex2 = it->second;

				LinkAniData linkData;
				linkData.node1Index = aniIndex1;
				linkData.node2Index = aniIndex2;
				rawAniData.memLinks.push_back(linkData);
			}
		}
	}

	for (uint i = 0; i < activeCellCount; i++) {
		for (uint j = 0; j < allocPara_m.maxIntnlNodePerCell; j++) {
			for (uint k = j + 1; k < allocPara_m.maxIntnlNodePerCell; k++) {
				index1 = i * maxNodePerCell + maxMemNodePerCell + j;
				index2 = i * maxNodePerCell + maxMemNodePerCell + k;
				if (hostIsActiveVec[index1] && hostIsActiveVec[index2]) {
					node1X = hostTmpVectorLocX[index1];
					node1Y = hostTmpVectorLocY[index1];
					node2X = hostTmpVectorLocX[index2];
					node2Y = hostTmpVectorLocY[index2];
					if (aniCri.isPairQualify_M(node1X, node1Y, node2X,
							node2Y)) {
						IndexMap::iterator it = locIndexToAniIndexMap.find(
								index1);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index1, curIndex));
							curIndex++;
							tmpPos = CVector(node1X, node1Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = cellColors[i];
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}
						it = locIndexToAniIndexMap.find(index2);
						if (it == locIndexToAniIndexMap.end()) {
							locIndexToAniIndexMap.insert(
									std::pair<uint, uint>(index2, curIndex));
							curIndex++;
							tmpPos = CVector(node2X, node2Y, 0);
							//aniVal = hostTmpVectorNodeType[index1];
							aniVal = cellColors[i];
							rawAniData.aniNodePosArr.push_back(tmpPos);
							rawAniData.aniNodeVal.push_back(aniVal);
						}

						it = locIndexToAniIndexMap.find(index1);
						uint aniIndex1 = it->second;
						it = locIndexToAniIndexMap.find(index2);
						uint aniIndex2 = it->second;

						LinkAniData linkData;
						linkData.node1Index = aniIndex1;
						linkData.node2Index = aniIndex2;
						rawAniData.internalLinks.push_back(linkData);
					}
				}
			}
		}
	}
	return rawAniData;
}

void SceCells::copyInitActiveNodeCount_M(
		std::vector<uint>& initMembrActiveNodeCounts,
		std::vector<uint>& initIntnlActiveNodeCounts,
		std::vector<double> &initGrowProgVec) {
	assert(
			initMembrActiveNodeCounts.size()
					== initIntnlActiveNodeCounts.size());
	totalNodeCountForActiveCells = initMembrActiveNodeCounts.size()
			* allocPara_m.maxAllNodePerCell;

	thrust::copy(initMembrActiveNodeCounts.begin(),
			initMembrActiveNodeCounts.end(),
			cellInfoVecs.activeMembrNodeCounts.begin());
	thrust::copy(initIntnlActiveNodeCounts.begin(),
			initIntnlActiveNodeCounts.end(),
			cellInfoVecs.activeIntnlNodeCounts.begin());
	thrust::copy(initGrowProgVec.begin(), initGrowProgVec.end(),
			cellInfoVecs.growthProgress.begin());
}

void SceCells::myDebugFunction() {

	uint maxActiveNodeCount = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxActiveCellCount = allocPara_m.currentActiveCellCount;
	std::cout << "totalNodeCountforActiveCells: "
			<< totalNodeCountForActiveCells << std::endl;
	std::cout << "maxAllNodePerCell: " << allocPara_m.maxAllNodePerCell
			<< std::endl;
	std::cout << "maxActiveCellCount: " << maxActiveCellCount << std::endl;
	std::cout << "bdryNodeCount: " << allocPara_m.bdryNodeCount << std::endl;

	std::cout << "grow threshold: " << miscPara.growThreshold << std::endl;

	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.growthProgress[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.isScheduledToGrow[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.lastCheckPoint[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeCount; i++) {
		if (nodes->getInfoVecs().nodeIsActive[i]
				&& nodes->getInfoVecs().nodeCellType[i] == CellIntnl) {
			std::cout << nodes->getInfoVecs().nodeVelX[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.activeIntnlNodeCounts[i] << " ";
	}
	std::cout << std::endl;

	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.expectedLength[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.smallestDistance[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.biggestDistance[i] << " ";
	}
	std::cout << std::endl;

	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.lengthDifference[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.centerCoordX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.centerCoordY[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.growthXDir[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveCellCount; i++) {
		std::cout << cellInfoVecs.growthYDir[i] << " ";
	}
	std::cout << std::endl;

	int jj;
	std::cin >> jj;
}

void SceCells::divDebug() {

	std::cout << "tmpIsActive_M: ";
	for (uint i = 0; i < divAuxData.tmpIsActive_M.size(); i++) {
		std::cout << divAuxData.tmpIsActive_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpNodePosX_M: ";
	for (uint i = 0; i < divAuxData.tmpNodePosX_M.size(); i++) {
		std::cout << divAuxData.tmpNodePosX_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpNodePosY_M : ";
	for (uint i = 0; i < divAuxData.tmpNodePosY_M.size(); i++) {
		std::cout << divAuxData.tmpNodePosY_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpCellRank_M : ";
	for (uint i = 0; i < divAuxData.tmpCellRank_M.size(); i++) {
		std::cout << divAuxData.tmpCellRank_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpDivDirX_M : ";
	for (uint i = 0; i < divAuxData.tmpDivDirX_M.size(); i++) {
		std::cout << divAuxData.tmpDivDirX_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpDivDirY_M : ";
	for (uint i = 0; i < divAuxData.tmpDivDirY_M.size(); i++) {
		std::cout << divAuxData.tmpDivDirY_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpCenterPosX_M : ";
	for (uint i = 0; i < divAuxData.tmpCenterPosX_M.size(); i++) {
		std::cout << divAuxData.tmpCenterPosX_M[i] << " ";
	}
	std::cout << std::endl;
	std::cout << "tmpCenterPosY_M : ";
	for (uint i = 0; i < divAuxData.tmpCenterPosY_M.size(); i++) {
		std::cout << divAuxData.tmpCenterPosY_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpIsActive1_M : ";
	for (uint i = 0; i < divAuxData.tmpIsActive1_M.size(); i++) {
		std::cout << divAuxData.tmpIsActive1_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpXPos1_M : ";
	for (uint i = 0; i < divAuxData.tmpXPos1_M.size(); i++) {
		std::cout << divAuxData.tmpXPos1_M[i] << " ";
		if (i > 0 && i < allocPara_m.maxMembrNodePerCell
				&& divAuxData.tmpIsActive1_M[i]
				&& divAuxData.tmpIsActive1_M[i - 1]
				&& fabs(divAuxData.tmpXPos1_M[i] - divAuxData.tmpXPos1_M[i - 1])
						> 0.1) {
			std::cout << "11111111111111111111111, " << i << std::endl;
			int jj;
			cin >> jj;
		}
	}
	std::cout << std::endl;
	std::cout << "XPos1_onDevice : ";
	for (uint i = 0; i < divAuxData.tmpCellRank_M.size(); i++) {
		for (uint j = 0; j < allocPara_m.maxAllNodePerCell; j++) {
			uint index = divAuxData.tmpCellRank_M[i]
					* allocPara_m.maxAllNodePerCell + j;
			std::cout << nodes->getInfoVecs().nodeLocX[index] << " ";
		}
	}
	std::cout << std::endl;

	std::cout << "tmpYPos1_M : ";
	for (uint i = 0; i < divAuxData.tmpYPos1_M.size(); i++) {
		std::cout << divAuxData.tmpYPos1_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpIsActive2_M: ";
	for (uint i = 0; i < divAuxData.tmpIsActive2_M.size(); i++) {
		std::cout << divAuxData.tmpIsActive2_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmpXPos2_M : ";
	for (uint i = 0; i < divAuxData.tmpXPos2_M.size(); i++) {
		std::cout << divAuxData.tmpXPos2_M[i] << " ";
		if (i > 0 && i < allocPara_m.maxMembrNodePerCell
				&& divAuxData.tmpIsActive2_M[i]
				&& divAuxData.tmpIsActive2_M[i - 1]
				&& fabs(divAuxData.tmpXPos2_M[i] - divAuxData.tmpXPos2_M[i - 1])
						> 0.1) {
			std::cout << "2222222222222222222, " << i << std::endl;
			int jj;
			cin >> jj;
		}
	}
	std::cout << std::endl;
	std::cout << "tmpYPos2_M : ";
	for (uint i = 0; i < divAuxData.tmpYPos2_M.size(); i++) {
		std::cout << divAuxData.tmpYPos2_M[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp1InternalActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp1InternalActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp1InternalActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp2InternalActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp2InternalActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp2InternalActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp1MemActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp1MemActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp1MemActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "tmp2MemActiveCounts: ";
	for (uint i = 0; i < divAuxData.tmp2MemActiveCounts.size(); i++) {
		std::cout << divAuxData.tmp2MemActiveCounts[i] << " ";
	}
	std::cout << std::endl;

	int jj;
	std::cin >> jj;
}

void SceCells::adjustGrowthInfo_M() {
	uint halfMax = allocPara_m.maxIntnlNodePerCell / 2;
	thrust::transform_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeIntnlNodeCounts.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.isScheduledToGrow.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isScheduledToGrow.begin(),
							cellInfoVecs.growthProgress.begin(),
							cellInfoVecs.lastCheckPoint.begin())),
			AdjustGrowth(halfMax), thrust::identity<bool>());
}

VtkAnimationData SceCells::outputVtkData(AniRawData& rawAniData,
		AnimationCriteria& aniCri) {
	VtkAnimationData vtkData;
	for (uint i = 0; i < rawAniData.aniNodePosArr.size(); i++) {
		PointAniData ptAniData;
		ptAniData.pos = rawAniData.aniNodePosArr[i];
		ptAniData.colorScale = rawAniData.aniNodeVal[i];
		vtkData.pointsAniData.push_back(ptAniData);
	}
	for (uint i = 0; i < rawAniData.internalLinks.size(); i++) {
		LinkAniData linkData = rawAniData.internalLinks[i];
		vtkData.linksAniData.push_back(linkData);
	}
	for (uint i = 0; i < rawAniData.memLinks.size(); i++) {
		LinkAniData linkData = rawAniData.memLinks[i];
		vtkData.linksAniData.push_back(linkData);
	}
	vtkData.isArrowIncluded = false;
	return vtkData;
}

void SceCells::copyToGPUConstMem() {
	double pI_CPU = acos(-1.0);
	double minLengthCPU =
			globalConfigVars.getConfigValue("MinLength").toDouble();
	cudaMemcpyToSymbol(minLength, &minLengthCPU, sizeof(double));
	double minDivisorCPU =
			globalConfigVars.getConfigValue("MinDivisor").toDouble();
	cudaMemcpyToSymbol(minDivisor, &minDivisorCPU, sizeof(double));
	cudaMemcpyToSymbol(membrEquLen, &membrPara.membrEquLenCPU, sizeof(double));
	cudaMemcpyToSymbol(membrStiff, &membrPara.membrStiffCPU, sizeof(double));
	cudaMemcpyToSymbol(pI, &pI_CPU, sizeof(double));

	cudaMemcpyToSymbol(bendCoeff, &membrPara.membrBendCoeff, sizeof(double));

	cudaMemcpyToSymbol(F_Ext_Incline_M2, &membrPara.F_Ext_Incline, sizeof(double)); //Ali
      
	uint maxAllNodePerCellCPU = globalConfigVars.getConfigValue(
			"MaxAllNodeCountPerCell").toInt();
	uint maxMembrNodePerCellCPU = globalConfigVars.getConfigValue(
			"MaxMembrNodeCountPerCell").toInt();
	uint maxIntnlNodePerCellCPU = globalConfigVars.getConfigValue(
			"MaxIntnlNodeCountPerCell").toInt();

	cudaMemcpyToSymbol(maxAllNodePerCell, &maxAllNodePerCellCPU, sizeof(uint));
	cudaMemcpyToSymbol(maxMembrPerCell, &maxMembrNodePerCellCPU, sizeof(uint));
	cudaMemcpyToSymbol(maxIntnlPerCell, &maxIntnlNodePerCellCPU, sizeof(uint));

	double sceIntnlBParaCPU_M[5];
	double sceIntraParaCPU_M[5];
	double sceIntraParaDivCPU_M[5];

	double U0_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_U0").toDouble();
	double V0_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_V0").toDouble();
	double k1_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_k1").toDouble();
	double k2_IntnlB =
			globalConfigVars.getConfigValue("SceIntnlB_k2").toDouble();
	double intnlBEffectiveRange = globalConfigVars.getConfigValue(
			"IntnlBEffectRange").toDouble();
	sceIntnlBParaCPU_M[0] = U0_IntnlB;
	sceIntnlBParaCPU_M[1] = V0_IntnlB;
	sceIntnlBParaCPU_M[2] = k1_IntnlB;
	sceIntnlBParaCPU_M[3] = k2_IntnlB;
	sceIntnlBParaCPU_M[4] = intnlBEffectiveRange;


        
 
	//////////////////////
	//// Block 3 /////////
	//////////////////////
	double U0_Intra =
			globalConfigVars.getConfigValue("IntraCell_U0").toDouble();
	double V0_Intra =
			globalConfigVars.getConfigValue("IntraCell_V0").toDouble();
	double k1_Intra =
			globalConfigVars.getConfigValue("IntraCell_k1").toDouble();
	double k2_Intra =
			globalConfigVars.getConfigValue("IntraCell_k2").toDouble();
	double intraLinkEffectiveRange = globalConfigVars.getConfigValue(
			"IntraEffectRange").toDouble();
	sceIntraParaCPU_M[0] = U0_Intra;
	sceIntraParaCPU_M[1] = V0_Intra;
	sceIntraParaCPU_M[2] = k1_Intra;
	sceIntraParaCPU_M[3] = k2_Intra;
	sceIntraParaCPU_M[4] = intraLinkEffectiveRange;

	//////////////////////
	//// Block 4 /////////
	//////////////////////
	double U0_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_U0_Div").toDouble();
	double V0_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_V0_Div").toDouble();
	double k1_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_k1_Div").toDouble();
	double k2_Intra_Div =
			globalConfigVars.getConfigValue("IntraCell_k2_Div").toDouble();
	double intraDivEffectiveRange = globalConfigVars.getConfigValue(
			"IntraDivEffectRange").toDouble();
	sceIntraParaDivCPU_M[0] = U0_Intra_Div;
	sceIntraParaDivCPU_M[1] = V0_Intra_Div;
	sceIntraParaDivCPU_M[2] = k1_Intra_Div;
	sceIntraParaDivCPU_M[3] = k2_Intra_Div;
	sceIntraParaDivCPU_M[4] = intraDivEffectiveRange;

	cudaMemcpyToSymbol(grthPrgrCriEnd_M, &growthAuxData.grthProgrEndCPU,
			sizeof(double));
	//cudaMemcpyToSymbol(grthPrgrCriVal_M, &growthPrgrCriVal, sizeof(double));
	cudaMemcpyToSymbol(sceIB_M, sceIntnlBParaCPU_M, 5 * sizeof(double));
	cudaMemcpyToSymbol(sceII_M, sceIntraParaCPU_M, 5 * sizeof(double));
	cudaMemcpyToSymbol(sceIIDiv_M, sceIntraParaDivCPU_M, 5 * sizeof(double));

	double IBDivHost[5];
	IBDivHost[0] =
			globalConfigVars.getConfigValue("SceIntnlB_U0_Div").toDouble();
	IBDivHost[1] =
			globalConfigVars.getConfigValue("SceIntnlB_V0_Div").toDouble();
	IBDivHost[2] =
			globalConfigVars.getConfigValue("SceIntnlB_k1_Div").toDouble();
	IBDivHost[3] =
			globalConfigVars.getConfigValue("SceIntnlB_k2_Div").toDouble();
	IBDivHost[4] =
			globalConfigVars.getConfigValue("IntnlBDivEffectRange").toDouble();
	cudaMemcpyToSymbol(sceIBDiv_M, IBDivHost, 5 * sizeof(double));
}

void SceCells::handleMembrGrowth_M() {
	// figure out membr growth speed
	calMembrGrowSpeed_M();
	// figure out which cells will add new point

	adjustMembrGrowSpeed_M();

	decideIfAddMembrNode_M();
// add membr nodes
	addMembrNodes_M();
	//membrDebug();
}

void SceCells::calMembrGrowSpeed_M() {
	membrPara.membrGrowCoeff = growthAuxData.prolifDecay
			* membrPara.membrGrowCoeff_Ori;
	membrPara.membrGrowLimit = growthAuxData.prolifDecay
			* membrPara.membrGrowLimit_Ori;
// reduce_by_key, find value of max tension and their index
	thrust::counting_iterator<uint> iBegin(0);
	uint maxNPerCell = allocPara_m.maxAllNodePerCell;
	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->getInfoVecs().membrTenMagRi.begin(),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxNPerCell)),
							nodes->getInfoVecs().membrLinkRiMidX.begin(),
							nodes->getInfoVecs().membrLinkRiMidY.begin())),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.maxTenRiVec.begin(),
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.maxTenRiMidXVec.begin(),
							cellInfoVecs.maxTenRiMidYVec.begin())),
			thrust::equal_to<uint>(), MaxWInfo());

	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxNPerCell))
					+ totalNodeCountForActiveCells,
			nodes->getInfoVecs().membrTensionMag.begin(),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.aveTension.begin(), thrust::equal_to<uint>(),
			thrust::plus<double>());

	thrust::transform(cellInfoVecs.aveTension.begin(),
			cellInfoVecs.aveTension.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.activeMembrNodeCounts.begin(),
			cellInfoVecs.aveTension.begin(), thrust::divides<double>());

	// linear relationship with highest tension; capped by a given value
	thrust::transform(cellInfoVecs.aveTension.begin(),
			cellInfoVecs.aveTension.begin()
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.membrGrowSpeed.begin(),
			MultiWithLimit(membrPara.membrGrowCoeff, membrPara.membrGrowLimit));
}

void SceCells::adjustMembrGrowSpeed_M() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.activeIntnlNodeCounts.begin()))
					+ allocPara_m.currentActiveCellCount,
			cellInfoVecs.membrGrowSpeed.begin(),
			AdjustMembrGrow(membrPara.growthConst_N, membrPara.initMembrCt_N,
					membrPara.initIntnlCt_N));
}

void SceCells::decideIfAddMembrNode_M() {
// decide if add membrane node given current active node count and
// membr growth progress
	uint curActCellCt = allocPara_m.currentActiveCellCount;
	thrust::transform(cellInfoVecs.membrGrowSpeed.begin(),
			cellInfoVecs.membrGrowSpeed.begin() + curActCellCt,
			cellInfoVecs.membrGrowProgress.begin(),
			cellInfoVecs.membrGrowProgress.begin(), SaxpyFunctor(dt));

	uint maxMembrNode = allocPara_m.maxMembrNodePerCell;
/**Ali	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.membrGrowProgress.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.membrGrowProgress.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin()))
					+ curActCellCt,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isMembrAddingNode.begin(),
							cellInfoVecs.membrGrowProgress.begin())),
			MemGrowFunc(maxMembrNode));
*/
         thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							   cellInfoVecs.membrGrowProgress.begin(),cellInfoVecs.maxTenRiVec.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.activeMembrNodeCounts.begin(),
							   cellInfoVecs.membrGrowProgress.begin(),cellInfoVecs.maxTenRiVec.begin()))
					+ curActCellCt,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellInfoVecs.isMembrAddingNode.begin(),
							cellInfoVecs.membrGrowProgress.begin())),
			MemGrowFunc(maxMembrNode));

}

/**
 * Add new membrane elements to cells.
 * This operation is relatively expensive because of memory rearrangement.
 */
void SceCells::addMembrNodes_M() {
	thrust::counting_iterator<uint> iBegin(0);
	uint curAcCCount = allocPara_m.currentActiveCellCount;
	uint maxNodePerCell = allocPara_m.maxAllNodePerCell;
	thrust::transform_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.maxTenRiMidXVec.begin(),
							cellInfoVecs.maxTenRiMidYVec.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(iBegin,
							cellInfoVecs.maxTenIndxVec.begin(),
							cellInfoVecs.activeMembrNodeCounts.begin(),
							cellInfoVecs.maxTenRiMidXVec.begin(),
							cellInfoVecs.maxTenRiMidYVec.begin()))
					+ curAcCCount, cellInfoVecs.isMembrAddingNode.begin(),
			cellInfoVecs.activeMembrNodeCounts.begin(),
			AddMemNode(maxNodePerCell, growthAuxData.nodeIsActiveAddress,
					growthAuxData.nodeXPosAddress,
					growthAuxData.nodeYPosAddress, growthAuxData.adhIndxAddr),
			thrust::identity<bool>());
}

void SceCells::membrDebug() {
	uint curAcCCount = allocPara_m.currentActiveCellCount;
	uint maxActiveNodeC = curAcCCount * allocPara_m.maxAllNodePerCell;
	uint maxNodePC = allocPara_m.maxAllNodePerCell;
	//uint tmp = 0;
	//for (uint i = 0; i < curAcCCount; i++) {
	//	tmp += cellInfoVecs.isMembrAddingNode[i];
	//}
	//if (tmp != 0) {
	//	tmpDebug = true;
	//}
	//if (!tmpDebug) {
	//	return;
	//}
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrTensionMag[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrTenMagRi[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrLinkRiMidX[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		if (i % maxNodePC == 0 || i % maxNodePC == 199
				|| i % maxNodePC == 200) {
			std::cout << nodes->getInfoVecs().membrLinkRiMidY[i] << " ";
		}
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendLeftX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendLeftY[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendRightX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < maxActiveNodeC; i++) {
		std::cout << nodes->getInfoVecs().membrBendRightX[i] << " ";
	}
	std::cout << std::endl;
	for (uint i = 0; i < curAcCCount; i++) {
		std::cout << "(" << cellInfoVecs.maxTenIndxVec[i] << ","
				<< cellInfoVecs.activeMembrNodeCounts[i] << ","
				<< cellInfoVecs.maxTenRiMidXVec[i] << ","
				<< cellInfoVecs.maxTenRiMidYVec[i] << ")" << std::endl;
	}
	int jj;
	std::cin >> jj;
}

void SceCells::assembleVecForTwoCells(uint i) {
	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	for (uint j = 0; j < membThreshold; j++) {
		index = i * maxAllNodePerCell + j;
		if (j < divAuxData.tmp1VecMem.size()) {
			divAuxData.tmpXPos1_M[index] = divAuxData.tmp1VecMem[j].x;
			divAuxData.tmpYPos1_M[index] = divAuxData.tmp1VecMem[j].y;
			divAuxData.tmpIsActive1_M[index] = true;
		} else {
			divAuxData.tmpIsActive1_M[index] = false;
		}
	}
	for (uint j = 0; j < membThreshold; j++) {
		index = i * maxAllNodePerCell + j;
		if (j < divAuxData.tmp2VecMem.size()) {
			divAuxData.tmpXPos2_M[index] = divAuxData.tmp2VecMem[j].x;
			divAuxData.tmpYPos2_M[index] = divAuxData.tmp2VecMem[j].y;
			divAuxData.tmpIsActive2_M[index] = true;
		} else {
			divAuxData.tmpIsActive2_M[index] = false;
		}
	}

	divAuxData.tmp1MemActiveCounts.push_back(divAuxData.tmp1VecMem.size());
	divAuxData.tmp2MemActiveCounts.push_back(divAuxData.tmp2VecMem.size());

	for (uint j = membThreshold; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		uint shift_j = j - membThreshold;
		if (shift_j < divAuxData.tmp1IntnlVec.size()) {
			divAuxData.tmpXPos1_M[index] = divAuxData.tmp1IntnlVec[shift_j].x;
			divAuxData.tmpYPos1_M[index] = divAuxData.tmp1IntnlVec[shift_j].y;
			divAuxData.tmpIsActive1_M[index] = true;
		} else {
			divAuxData.tmpIsActive1_M[index] = false;
		}
		if (shift_j < divAuxData.tmp2IntnlVec.size()) {
			divAuxData.tmpXPos2_M[index] = divAuxData.tmp2IntnlVec[shift_j].x;
			divAuxData.tmpYPos2_M[index] = divAuxData.tmp2IntnlVec[shift_j].y;
			divAuxData.tmpIsActive2_M[index] = true;
		} else {
			divAuxData.tmpIsActive2_M[index] = false;
		}
	}
	divAuxData.tmp1InternalActiveCounts.push_back(
			divAuxData.tmp1IntnlVec.size());
	divAuxData.tmp2InternalActiveCounts.push_back(
			divAuxData.tmp2IntnlVec.size());
}

void SceCells::shiftIntnlNodesByCellCenter(CVector cell1Center,
		CVector cell2Center) {
	CVector tmpCell1Center(0, 0, 0);
	for (uint j = 0; j < divAuxData.tmp1IntnlVec.size(); j++) {
		tmpCell1Center = tmpCell1Center + divAuxData.tmp1IntnlVec[j];
	}
	tmpCell1Center = tmpCell1Center / divAuxData.tmp1IntnlVec.size();
	CVector shiftVec1 = cell1Center - tmpCell1Center;
	for (uint j = 0; j < divAuxData.tmp1IntnlVec.size(); j++) {
		divAuxData.tmp1IntnlVec[j] = divAuxData.tmp1IntnlVec[j] + shiftVec1;
	}

	CVector tmpCell2Center(0, 0, 0);
	for (uint j = 0; j < divAuxData.tmp2IntnlVec.size(); j++) {
		tmpCell2Center = tmpCell2Center + divAuxData.tmp2IntnlVec[j];
	}
	tmpCell2Center = tmpCell2Center / divAuxData.tmp2IntnlVec.size();
	CVector shiftVec2 = cell2Center - tmpCell2Center;
	for (uint j = 0; j < divAuxData.tmp2IntnlVec.size(); j++) {
		divAuxData.tmp2IntnlVec[j] = divAuxData.tmp2IntnlVec[j] + shiftVec2;
	}
}

void SceCells::processMemVec(std::vector<VecVal>& tmp1,
		std::vector<VecVal>& tmp2) {
	divAuxData.tmp1VecMem.clear();
	divAuxData.tmp2VecMem.clear();

	uint membThreshold = allocPara_m.maxMembrNodePerCell;

	std::sort(tmp1.begin(), tmp1.end());
	std::sort(tmp2.begin(), tmp2.end());

	//assert(tmp1.size() < allocPara_m.maxMembrNodePerCell);
	//assert(tmp2.size() < allocPara_m.maxMembrNodePerCell);

	uint maxDivMembrNodeCount1 = allocPara_m.maxMembrNodePerCell - tmp1.size();
	uint maxDivMembrNodeCount2 = allocPara_m.maxMembrNodePerCell - tmp2.size();

	std::vector<CVector> ptsBetween1, ptsBetween2;

	// if size is less than 1, the situation would have already been very bad.
	// Just keep this statement so no seg fault would happen.
	if (tmp1.size() >= 1) {
		ptsBetween1 = obtainPtsBetween(tmp1[tmp1.size() - 1].vec, tmp1[0].vec,
				memNewSpacing, maxDivMembrNodeCount1);
	}
	// if size is less than 1, the situation would have already been very bad.
	// Just keep this statement so no seg fault would happen.
	if (tmp2.size() >= 1) {
		ptsBetween2 = obtainPtsBetween(tmp2[tmp2.size() - 1].vec, tmp2[0].vec,
				memNewSpacing, maxDivMembrNodeCount2);
	}

	for (uint j = 0; j < tmp1.size(); j++) {
		divAuxData.tmp1VecMem.push_back(tmp1[j].vec);
	}
	for (uint j = 0; j < tmp2.size(); j++) {
		divAuxData.tmp2VecMem.push_back(tmp2[j].vec);
	}
	for (uint j = 0; j < ptsBetween1.size(); j++) {
		divAuxData.tmp1VecMem.push_back(ptsBetween1[j]);
	}
	for (uint j = 0; j < ptsBetween2.size(); j++) {
		divAuxData.tmp2VecMem.push_back(ptsBetween2[j]);
	}

	assert(divAuxData.tmp1VecMem.size() <= membThreshold);
	assert(divAuxData.tmp2VecMem.size() <= membThreshold);
}

void SceCells::obtainMembrAndIntnlNodes(uint i, vector<CVector>& membrNodes,
		vector<CVector>& intnlNodes) {
	membrNodes.clear();
	intnlNodes.clear();

	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (divAuxData.tmpIsActive_M[index] != true) {
			continue;
		}
		double posX = divAuxData.tmpNodePosX_M[index];
		double posY = divAuxData.tmpNodePosY_M[index];
		if (j < membThreshold) {
			// means node type is membrane
			CVector memPos(posX, posY, 0);
			membrNodes.push_back(memPos);
		} else {
			CVector intnlPos(posX, posY, 0);
			intnlNodes.push_back(intnlPos);
		}
	}
}

CVector SceCells::obtainCenter(uint i) {
	double oldCenterX = divAuxData.tmpCenterPosX_M[i];
	double oldCenterY = divAuxData.tmpCenterPosY_M[i];
	CVector centerPos(oldCenterX, oldCenterY, 0);
	return centerPos;
}

CVector SceCells::calDivDir_MajorAxis(CVector center,
		vector<CVector>& membrNodes, double& lenAlongMajorAxis) {
// not the optimal algorithm but easy to code
	double maxDiff = 0;
	CVector majorAxisDir;
	for (uint i = 0; i < membrNodes.size(); i++) {
		CVector tmpDir = membrNodes[i] - center;
		CVector tmpUnitDir = tmpDir.getUnitVector();
		double min = 0, max = 0;
		for (uint j = 0; j < membrNodes.size(); j++) {
			CVector tmpDir2 = membrNodes[j] - center;
			double tmpVecProduct = tmpDir2 * tmpUnitDir;
			if (tmpVecProduct < min) {
				min = tmpVecProduct;
			}
			if (tmpVecProduct > max) {
				max = tmpVecProduct;
			}
		}
		double diff = max - min;
		if (diff > maxDiff) {
			maxDiff = diff;
			majorAxisDir = tmpUnitDir;
		}
	}
	lenAlongMajorAxis = maxDiff;
	return majorAxisDir;
}

void SceCells::obtainTwoNewCenters(CVector& oldCenter, CVector& divDir,
		double len_MajorAxis, CVector& centerNew1, CVector& centerNew2) {

	CVector divDirUnit = divDir.getUnitVector();
	double lenChange = len_MajorAxis / 2.0 * centerShiftRatio;
	centerNew1 = oldCenter + lenChange * divDirUnit;
	centerNew2 = oldCenter - lenChange * divDirUnit;
}

void SceCells::prepareTmpVec(uint i, CVector divDir, CVector oldCenter,
		std::vector<VecVal>& tmp1, std::vector<VecVal>& tmp2) {
	tmp1.clear();
	tmp2.clear();
	uint membThreshold = allocPara_m.maxMembrNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint index;
	VecVal tmpData;
	CVector splitDir = divDir.rotateNintyDeg_XY_CC();
	for (uint j = 0; j < maxAllNodePerCell; j++) {
		index = i * maxAllNodePerCell + j;
		if (j < membThreshold) {
			// means node type is membrane
			if (divAuxData.tmpIsActive_M[index] == true) {
				CVector memPos(divAuxData.tmpNodePosX_M[index],
						divAuxData.tmpNodePosY_M[index], 0);
				CVector centerToPosDir = memPos - oldCenter;
				CVector centerToPosUnit = centerToPosDir.getUnitVector();
				CVector crossProduct = Cross(centerToPosDir, splitDir);
				double dotProduct = centerToPosUnit * splitDir;
				tmpData.val = dotProduct;
				tmpData.vec = memPos;
				if (crossProduct.z >= 0) {
					// counter-cloce wise
					tmp1.push_back(tmpData);
				} else {
					// cloce wise
					tmp2.push_back(tmpData);
				}
			}
		} else {
			if (divAuxData.tmpIsActive_M[index] == true) {
				CVector internalPos(divAuxData.tmpNodePosX_M[index],
						divAuxData.tmpNodePosY_M[index], 0);
				CVector centerToPosDir = internalPos - oldCenter;
				CVector shrinkedPos = centerToPosDir * shrinkRatio + oldCenter;
				double dotProduct = centerToPosDir * divDir;
				if (dotProduct > 0) {
					divAuxData.tmp1IntnlVec.push_back(shrinkedPos);
				} else {
					divAuxData.tmp2IntnlVec.push_back(shrinkedPos);
				}
			}
		}
	}
}

void SceCells::calCellArea() {
	thrust::counting_iterator<uint> iBegin(0), iBegin2(0);
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));

	thrust::reduce_by_key(
			make_transform_iterator(iBegin, DivideFunctor(maxAllNodePerCell)),
			make_transform_iterator(iBegin, DivideFunctor(maxAllNodePerCell))
					+ totalNodeCountForActiveCells,
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(
									thrust::make_permutation_iterator(
											cellInfoVecs.activeMembrNodeCounts.begin(),
											make_transform_iterator(iBegin,
													DivideFunctor(
															maxAllNodePerCell))),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell)),
									make_transform_iterator(iBegin,
											ModuloFunctor(maxAllNodePerCell)),
									make_permutation_iterator(
											cellInfoVecs.centerCoordX.begin(),
											make_transform_iterator(iBegin,
													DivideFunctor(
															maxAllNodePerCell))),
									make_permutation_iterator(
											cellInfoVecs.centerCoordY.begin(),
											make_transform_iterator(iBegin,
													DivideFunctor(
															maxAllNodePerCell))))),
					CalTriArea(maxAllNodePerCell, nodeIsActiveAddr,
							nodeLocXAddr, nodeLocYAddr)),
			cellInfoVecs.cellRanksTmpStorage.begin(),
			cellInfoVecs.cellAreaVec.begin(), thrust::equal_to<uint>(),
			thrust::plus<double>());
}

CellsStatsData SceCells::outputPolyCountData() {
       
        cout << " I am at begining of outpolycount"<< std::flush  ; 
	std::cout.flush();
       double sumX,sumY,cntr_X_Domain,cntr_Y_Domain ; 
       int BdryApproach ; 
       BdryApproach=2 ;  
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
        cout << " I am before cells area"<< endl ; 
	calCellArea();
        cout << " I am after cells area" << endl ; 
	CellsStatsData result;

        cout << " I am after result" << endl ; 
	uint bdryCriteria =
			globalConfigVars.getConfigValue("BdryCellCriteria").toInt();
	// already on host; no need to call thrust::copy
	thrust::host_vector<int> adhIndxHost =
			nodes->getInfoVecs().nodeAdhIndxHostCopy;

	thrust::host_vector<double> growthProVecHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.growthProgress.begin(),
			cellInfoVecs.growthProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			growthProVecHost.begin());
	thrust::host_vector<double> growthProMembrVecHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.membrGrowProgress.begin(),
			cellInfoVecs.membrGrowProgress.begin()
					+ allocPara_m.currentActiveCellCount,
			growthProMembrVecHost.begin());
	thrust::host_vector<uint> activeMembrNodeCountHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.activeMembrNodeCounts.begin(),
			cellInfoVecs.activeMembrNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount,
			activeMembrNodeCountHost.begin());
	thrust::host_vector<uint> activeIntnlNodeCountHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.activeIntnlNodeCounts.begin(),
			cellInfoVecs.activeIntnlNodeCounts.begin()
					+ allocPara_m.currentActiveCellCount,
			activeIntnlNodeCountHost.begin());
	thrust::host_vector<double> centerCoordXHost(
			allocPara_m.currentActiveCellCount);
	thrust::host_vector<double> centerCoordYHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.centerCoordX.begin(),
			cellInfoVecs.centerCoordX.begin()
					+ allocPara_m.currentActiveCellCount,
			centerCoordXHost.begin());
	thrust::copy(cellInfoVecs.centerCoordY.begin(),
			cellInfoVecs.centerCoordY.begin()
					+ allocPara_m.currentActiveCellCount,
			centerCoordYHost.begin());

	thrust::host_vector<double> cellAreaHost(
			allocPara_m.currentActiveCellCount);
	thrust::copy(cellInfoVecs.cellAreaVec.begin(),
			cellInfoVecs.cellAreaVec.begin()
					+ allocPara_m.currentActiveCellCount, cellAreaHost.begin());
        sumX=0 ; 
        sumY=0 ; 
	for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
		CellStats cellStatsData;
		cellStatsData.cellGrowthProgress = growthProVecHost[i];
		cellStatsData.cellRank = i;
		bool isBdry = false;
		std::set<int> neighbors;
		int continousNoAdh = 0;
		//std::cout << "printing adhesion indicies ";
		for (uint j = 0; j < activeMembrNodeCountHost[i]; j++) {
			uint index = i * allocPara_m.maxAllNodePerCell + j;
			//std::cout << adhIndxHost[index] << ",";
			if (adhIndxHost[index] != -1) {
				uint adhCellRank = adhIndxHost[index]
						/ allocPara_m.maxAllNodePerCell;
				//std::cout << adhCellRank << " ";
				neighbors.insert(adhCellRank);
				continousNoAdh = 0;
			} else {
				continousNoAdh = continousNoAdh + 1;
				if (continousNoAdh > bdryCriteria) {
					isBdry = true;
				}
			}
			if (j == activeMembrNodeCountHost[i] - 1
					&& adhIndxHost[index] == -1) {
				int k = 0;
				uint indexNew;
				while (k < activeMembrNodeCountHost[i] - 1) {
					indexNew = i * allocPara_m.maxAllNodePerCell + k;
					if (adhIndxHost[indexNew] == -1) {
						continousNoAdh = continousNoAdh + 1;
						if (continousNoAdh > bdryCriteria) {
							isBdry = true;
						}
						k++;
					} else {
						break;
					}
				}
			}

		}
                 
                
		cellStatsData.isBdryCell = isBdry;
		cellStatsData.numNeighbors = neighbors.size();
		cellStatsData.currentActiveMembrNodes = activeMembrNodeCountHost[i];
		cellStatsData.currentActiveIntnlNodes = activeIntnlNodeCountHost[i];
		cellStatsData.neighborVec = neighbors;
		cellStatsData.membrGrowthProgress = growthProMembrVecHost[i];
		cellStatsData.cellCenter = CVector(centerCoordXHost[i],
				centerCoordYHost[i], 0);
		cellStatsData.cellArea = cellAreaHost[i];
		result.cellsStats.push_back(cellStatsData);
                sumX=sumX+cellStatsData.cellCenter.x ; 
                sumY=sumY+cellStatsData.cellCenter.y ;
                
	}
//Ali
        if (BdryApproach==2) {  
          cout << "sumX=" << sumX << endl ; 
          cout << "sumY=" << sumY << endl ; 
          cntr_X_Domain=sumX/result.cellsStats.size() ; 
          cntr_Y_Domain=sumY/result.cellsStats.size() ;  
          cout << "cntr_X=" << cntr_X_Domain << endl ; 
          cout << "cntr_Y=" << cntr_Y_Domain << endl ;

          double R_Max ;
          double Distance ;
          R_Max=0 ;  
	  for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
            Distance=sqrt( pow(centerCoordXHost[i]-cntr_X_Domain,2) +pow(centerCoordYHost[i]-cntr_Y_Domain,2) ) ; 
            if (Distance > R_Max) {
              R_Max=Distance ; 
            }
          }
        
          cout << "R_Max=" << R_Max << endl ;

	  for (uint i = 0; i < allocPara_m.currentActiveCellCount; i++) {
            Distance=sqrt( pow(centerCoordXHost[i]-cntr_X_Domain,2) +pow(centerCoordYHost[i]-cntr_Y_Domain,2) ) ; 
            if (Distance > 0.9* R_Max) {
	      result.cellsStats[i].isBdryCell = true;
              cout << "isBdryCell"<< i<< endl ; 
            }
            else {
	      result.cellsStats[i].isBdryCell = false;
              cout << "isNormalCell"<< i << endl ; 
            }
          }
        }
        //Ali
        cout << "I want to write data" << endl ;  
       // ofstream  Stress_Strain_Single ; 
        //Stress_Strain_Single.open("Stress_Strain_Single.txt"); 
        //Stress_Strain_Single.close() ;
       //Ali
        result.MaxDistanceX=abs(centerCoordXHost[1]-centerCoordXHost[0]); //Ali
        result.Cells_Extrem_Loc[0]=MinX; 
        result.Cells_Extrem_Loc[1]=MaxX; 
        result.Cells_Extrem_Loc[2]=MinY;
        result.Cells_Extrem_Loc[3]=MaxY ;
        result.F_Ext_Out=membrPara.F_Ext_Incline*curTime ; 
        //if (dt==curTime) { 
        //result.Init_Displace=MaxX-MinX ; 
       // }
       //Ali
	return result;
}

__device__ bool bigEnough(double& num) {
	if (num > minDivisor) {
		return true;
	} else {
		return false;
	}
}

__device__ double cross_Z(double vecA_X, double vecA_Y, double vecB_X,
		double vecB_Y) {
	return vecA_X * vecB_Y - vecA_Y * vecB_X;
}

__device__ double calBendMulti(double& angle, uint activeMembrCt) {
	double equAngle = PI - PI / activeMembrCt;
	return bendCoeff * (angle - equAngle);
}

void SceCells::applySceCellDisc_M() {
	totalNodeCountForActiveCells = allocPara_m.currentActiveCellCount
			* allocPara_m.maxAllNodePerCell;
	uint maxAllNodePerCell = allocPara_m.maxAllNodePerCell;
	uint maxMemNodePerCell = allocPara_m.maxMembrNodePerCell;
	thrust::counting_iterator<uint> iBegin(0);

	double* nodeLocXAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocX[0]));
	double* nodeLocYAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeLocY[0]));
	bool* nodeIsActiveAddr = thrust::raw_pointer_cast(
			&(nodes->getInfoVecs().nodeIsActive[0]));

	double grthPrgrCriVal_M = growthAuxData.grthProgrEndCPU
			- growthAuxData.prolifDecay
					* (growthAuxData.grthProgrEndCPU
							- growthAuxData.grthPrgrCriVal_M_Ori);

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeIntnlNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							thrust::make_permutation_iterator(
									cellInfoVecs.activeMembrNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							thrust::make_permutation_iterator(
									cellInfoVecs.activeIntnlNodeCounts.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							make_transform_iterator(iBegin,
									DivideFunctor(maxAllNodePerCell)),
							make_transform_iterator(iBegin,
									ModuloFunctor(maxAllNodePerCell)),
							thrust::make_permutation_iterator(
									cellInfoVecs.growthProgress.begin(),
									make_transform_iterator(iBegin,
											DivideFunctor(maxAllNodePerCell))),
							nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin()))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->getInfoVecs().nodeVelX.begin(),
							nodes->getInfoVecs().nodeVelY.begin())),
			AddSceCellForce(maxAllNodePerCell, maxMemNodePerCell, nodeLocXAddr,
					nodeLocYAddr, nodeIsActiveAddr, grthPrgrCriVal_M));
}

__device__
void calAndAddIB_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& xRes, double& yRes, double grthPrgrCriVal_M) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	if (growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceIBDiv_M[4]) {
			forceValue = -sceIBDiv_M[0] / sceIBDiv_M[2]
					* exp(-linkLength / sceIBDiv_M[2])
					+ sceIBDiv_M[1] / sceIBDiv_M[3]
							* exp(-linkLength / sceIBDiv_M[3]);
		}
	} else if (growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIBDiv_M[4])
				+ (1.0 - percent) * sceIB_M[4];
		if (linkLength < lenLimit) {
			double intnlBPara0 = percent * (sceIBDiv_M[0])
					+ (1.0 - percent) * sceIB_M[0];
			double intnlBPara1 = percent * (sceIBDiv_M[1])
					+ (1.0 - percent) * sceIB_M[1];
			double intnlBPara2 = percent * (sceIBDiv_M[2])
					+ (1.0 - percent) * sceIB_M[2];
			double intnlBPara3 = percent * (sceIBDiv_M[3])
					+ (1.0 - percent) * sceIB_M[3];
			forceValue = -intnlBPara0 / intnlBPara2
					* exp(-linkLength / intnlBPara2)
					+ intnlBPara1 / intnlBPara3
							* exp(-linkLength / intnlBPara3);
		}
	} else {
		if (linkLength < sceIB_M[4]) {
			forceValue = -sceIB_M[0] / sceIB_M[2]
					* exp(-linkLength / sceIB_M[2])
					+ sceIB_M[1] / sceIB_M[3] * exp(-linkLength / sceIB_M[3]);
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}

__device__
void calAndAddII_M(double& xPos, double& yPos, double& xPos2, double& yPos2,
		double& growPro, double& xRes, double& yRes, double grthPrgrCriVal_M) {
	double linkLength = compDist2D(xPos, yPos, xPos2, yPos2);

	double forceValue = 0;
	if (growPro > grthPrgrCriEnd_M) {
		if (linkLength < sceIIDiv_M[4]) {
			forceValue = -sceIIDiv_M[0] / sceIIDiv_M[2]
					* exp(-linkLength / sceIIDiv_M[2])
					+ sceIIDiv_M[1] / sceIIDiv_M[3]
							* exp(-linkLength / sceIIDiv_M[3]);
		}
	} else if (growPro > grthPrgrCriVal_M) {
		double percent = (growPro - grthPrgrCriVal_M)
				/ (grthPrgrCriEnd_M - grthPrgrCriVal_M);
		double lenLimit = percent * (sceIIDiv_M[4])
				+ (1.0 - percent) * sceII_M[4];
		if (linkLength < lenLimit) {
			double intraPara0 = percent * (sceIIDiv_M[0])
					+ (1.0 - percent) * sceII_M[0];
			double intraPara1 = percent * (sceIIDiv_M[1])
					+ (1.0 - percent) * sceII_M[1];
			double intraPara2 = percent * (sceIIDiv_M[2])
					+ (1.0 - percent) * sceII_M[2];
			double intraPara3 = percent * (sceIIDiv_M[3])
					+ (1.0 - percent) * sceII_M[3];
			forceValue = -intraPara0 / intraPara2
					* exp(-linkLength / intraPara2)
					+ intraPara1 / intraPara3 * exp(-linkLength / intraPara3);
		}
	} else {
		if (linkLength < sceII_M[4]) {
			forceValue = -sceII_M[0] / sceII_M[2]
					* exp(-linkLength / sceII_M[2])
					+ sceII_M[1] / sceII_M[3] * exp(-linkLength / sceII_M[3]);
		}
	}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
}
