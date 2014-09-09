#include "SceCells.h"

__constant__ uint GridDimension[2];
__constant__ double gridSpacing;

double epsilon = 1.0e-12;

void SceCells_M::distributeBdryIsActiveInfo() {
	thrust::fill(nodes->nodeIsActive.begin(),
			nodes->nodeIsActive.begin() + beginPosOfEpiNode, true);
}

void SceCells_M::distributeProfileIsActiveInfo() {
	thrust::fill(nodes->nodeIsActive.begin() + beginPosOfEpiNode,
			nodes->nodeIsActive.begin() + beginPosOfEpiNode
					+ nodes->currentActiveProfileNodeCount, true);
}

void SceCells_M::distributeECMIsActiveInfo() {
	uint totalNodeCountForActiveECM = currentActiveECMCount * maxNodeOfECM;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveECM);
	thrust::fill(nodes->nodeIsActive.begin() + beginPosOfECMNode,
			nodes->nodeIsActive.begin() + totalNodeCountForActiveECM
					+ beginPosOfECMNode, true);
}

void SceCells_M::distributeCellIsActiveInfo() {
	//totalNodeCountForActiveCells = currentActiveCellCount * maxNodeOfOneCell;
	//std::cout << "currentActive cell count" << currentActiveCellCount
	//		<< std::endl;
	//thrust::host_vector<uint> activeCountsFromGPU = activeNodeCountOfThisCell;
	//int count = 0;
	//for (int i = 0; i < activeCountsFromGPU.size(); i++) {
	//	if (activeCountsFromGPU[i] != 0) {
	//		count++;
	//	}
	//}
	//std::cout << "count of non zero init counts:" << count << std::endl;
	//int jj;
	//std::cin >> jj;
	totalNodeCountForActiveCells = currentActiveCellCount
				* maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);
	std::cout << "started distri cell active" << std::endl;
	std::cout << "active node count = " << totalNodeCountForActiveCells
			<< std::endl;
	thrust::transform(
			thrust::make_transform_iterator(countingBegin,
					ModuloFunctor(maxNodeOfOneCell)),
			thrust::make_transform_iterator(countingEnd,
					ModuloFunctor(maxNodeOfOneCell)),
			thrust::make_permutation_iterator(activeNodeCountOfThisCell.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(maxNodeOfOneCell))),
			nodes->nodeIsActive.begin() + beginPosOfCellsNode,
			thrust::less<uint>());
	/*
	 thrust::host_vector<bool> hostIsActive = nodes->nodeIsActive;
	 count = 0;
	 for (int i = 0; i < (nodes->startPosCells + totalNodeCountForActiveCells);
	 i++) {
	 if (hostIsActive[i] == true) {
	 count++;
	 }
	 }
	 */
	//std::cout << "total number of active nodes:" << count << std::endl;
	//std::cin >> jj;
}

/**
 * constructor for SceCells_M.
 * takes SceNodes, which is a pre-allocated multi-array, as input argument.
 * This might be strange from a design perspective but has a better performance
 * while running on parallel.
 */
SceCells_M::SceCells_M(SceNodes* nodesInput) :
		maxNodeOfOneCell(nodesInput->getMaxNodeOfOneCell()), countingBegin(0), initCellCount(
				maxNodeOfOneCell / 2), initGrowthProgress(0.0) {

	addNodeDistance =
			globalConfigVars.getConfigValue("DistanceForAddingNode").toDouble();
	minDistanceToOtherNode = globalConfigVars.getConfigValue(
			"MinDistanceToOtherNode").toDouble();
	cellInitLength =
			globalConfigVars.getConfigValue("CellInitLength").toDouble();
	cellFinalLength =
			globalConfigVars.getConfigValue("CellFinalLength").toDouble();
	elongationCoefficient = globalConfigVars.getConfigValue(
			"ElongateCoefficient").toDouble();
	chemoCoefficient =
			globalConfigVars.getConfigValue("ChemoCoefficient").toDouble();
	isDivideCriticalRatio = globalConfigVars.getConfigValue(
			"IsDivideCrticalRatio").toDouble();

	beginPosOfBdry = 0;
//maxNodeOfBdry = nodesInput->cellSpaceForBdry;

	beginPosOfEpi = 1;
	maxNodeOfEpi = nodesInput->maxProfileNodeCount; // represents maximum number of nodes of epithilum layer.
	beginPosOfEpiNode = nodesInput->startPosProfile;

	maxNodeOfECM = nodesInput->maxNodePerECM; // represents maximum number of nodes per ECM
// represents begining position of ECM (in node perspective) value is 2
	beginPosOfECM = 2;
// represents begining position of ECM (in cell perspective)
	beginPosOfECMNode = nodesInput->startPosECM;
	maxECMCount = nodesInput->maxECMCount;  // represents maximum number of ECM.

	beginPosOfCells = beginPosOfECM + maxECMCount; // represents begining position of cells (in cell perspective)
	beginPosOfCellsNode = nodesInput->startPosCells; // represents begining position of cells (in node perspective)
// maxCellCount only take FNM and MX cells into consideration.
	maxNodeOfOneCell = nodesInput->getMaxNodeOfOneCell();
	maxCellCount = nodesInput->getMaxCellCount();
// maxCellCountAll counts all cells include pseduo cells
	maxCellCountAll = nodesInput->getMaxCellCountAll();
	maxTotalNodeCountCellOnly = nodesInput->getMaxTotalCellNodeCount();
	currentActiveCellCount = nodesInput->getCurrentActiveCellCount();

// cellSpaceForBdry = nodesInput->getCellSpaceForBdry();

	nodes = nodesInput;
	growthProgress.resize(maxCellCount, 0.0);
	expectedLength.resize(maxCellCount, cellInitLength);
	lengthDifference.resize(maxCellCount, 0.0);
	smallestDistance.resize(maxCellCount);
	biggestDistance.resize(maxCellCount);
	activeNodeCountOfThisCell.resize(maxCellCount);
	lastCheckPoint.resize(maxCellCount, 0.0);
	isDivided.resize(maxCellCount);
	cellTypes.resize(maxCellCount, MX);
	isScheduledToGrow.resize(maxCellCount, false);
	centerCoordX.resize(maxCellCount);
	centerCoordY.resize(maxCellCount);
	centerCoordZ.resize(maxCellCount);
	cellRanksTmpStorage.resize(maxCellCount);
	growthSpeed.resize(maxCellCount, 0.0);
	growthXDir.resize(maxCellCount);
	growthYDir.resize(maxCellCount);

	//xCoordTmp.resize(maxTotalCellNodeCount);
	//yCoordTmp.resize(maxTotalCellNodeCount);
	//zCoordTmp.resize(maxTotalCellNodeCount);
	cellRanks.resize(maxTotalNodeCountCellOnly);
	activeXPoss.resize(maxTotalNodeCountCellOnly);
	activeYPoss.resize(maxTotalNodeCountCellOnly);
	activeZPoss.resize(maxTotalNodeCountCellOnly);
	distToCenterAlongGrowDir.resize(maxTotalNodeCountCellOnly);

// reason for adding a small term here is to avoid scenario when checkpoint might add many times
// up to 0.99999999 which is theoretically 1.0 but not in computer memory. If we don't include
// this small term we might risk adding one more node.
	growThreshold = 1.0 / (maxNodeOfOneCell - maxNodeOfOneCell / 2) + epsilon;

	nodeIsActiveAddress = thrust::raw_pointer_cast(
			&(nodes->nodeIsActive[beginPosOfCellsNode]));
	nodeXPosAddress = thrust::raw_pointer_cast(
			&(nodes->nodeLocX[beginPosOfCellsNode]));
	nodeYPosAddress = thrust::raw_pointer_cast(
			&(nodes->nodeLocY[beginPosOfCellsNode]));


	distributeIsCellRank();
}

/**
 * This is a method for cell growth. growth is influened by two chemical fields.
 * first step: assign the growth magnitude and direction info that was calculated outside
 *     to internal values
 *     please note that a cell should not grow if its type is boundary.
 *
 * second step: use the growth magnitude and dt to update growthProgress
 *
 * third step: use lastCheckPoint and growthProgress to decide whether add point or not
 *
 * fourth step: use growthProgress and growthXDir&growthYDir to compute
 *     expected length along the growth direction.
 *
 * fifth step:  reducing the smallest value and biggest value
 *     a cell's node to its center point
 *
 * sixth step: compute the current length and then
 *     compute its difference with expected length
 *
 * seventh step: use the difference that just computed and growthXDir&growthYDir
 *     to apply stretching force (velocity) on nodes of all cells
 *
 * eighth step: cell move according to the velocity computed
 *
 * ninth step: also add a point if scheduled to grow.
 *     This step does not guarantee success ; If adding new point failed, it will not change
 *     isScheduleToGrow and activeNodeCount;
 */
void SceCells_M::grow2DTwoRegions(double d_t, GrowthDistriMap &region1,
		GrowthDistriMap &region2) {

	dt = d_t;

// obtain pointer address for first region
	growthFactorMagAddress = thrust::raw_pointer_cast(
			&(region1.growthFactorMag[0]));
	growthFactorDirXAddress = thrust::raw_pointer_cast(
			&(region1.growthFactorDirXComp[0]));
	growthFactorDirYAddress = thrust::raw_pointer_cast(
			&(region1.growthFactorDirYComp[0]));

// obtain pointer address for second region
	growthFactorMagAddress2 = thrust::raw_pointer_cast(
			&(region2.growthFactorMag[0]));
	growthFactorDirXAddress2 = thrust::raw_pointer_cast(
			&(region2.growthFactorDirXComp[0]));
	growthFactorDirYAddress2 = thrust::raw_pointer_cast(
			&(region2.growthFactorDirYComp[0]));

	totalNodeCountForActiveCells = currentActiveCellCount * maxNodeOfOneCell;
	std::cout << "totalNodeCount = " << totalNodeCountForActiveCells
			<< std::endl;
	std::cout << "before all functions start" << std::endl;
	copyGrowInfoFromGridToCells(region1, region2);
	std::cout << "after copy grow info" << std::endl;
	updateGrowthProgress();
	std::cout << "after update growth progress" << std::endl;
	decideIsScheduleToGrow();
	std::cout << "after decode os schedule to grow" << std::endl;
	computeCellTargetLength();
	std::cout << "after compute cell target length" << std::endl;
	computeDistToCellCenter();
	std::cout << "after compute dist to center" << std::endl;
	findMinAndMaxDistToCenter();
	std::cout << "after find min and max dist" << std::endl;
	computeLenDiffExpCur();
	std::cout << "after compute diff " << std::endl;
	stretchCellGivenLenDiff();
	std::cout << "after apply stretch force" << std::endl;
	cellChemotaxis();
	std::cout << "after apply cell chemotaxis" << std::endl;
	addPointIfScheduledToGrow();
	std::cout << "after adding node" << std::endl;
}

/**
 * we need to copy the growth information from grid for chemical to cell nodes.
 *
 * checked.
 */
void SceCells_M::copyGrowInfoFromGridToCells(GrowthDistriMap &region1,
		GrowthDistriMap &region2) {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(centerCoordX.begin(),
							centerCoordY.begin(), cellTypes.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(centerCoordX.begin(),
							centerCoordY.begin(), cellTypes.begin()))
					+ currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(growthSpeed.begin(), growthXDir.begin(),
							growthYDir.begin())),
			LoadChemDataToNode(region1.gridDimensionX, region1.gridDimensionY,
					region1.gridSpacing, growthFactorMagAddress,
					growthFactorDirXAddress, growthFactorDirYAddress,
					region2.gridDimensionX, region2.gridDimensionY,
					region2.gridSpacing, growthFactorMagAddress2,
					growthFactorDirXAddress2, growthFactorDirYAddress2));
}

/**
 * Use the growth magnitude and dt to update growthProgress.
 */
void SceCells_M::updateGrowthProgress() {
	thrust::transform(growthSpeed.begin(),
			growthSpeed.begin() + currentActiveCellCount,
			growthProgress.begin(), growthProgress.begin(),
			SaxpyFunctorWithMaxOfOne(dt));
}

/**
 * Decide if the cells are going to add a node or not.
 * Use lastCheckPoint and growthProgress to decide whether add point or not
 */
void SceCells_M::decideIsScheduleToGrow() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(growthProgress.begin(),
							lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(growthProgress.begin(),
							lastCheckPoint.begin())) + currentActiveCellCount,
			isScheduledToGrow.begin(), PtCondiOp(growThreshold));
}

/**
 * Calculate target length of cell given the cell growth progress.
 * length is along the growth direction.
 */
void SceCells_M::computeCellTargetLength() {
	thrust::transform(growthProgress.begin(),
			growthProgress.begin() + currentActiveCellCount,
			expectedLength.begin(),
			CompuTarLen(cellInitLength, cellFinalLength));
}

/**
 * Compute distance of each node to its corresponding cell center.
 * The distantce could be either positive or negative, depending on the pre-defined
 * growth direction.
 */
void SceCells_M::computeDistToCellCenter() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(centerCoordX.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(centerCoordY.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeIsActive.begin() + beginPosOfCellsNode)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(centerCoordX.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(centerCoordY.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeIsActive.begin() + beginPosOfCellsNode))
					+ totalNodeCountForActiveCells,
			distToCenterAlongGrowDir.begin(), CompuDist());
}

/**
 * For nodes of each cell, find the maximum and minimum distance to the center.
 * We will then calculate the current length of a cell along its growth direction
 * using max and min distance to the center.
 */
void SceCells_M::findMinAndMaxDistToCenter() {
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(maxNodeOfOneCell))
					+ totalNodeCountForActiveCells,
			distToCenterAlongGrowDir.begin(), cellRanksTmpStorage.begin(),
			smallestDistance.begin(), thrust::equal_to<uint>(),
			thrust::minimum<double>());
// for nodes of each cell, find the maximum distance from the node to the corresponding
// cell center along the pre-defined growth direction.
	thrust::reduce_by_key(
			make_transform_iterator(countingBegin,
					DivideFunctor(maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					DivideFunctor(maxNodeOfOneCell))
					+ totalNodeCountForActiveCells,
			distToCenterAlongGrowDir.begin(), cellRanksTmpStorage.begin(),
			biggestDistance.begin(), thrust::equal_to<uint>(),
			thrust::maximum<double>());
}

/**
 * Compute the difference for cells between their expected length and current length.
 */
void SceCells_M::computeLenDiffExpCur() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(expectedLength.begin(),
							smallestDistance.begin(), biggestDistance.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(expectedLength.begin(),
							smallestDistance.begin(), biggestDistance.begin()))
					+ currentActiveCellCount, lengthDifference.begin(),
			CompuDiff());
}

/**
 * Use the difference that just computed and growthXDir&growthYDir
 * to apply stretching force (velocity) on nodes of all cells
 */
void SceCells_M::stretchCellGivenLenDiff() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(lengthDifference.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							nodes->nodeVelX.begin() + beginPosOfCellsNode,
							nodes->nodeVelY.begin() + beginPosOfCellsNode)),
			thrust::make_zip_iterator(
					thrust::make_tuple(distToCenterAlongGrowDir.begin(),
							make_permutation_iterator(lengthDifference.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							nodes->nodeVelX.begin() + beginPosOfCellsNode,
							nodes->nodeVelY.begin() + beginPosOfCellsNode))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->nodeVelX.begin() + beginPosOfCellsNode,
							nodes->nodeVelY.begin() + beginPosOfCellsNode)),
			ApplyStretchForce(elongationCoefficient));
}
/**
 * This is just an attempt. Cells move according to chemicals.
 */
void SceCells_M::cellChemotaxis() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(growthSpeed.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							nodes->nodeVelX.begin() + beginPosOfCellsNode,
							nodes->nodeVelY.begin() + beginPosOfCellsNode)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_permutation_iterator(growthSpeed.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthXDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							make_permutation_iterator(growthYDir.begin(),
									make_transform_iterator(countingBegin,
											DivideFunctor(maxNodeOfOneCell))),
							nodes->nodeVelX.begin() + beginPosOfCellsNode,
							nodes->nodeVelY.begin() + beginPosOfCellsNode))
					+ totalNodeCountForActiveCells,
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->nodeVelX.begin() + beginPosOfCellsNode,
							nodes->nodeVelY.begin() + beginPosOfCellsNode)),
			ApplyChemoVel(chemoCoefficient));
}
/**
 * Adjust the velocities of nodes.
 * For example, velocity of boundary nodes must be zero.
 */
void SceCells_M::adjustNodeVel() {
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(
			totalNodeCountForActiveCells + beginPosOfCellsNode);

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeVelX.begin(),
							nodes->nodeVelY.begin(),
							nodes->nodeIsActive.begin(),
							nodes->nodeCellType.begin(), countingIterBegin)),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeVelX.begin(),
							nodes->nodeVelY.begin(),
							nodes->nodeIsActive.begin(),
							nodes->nodeCellType.begin(), countingIterBegin))
					+ totalNodeCountForActiveCells + beginPosOfCellsNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeVelX.begin(),
							nodes->nodeVelY.begin())),
			VelocityModifier(beginPosOfEpiNode, currectActiveEpiNode));
}
/**
 * Move nodes according to the velocity we just adjusted.
 */
void SceCells_M::moveNodes() {
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeVelX.begin(),
							nodes->nodeVelY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeVelX.begin(),
							nodes->nodeVelY.begin()))
					+ totalNodeCountForActiveCells + beginPosOfCellsNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeLocX.begin(),
							nodes->nodeLocY.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes->nodeLocX.begin(),
							nodes->nodeLocY.begin())), SaxpyFunctorDim2(dt));
}
/**
 * Add a point to a cell if it is scheduled to grow.
 * This step does not guarantee success ; If adding new point failed, it will not change
 * isScheduleToGrow and activeNodeCount;
 */
void SceCells_M::addPointIfScheduledToGrow() {

	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(isScheduledToGrow.begin(),
							activeNodeCountOfThisCell.begin(),
							centerCoordX.begin(), centerCoordY.begin(),
							countingBegin, lastCheckPoint.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(isScheduledToGrow.begin(),
							activeNodeCountOfThisCell.begin(),
							centerCoordX.begin(), centerCoordY.begin(),
							countingBegin, lastCheckPoint.begin()))
					+ currentActiveCellCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(isScheduledToGrow.begin(),
							activeNodeCountOfThisCell.begin(),
							lastCheckPoint.begin())),
			AddPtOp(maxNodeOfOneCell, addNodeDistance, minDistanceToOtherNode,
					nodeIsActiveAddress, nodeXPosAddress, nodeYPosAddress,
					time(NULL), growThreshold));
}

/**
 * To run all the cell level logics.
 * First step we got center positions of cells.
 * Grow.
 */
void SceCells_M::runAllCellLevelLogics(double dt, GrowthDistriMap &region1,
		GrowthDistriMap &region2) {
	std::cerr << "enter run all cell level logics" << std::endl;
	computeCenterPos();
	std::cerr << "after compute center position." << std::endl;
	grow2DTwoRegions(dt, region1, region2);
	std::cerr << "after grow cells" << std::endl;
	distributeIsActiveInfo();
	std::cerr << "after distribute is active info." << std::endl;
	divide2DSimplified();
	std::cerr << "after divide 2D simplified." << std::endl;
	distributeIsActiveInfo();
	std::cerr << "after distribute is active info." << std::endl;
	allComponentsMove();
	std::cerr << "after all components move." << std::endl;
}

void SceCells_M::allComponentsMove() {
	adjustNodeVel();
	moveNodes();
}

/**
 * Mark cell node as either activdistributeIsActiveInfo()e or inactive.
 * left part of the node array will be active and right part will be inactive.
 * the threshold is defined by array @activeNodeCountOfThisCell.
 * e.g. activeNodeCountOfThisCell = {2,3} and  maxNodeOfOneCell = 5,
 *
 * @Checked.
 */
void SceCells_M::distributeIsActiveInfo() {
	std::cout << "before distribute bdry isActive" << std::endl;
	distributeBdryIsActiveInfo();
	std::cout << "before distribute profile isActive" << std::endl;
	distributeProfileIsActiveInfo();
	std::cout << "before distribute ecm isActive" << std::endl;
	distributeECMIsActiveInfo();
	std::cout << "before distribute cells isActive" << std::endl;
	distributeCellIsActiveInfo();
}

void SceCells_M::distributeIsCellRank() {
	uint totalNodeCountForActiveCells = currentActiveCellCount
			* maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingCellEnd(
			totalNodeCountForActiveCells);
	std::cerr << "totalNodeCount for active cells "
			<< totalNodeCountForActiveCells << std::endl;

//thrust::counting_iterator<uint> countingECMEnd(countingECMEnd);

// only computes the cell ranks of cells. the rest remain unchanged.
	thrust::transform(countingBegin, countingCellEnd,
			nodes->nodeCellRank.begin() + beginPosOfCellsNode,
			DivideFunctor(maxNodeOfOneCell));
	std::cerr << "finished cellRank transformation" << std::endl;
}

/**
 * This method computes center of all cells.
 * more efficient then simply iterating the cell because of parallel reducing.
 *
 * @Checked.
 */
void SceCells_M::computeCenterPos() {
	uint totalNodeCountForActiveCells = currentActiveCellCount
			* maxNodeOfOneCell;
	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd(totalNodeCountForActiveCells);
	uint totalNumberOfActiveNodes = thrust::reduce(
			activeNodeCountOfThisCell.begin(),
			activeNodeCountOfThisCell.begin() + currentActiveCellCount);

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(maxNodeOfOneCell)),
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeLocZ.begin() + beginPosOfCellsNode)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(maxNodeOfOneCell)),
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeLocZ.begin() + beginPosOfCellsNode))
					+ totalNodeCountForActiveCells,
			nodes->nodeIsActive.begin() + beginPosOfCellsNode,
			thrust::make_zip_iterator(
					thrust::make_tuple(cellRanks.begin(), activeXPoss.begin(),
							activeYPoss.begin(), activeZPoss.begin())),
			isTrue());

	thrust::reduce_by_key(cellRanks.begin(),
			cellRanks.begin() + totalNumberOfActiveNodes,
			thrust::make_zip_iterator(
					thrust::make_tuple(activeXPoss.begin(), activeYPoss.begin(),
							activeZPoss.begin())), cellRanksTmpStorage.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(centerCoordX.begin(),
							centerCoordY.begin(), centerCoordZ.begin())),
			thrust::equal_to<uint>(), CVec3Add());
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(centerCoordX.begin(),
							centerCoordY.begin(), centerCoordZ.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(centerCoordX.begin(),
							centerCoordY.begin(), centerCoordZ.begin()))
					+ currentActiveCellCount, activeNodeCountOfThisCell.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(centerCoordX.begin(),
							centerCoordY.begin(), centerCoordZ.begin())),
			CVec3Divide());
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

//TODO: also pay attention to number of active nodes per cell. This seems to be omitted.
void SceCells_M::divide2DSimplified() {
	decideIfGoingToDivide();
	copyCellsPreDivision();
	sortNodesAccordingToDist();
	copyLeftAndRightToSeperateArrays();
	transformIsActiveArrayOfBothArrays();
	addSecondArrayToCellArray();
	copyFirstArrayToPreviousPos();
	updateActiveCellCount();
	markIsDivideFalse();
}

void SceCells_M::decideIfGoingToDivide() {
// step 1
	thrust::transform(
			thrust::make_zip_iterator(
					thrust::make_tuple(lengthDifference.begin(),
							expectedLength.begin(), growthProgress.begin(),
							activeNodeCountOfThisCell.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(lengthDifference.begin(),
							expectedLength.begin(), growthProgress.begin(),
							activeNodeCountOfThisCell.begin()))
					+ currentActiveCellCount, isDivided.begin(),
			CompuIsDivide(isDivideCriticalRatio, maxNodeOfOneCell));
}

void SceCells_M::copyCellsPreDivision() {
// step 2 : copy all cell rank and distance to its corresponding center with divide flag = 1
	totalNodeCountForActiveCells = currentActiveCellCount * maxNodeOfOneCell;
// sum all bool values which indicate whether the cell is going to divide.
// toBeDivideCount is the total number of cells going to divide.
	toBeDivideCount = thrust::reduce(isDivided.begin(),
			isDivided.begin() + currentActiveCellCount, (uint) (0));
	nodeStorageCount = toBeDivideCount * maxNodeOfOneCell;
	tmpIsActiveHold1 = thrust::device_vector<bool>(nodeStorageCount, true);
	tmpDistToCenter1 = thrust::device_vector<double>(nodeStorageCount, 0.0);
	tmpCellRankHold1 = thrust::device_vector<uint>(nodeStorageCount, 0.0);
	tmpXValueHold1 = thrust::device_vector<double>(nodeStorageCount, 0.0);
	tmpYValueHold1 = thrust::device_vector<double>(nodeStorageCount, 0.0);
	tmpZValueHold1 = thrust::device_vector<double>(nodeStorageCount, 0.0);

	tmpCellTypes = thrust::device_vector<SceNodeType>(nodeStorageCount);

	tmpIsActiveHold2 = thrust::device_vector<bool>(nodeStorageCount, false);
	tmpDistToCenter2 = thrust::device_vector<double>(nodeStorageCount, 0.0);
	tmpXValueHold2 = thrust::device_vector<double>(nodeStorageCount, 0.0);
	tmpYValueHold2 = thrust::device_vector<double>(nodeStorageCount, 0.0);
	tmpZValueHold2 = thrust::device_vector<double>(nodeStorageCount, 0.0);

// step 2 , continued
	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(maxNodeOfOneCell)),
							distToCenterAlongGrowDir.begin(),
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeLocZ.begin() + beginPosOfCellsNode,
							nodes->nodeCellType.begin() + beginPosOfCellsNode)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							make_transform_iterator(countingBegin,
									DivideFunctor(maxNodeOfOneCell)),
							distToCenterAlongGrowDir.begin(),
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeLocZ.begin() + beginPosOfCellsNode,
							nodes->nodeCellType.begin() + beginPosOfCellsNode))
					+ totalNodeCountForActiveCells,
			thrust::make_permutation_iterator(isDivided.begin(),
					make_transform_iterator(countingBegin,
							DivideFunctor(maxNodeOfOneCell))),
			thrust::make_zip_iterator(
					thrust::make_tuple(tmpCellRankHold1.begin(),
							tmpDistToCenter1.begin(), tmpXValueHold1.begin(),
							tmpYValueHold1.begin(), tmpZValueHold1.begin(),
							tmpCellTypes.begin())), isTrue());
}

/**
 * performance wise, this implementation is not the best because I can use only one sort_by_key
 * with speciialized comparision operator. However, This implementation is more robust and won't
 * compromise performance too much.
 */
void SceCells_M::sortNodesAccordingToDist() {
//step 3
	for (uint i = 0; i < toBeDivideCount; i++) {
		thrust::sort_by_key(tmpDistToCenter1.begin() + i * maxNodeOfOneCell,
				tmpDistToCenter1.begin() + (i + 1) * maxNodeOfOneCell,
				thrust::make_zip_iterator(
						thrust::make_tuple(
								tmpXValueHold1.begin() + i * maxNodeOfOneCell,
								tmpYValueHold1.begin() + i * maxNodeOfOneCell,
								tmpZValueHold1.begin()
										+ i * maxNodeOfOneCell)));
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
void SceCells_M::copyLeftAndRightToSeperateArrays() {
//step 4.
	thrust::scatter_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(tmpXValueHold1.begin(),
							tmpYValueHold1.begin(), tmpZValueHold1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(tmpXValueHold1.end(),
							tmpYValueHold1.end(), tmpZValueHold1.end())),
			make_transform_iterator(countingBegin,
					LeftShiftFunctor(maxNodeOfOneCell)),
			make_transform_iterator(countingBegin,
					IsRightSide(maxNodeOfOneCell)),
			thrust::make_zip_iterator(
					thrust::make_tuple(tmpXValueHold2.begin(),
							tmpYValueHold2.begin(), tmpZValueHold2.begin())));
}

void SceCells_M::transformIsActiveArrayOfBothArrays() {
	thrust::transform(countingBegin, countingBegin + nodeStorageCount,
			tmpIsActiveHold1.begin(), IsLeftSide(maxNodeOfOneCell));
	thrust::transform(countingBegin, countingBegin + nodeStorageCount,
			tmpIsActiveHold2.begin(), IsLeftSide(maxNodeOfOneCell));
	if (toBeDivideCount != 0) {
		std::cout << "before insert, active cell count in nodes:"
				<< nodes->getCurrentActiveCellCount() << std::endl;
	}
}

void SceCells_M::addSecondArrayToCellArray() {
/// step 6. call SceNodes function to add newly divided cells
	nodes->addNewlyDividedCells(tmpXValueHold2, tmpYValueHold2, tmpZValueHold2,
			tmpIsActiveHold2, tmpCellTypes);
}

void SceCells_M::copyFirstArrayToPreviousPos() {
	thrust::scatter(
			thrust::make_zip_iterator(
					thrust::make_tuple(tmpIsActiveHold1.begin(),
							tmpXValueHold1.begin(), tmpYValueHold1.begin(),
							tmpZValueHold1.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(tmpIsActiveHold1.end(),
							tmpXValueHold1.end(), tmpYValueHold1.end(),
							tmpZValueHold1.end())),
			thrust::make_transform_iterator(
					thrust::make_zip_iterator(
							thrust::make_tuple(countingBegin,
									tmpCellRankHold1.begin())),
					CompuPos(maxNodeOfOneCell)),
			thrust::make_zip_iterator(
					thrust::make_tuple(
							nodes->nodeIsActive.begin() + beginPosOfCellsNode,
							nodes->nodeLocX.begin() + beginPosOfCellsNode,
							nodes->nodeLocY.begin() + beginPosOfCellsNode,
							nodes->nodeLocZ.begin() + beginPosOfCellsNode)));

	thrust::scatter(initCellCount, initCellCount + toBeDivideCount,
			tmpCellRankHold1.begin(), activeNodeCountOfThisCell.begin());

	thrust::scatter(initGrowthProgress, initGrowthProgress + toBeDivideCount,
			tmpCellRankHold1.begin(), growthProgress.begin());
	thrust::scatter(initGrowthProgress, initGrowthProgress + toBeDivideCount,
			tmpCellRankHold1.begin(), lastCheckPoint.begin());

	thrust::fill(activeNodeCountOfThisCell.begin() + currentActiveCellCount,
			activeNodeCountOfThisCell.begin() + currentActiveCellCount
					+ toBeDivideCount, maxNodeOfOneCell / 2);

}

void SceCells_M::updateActiveCellCount() {
	currentActiveCellCount = currentActiveCellCount + toBeDivideCount;
	nodes->setCurrentActiveCellCount(currentActiveCellCount);
}

void SceCells_M::markIsDivideFalse() {
	thrust::fill(isDivided.begin(), isDivided.begin() + currentActiveCellCount,
			false);
}
