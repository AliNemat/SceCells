/**
 * @file SimulationDomainGPU.cpp
 * @brief this file contains domain level logic.
 * @author Wenzhao Sun wsun2@nd.edu
 * @bug no know bugs
 */

#include "SimulationDomainGPU.h"

using namespace std;

/**
 * Constructor.
 * Focus on initializing domain-level parameter values.
 */
SimulationDomainGPU::SimulationDomainGPU() {
	cout << "before allocation memory" << endl;
	//thrust::host_vector<int> aa;
	//aa.resize(50000);
	//thrust::device_vector<int> bb = aa;
	//thrust::device_vector<double> cc(5000);
	cout << "after allocate memory" << endl;

	cout << "start to create simulatonDomainGPU object" << endl;

	maxCellInDomain =
			globalConfigVars.getConfigValue("MaxCellInDomain").toInt();
	maxNodePerCell =
			globalConfigVars.getConfigValue("MaxNodePerCell").toDouble();
	maxECMInDomain =
			globalConfigVars.getConfigValue("MaxECMInDomain").toDouble();
	maxNodePerECM = globalConfigVars.getConfigValue("MaxNodePerECM").toDouble();
	//initECMCount = globalConfigVars.getConfigValue("InitECMCount").toInt();
	FinalToInitProfileNodeCountRatio = globalConfigVars.getConfigValue(
			"FinalToInitProfileNodeCountRatio").toDouble();

	minX = globalConfigVars.getConfigValue("DOMAIN_XMIN").toDouble();
	maxX = globalConfigVars.getConfigValue("DOMAIN_XMAX").toDouble();
	minY = globalConfigVars.getConfigValue("DOMAIN_YMIN").toDouble();
	maxY = globalConfigVars.getConfigValue("DOMAIN_YMAX").toDouble();

	gridSpacing =
			globalConfigVars.getConfigValue("Cell_Center_Interval").toDouble();

	growthGridXDim = globalConfigVars.getConfigValue("GrowthGridXDim").toInt();
	growthGridYDim = globalConfigVars.getConfigValue("GrowthGridYDim").toInt();
	growthGridSpacing =
			globalConfigVars.getConfigValue("GrowthGridSpacing").toDouble();
	growthGridLowerLeftPtX = globalConfigVars.getConfigValue(
			"GrowthGridLowerLeftPtX").toDouble();
	growthGridLowerLeftPtY = globalConfigVars.getConfigValue(
			"GrowthGridLowerLeftPtY").toDouble();

	growthMorCenterXCoord = globalConfigVars.getConfigValue(
			"GrowthMorCenterXCoord").toDouble();
	growthMorCenterYCoord = globalConfigVars.getConfigValue(
			"GrowthMorCenterYCoord").toDouble();
	growthMorHighConcen =
			globalConfigVars.getConfigValue("GrowthMorHighConcen").toDouble();
	growthMorLowConcen =
			globalConfigVars.getConfigValue("GrowthMorLowConcen").toDouble();
	growthMorDiffSlope =
			globalConfigVars.getConfigValue("GrowthMorDiffSlope").toDouble();

	growthMorCenterXCoordMX = globalConfigVars.getConfigValue(
			"GrowthMorCenterXCoordMX").toDouble();
	growthMorCenterYCoordMX = globalConfigVars.getConfigValue(
			"GrowthMorCenterYCoordMX").toDouble();
	growthMorHighConcenMX = globalConfigVars.getConfigValue(
			"GrowthMorHighConcenMX").toDouble();
	growthMorLowConcenMX = globalConfigVars.getConfigValue(
			"GrowthMorLowConcenMX").toDouble();
	growthMorDiffSlopeMX = globalConfigVars.getConfigValue(
			"GrowthMorDiffSlopeMX").toDouble();

	intraLinkDisplayRange = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();

	cout << "after reading values from config file" << endl;
	cout << "key parameters are : maxCellInDomain = " << maxCellInDomain
			<< "maxNodePerCell = " << maxNodePerCell << "maxECMInDomain = "
			<< maxECMInDomain << "maxNodePerECM = " << maxNodePerECM << endl;
	cout << "after created nodes object" << endl;
	//cells = SceCells(&nodes);

	cout << "after created nodes and cells object" << endl;
	growthMap = GrowthDistriMap(growthGridXDim, growthGridYDim,
			growthGridSpacing);
	growthMap.initialize(growthGridLowerLeftPtX, growthGridLowerLeftPtY,
			growthMorCenterXCoord, growthMorCenterYCoord, growthMorHighConcen,
			growthMorLowConcen, growthMorDiffSlope);

	cout << "after created growthMap1" << endl;
	growthMap2 = GrowthDistriMap(growthGridXDim, growthGridYDim,
			growthGridSpacing);
	growthMap2.initialize(growthGridLowerLeftPtX, growthGridLowerLeftPtY,
			growthMorCenterXCoordMX, growthMorCenterYCoordMX,
			growthMorHighConcenMX, growthMorLowConcenMX, growthMorDiffSlopeMX);
	cout << "after created growthMap2" << endl;
}

/**
 * Initialize five different cell nodes.
 * We have to initialize five types of cells:
 * first is Boundary (B), fixed nodes on the boundary;
 * second is Profile (P), Epithilum cells;
 * third is ECM (E), extra-cellular matrix;
 * fourth is FNM (F), front nasal mass;
 * fifth is MX (M) maxillary cells.
 *
 * like this:
 * B-B-B-B-B-B-B-B-B-P-P-P-P-E-E-E-E-F-F-F-F-F-F-F-F-M-M-M-M-M-M-M-M
 * B, P and E is fixed. F and M will grow.
 * Rules:
 * 1a, Number of boundary nodes is fixed.
 * 1b, Profile nodes may or may not increase. Still testing the model.
 *     however the space is always reserved for it, though some spaces may not be active.
 * 1c, Extra-cellular matrix will grow.
 *     however the space is always reserved for it, but some spaces are not active.
 * 1d, In F part, each input vector must be divided exactly by (max node per cell)
 * 1e, In M part, each input vector must be divided exactly by (max node per cell)
 * 2a, Sum of number of cells from init FNM and init MX input vectors must be size of cellTypes
 *     so that all cells will have its own type
 * 2b, Read number of node per cell, etc, from config file.
 * 3a, First part of this function is error checking.
 * 3b, Second part of this function is the actual initialization
 */
void SimulationDomainGPU::initialCellsOfFiveTypes(
		std::vector<CellType> &cellTypes,
		std::vector<uint> &numOfInitActiveNodesOfCells,
		std::vector<double> &initBdryCellNodePosX,
		std::vector<double> &initBdryCellNodePosY,
		std::vector<double> &initProfileNodePosX,
		std::vector<double> &initProfileNodePosY,
		std::vector<double> &initECMNodePosX,
		std::vector<double> &initECMNodePosY,
		std::vector<double> &initFNMCellNodePosX,
		std::vector<double> &initFNMCellNodePosY,
		std::vector<double> &initMXCellNodePosX,
		std::vector<double> &initMXCellNodePosY) {

	/*
	 * zero step: redefine nodes.
	 */
	/*
	 * number of boundary nodes is always fixed.
	 */
	uint bdryNodeCount = initBdryCellNodePosX.size();
	// total potential number of profile nodes could be more than initial number provided,
	// because the total length of profile might increase.
	// FinalToInitProfileNodeCountRatio is defined in Config file.
	uint maxProfileNodeCount = initProfileNodePosX.size()
			* FinalToInitProfileNodeCountRatio;

	/*
	 * Initialize SceNodes by constructor. first two parameters come from input parameters
	 * while the last four parameters come from Config file.
	 */
	nodes = SceNodes(bdryNodeCount, maxProfileNodeCount, maxECMInDomain,
			maxNodePerECM, maxCellInDomain, maxNodePerCell);
	/*
	 * Initialize SceCells_M ( M means modified) by nodes information.
	 */
	cells_m = SceCells_M(&nodes);

	/*
	 * first step: error checking.
	 * we need to first check if inputs are valid
	 */
	cout << "begin init cells of five types" << endl;
	// get max node per cell. should be defined previously.
	uint maxNodePerCell = nodes.maxNodeOfOneCell;
	// get max node per ECM. Should be defined previously.
	uint maxNodePerECM = nodes.getMaxNodePerEcm();
	// check if we successfully loaded maxNodePerCell
	assert(maxNodePerCell != 0);

	// obtain sizes of the input arrays
	uint bdryNodeCountX = initBdryCellNodePosX.size();
	uint bdryNodeCountY = initBdryCellNodePosY.size();
	uint ProfileNodeCountX = initProfileNodePosX.size();
	uint ProfileNodeCountY = initProfileNodePosY.size();
	uint ECMNodeCountX = initECMNodePosX.size();
	uint ECMNodeCountY = initECMNodePosY.size();
	uint FNMNodeCountX = initFNMCellNodePosX.size();
	uint FNMNodeCountY = initFNMCellNodePosY.size();
	uint MXNodeCountX = initMXCellNodePosX.size();
	uint MXNodeCountY = initMXCellNodePosY.size();

	cout << "size of all node vectors: Boundary = " << bdryNodeCountX
			<< ", init profile node count =" << "," << ProfileNodeCountX
			<< "max profile node count = " << maxProfileNodeCount
			<< ", initial ECM node count = " << ECMNodeCountX
			<< ", init FNM node count = " << FNMNodeCountX
			<< ", init MX node count = " << MXNodeCountX << endl;

	// array size of cell type array
	uint cellTypeSize = cellTypes.size();
	// array size of initial active node count of cells array.
	uint initNodeCountSize = numOfInitActiveNodesOfCells.size();
	// two sizes must match.
	assert(cellTypeSize == initNodeCountSize);
	// size of X and Y must match.
	assert(bdryNodeCountX == bdryNodeCountY);
	assert(ECMNodeCountX == ECMNodeCountY);
	assert(ProfileNodeCountX == ProfileNodeCountY);
	assert(FNMNodeCountX == FNMNodeCountY);
	assert(MXNodeCountX == MXNodeCountY);

	cout << "passed init checks" << endl;

	// size of inputs must be divided exactly by max node per cell.
	// uint bdryRemainder = bdryNodeCountX % maxNodePerCell;
	uint ecmRemainder = ECMNodeCountX % maxNodePerECM;
	uint fnmRemainder = FNMNodeCountX % maxNodePerCell;
	uint mxRemainder = MXNodeCountX % maxNodePerCell;

	// uint bdryQuotient = bdryNodeCountX / maxNodePerCell;
	uint ecmQuotient = ECMNodeCountX / maxNodePerECM;
	uint fnmQuotient = FNMNodeCountX / maxNodePerCell;
	uint mxQuotient = MXNodeCountX / maxNodePerCell;

	// for now we try to make boundary cells one complete part so ....
	uint bdryRemainder = 0;
	uint bdryQuotient = 1;

	// for now we try to make profile nodes one complete part soremdiner = 0 and quotient = 1
	uint profileRemainder = 0;
	uint profileQuotient = 1;

	// remainder must be zero.
	assert((fnmRemainder == 0) && (mxRemainder == 0) && (ecmRemainder == 0));
	// size of cellType array and sum of all cell types must match.
	assert(fnmQuotient + mxQuotient == cellTypeSize);

	cerr << "passed size assertion" << endl;

	// make sure the cell types follow format requirement.
	// must follow sequence : B - P - E - F - M
	int counter = 0;
	CellType cellTypesForEachLevel[5] = { Boundary, Profile, ECM, FNM, MX };
	int bounds[5];
	bounds[0] = bdryQuotient;
	bounds[1] = bounds[0] + profileQuotient;
	bounds[2] = bounds[1] + ecmQuotient;
	bounds[3] = bounds[2] + fnmQuotient;
	bounds[4] = bounds[3] + mxQuotient;
	int level = 0;
	while (counter < cellTypeSize) {
		// if count is already beyond the bound, we need to increase the current level.
		if (counter == bounds[level]) {
			level++;
		}
		// make sure that the input cell types array fits the calculated result.
		// depreciated -- requirement changed.
		// assert(cellTypes[counter] == cellTypesForEachLevel[level]);
		counter++;
	}
	cerr << "before set parameters" << endl;
	/*
	 * second part: actual initialization
	 * copy data from main system memory to GPU memory
	 */

	nodes.setCurrentActiveCellCount(fnmQuotient + mxQuotient);
	cells_m.currentActiveCellCount = fnmQuotient + mxQuotient;
	nodes.setCurrentActiveEcm(ecmQuotient);
	cells_m.currentActiveECMCount = ecmQuotient;

	nodes.currentActiveProfileNodeCount = ProfileNodeCountX;

	cerr << "after set parameters" << endl;
	// copy input of initial active node of cells to our actual data location
	//thrust::copy(numOfInitActiveNodesOfCells.begin(),
	//		numOfInitActiveNodesOfCells.end(),
	//		cells_m.activeNodeCountOfThisCell.begin());

	// find the begining position of Profile.

	assert(nodes.startPosProfile == bdryNodeCountX);
	uint beginAddressOfProfile = nodes.startPosProfile;
	// find the begining position of ECM.
	uint beginAddressOfECM = nodes.startPosECM;
	// find the begining position of FNM cells.
	uint beginAddressOfFNM = nodes.startPosCells;
	// find the begining position of MX cells.
	uint beginAddressOfMX = beginAddressOfFNM + FNMNodeCountX;

	// initialize the cell ranks and types

	int numOfCellEachLevel[] = { bdryQuotient, profileQuotient, ecmQuotient,
			fnmQuotient, mxQuotient };
	int nodesPerCellEachLevel[] = { bdryNodeCount, maxProfileNodeCount,
			maxNodePerECM, maxNodePerCell, maxNodePerCell };
	uint totalSize = nodes.nodeLocX.size();
	//thrust::host_vector<CellType> allNodeTypes;
	//thrust::host_vector<int> cellRanks;
	//allNodeTypes.resize(totalSize);
	//cellRanks.resize(totalSize);
	//int currentRank = 0;
	//level = 0;
	//counter = 0;
	//while (counter < cellTypeSize) {
	// if count is already beyond the bound, we need to increase the current level.
	//	if (counter == bounds[level]) {
	//		level++;
	//	}
	//	allNodeTypes[counter] = cellTypesForEachLevel[level];
	//	if (level == 0) {
	//		cellRanks[counter] = counter / nodesPerCellEachLevel[0];
	//		//currentRank = cellRanks[counter];
	//	} else {
	//		cellRanks[counter] = (counter - bounds[level - 1])
	//				/ nodesPerCellEachLevel[level];
	//	}
	//	counter++;
	//}
	// copy node rank and node type information to GPU
	//thrust::copy(allNodeTypes.begin(), allNodeTypes.end(),
	//		nodes.nodeCellType.begin());
	//thrust::copy(cellRanks.begin(), cellRanks.end(),
	//		nodes.nodeCellRank.begin());

	// copy x and y position of nodes of boundary cells to actual node position.

	// set cell types
	thrust::device_vector<CellType> cellTypesToPass = cellTypes;

	cout << "Break point 1111" << endl;

	cout << "size1 = " << numOfInitActiveNodesOfCells.size() << endl;
	cout << "size2 = " << cells_m.activeNodeCountOfThisCell.size() << endl;

	// copy initial active node count info to GPU
	thrust::copy(numOfInitActiveNodesOfCells.begin(),
			numOfInitActiveNodesOfCells.end(),
			cells_m.activeNodeCountOfThisCell.begin());

	cout << "Break point 2222" << endl;

	cout << "number of initial cells: " << numOfInitActiveNodesOfCells.size();

	nodes.initValues(initBdryCellNodePosX, initBdryCellNodePosY,
			initProfileNodePosX, initProfileNodePosY, initECMNodePosX,
			initECMNodePosY, initFNMCellNodePosX, initFNMCellNodePosY,
			initMXCellNodePosX, initMXCellNodePosY);

	//for (int i = 0; i < numOfInitActiveNodesOfCells.size(); i++) {
	//	cout << cells_m.activeNodeCountOfThisCell[i] << " " << endl;
	//}
	// set isActiveInfo
	// allocate space for isActive info
	// uint sizeOfTmpVector = maxNodePerCell * initNodeCountSize;
	// thrust::host_vector<bool> isActive(sizeOfTmpVector, false);

	//thrust::host_vector<bool> isActive(totalSize, false);

	//for (int i = 0; i < initNodeCountSize; i++) {
	//	int j = 0;
	//	int index;
	//	while (j < numOfInitActiveNodesOfCells[i]) {
	//		index = i * maxNodePerCell + j;
	//		isActive[index] = true;
	//		j++;
	//	}
	//}
	/*
	 for (int i = 0; i < totalSize; i++) {
	 if (i < beginAddressOfProfile) {
	 isActive[i] = true;
	 } else if (i < beginAddressOfECM) {
	 if (i - beginAddressOfProfile < ProfileNodeCountX) {
	 isActive[i] = true;
	 } else {
	 isActive[i] = false;
	 }
	 } else if (i < beginAddressOfFNM) {
	 if (i - beginAddressOfECM < initECMCount * maxNodePerECM) {
	 isActive[i] = true;
	 } else {
	 isActive[i] = false;
	 }

	 } else {
	 // else : initially we don't need to set active for FNM and MX because
	 // this information will be updated while cell logic.
	 isActive[i] = false;
	 }
	 }
	 // copy is active info to GPU
	 thrust::copy(isActive.begin(), isActive.end(), nodes.nodeIsActive.begin());

	 */
	// set cell types
	cells_m.setCellTypes(cellTypesToPass);

	cout << "finished initialize cells of five types" << endl;
	cout << "press any key to continue" << endl;
	//int jj;
	//cin >> jj;
}

void SimulationDomainGPU::initialize_V2(SimulationInitData &initData) {
	initialCellsOfFiveTypes(initData.cellTypes,
			initData.numOfInitActiveNodesOfCells, initData.initBdryCellNodePosX,
			initData.initBdryCellNodePosY, initData.initProfileNodePosX,
			initData.initProfileNodePosY, initData.initECMNodePosX,
			initData.initECMNodePosY, initData.initFNMCellNodePosX,
			initData.initFNMCellNodePosY, initData.initMXCellNodePosX,
			initData.initMXCellNodePosY);
	cout << "finished init cells of five types" << endl;
	nodes.initDimension(minX, maxX, minY, maxY, gridSpacing);
}

/**
 * Highest level logic of domain.
 *
 */
void SimulationDomainGPU::runAllLogic(double dt) {

	nodes.calculateAndApplySceForces();
	//nodes.move(dt);
	//cells.growAndDivide(dt, growthMap.growthFactorMag,
	//		growthMap.growthFactorDirXComp, growthMap.growthFactorDirYComp,
	//		growthMap.gridDimensionX, growthMap.gridDimensionY,
	//		growthMap.gridSpacing);
	// cells.growAndDivide(dt, growthMap, growthMap2);
	cells_m.runAllCellLevelLogics(dt, growthMap, growthMap2);
	// cells.runAllCellLevelLogics(dt,growthMap,growthMap2);
}

/**
 * Depreciated.
 *
 */
void SimulationDomainGPU::outputVtkFilesWithColor(std::string scriptNameBase,
		int rank) {

	uint activeTotalNodeCount = nodes.currentActiveCellCount
			* nodes.maxNodeOfOneCell;

	thrust::host_vector<uint> hostActiveCountOfCells(
			nodes.currentActiveCellCount);

	thrust::host_vector<double> hostTmpVectorLocX(activeTotalNodeCount);
	thrust::host_vector<double> hostTmpVectorLocY(activeTotalNodeCount);
	thrust::host_vector<double> hostTmpVectorLocZ(activeTotalNodeCount);
	thrust::host_vector<bool> hostTmpVectorIsActive(activeTotalNodeCount);

	std::vector<uint> prefixSum(nodes.currentActiveCellCount);
	std::vector<uint> prefixSumLinks(nodes.currentActiveCellCount);

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes.nodeLocX.begin(),
							nodes.nodeLocY.begin(), nodes.nodeLocZ.begin(),
							nodes.nodeIsActive.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes.nodeLocX.begin(),
							nodes.nodeLocY.begin(), nodes.nodeLocZ.begin(),
							nodes.nodeIsActive.begin())) + activeTotalNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(hostTmpVectorLocX.begin(),
							hostTmpVectorLocY.begin(),
							hostTmpVectorLocZ.begin(),
							hostTmpVectorIsActive.begin())));

	thrust::copy(cells_m.activeNodeCountOfThisCell.begin(),
			cells_m.activeNodeCountOfThisCell.begin()
					+ cells_m.currentActiveCellCount,
			hostActiveCountOfCells.begin());

	int i, j, k;
	int totalNNum = thrust::reduce(hostActiveCountOfCells.begin(),
			hostActiveCountOfCells.end());
	std::vector<std::pair<uint, uint> > links;
	uint tmpRes = 0;
	for (i = 0; i < nodes.currentActiveCellCount; i++) {
		prefixSum[i] = tmpRes;
		tmpRes = tmpRes + hostActiveCountOfCells[i];
	}

	thrust::host_vector<CellType> cellTypesFromGPU = cells_m.getCellTypes();
	cout << "size of cellTypes from GPU is " << cellTypesFromGPU.size()
			<< ", size of currentActiveCellCount is"
			<< nodes.currentActiveCellCount << endl;
//assert(cellTypesFromGPU.size() == nodes.currentActiveCellCount);
// using string stream is probably not the best solution,
// but I can't use c++ 11 features for backward compatibility
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string scriptNameRank = ss.str();
	std::string vtkFileName = scriptNameBase + scriptNameRank + ".vtk";
	std::cout << "start to create vtk file" << vtkFileName << std::endl;
	std::ofstream fs;
	fs.open(vtkFileName.c_str());

//int totalNNum = getTotalNodeCount();
//int LNum = 0;
//int NNum;
	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << "Lines and points representing subcelluar element cells "
			<< std::endl;
	fs << "ASCII" << std::endl;
	fs << std::endl;
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << totalNNum << " float" << std::endl;

	uint counterForLink = 0;
	for (i = 0; i < nodes.currentActiveCellCount; i++) {
		uint activeNodeCount = hostActiveCountOfCells[i];
		for (j = 0; j < activeNodeCount; j++) {
			uint pos = i * nodes.maxNodeOfOneCell + j;
			fs << hostTmpVectorLocX[pos] << " " << hostTmpVectorLocY[pos] << " "
					<< hostTmpVectorLocZ[pos] << std::endl;
		}
		uint tmpCount = 0;
		for (j = 0; j < activeNodeCount; j++) {
			for (k = j; k < activeNodeCount; k++) {
				if (j == k) {
					continue;
				} else {
					uint pos1 = prefixSum[i] + j;
					uint pos2 = prefixSum[i] + k;
					uint pos1InVec = i * nodes.maxNodeOfOneCell + j;
					uint pos2InVec = i * nodes.maxNodeOfOneCell + k;
					if (compuDist(hostTmpVectorLocX[pos1InVec],
							hostTmpVectorLocY[pos1InVec],
							hostTmpVectorLocZ[pos1InVec],
							hostTmpVectorLocX[pos2InVec],
							hostTmpVectorLocY[pos2InVec],
							hostTmpVectorLocZ[pos2InVec])
							<= intraLinkDisplayRange) {
						links.push_back(std::make_pair<uint, uint>(pos1, pos2));
						tmpCount++;
					} else {
						continue;
					}
				}
			}
		}
		prefixSumLinks[i] = counterForLink;
		counterForLink = counterForLink + tmpCount;
	}
	fs << std::endl;
	fs << "CELLS " << counterForLink << " " << 3 * counterForLink << std::endl;
	uint linkSize = links.size();
	for (uint i = 0; i < linkSize; i++) {
		fs << 2 << " " << links[i].first << " " << links[i].second << std::endl;
	}
	uint LNum = links.size();
	fs << "CELL_TYPES " << LNum << endl;
	for (i = 0; i < LNum; i++) {
		fs << "3" << endl;
	}
	fs << "POINT_DATA " << totalNNum << endl;
	fs << "SCALARS point_scalars float" << endl;
	fs << "LOOKUP_TABLE default" << endl;

//for (i = 0; i < nodes.currentActiveCellCount; i++) {
//	uint activeNodeCount = hostActiveCountOfCells[i];
//	for (j = 0; j < activeNodeCount; j++) {
//		fs << i << endl;
//	}
//}

	for (i = 0; i < nodes.currentActiveCellCount; i++) {
		uint activeNodeCount = hostActiveCountOfCells[i];
		if (cellTypesFromGPU[i] == Boundary) {
			for (j = 0; j < activeNodeCount; j++) {
				fs << 1 << endl;
			}
		} else if (cellTypesFromGPU[i] == FNM) {
			for (j = 0; j < activeNodeCount; j++) {
				fs << 2 << endl;
			}
		} else if (cellTypesFromGPU[i] == MX) {
			for (j = 0; j < activeNodeCount; j++) {
				fs << 3 << endl;
			}
		}

	}

	fs.flush();
	fs.close();
}

/**
 * This is the second version of visualization.
 * Inter-links are not colored for a clearer representation
 * Fixed Boundary: Black
 * FNM cells: Red
 * MX cells: Blue
 * ECM: Yellow
 * Moving Epithilum layer: Green
 */
void SimulationDomainGPU::outputVtkFilesWithColor_v2(std::string scriptNameBase,
		int rank) {
	cerr << "start to output animation " << rank << endl;
	//cells_m.distributeIsCellRank();
	cells_m.distributeIsActiveInfo();
	thrust::host_vector<bool> isActiveHost = nodes.nodeIsActive;
	cout << "finished distribute cell is active" << endl;

	int num1 = nodes.startPosCells + cells_m.totalNodeCountForActiveCells;
	int num2 = cells_m.beginPosOfCellsNode
			+ nodes.currentActiveCellCount * nodes.maxNodeOfOneCell;
	cout << "Total active nodes estimation1 = " << num1 << ", part 1 = "
			<< nodes.startPosCells << ", part 2 = "
			<< cells_m.totalNodeCountForActiveCells << endl;
	cout << "Total active nodes estimation2 = " << num2 << ", part 1 = "
			<< cells_m.beginPosOfCellsNode << ", part 2 = "
			<< nodes.currentActiveCellCount * nodes.maxNodeOfOneCell << endl;

	int count = 0;
	for (int i = 0;
			i < (nodes.startPosCells + cells_m.totalNodeCountForActiveCells);
			i++) {
		if (isActiveHost[i] == true) {
			count++;
		}
	}
	cout << "before printing, number of active nodes in domain: " << count
			<< endl;
	//int jj;
	//cin >> jj;

	uint activeTotalNodeCount = cells_m.beginPosOfCellsNode
			+ nodes.currentActiveCellCount * nodes.maxNodeOfOneCell;

	cout << "number of nodes = " << activeTotalNodeCount << endl;
	cout << "total nodes in the bool array = " << nodes.nodeIsActive.size()
			<< endl;

	uint totalActiveCount = thrust::reduce(nodes.nodeIsActive.begin(),
			nodes.nodeIsActive.begin() + activeTotalNodeCount, (int) (0));

	//uint totalActiveCount = thrust::reduce(nodes.nodeIsActive.begin(),
	//		nodes.nodeIsActive.begin() + 1);
	cerr << "number of active nodes = " << totalActiveCount << endl;
	cerr << "number of total active nodes possible= " << activeTotalNodeCount
			<< endl;

	cout.flush();

	thrust::device_vector<double> deviceTmpVectorLocX(totalActiveCount);
	thrust::device_vector<double> deviceTmpVectorLocY(totalActiveCount);
	thrust::device_vector<double> deviceTmpVectorLocZ(totalActiveCount);
	thrust::device_vector<bool> deviceTmpVectorIsActive(totalActiveCount);
	thrust::device_vector<CellType> deviceTmpVectorNodeType(totalActiveCount);
	thrust::device_vector<uint> deviceTmpVectorNodeRank(totalActiveCount);

	thrust::host_vector<double> hostTmpVectorLocX(totalActiveCount);
	thrust::host_vector<double> hostTmpVectorLocY(totalActiveCount);
	thrust::host_vector<double> hostTmpVectorLocZ(totalActiveCount);
	thrust::host_vector<bool> hostTmpVectorIsActive(totalActiveCount);
	thrust::host_vector<CellType> hostTmpVectorNodeType(totalActiveCount);
	thrust::host_vector<uint> hostTmpVectorNodeRank(totalActiveCount);

	cout << "finished initialization space for GPU and CPU" << endl;

	cout << "During output animation, active total node count is "
			<< activeTotalNodeCount << endl;

	thrust::copy_if(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes.nodeLocX.begin(),
							nodes.nodeLocY.begin(), nodes.nodeLocZ.begin(),
							nodes.nodeCellType.begin(),
							nodes.nodeCellRank.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes.nodeLocX.begin(),
							nodes.nodeLocY.begin(), nodes.nodeLocZ.begin(),
							nodes.nodeCellType.begin(),
							nodes.nodeCellRank.begin())) + activeTotalNodeCount,
			nodes.nodeIsActive.begin(),
			thrust::make_zip_iterator(
					thrust::make_tuple(deviceTmpVectorLocX.begin(),
							deviceTmpVectorLocY.begin(),
							deviceTmpVectorLocZ.begin(),
							deviceTmpVectorNodeType.begin(),
							deviceTmpVectorNodeRank.begin())), isTrue());

	cout << "finished cpu data from GPU to CPU" << endl;

	hostTmpVectorLocX = deviceTmpVectorLocX;
	hostTmpVectorLocY = deviceTmpVectorLocY;
	hostTmpVectorLocZ = deviceTmpVectorLocZ;
	hostTmpVectorIsActive = deviceTmpVectorIsActive;
	hostTmpVectorNodeType = deviceTmpVectorNodeType;
	hostTmpVectorNodeRank = deviceTmpVectorNodeRank;

	int i, j;
	std::vector<std::pair<uint, uint> > links;

	//assert(cellTypesFromGPU.size() == nodes.currentActiveCellCount);
	// using string stream is probably not the best solution,
	// but I can't use c++ 11 features for backward compatibility
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string scriptNameRank = ss.str();
	std::string vtkFileName = scriptNameBase + scriptNameRank + ".vtk";
	std::cout << "start to create vtk file" << vtkFileName << std::endl;
	std::ofstream fs;
	fs.open(vtkFileName.c_str());

	//int totalNNum = getTotalNodeCount();
	//int LNum = 0;
	//int NNum;
	assert(hostTmpVectorLocX.size() == totalActiveCount);

	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << "Lines and points representing subcelluar element cells "
			<< std::endl;
	fs << "ASCII" << std::endl;
	fs << std::endl;
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << totalActiveCount << " float" << std::endl;

	uint counterForLink = 0;
	for (i = 0; i < totalActiveCount; i++) {
		fs << hostTmpVectorLocX[i] << " " << hostTmpVectorLocY[i] << " "
				<< hostTmpVectorLocZ[i] << std::endl;
		for (j = 0; j < totalActiveCount; j++) {
			if (i == j) {
				continue;
			}
			if (compuDist(hostTmpVectorLocX[i], hostTmpVectorLocY[i],
					hostTmpVectorLocZ[i], hostTmpVectorLocX[j],
					hostTmpVectorLocY[j], hostTmpVectorLocZ[j])
					<= intraLinkDisplayRange) {
				// have this extra if because we don't want to include visualization
				// for inter-cell interactions.
				// to achieve that, we don't include links between different types
				// and also avoid display links between different cells of same type
				if (hostTmpVectorNodeType[i] == hostTmpVectorNodeType[j]
						&& hostTmpVectorNodeRank[i]
								== hostTmpVectorNodeRank[j]) {
					if (hostTmpVectorNodeType[i] == Boundary) {
						if (j - i == 1) {
							links.push_back(std::make_pair<uint, uint>(i, j));
							counterForLink++;
						}
					} else if (hostTmpVectorNodeType[i] == Profile) {
						if (j - i == 1) {
							links.push_back(std::make_pair<uint, uint>(i, j));
							counterForLink++;
						}
					} else if (hostTmpVectorNodeType[i] == ECM) {
						if (j - i == 1) {
							links.push_back(std::make_pair<uint, uint>(i, j));
							counterForLink++;
						}
					} else if (hostTmpVectorNodeType[i] == MX
							|| hostTmpVectorNodeType[i] == FNM) {
						links.push_back(std::make_pair<uint, uint>(i, j));
						counterForLink++;
					}
				}
			} else {
				//if (hostTmpVectorNodeType[i] == Profile) {
				//	if (j - i == 1) {
				//		links.push_back(std::make_pair<uint, uint>(i, j));
				//		counterForLink++;
				//	}
				//}
			}
		}
		//std::cout << "link count = " << counterForLink << std::endl;
	}
	fs << std::endl;
	fs << "CELLS " << counterForLink << " " << 3 * counterForLink << std::endl;
	uint linkSize = links.size();
	for (uint i = 0; i < linkSize; i++) {
		fs << 2 << " " << links[i].first << " " << links[i].second << std::endl;
	}
	uint LNum = links.size();
	fs << "CELL_TYPES " << LNum << endl;
	for (i = 0; i < LNum; i++) {
		fs << "3" << endl;
	}
	fs << "POINT_DATA " << totalActiveCount << endl;
	fs << "SCALARS point_scalars float" << endl;
	fs << "LOOKUP_TABLE default" << endl;

	for (i = 0; i < totalActiveCount; i++) {
		if (hostTmpVectorNodeType[i] == Boundary) {
			fs << 1 << endl;
		} else if (hostTmpVectorNodeType[i] == Profile) {
			fs << 2 << endl;
		} else if (hostTmpVectorNodeType[i] == ECM) {
			fs << 3 << endl;
		} else if (hostTmpVectorNodeType[i] == FNM) {
			fs << 4 << endl;
		} else if (hostTmpVectorNodeType[i] == MX) {
			fs << 5 << endl;
		}
	}

	//fs.flush();
	fs.close();
}

void SimulationDomainGPU::outputVtkFilesWithColor_v3(std::string scriptNameBase,
		int rank) {
	cout << "before calling obtaining neighbors" << endl;
	//std::vector<std::pair<uint, uint> > pairs = nodes.obtainNeighborPairs();
	uint activeTotalNodeCount = cells_m.beginPosOfCellsNode
			+ nodes.currentActiveCellCount * nodes.maxNodeOfOneCell;

	uint activeNodeCount = thrust::reduce(nodes.nodeIsActive.begin(),
			nodes.nodeIsActive.begin() + activeTotalNodeCount, (int) (0));

	cout << "before creating vectors" << endl;

	thrust::device_vector<double> deviceTmpVectorLocX(activeTotalNodeCount);
	thrust::device_vector<double> deviceTmpVectorLocY(activeTotalNodeCount);
	thrust::device_vector<double> deviceTmpVectorLocZ(activeTotalNodeCount);
	thrust::device_vector<CellType> deviceTmpVectorNodeType(
			activeTotalNodeCount);
	thrust::device_vector<uint> deviceTmpVectorNodeRank(activeTotalNodeCount);
	thrust::device_vector<bool> deviceTmpVectorIsActive(activeTotalNodeCount);

	thrust::host_vector<double> hostTmpVectorLocX(activeTotalNodeCount);
	thrust::host_vector<double> hostTmpVectorLocY(activeTotalNodeCount);
	thrust::host_vector<double> hostTmpVectorLocZ(activeTotalNodeCount);
	thrust::host_vector<CellType> hostTmpVectorNodeType(activeTotalNodeCount);
	thrust::host_vector<uint> hostTmpVectorNodeRank(activeTotalNodeCount);
	thrust::host_vector<bool> hostTmpVectorIsActive(activeTotalNodeCount);

	cout << "finished initialization space for GPU and CPU" << endl;

	cout << "During output animation, active total node count is "
			<< activeTotalNodeCount << endl;

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes.nodeLocX.begin(),
							nodes.nodeLocY.begin(), nodes.nodeLocZ.begin(),
							nodes.nodeCellType.begin(),
							nodes.nodeCellRank.begin(),
							nodes.nodeIsActive.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodes.nodeLocX.begin(),
							nodes.nodeLocY.begin(), nodes.nodeLocZ.begin(),
							nodes.nodeCellType.begin(),
							nodes.nodeCellRank.begin(),
							nodes.nodeIsActive.begin())) + activeTotalNodeCount,
			thrust::make_zip_iterator(
					thrust::make_tuple(deviceTmpVectorLocX.begin(),
							deviceTmpVectorLocY.begin(),
							deviceTmpVectorLocZ.begin(),
							deviceTmpVectorNodeType.begin(),
							deviceTmpVectorNodeRank.begin(),
							deviceTmpVectorIsActive.begin())));
	hostTmpVectorLocX = deviceTmpVectorLocX;
	hostTmpVectorLocY = deviceTmpVectorLocY;
	hostTmpVectorLocZ = deviceTmpVectorLocZ;
	hostTmpVectorNodeType = deviceTmpVectorNodeType;
	hostTmpVectorNodeRank = deviceTmpVectorNodeRank;
	hostTmpVectorIsActive = deviceTmpVectorIsActive;

	int i, j;

	//assert(cellTypesFromGPU.size() == nodes.currentActiveCellCount);
	// using string stream is probably not the best solution,
	// but I can't use c++ 11 features for backward compatibility
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string scriptNameRank = ss.str();
	std::string vtkFileName = scriptNameBase + scriptNameRank + ".vtk";
	std::cout << "start to create vtk file" << vtkFileName << std::endl;
	std::ofstream fs;
	fs.open(vtkFileName.c_str());

	//int totalNNum = getTotalNodeCount();
	//int LNum = 0;
	//int NNum;
	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << "Lines and points representing subcelluar element cells "
			<< std::endl;
	fs << "ASCII" << std::endl;
	fs << std::endl;
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << activeNodeCount << " float" << std::endl;

	for (i = 0; i < activeTotalNodeCount; i++) {
		if (hostTmpVectorIsActive[i] == true) {
			fs << hostTmpVectorLocX[i] << " " << hostTmpVectorLocY[i] << " "
					<< hostTmpVectorLocZ[i] << std::endl;
			assert(!isnan(hostTmpVectorLocX[i]));
		}
	}
	std::vector<std::pair<uint, uint> > links = nodes.obtainNeighborPairs();
	uint pairSize = links.size();
	std::vector<std::pair<uint, uint> > processedLinks;

	for (i = 0; i < pairSize; i++) {
		int n1 = links[i].first;
		int n2 = links[i].second;
		if (compuDist(hostTmpVectorLocX[n1], hostTmpVectorLocY[n1],
				hostTmpVectorLocZ[n1], hostTmpVectorLocX[n2],
				hostTmpVectorLocY[n2], hostTmpVectorLocZ[n2])
				<= intraLinkDisplayRange) {
			// have this extra if because we don't want to include visualization
			// for inter-cell interactions.
			// to achieve that, we don't include links between different types
			// and also avoid display links between different cells of same type
			if (hostTmpVectorNodeType[n1] == hostTmpVectorNodeType[n2]
					&& hostTmpVectorNodeRank[n1] == hostTmpVectorNodeRank[n2]) {
				if (hostTmpVectorNodeType[n1] == Boundary) {
					if (n2 - n1 == 1) {
						processedLinks.push_back(
								std::make_pair<uint, uint>(n1, n2));
					}
				} else if (hostTmpVectorNodeType[n1] == Profile) {
					if (n2 - n1 == 1) {
						processedLinks.push_back(
								std::make_pair<uint, uint>(n1, n2));
					}
				} else if (hostTmpVectorNodeType[n1] == ECM) {
					if (n2 - n1 == 1) {
						processedLinks.push_back(
								std::make_pair<uint, uint>(n1, n2));
					}
				} else if (hostTmpVectorNodeType[n1] == MX
						|| hostTmpVectorNodeType[n1] == FNM) {
					processedLinks.push_back(
							std::make_pair<uint, uint>(n1, n2));
				}
			}
		}

		//std::cout << "link count = " << counterForLink << std::endl;
	}
	fs << std::endl;
	fs << "CELLS " << processedLinks.size() << " " << 3 * processedLinks.size()
			<< std::endl;
	uint linkSize = processedLinks.size();
	for (uint i = 0; i < linkSize; i++) {
		fs << 2 << " " << processedLinks[i].first << " "
				<< processedLinks[i].second << std::endl;
	}
	uint LNum = processedLinks.size();
	fs << "CELL_TYPES " << LNum << endl;
	for (i = 0; i < LNum; i++) {
		fs << "3" << endl;
	}
	fs << "POINT_DATA " << activeNodeCount << endl;
	fs << "SCALARS point_scalars float" << endl;
	fs << "LOOKUP_TABLE default" << endl;

//for (i = 0; i < nodes.currentActiveCellCount; i++) {
//	uint activeNodeCount = hostActiveCountOfCells[i];
//	for (j = 0; j < activeNodeCount; j++) {
//		fs << i << endl;
//	}
//}

	for (i = 0; i < activeNodeCount; i++) {
		if (hostTmpVectorNodeType[i] == Boundary) {
			fs << 1 << endl;
		} else if (hostTmpVectorNodeType[i] == Profile) {
			fs << 2 << endl;
		} else if (hostTmpVectorNodeType[i] == ECM) {
			fs << 3 << endl;
		} else if (hostTmpVectorNodeType[i] == FNM) {
			fs << 4 << endl;
		} else if (hostTmpVectorNodeType[i] == MX) {
			fs << 5 << endl;
		}
	}

	fs.flush();
	fs.close();
}

void SimulationDomainGPU::checkIfAllDataFieldsValid() {
	cout << "Begin output information about nodes:" << endl;
	cout << "size of isActive:" << nodes.nodeIsActive.size() << endl;
	cout << "size of nodeLocX:" << nodes.nodeLocX.size() << endl;
	cout << "size of nodeLocY:" << nodes.nodeLocY.size() << endl;
	cout << "size of nodeLocZ:" << nodes.nodeLocZ.size() << endl;
	cout << "size of nodeVelX:" << nodes.nodeVelX.size() << endl;
	cout << "size of nodeVelY:" << nodes.nodeVelY.size() << endl;
	cout << "size of nodeVelZ:" << nodes.nodeVelZ.size() << endl;
	cout << "size of CellType:" << nodes.nodeCellType.size() << endl;
	cout << "size of nodeCellRank:" << nodes.nodeCellRank.size() << endl;

	cout << "start position of Profile is " << nodes.startPosProfile << endl;
	cout << "start position of ECM is " << nodes.startPosECM << endl;
	cout << "start position of Cells is " << nodes.startPosCells << endl;

	cout << "max node of one cell is " << nodes.maxNodeOfOneCell << endl;
	cout << "max number of cells is " << nodes.maxCellCount << endl;
	cout << "max total cell node count is " << nodes.maxTotalCellNodeCount
			<< endl;
	cout << "current active cell count is " << nodes.currentActiveCellCount
			<< endl;

	cout << "max node of one ECM is " << nodes.maxNodePerECM << endl;
	cout << "max number of ECm is " << nodes.maxECMCount << endl;
	cout << "max total ECM node count is " << nodes.maxTotalECMNodeCount
			<< endl;
	cout << "current active ECM count is " << nodes.currentActiveECM << endl;

	cout << "max profile node count is " << nodes.maxProfileNodeCount << endl;
	cout << "current active profile node count is "
			<< nodes.currentActiveProfileNodeCount << endl;
	//int jj;
	//cin >> jj;
}

