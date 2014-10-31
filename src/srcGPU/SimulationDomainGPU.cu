/**
 * @file SimulationDomainGPU.cu
 * @brief this file contains domain level logic.
 * @author Wenzhao Sun wsun2@nd.edu
 * @bug no know bugs
 */

#include "SimulationDomainGPU.h"

using namespace std;

/**
 * Constructor.
 * reads values from config file.
 */
SimulationDomainGPU::SimulationDomainGPU() {
	readAllParameters();
	initializeGrowthMap();
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
		std::vector<SceNodeType> &cellTypes,
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
			* memPara.FinalToInitProfileNodeCountRatio;

	/*
	 * Initialize SceNodes by constructor. first two parameters come from input parameters
	 * while the last four parameters come from Config file.
	 */
	nodes = SceNodes(bdryNodeCount, maxProfileNodeCount, memPara.maxECMInDomain,
			memPara.maxNodePerECM, memPara.maxCellInDomain,
			memPara.maxNodePerCell);

	/*
	 * first step: error checking.
	 * we need to first check if inputs are valid
	 */
	cout << "begin init cells of five types" << endl;
	// get max node per cell. should be defined previously.
	uint maxNodePerCell = nodes.getAllocPara().maxNodeOfOneCell;
	// get max node per ECM. Should be defined previously.
	uint maxNodePerECM = nodes.getAllocPara().maxNodePerECM;
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

	//int jj;
	//cin >> jj;

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
	//uint bdryRemainder = 0;
	uint bdryQuotient = 1;

	// for now we try to make profile nodes one complete part soremdiner = 0 and quotient = 1
	//uint profileRemainder = 0;
	uint profileQuotient = 1;

	// remainder must be zero.
	assert((fnmRemainder == 0) && (mxRemainder == 0) && (ecmRemainder == 0));
	// size of cellType array and sum of all cell types must match.
	assert(fnmQuotient + mxQuotient == cellTypeSize);

	cerr << "passed size assertion" << endl;

	// make sure the cell types follow format requirement.
	// must follow sequence : B - P - E - F - M
	int counter = 0;
	//SceNodeType cellTypesForEachLevel[5] = { Boundary, Profile, ECM, FNM, MX };
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

	NodeAllocPara para = nodes.getAllocPara();
	para.currentActiveCellCount = fnmQuotient + mxQuotient;
	para.currentActiveECM = ecmQuotient;
	para.currentActiveProfileNodeCount = ProfileNodeCountX;
	nodes.setAllocPara(para);

	//NodeAllocPara nodePara = nodes.getAllocPara();

	assert(nodes.getAllocPara().startPosProfile == bdryNodeCountX);

	uint totalSize = nodes.getInfoVecs().nodeLocX.size();

	// set cell types
	//thrust::device_vector<SceNodeType> cellTypesToPass = cellTypes;

	nodes.initValues(initBdryCellNodePosX, initBdryCellNodePosY,
			initProfileNodePosX, initProfileNodePosY, initECMNodePosX,
			initECMNodePosY, initFNMCellNodePosX, initFNMCellNodePosY,
			initMXCellNodePosX, initMXCellNodePosY);

	/*
	 * Initialize SceCells_M ( M means modified) by nodes information.
	 */
	//cells = SceCells(&nodes);
	// copy initial active node count info to GPU
	//thrust::copy(numOfInitActiveNodesOfCells.begin(),
	//		numOfInitActiveNodesOfCells.end(),
	//		cells.activeNodeCountOfThisCell.begin());
	// set cell types
	//cells.setCellTypes(cellTypesToPass);
	//cells.distributeIsActiveInfo();
	cells = SceCells(&nodes, numOfInitActiveNodesOfCells, cellTypes);
}

void SimulationDomainGPU::initializeNodes(CartPara &cartPara,
		std::vector<SceNodeType>& cellTypes,
		std::vector<uint>& numOfInitActiveNodesOfCells,
		std::vector<CVector>& initBdryNodeVec,
		std::vector<CVector>& initProfileNodeVec,
		std::vector<CVector>& initCartNodeVec,
		std::vector<CVector>& initECMNodeVec,
		std::vector<CVector>& initFNMNodeVec,
		std::vector<CVector>& initMXNodeVec) {

	/*
	 * number of boundary nodes is always fixed.
	 */
	uint bdryNodeCount = initBdryNodeVec.size();
	/*
	 * total potential number of profile nodes could be more than initial number provided,
	 * because the total length of profile might increase.
	 * FinalToInitProfileNodeCountRatio is defined in Config file.
	 */
	uint maxProfileNodeCount = initProfileNodeVec.size()
			* memPara.FinalToInitProfileNodeCountRatio;

	/**
	 * different from profile nodes, the cartilage nodes should already
	 * have buffer before executing this method.
	 */
	uint maxCartNodeCount = initCartNodeVec.size();

	/*
	 * Initialize SceNodes by constructor. first two parameters come from input parameters
	 * while the last four parameters come from Config file.
	 */
	nodes = SceNodes(bdryNodeCount, maxProfileNodeCount, maxCartNodeCount,
			memPara.maxECMInDomain, memPara.maxNodePerECM,
			memPara.maxCellInDomain, memPara.maxNodePerCell);

	//cartilage.distributeIsActive();
	//cartilage.initializeNodes(initCartNodeVec);

	/*
	 * first step: error checking.
	 * we need to first check if inputs are valid
	 */
	// get max node per cell. should be defined previously.
	uint maxNodePerCell = nodes.getAllocPara().maxNodeOfOneCell;
	// get max node per ECM. Should be defined previously.
	uint maxNodePerECM = nodes.getAllocPara().maxNodePerECM;
	// check if we successfully loaded maxNodePerCell
	assert(maxNodePerCell != 0);

	// obtain sizes of the input arrays
	//uint bdryNodeCount = initBdryNodeVec.size();
	uint ProfileNodeCount = initProfileNodeVec.size();
	uint CartNodeCount = initCartNodeVec.size();
	uint ECMNodeCount = initECMNodeVec.size();
	uint FNMNodeCount = initFNMNodeVec.size();
	uint MXNodeCount = initMXNodeVec.size();

	// array size of cell type array
	uint cellTypeSize = cellTypes.size();
	// array size of initial active node count of cells array.
	uint initNodeCountSize = numOfInitActiveNodesOfCells.size();
	// two sizes must match.
	assert(cellTypeSize == initNodeCountSize);

	// size of inputs must be divided exactly by max node per cell.
	// uint bdryRemainder = bdryNodeCountX % maxNodePerCell;
	uint ecmRemainder = ECMNodeCount % maxNodePerECM;
	uint fnmRemainder = FNMNodeCount % maxNodePerCell;
	uint mxRemainder = MXNodeCount % maxNodePerCell;

	// uint bdryQuotient = bdryNodeCountX / maxNodePerCell;
	uint ecmQuotient = ECMNodeCount / maxNodePerECM;
	uint fnmQuotient = FNMNodeCount / maxNodePerCell;
	uint mxQuotient = MXNodeCount / maxNodePerCell;

	// remainder must be zero.
	assert((fnmRemainder == 0) && (mxRemainder == 0) && (ecmRemainder == 0));
	// size of cellType array and sum of all cell types must match.
	assert(fnmQuotient + mxQuotient == cellTypeSize);

	for (uint i = 0; i < cellTypeSize; i++) {
		if (i < fnmQuotient) {
			assert(cellTypes[i] == FNM);
		} else {
			assert(cellTypes[i] == MX);
		}
	}
	/*
	 * second part: actual initialization
	 * copy data from main system memory to GPU memory
	 */
	NodeAllocPara para = nodes.getAllocPara();
	para.currentActiveCellCount = fnmQuotient + mxQuotient;
	para.currentActiveECM = ecmQuotient;
	para.currentActiveProfileNodeCount = ProfileNodeCount;
	nodes.setAllocPara(para);

	assert(nodes.getAllocPara().startPosProfile == bdryNodeCount);

	nodes.initValues_v2(initBdryNodeVec, initProfileNodeVec, initCartNodeVec,
			initECMNodeVec, initFNMNodeVec, initMXNodeVec);

	/**
	 * setting the cartilage related parameters in the simulation domain.
	 */
	if (memPara.simuType == Beak) {
		cartilage.setCartPara(cartPara);
		cartilage.initializeMem(&nodes);
	}

	cells = SceCells(&nodes, numOfInitActiveNodesOfCells, cellTypes);
}

void SimulationDomainGPU::initialize(SimulationInitData &initData) {
	initialCellsOfFiveTypes(initData.cellTypes,
			initData.numOfInitActiveNodesOfCells, initData.initBdryCellNodePosX,
			initData.initBdryCellNodePosY, initData.initProfileNodePosX,
			initData.initProfileNodePosY, initData.initECMNodePosX,
			initData.initECMNodePosY, initData.initFNMCellNodePosX,
			initData.initFNMCellNodePosY, initData.initMXCellNodePosX,
			initData.initMXCellNodePosY);
	std::cout << "finished init cells of five types" << std::endl;
	nodes.initDimension(domainPara.minX, domainPara.maxX, domainPara.minY,
			domainPara.maxY, domainPara.gridSpacing);
	std::cout << "finished init nodes dimension" << std::endl;
}

void SimulationDomainGPU::initialize_v2(SimulationInitData_V2& initData) {
	initializeNodes(initData.cartPara, initData.cellTypes,
			initData.numOfInitActiveNodesOfCells, initData.initBdryNodeVec,
			initData.initProfileNodeVec, initData.initCartNodeVec,
			initData.initECMNodeVec, initData.initFNMNodeVec,
			initData.initMXNodeVec);
	nodes.initDimension(domainPara.minX, domainPara.maxX, domainPara.minY,
			domainPara.maxY, domainPara.gridSpacing);
}

/**
 * Highest level logic of domain.
 *
 */
void SimulationDomainGPU::runAllLogic(double dt) {
	if (memPara.simuType == Beak) {
		nodes.processCartGrowthDir(cartilage.getCartPara().growthDir);
		//std::cout << "growth direction is ";
		cartilage.getCartPara().growthDir.Print();
	}

	// This function only calculates velocity.
	nodes.calculateAndApplySceForces();
	// Only beak simulation need to take care of cartilage.
	if (memPara.simuType == Beak) {
		// cartilage logics must come before cell logics, because node velocities will be modified
		// in cell logic and consequently we won't be able to compute cartilage data.
		// also responsible for handling interaction between epithelium layer and carilage.
		cartilage.runAllLogics(dt);
	}
	// This function applies velocity so nodes actually move inside this function.
	cells.runAllCellLevelLogics(dt, growthMap, growthMap2);

}

void SimulationDomainGPU::readMemPara() {
	int simuTypeConfigValue =
			globalConfigVars.getConfigValue("SimulationType").toInt();
	if (simuTypeConfigValue == 0) {
		memPara.simuType = Beak;
	} else if (simuTypeConfigValue == 1) {
		memPara.simuType = Disc;
	} else {
		throw SceException("Simulation Type in config file is not recognized!",
				ConfigValueException);
	}

	memPara.maxCellInDomain =
			globalConfigVars.getConfigValue("MaxCellInDomain").toInt();
	memPara.maxNodePerCell =
			globalConfigVars.getConfigValue("MaxNodePerCell").toInt();
	memPara.maxECMInDomain =
			globalConfigVars.getConfigValue("MaxECMInDomain").toInt();
	memPara.maxNodePerECM =
			globalConfigVars.getConfigValue("MaxNodePerECM").toInt();
	memPara.FinalToInitProfileNodeCountRatio = globalConfigVars.getConfigValue(
			"FinalToInitProfileNodeCountRatio").toDouble();
	//memPara.FinalToInitCartNodeCountRatio = globalConfigVars.getConfigValue(
	//		"FinalToInitCartNodeCountRatio").toDouble();
}

void SimulationDomainGPU::readDomainPara() {
	domainPara.minX = globalConfigVars.getConfigValue("DOMAIN_XMIN").toDouble();
	domainPara.maxX = globalConfigVars.getConfigValue("DOMAIN_XMAX").toDouble();
	domainPara.minY = globalConfigVars.getConfigValue("DOMAIN_YMIN").toDouble();
	domainPara.maxY = globalConfigVars.getConfigValue("DOMAIN_YMAX").toDouble();
	domainPara.minZ = globalConfigVars.getConfigValue("DOMAIN_ZMIN").toDouble();
	domainPara.maxZ = globalConfigVars.getConfigValue("DOMAIN_ZMAX").toDouble();
	domainPara.gridSpacing = globalConfigVars.getConfigValue(
			"DOMAIN_GRID_SPACING").toDouble();
	domainPara.numOfBucketsInXDim = (domainPara.maxX - domainPara.minX)
			/ domainPara.gridSpacing + 1;
	domainPara.numOfBucketsInYDim = (domainPara.maxY - domainPara.minY)
			/ domainPara.gridSpacing + 1;
}

void SimulationDomainGPU::readChemPara() {
	chemPara.growthGridXDim =
			globalConfigVars.getConfigValue("GrowthGridXDim").toInt();
	chemPara.growthGridYDim =
			globalConfigVars.getConfigValue("GrowthGridYDim").toInt();
	chemPara.growthGridSpacing = globalConfigVars.getConfigValue(
			"GrowthGridSpacing").toDouble();
	chemPara.growthGridLowerLeftPtX = globalConfigVars.getConfigValue(
			"GrowthGridLowerLeftPtX").toDouble();
	chemPara.growthGridLowerLeftPtY = globalConfigVars.getConfigValue(
			"GrowthGridLowerLeftPtY").toDouble();

	chemPara.growthMorCenterXCoord = globalConfigVars.getConfigValue(
			"GrowthMorCenterXCoord").toDouble();
	chemPara.growthMorCenterYCoord = globalConfigVars.getConfigValue(
			"GrowthMorCenterYCoord").toDouble();
	chemPara.growthMorHighConcen = globalConfigVars.getConfigValue(
			"GrowthMorHighConcen").toDouble();
	chemPara.growthMorLowConcen = globalConfigVars.getConfigValue(
			"GrowthMorLowConcen").toDouble();
	chemPara.growthMorDiffSlope = globalConfigVars.getConfigValue(
			"GrowthMorDiffSlope").toDouble();

	chemPara.growthMorCenterXCoordMX = globalConfigVars.getConfigValue(
			"GrowthMorCenterXCoordMX").toDouble();
	chemPara.growthMorCenterYCoordMX = globalConfigVars.getConfigValue(
			"GrowthMorCenterYCoordMX").toDouble();
	chemPara.growthMorHighConcenMX = globalConfigVars.getConfigValue(
			"GrowthMorHighConcenMX").toDouble();
	chemPara.growthMorLowConcenMX = globalConfigVars.getConfigValue(
			"GrowthMorLowConcenMX").toDouble();
	chemPara.growthMorDiffSlopeMX = globalConfigVars.getConfigValue(
			"GrowthMorDiffSlopeMX").toDouble();
}

void SimulationDomainGPU::readAllParameters() {
	readMemPara();
	readDomainPara();
	readChemPara();
}

void SimulationDomainGPU::initializeGrowthMap() {
	growthMap = GrowthDistriMap(chemPara.growthGridXDim,
			chemPara.growthGridYDim, chemPara.growthGridSpacing);
	growthMap.initialize(chemPara.growthGridLowerLeftPtX,
			chemPara.growthGridLowerLeftPtY, chemPara.growthMorCenterXCoord,
			chemPara.growthMorCenterYCoord, chemPara.growthMorHighConcen,
			chemPara.growthMorLowConcen, chemPara.growthMorDiffSlope);

	//cout << "after created growthMap1" << endl;
	growthMap2 = GrowthDistriMap(chemPara.growthGridXDim,
			chemPara.growthGridYDim, chemPara.growthGridSpacing);
	growthMap2.initialize(chemPara.growthGridLowerLeftPtX,
			chemPara.growthGridLowerLeftPtY, chemPara.growthMorCenterXCoordMX,
			chemPara.growthMorCenterYCoordMX, chemPara.growthMorHighConcenMX,
			chemPara.growthMorLowConcenMX, chemPara.growthMorDiffSlopeMX);
	//cout << "after created growthMap2" << endl;
}

std::vector<CVector> SimulationDomainGPU::stablizeCellCenters(
		SimulationInitData& initData) {

	std::vector<CVector> result;

	stabPara.outputFrameCount = globalConfigVars.getConfigValue(
			"StabFrameCount").toInt();
	stabPara.totalIterCount = globalConfigVars.getConfigValue(
			"StabTotalIterCount").toInt();
	stabPara.bdrySpacingRatio = globalConfigVars.getConfigValue(
			"StabBdrySpacingRatio").toDouble();
	stabPara.dt = globalConfigVars.getConfigValue("StabDt").toDouble();
	stabPara.outputAniName =
			globalConfigVars.getConfigValue("StabAniName").toString();

	initialize(initData);

	int aniAuxPara = (double) (stabPara.totalIterCount)
			/ stabPara.outputFrameCount;

	AnimationCriteria aniCri;
	aniCri.defaultEffectiveDistance = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();
	aniCri.isStressMap = false;

	for (int i = 0; i < stabPara.totalIterCount; i++) {
		//std::cout << "in stablizing, before run all logics" << std::endl;
		runAllLogic(stabPara.dt);
		if (i % aniAuxPara == 0) {
			outputVtkFilesWithColor(stabPara.outputAniName, i, aniCri);
		}
	}

	result = cells.getAllCellCenters();

	cout << "finished stablizeCellCenters" << endl;
	cout.flush();
	return result;
}

void SimulationDomainGPU::outputVtkFilesWithColor(std::string scriptNameBase,
		int rank, AnimationCriteria aniCri) {
	nodes.prepareSceForceComputation();
	VtkAnimationData aniData = nodes.obtainAnimationData(aniCri);
	aniData.outputVtkAni(scriptNameBase, rank);
}

void SimulationDomainGPU::checkIfAllDataFieldsValid() {
	cout << "Begin output information about nodes:" << endl;
	cout << "size of isActive:" << nodes.getInfoVecs().nodeIsActive.size()
			<< endl;
	cout << "size of nodeLocX:" << nodes.getInfoVecs().nodeLocX.size() << endl;
	cout << "size of nodeLocY:" << nodes.getInfoVecs().nodeLocY.size() << endl;
	cout << "size of nodeLocZ:" << nodes.getInfoVecs().nodeLocZ.size() << endl;
	cout << "size of nodeVelX:" << nodes.getInfoVecs().nodeVelX.size() << endl;
	cout << "size of nodeVelY:" << nodes.getInfoVecs().nodeVelY.size() << endl;
	cout << "size of nodeVelZ:" << nodes.getInfoVecs().nodeVelZ.size() << endl;
	cout << "size of CellType:" << nodes.getInfoVecs().nodeCellType.size()
			<< endl;
	cout << "size of nodeCellRank:" << nodes.getInfoVecs().nodeCellRank.size()
			<< endl;

	cout << "start position of Profile is "
			<< nodes.getAllocPara().startPosProfile << endl;
	cout << "start position of ECM is " << nodes.getAllocPara().startPosECM
			<< endl;
	cout << "start position of Cells is " << nodes.getAllocPara().startPosCells
			<< endl;

	cout << "max node of one cell is " << nodes.getAllocPara().maxNodeOfOneCell
			<< endl;
	cout << "max number of cells is " << nodes.getAllocPara().maxCellCount
			<< endl;
	cout << "max total cell node count is "
			<< nodes.getAllocPara().maxTotalCellNodeCount << endl;
	cout << "current active cell count is "
			<< nodes.getAllocPara().currentActiveCellCount << endl;

	cout << "max node of one ECM is " << nodes.getAllocPara().maxNodePerECM
			<< endl;
	cout << "max number of ECm is " << nodes.getAllocPara().maxECMCount << endl;
	cout << "max total ECM node count is "
			<< nodes.getAllocPara().maxTotalECMNodeCount << endl;
	cout << "current active ECM count is "
			<< nodes.getAllocPara().currentActiveECM << endl;

	cout << "max profile node count is "
			<< nodes.getAllocPara().maxProfileNodeCount << endl;
	cout << "current active profile node count is "
			<< nodes.getAllocPara().currentActiveProfileNodeCount << endl;
//int jj;
//cin >> jj;
}


void SimulationDomainGPU::outputLabelMatrix(std::string resultNameBase,
		int rank, PixelizePara& pixelPara) {
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string resultNameRank = ss.str();
	std::string matrixFileName = resultNameBase + resultNameRank + ".dat";
	vector<vector<int> > matrix = nodes.obtainLabelMatrix(pixelPara);
	printMatrixToFile(matrix, matrixFileName);
}
