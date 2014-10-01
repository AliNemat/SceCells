#include "SceNodes.h"

__constant__ double sceInterPara[5];
__constant__ double sceIntraPara[4];
__constant__ double sceInterDiffPara[5];
__constant__ double sceProfilePara[7];
__constant__ double sceECMPara[5];
__constant__ double sceDiffPara[5];
__constant__ uint ProfilebeginPos;
__constant__ uint ECMbeginPos;
__constant__ uint cellNodeBeginPos;
__constant__ uint nodeCountPerECM;
__constant__ uint nodeCountPerCell;

// This template method expands an input sequence by
// replicating each element a variable number of times. For example,
//
//   expand([2,2,2],[A,B,C]) -> [A,A,B,B,C,C]
//   expand([3,0,1],[A,B,C]) -> [A,A,A,C]
//   expand([1,3,2],[A,B,C]) -> [A,B,B,B,C,C]
//
// The element counts are assumed to be non-negative integers
template<typename InputIterator1, typename InputIterator2,
		typename OutputIterator>
OutputIterator expand(InputIterator1 first1, InputIterator1 last1,
		InputIterator2 first2, OutputIterator output) {
	typedef typename thrust::iterator_difference<InputIterator1>::type difference_type;

	difference_type input_size = thrust::distance(first1, last1);
	difference_type output_size = thrust::reduce(first1, last1);

	// scan the counts to obtain output offsets for each input element
	thrust::device_vector<difference_type> output_offsets(input_size, 0);
	thrust::exclusive_scan(first1, last1, output_offsets.begin());

	// scatter the nonzero counts into their corresponding output positions
	thrust::device_vector<difference_type> output_indices(output_size, 0);
	thrust::scatter_if(thrust::counting_iterator<difference_type>(0),
			thrust::counting_iterator<difference_type>(input_size),
			output_offsets.begin(), first1, output_indices.begin());

	// compute max-scan over the output indices, filling in the holes
	thrust::inclusive_scan(output_indices.begin(), output_indices.end(),
			output_indices.begin(), thrust::maximum<difference_type>());

	// gather input values according to index array (output = first2[output_indices])
	OutputIterator output_end = output;
	thrust::advance(output_end, output_size);
	thrust::gather(output_indices.begin(), output_indices.end(), first2,
			output);

	// return output + output_size
	thrust::advance(output, output_size);
	return output;
}

SceNodes::SceNodes() {
	readDomainPara();
}

void SceNodes::readDomainPara() {
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

void SceNodes::readMechPara() {
	static const double U0 =
			globalConfigVars.getConfigValue("InterCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_U0_DivFactor").toDouble();
	static const double V0 =
			globalConfigVars.getConfigValue("InterCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_V0_DivFactor").toDouble();
	static const double k1 =
			globalConfigVars.getConfigValue("InterCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_k1_DivFactor").toDouble();
	static const double k2 =
			globalConfigVars.getConfigValue("InterCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue("InterCell_k2_DivFactor").toDouble();
	static const double interLinkEffectiveRange =
			globalConfigVars.getConfigValue("InterCellLinkBreakRange").toDouble();

	mechPara.sceInterParaCPU[0] = U0;
	mechPara.sceInterParaCPU[1] = V0;
	mechPara.sceInterParaCPU[2] = k1;
	mechPara.sceInterParaCPU[3] = k2;
	mechPara.sceInterParaCPU[4] = interLinkEffectiveRange;

	static const double U0_Intra =
			globalConfigVars.getConfigValue("IntraCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_U0_DivFactor").toDouble();
	static const double V0_Intra =
			globalConfigVars.getConfigValue("IntraCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_V0_DivFactor").toDouble();
	static const double k1_Intra =
			globalConfigVars.getConfigValue("IntraCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_k1_DivFactor").toDouble();
	static const double k2_Intra =
			globalConfigVars.getConfigValue("IntraCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue("IntraCell_k2_DivFactor").toDouble();

	mechPara.sceIntraParaCPU[0] = U0_Intra;
	mechPara.sceIntraParaCPU[1] = V0_Intra;
	mechPara.sceIntraParaCPU[2] = k1_Intra;
	mechPara.sceIntraParaCPU[3] = k2_Intra;

	static const double U0_Diff =
			globalConfigVars.getConfigValue("InterCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Diff_U0_DivFactor").toDouble();
	static const double V0_Diff =
			globalConfigVars.getConfigValue("InterCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Diff_V0_DivFactor").toDouble();
	static const double k1_Diff =
			globalConfigVars.getConfigValue("InterCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Diff_k1_DivFactor").toDouble();
	static const double k2_Diff =
			globalConfigVars.getConfigValue("InterCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Diff_k2_DivFactor").toDouble();

	mechPara.sceInterDiffParaCPU[0] = U0_Diff;
	mechPara.sceInterDiffParaCPU[1] = V0_Diff;
	mechPara.sceInterDiffParaCPU[2] = k1_Diff;
	mechPara.sceInterDiffParaCPU[3] = k2_Diff;
	mechPara.sceInterDiffParaCPU[4] = interLinkEffectiveRange;

	static const double U0_Bdry =
			globalConfigVars.getConfigValue("InterCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Bdry_U0_DivFactor").toDouble();
	static const double V0_Bdry =
			globalConfigVars.getConfigValue("InterCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Bdry_V0_DivFactor").toDouble();
	static const double k1_Bdry =
			globalConfigVars.getConfigValue("InterCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Bdry_k1_DivFactor").toDouble();
	static const double k2_Bdry =
			globalConfigVars.getConfigValue("InterCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Bdry_k2_DivFactor").toDouble();
	// 1.8 comes from standard
	static const double neutralLength =
			globalConfigVars.getConfigValue("Bdry_base_neutral_dist").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_Bdry_k2_DivFactor").toDouble()
					* globalConfigVars.getConfigValue("baseline_k_value").toDouble();

	static const double linearParameter = globalConfigVars.getConfigValue(
			"Profile_linear_parameter").toDouble();

	mechPara.sceProfileParaCPU[0] = U0_Bdry;
	mechPara.sceProfileParaCPU[1] = V0_Bdry;
	mechPara.sceProfileParaCPU[2] = k1_Bdry;
	mechPara.sceProfileParaCPU[3] = k2_Bdry;
	mechPara.sceProfileParaCPU[4] = interLinkEffectiveRange;
	mechPara.sceProfileParaCPU[5] = linearParameter;
	mechPara.sceProfileParaCPU[6] = neutralLength;

	static const double U0_ECM =
			globalConfigVars.getConfigValue("InterCell_U0_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_ECM_U0_DivFactor").toDouble();
	static const double V0_ECM =
			globalConfigVars.getConfigValue("InterCell_V0_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_ECM_V0_DivFactor").toDouble();
	static const double k1_ECM =
			globalConfigVars.getConfigValue("InterCell_k1_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_ECM_k1_DivFactor").toDouble();
	static const double k2_ECM =
			globalConfigVars.getConfigValue("InterCell_k2_Original").toDouble()
					/ globalConfigVars.getConfigValue(
							"InterCell_ECM_k2_DivFactor").toDouble();
	mechPara.sceECMParaCPU[0] = U0_ECM;
	mechPara.sceECMParaCPU[1] = V0_ECM;
	mechPara.sceECMParaCPU[2] = k1_ECM;
	mechPara.sceECMParaCPU[3] = k2_ECM;
	mechPara.sceECMParaCPU[4] = interLinkEffectiveRange;
}

SceNodes::SceNodes(uint totalBdryNodeCount, uint maxProfileNodeCount,
		uint maxTotalECMCount, uint maxNodeInECM, uint maxTotalCellCount,
		uint maxNodeInCell) {
	readDomainPara();
	initNodeAllocPara(totalBdryNodeCount, maxProfileNodeCount, maxTotalECMCount,
			maxNodeInECM, maxTotalCellCount, maxNodeInCell);
	uint maxTotalNodeCount = totalBdryNodeCount + maxProfileNodeCount
			+ allocPara.maxTotalECMNodeCount + allocPara.maxTotalCellNodeCount;
	allocSpaceForNodes(maxTotalNodeCount);
	thrust::host_vector<SceNodeType> hostTmpVector(maxTotalNodeCount);
	thrust::host_vector<bool> hostTmpVector2(maxTotalNodeCount);
	thrust::host_vector<int> hostTmpVector3(maxTotalNodeCount);
	for (int i = 0; i < maxTotalNodeCount; i++) {
		if (i < allocPara.startPosProfile) {
			hostTmpVector[i] = Boundary;
			hostTmpVector3[i] = 0;
		} else if (i < allocPara.startPosECM) {
			hostTmpVector[i] = Profile;
			hostTmpVector3[i] = 0;
		} else if (i < allocPara.startPosCells) {
			hostTmpVector[i] = ECM;
			hostTmpVector3[i] = (i - allocPara.startPosECM)
					/ allocPara.maxNodePerECM;
		} else {
			// all initialized as FNM
			hostTmpVector[i] = FNM;
			hostTmpVector3[i] = (i - allocPara.startPosCells)
					/ allocPara.maxNodeOfOneCell;
		}
		hostTmpVector2[i] = false;
	}
	infoVecs.nodeCellType = hostTmpVector;
	infoVecs.nodeIsActive = hostTmpVector2;
	infoVecs.nodeCellRank = hostTmpVector3;
	copyParaToGPUConstMem();
}

SceNodes::SceNodes(uint totalBdryNodeCount, uint maxProfileNodeCount,
		uint maxCartNodeCount, uint maxTotalECMCount, uint maxNodeInECM,
		uint maxTotalCellCount, uint maxNodeInCell) {
	readDomainPara();
	initNodeAllocPara_v2(totalBdryNodeCount, maxProfileNodeCount,
			maxCartNodeCount, maxTotalECMCount, maxNodeInECM, maxTotalCellCount,
			maxNodeInCell);
	uint maxTotalNodeCount = totalBdryNodeCount + maxProfileNodeCount
			+ maxCartNodeCount + allocPara.maxTotalECMNodeCount
			+ allocPara.maxTotalCellNodeCount;
	allocSpaceForNodes(maxTotalNodeCount);
	thrust::host_vector<SceNodeType> hostTmpVector(maxTotalNodeCount);
	thrust::host_vector<bool> hostTmpVector2(maxTotalNodeCount);
	thrust::host_vector<int> hostTmpVector3(maxTotalNodeCount);
	for (int i = 0; i < maxTotalNodeCount; i++) {
		if (i < allocPara.startPosProfile) {
			hostTmpVector[i] = Boundary;
			hostTmpVector3[i] = 0;
		} else if (i < allocPara.startPosCart) {
			hostTmpVector[i] = Profile;
			hostTmpVector3[i] = 0;
		} else if (i < allocPara.startPosECM) {
			hostTmpVector[i] = Cart;
			hostTmpVector3[i] = 0;
		} else if (i < allocPara.startPosCells) {
			hostTmpVector[i] = ECM;
			hostTmpVector3[i] = (i - allocPara.startPosECM)
					/ allocPara.maxNodePerECM;
		} else {
			// all initialized as FNM
			hostTmpVector[i] = FNM;
			hostTmpVector3[i] = (i - allocPara.startPosCells)
					/ allocPara.maxNodeOfOneCell;
		}
		hostTmpVector2[i] = false;
	}
	infoVecs.nodeCellType = hostTmpVector;
	infoVecs.nodeIsActive = hostTmpVector2;
	infoVecs.nodeCellRank = hostTmpVector3;
	copyParaToGPUConstMem();
}

void SceNodes::copyParaToGPUConstMem() {

	readMechPara();

	cudaMemcpyToSymbol(sceInterPara, mechPara.sceInterParaCPU,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceIntraPara, mechPara.sceIntraParaCPU,
			4 * sizeof(double));
	cudaMemcpyToSymbol(ProfilebeginPos, &allocPara.startPosProfile,
			sizeof(uint));
	cudaMemcpyToSymbol(ECMbeginPos, &allocPara.startPosECM, sizeof(uint));
	cudaMemcpyToSymbol(cellNodeBeginPos, &allocPara.startPosCells,
			sizeof(uint));
	cudaMemcpyToSymbol(nodeCountPerECM, &allocPara.maxNodePerECM, sizeof(uint));
	cudaMemcpyToSymbol(nodeCountPerCell, &allocPara.maxNodeOfOneCell,
			sizeof(uint));
	cudaMemcpyToSymbol(sceProfilePara, mechPara.sceProfileParaCPU,
			7 * sizeof(double));
	cudaMemcpyToSymbol(sceInterDiffPara, mechPara.sceInterDiffParaCPU,
			5 * sizeof(double));
	cudaMemcpyToSymbol(sceECMPara, mechPara.sceECMParaCPU, 5 * sizeof(double));
	//std::cout << "finished SceNodes:" << std::endl;
}

void SceNodes::initDimension(double domainMinX, double domainMaxX,
		double domainMinY, double domainMaxY, double domainBucketSize) {
	domainPara.minX = domainMinX;
	domainPara.maxX = domainMaxX;
	domainPara.minY = domainMinY;
	domainPara.maxY = domainMaxY;
	domainPara.gridSpacing = domainBucketSize;
	domainPara.numOfBucketsInXDim = (domainPara.maxX - domainPara.minX)
			/ domainPara.gridSpacing + 1;
	domainPara.numOfBucketsInYDim = (domainPara.maxY - domainPara.minY)
			/ domainPara.gridSpacing + 1;
	domainPara.totalBucketCount = domainPara.numOfBucketsInXDim
			* domainPara.numOfBucketsInYDim;

	auxVecs.keyBegin.resize(domainPara.totalBucketCount);
	auxVecs.keyEnd.resize(domainPara.totalBucketCount);
}

std::vector<std::pair<uint, uint> > SceNodes::obtainPossibleNeighborPairs() {

	std::vector<std::pair<uint, uint> > result;
	thrust::host_vector<uint> keyBeginCPU = auxVecs.keyBegin;
	thrust::host_vector<uint> keyEndCPU = auxVecs.keyEnd;
	thrust::host_vector<uint> bucketKeysCPU = auxVecs.bucketKeys;
	thrust::host_vector<uint> bucketValuesCPU = auxVecs.bucketValues;
	thrust::host_vector<uint> bucketValuesExtendedCPU =
			auxVecs.bucketValuesIncludingNeighbor;

	//cout << "sizes = " << keyBeginCPU.size() << " " << keyEndCPU.size() << " "
	//		<< bucketKeysCPU.size() << " " << bucketValuesCPU.size() << " "
	//		<< bucketValuesExtendedCPU.size() << endl;

	uint iterationCounter = 0;

	int size = bucketKeysCPU.size();
	for (int i = 0; i < size; i++) {
		for (int j = keyBeginCPU[bucketKeysCPU[i]];
				j < keyEndCPU[bucketKeysCPU[i]]; j++) {
			//std::cout << "pair node 1: " << bucketValues[i] << ",pair node2: "
			//		<< bucketValuesIncludingNeighbor[j] << std::endl;
			int node1 = bucketValuesCPU[i];
			int node2 = bucketValuesExtendedCPU[j];
			//cout << "node1 = " << node1 << "node2 = " << node2 << endl;
			if (node1 >= node2) {
				//cout << "not inserted, continue" << endl;
				continue;
			} else {
				//cout << "inserted!" << endl;
				result.push_back(std::make_pair<uint, uint>(node1, node2));
			}
			iterationCounter++;
		}
	}
	//cout << "iter counter = " << iterationCounter << endl;
	return result;
}

void SceNodes::initValues(std::vector<double>& initBdryCellNodePosX,
		std::vector<double>& initBdryCellNodePosY,
		std::vector<double>& initProfileNodePosX,
		std::vector<double>& initProfileNodePosY,
		std::vector<double>& initECMNodePosX,
		std::vector<double>& initECMNodePosY,
		std::vector<double>& initFNMCellNodePosX,
		std::vector<double>& initFNMCellNodePosY,
		std::vector<double>& initMXCellNodePosX,
		std::vector<double>& initMXCellNodePosY) {

	uint FNMNodeCountX = initFNMCellNodePosX.size();
	uint MXNodeCountX = initMXCellNodePosX.size();

	uint beginAddressOfProfile = allocPara.startPosProfile;
// find the begining position of ECM.
	uint beginAddressOfECM = allocPara.startPosECM;
// find the begining position of FNM cells.
	uint beginAddressOfFNM = allocPara.startPosCells;
// find the begining position of MX cells.
	uint beginAddressOfMX = beginAddressOfFNM + FNMNodeCountX;

//std::cerr << "before copying arrays" << endl;

	thrust::copy(initBdryCellNodePosX.begin(), initBdryCellNodePosX.end(),
			infoVecs.nodeLocX.begin());
	thrust::copy(initBdryCellNodePosY.begin(), initBdryCellNodePosY.end(),
			infoVecs.nodeLocY.begin());

//std::cerr << "copy 1" << endl;

// copy x and y position of nodes of Profile to actual node position.
	thrust::copy(initProfileNodePosX.begin(), initProfileNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfProfile);
	thrust::copy(initProfileNodePosY.begin(), initProfileNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfProfile);

//std::cerr << "copy 2" << endl;

// copy x and y position of nodes of ECM to actual node position.
	thrust::copy(initECMNodePosX.begin(), initECMNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfECM);
	thrust::copy(initECMNodePosY.begin(), initECMNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfECM);

// debug
	for (int i = 0; i < initECMNodePosX.size(); i++) {
		std::cout << "i + beginAddressOfECM = " << (i + beginAddressOfECM)
				<< "nodeLocX =" << infoVecs.nodeLocX[i + beginAddressOfECM]
				<< std::endl;
		assert(infoVecs.nodeLocX[i + beginAddressOfECM] == initECMNodePosX[i]);
		assert(!isnan(initECMNodePosX[i]));
	}

// std::cerr << "copy 3" << endl;

// copy x and y position of nodes of FNM cells to actual node position.
	thrust::copy(initFNMCellNodePosX.begin(), initFNMCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfFNM);
	thrust::copy(initFNMCellNodePosY.begin(), initFNMCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfFNM);

// std::cerr << "copy 4" << endl;

	thrust::fill(infoVecs.nodeCellType.begin() + beginAddressOfFNM,
			infoVecs.nodeCellType.begin() + beginAddressOfMX, FNM);

// copy x and y position of nodes of MX cells to actual node position.
	thrust::copy(initMXCellNodePosX.begin(), initMXCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfMX);
	thrust::copy(initMXCellNodePosY.begin(), initMXCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfMX);

//std::cerr << "after copying arrays" << endl;

	thrust::fill(infoVecs.nodeCellType.begin() + beginAddressOfMX,
			infoVecs.nodeCellType.begin() + beginAddressOfMX + MXNodeCountX,
			MX);

//std::cout << "initial MX cell numbers: " << mxQuotient << std::endl;
}

void SceNodes::initValues_v2(std::vector<CVector>& initBdryCellNodePos,
		std::vector<CVector>& initProfileNodePos,
		std::vector<CVector>& initCartNodePos,
		std::vector<CVector>& initECMNodePos,
		std::vector<CVector>& initFNMCellNodePos,
		std::vector<CVector>& initMXCellNodePos) {

	uint FNMNodeCount = initFNMCellNodePos.size();
	uint MXNodeCount = initMXCellNodePos.size();

	uint beginAddressOfProfile = allocPara.startPosProfile;
	uint beginAddressOfCart = allocPara.startPosCart;
	// find the begining position of ECM.
	uint beginAddressOfECM = allocPara.startPosECM;
	// find the begining position of FNM cells.
	uint beginAddressOfFNM = allocPara.startPosCells;
	// find the begining position of MX cells.
	uint beginAddressOfMX = beginAddressOfFNM + FNMNodeCount;

	//std::cerr << "before copying arrays" << endl;
	std::vector<double> initBdryCellNodePosX = getArrayXComp(
			initBdryCellNodePos);
	thrust::copy(initBdryCellNodePosX.begin(), initBdryCellNodePosX.end(),
			infoVecs.nodeLocX.begin());
	std::vector<double> initBdryCellNodePosY = getArrayYComp(
			initBdryCellNodePos);
	thrust::copy(initBdryCellNodePosY.begin(), initBdryCellNodePosY.end(),
			infoVecs.nodeLocY.begin());

	// copy x and y position of nodes of Profile to actual node position.
	std::vector<double> initProfileNodePosX = getArrayXComp(initProfileNodePos);
	thrust::copy(initProfileNodePosX.begin(), initProfileNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfProfile);
	std::vector<double> initProfileNodePosY = getArrayYComp(initProfileNodePos);
	thrust::copy(initProfileNodePosY.begin(), initProfileNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfProfile);

	// copy x and y position of nodes of Profile to actual node position.
	std::vector<double> initCartNodePosX = getArrayXComp(initCartNodePos);
	thrust::copy(initCartNodePosX.begin(), initCartNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfCart);
	std::vector<double> initCartNodePosY = getArrayYComp(initCartNodePos);
	thrust::copy(initCartNodePosY.begin(), initCartNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfCart);

	// copy x and y position of nodes of ECM to actual node position.
	std::vector<double> initECMNodePosX = getArrayXComp(initECMNodePos);
	thrust::copy(initECMNodePosX.begin(), initECMNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfECM);
	std::vector<double> initECMNodePosY = getArrayYComp(initECMNodePos);
	thrust::copy(initECMNodePosY.begin(), initECMNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfECM);

	// debug
	for (int i = 0; i < initECMNodePosX.size(); i++) {
		std::cout << "i + beginAddressOfECM = " << (i + beginAddressOfECM)
				<< "nodeLocX =" << infoVecs.nodeLocX[i + beginAddressOfECM]
				<< std::endl;
		assert(infoVecs.nodeLocX[i + beginAddressOfECM] == initECMNodePosX[i]);
		assert(!isnan(initECMNodePosX[i]));
	}

	// std::cerr << "copy 3" << endl;

	// copy x and y position of nodes of FNM cells to actual node position.
	std::vector<double> initFNMCellNodePosX = getArrayXComp(initFNMCellNodePos);
	thrust::copy(initFNMCellNodePosX.begin(), initFNMCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfFNM);
	std::vector<double> initFNMCellNodePosY = getArrayYComp(initFNMCellNodePos);
	thrust::copy(initFNMCellNodePosY.begin(), initFNMCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfFNM);

	// std::cerr << "copy 4" << endl;

	thrust::fill(infoVecs.nodeCellType.begin() + beginAddressOfFNM,
			infoVecs.nodeCellType.begin() + beginAddressOfMX, FNM);

	// copy x and y position of nodes of MX cells to actual node position.
	std::vector<double> initMXCellNodePosX = getArrayXComp(initMXCellNodePos);
	thrust::copy(initMXCellNodePosX.begin(), initMXCellNodePosX.end(),
			infoVecs.nodeLocX.begin() + beginAddressOfMX);
	std::vector<double> initMXCellNodePosY = getArrayYComp(initMXCellNodePos);
	thrust::copy(initMXCellNodePosY.begin(), initMXCellNodePosY.end(),
			infoVecs.nodeLocY.begin() + beginAddressOfMX);

	//std::cerr << "after copying arrays" << endl;
	thrust::fill(infoVecs.nodeCellType.begin() + beginAddressOfMX,
			infoVecs.nodeCellType.begin() + beginAddressOfMX + MXNodeCount, MX);
}

void SceNodes::applyProfileForces() {
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(
			allocPara.currentActiveProfileNodeCount);

	double* nodeLocXAddressEpiBegin = thrust::raw_pointer_cast(
			&infoVecs.nodeLocX[allocPara.startPosProfile]);
	double* nodeLocYAddressEpiBegin = thrust::raw_pointer_cast(
			&infoVecs.nodeLocY[allocPara.startPosProfile]);
	double* nodeLocZAddressEpiBegin = thrust::raw_pointer_cast(
			&infoVecs.nodeLocZ[allocPara.startPosProfile]);

	double* nodeVelXAddressEpiBegin = thrust::raw_pointer_cast(
			&infoVecs.nodeVelX[allocPara.startPosProfile]);
	double* nodeVelYAddressEpiBegin = thrust::raw_pointer_cast(
			&infoVecs.nodeVelY[allocPara.startPosProfile]);
	double* nodeVelZAddressEpiBegin = thrust::raw_pointer_cast(
			&infoVecs.nodeVelZ[allocPara.startPosProfile]);

	thrust::transform(countingIterBegin, countingIterEnd,
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeVelX.begin(),
							infoVecs.nodeVelY.begin(),
							infoVecs.nodeVelZ.begin()))
					+ allocPara.startPosProfile,
			AddLinkForces(nodeLocXAddressEpiBegin, nodeLocYAddressEpiBegin,
					nodeLocZAddressEpiBegin, nodeVelXAddressEpiBegin,
					nodeVelYAddressEpiBegin, nodeVelZAddressEpiBegin,
					allocPara.currentActiveProfileNodeCount));
}

VtkAnimationData SceNodes::obtainAnimationData(AnimationCriteria aniCri) {
	VtkAnimationData vtkData;
	std::vector<std::pair<uint, uint> > pairs = obtainPossibleNeighborPairs();
	cout << "size of potential pairs = " << pairs.size() << endl;
	std::vector<std::pair<uint, uint> > pairsTobeAnimated;

	// unordered_map is more efficient than map, but it is a c++ 11 feature
	// and c++ 11 seems to be incompatible with Thrust.
	IndexMap locIndexToAniIndexMap;

	// Doesn't have to copy the entire nodeLocX array.
	// Only copy the first half will be sufficient
	thrust::host_vector<double> hostTmpVectorLocX = infoVecs.nodeLocX;
	thrust::host_vector<double> hostTmpVectorLocY = infoVecs.nodeLocY;
	thrust::host_vector<double> hostTmpVectorLocZ = infoVecs.nodeLocZ;
	thrust::host_vector<SceNodeType> hostTmpVectorNodeType =
			infoVecs.nodeCellType;
	thrust::host_vector<uint> hostTmpVectorNodeRank = infoVecs.nodeCellRank;
	thrust::host_vector<double> hostTmpVectorNodeStress;
	if (aniCri.isStressMap) {
		hostTmpVectorNodeStress = infoVecs.nodeMaxForce;
	}

	uint curIndex = 0;
	for (uint i = 0; i < pairs.size(); i++) {
		uint node1Index = pairs[i].first;
		uint node2Index = pairs[i].second;
		double node1X = hostTmpVectorLocX[node1Index];
		double node1Y = hostTmpVectorLocY[node1Index];
		double node1Z = hostTmpVectorLocZ[node1Index];
		SceNodeType node1T = hostTmpVectorNodeType[node1Index];
		uint node1R = hostTmpVectorNodeRank[node1Index];
		double node2X = hostTmpVectorLocX[node2Index];
		double node2Y = hostTmpVectorLocY[node2Index];
		double node2Z = hostTmpVectorLocZ[node2Index];
		SceNodeType node2T = hostTmpVectorNodeType[node2Index];
		uint node2R = hostTmpVectorNodeRank[node2Index];

		if (aniCri.isPairQualify(node1Index, node2Index, node1X, node1Y, node1Z,
				node1T, node1R, node2X, node2Y, node2Z, node2T, node2R)) {
			pairsTobeAnimated.push_back(pairs[i]);
			IndexMap::iterator it = locIndexToAniIndexMap.find(pairs[i].first);
			if (it == locIndexToAniIndexMap.end()) {
				locIndexToAniIndexMap.insert(
						std::pair<uint, uint>(pairs[i].first, curIndex));
				curIndex++;
				PointAniData ptAniData;
				if (aniCri.isStressMap) {
					ptAniData.colorScale = hostTmpVectorNodeStress[node1Index];
				} else {
					ptAniData.colorScale = nodeTypeToScale(node1T);
				}
				ptAniData.pos = CVector(node1X, node1Y, node1Z);
				vtkData.pointsAniData.push_back(ptAniData);
			}
			it = locIndexToAniIndexMap.find(pairs[i].second);
			if (it == locIndexToAniIndexMap.end()) {
				locIndexToAniIndexMap.insert(
						std::pair<uint, uint>(pairs[i].second, curIndex));
				curIndex++;
				PointAniData ptAniData;
				if (aniCri.isStressMap) {
					ptAniData.colorScale = hostTmpVectorNodeStress[node2Index];
				} else {
					ptAniData.colorScale = nodeTypeToScale(node2T);
				}
				ptAniData.pos = CVector(node2X, node2Y, node2Z);
				vtkData.pointsAniData.push_back(ptAniData);
			}

			it = locIndexToAniIndexMap.find(pairs[i].first);
			uint aniIndex1 = it->second;
			it = locIndexToAniIndexMap.find(pairs[i].second);
			uint aniIndex2 = it->second;

			LinkAniData linkData;
			linkData.node1Index = aniIndex1;
			linkData.node2Index = aniIndex2;
			vtkData.linksAniData.push_back(linkData);
		}
	}

	uint profileStartIndex = allocPara.startPosProfile;
	uint profileEndIndex = profileStartIndex
			+ allocPara.currentActiveProfileNodeCount;
	cout << "start index is : " << profileStartIndex << endl;
	cout << "end index is : " << profileEndIndex << endl;
	for (uint i = profileStartIndex; i < profileEndIndex; i++) {
		PointAniData ptAniData;
		ptAniData.pos = CVector(hostTmpVectorLocX[i], hostTmpVectorLocY[i],
				hostTmpVectorLocZ[i]);
		if (aniCri.isStressMap) {
			ptAniData.colorScale = hostTmpVectorNodeStress[i];
		} else {
			ptAniData.colorScale = nodeTypeToScale(hostTmpVectorNodeType[i]);
		}
		vtkData.pointsAniData.push_back(ptAniData);
		LinkAniData linkData;
		linkData.node1Index = curIndex;
		linkData.node2Index = curIndex + 1;
		if (i != profileEndIndex - 1) {
			vtkData.linksAniData.push_back(linkData);
		}
		curIndex++;
	}

	return vtkData;
}

void SceNodes::findBucketBounds() {
	thrust::counting_iterator<unsigned int> search_begin(0);
	thrust::lower_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(), search_begin,
			search_begin + domainPara.totalBucketCount,
			auxVecs.keyBegin.begin());
	thrust::upper_bound(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(), search_begin,
			search_begin + domainPara.totalBucketCount, auxVecs.keyEnd.begin());
}

void SceNodes::prepareSceForceComputation() {
	buildBuckets2D();
	extendBuckets2D();
	findBucketBounds();
}

void SceNodes::addNewlyDividedCells(
		thrust::device_vector<double> &nodeLocXNewCell,
		thrust::device_vector<double> &nodeLocYNewCell,
		thrust::device_vector<double> &nodeLocZNewCell,
		thrust::device_vector<bool> &nodeIsActiveNewCell,
		thrust::device_vector<SceNodeType> &nodeCellTypeNewCell) {

// data validation
	uint nodesSize = nodeLocXNewCell.size();
	assert(nodesSize % allocPara.maxNodeOfOneCell == 0);
	uint addCellCount = nodesSize / allocPara.maxNodeOfOneCell;

// position that we will add newly divided cells.
	uint shiftStartPosNewCell = allocPara.startPosCells
			+ allocPara.currentActiveCellCount * allocPara.maxNodeOfOneCell;

	thrust::copy(
			thrust::make_zip_iterator(
					thrust::make_tuple(nodeLocXNewCell.begin(),
							nodeLocYNewCell.begin(), nodeLocZNewCell.begin(),
							nodeIsActiveNewCell.begin(),
							nodeCellTypeNewCell.begin())),
			thrust::make_zip_iterator(
					thrust::make_tuple(nodeLocXNewCell.end(),
							nodeLocYNewCell.end(), nodeLocZNewCell.end(),
							nodeIsActiveNewCell.end(),
							nodeCellTypeNewCell.end())),
			thrust::make_zip_iterator(
					thrust::make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(),
							infoVecs.nodeCellType.begin()))
					+ shiftStartPosNewCell);

// total number of cells has increased.
	allocPara.currentActiveCellCount = allocPara.currentActiveCellCount
			+ addCellCount;
}

void SceNodes::buildBuckets2D() {
	int totalActiveNodes = allocPara.startPosCells
			+ allocPara.currentActiveCellCount * allocPara.maxNodeOfOneCell;

	auxVecs.bucketKeys.resize(totalActiveNodes);
	auxVecs.bucketValues.resize(totalActiveNodes);
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(totalActiveNodes);

	// takes counting iterator and coordinates
	// return tuple of keys and values
	// transform the points to their bucket indices
	thrust::transform(
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), countingIterBegin)),
			make_zip_iterator(
					make_tuple(infoVecs.nodeLocX.begin(),
							infoVecs.nodeLocY.begin(),
							infoVecs.nodeLocZ.begin(),
							infoVecs.nodeIsActive.begin(), countingIterBegin))
					+ totalActiveNodes,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			pointToBucketIndex2D(domainPara.minX, domainPara.maxX,
					domainPara.minY, domainPara.maxY, domainPara.gridSpacing));

	// sort the points by their bucket index
	thrust::sort_by_key(auxVecs.bucketKeys.begin(), auxVecs.bucketKeys.end(),
			auxVecs.bucketValues.begin());
	// for those nodes that are inactive, key value of UINT_MAX will be returned.
	// we need to removed those keys along with their values.
	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeys.begin(),
			auxVecs.bucketKeys.end(), UINT_MAX);

	auxVecs.bucketKeys.erase(auxVecs.bucketKeys.end() - numberOfOutOfRange,
			auxVecs.bucketKeys.end());
	auxVecs.bucketValues.erase(auxVecs.bucketValues.end() - numberOfOutOfRange,
			auxVecs.bucketValues.end());
}
__device__
double computeDist(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2)
					+ (zPos - zPos2) * (zPos - zPos2));
}

__device__
void calculateAndAddECMForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {

	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceECMPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceECMPara[0] / sceECMPara[2]
				* exp(-linkLength / sceECMPara[2])
				+ sceECMPara[1] / sceECMPara[3]
						* exp(-linkLength / sceECMPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.3;
		}
	}
	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}

}
__device__
void calculateAndAddProfileForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	forceValue = -sceProfilePara[5] * (linkLength - sceProfilePara[6]);

	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}
}
__device__
void calculateAndAddInterForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterPara[4]) {
		//forceValue = 0;
	} else {
		forceValue = -sceInterPara[0] / sceInterPara[2]
				* exp(-linkLength / sceInterPara[2])
				+ sceInterPara[1] / sceInterPara[3]
						* exp(-linkLength / sceInterPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.5;
		}
	}
	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}
}
__device__
void calculateAndAddDiffInterCellForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterDiffPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterDiffPara[0] / sceInterDiffPara[2]
				* exp(-linkLength / sceInterDiffPara[2])
				+ sceInterDiffPara[1] / sceInterDiffPara[3]
						* exp(-linkLength / sceInterDiffPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.2;
		}
	}
	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}
}
__device__
void calculateAndAddInterForceDiffType(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = 0;
	if (linkLength > sceInterPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterPara[0] / sceInterPara[2]
				* exp(-linkLength / sceInterPara[2])
				+ sceInterPara[1] / sceInterPara[3]
						* exp(-linkLength / sceInterPara[3]);
		if (forceValue > 0) {
			//forceValue = 0;
			forceValue = forceValue * 0.3;
		}
	}
	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}
}
__device__
void calculateAndAddIntraForce(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes) {
	double linkLength = computeDist(xPos, yPos, zPos, xPos2, yPos2, zPos2);
	double forceValue = -sceIntraPara[0] / sceIntraPara[2]
			* exp(-linkLength / sceIntraPara[2])
			+ sceIntraPara[1] / sceIntraPara[3]
					* exp(-linkLength / sceIntraPara[3]);
	if (linkLength > 1.0e-12) {
		xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
		yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
		zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
	}
}
__device__ bool bothNodesCellNode(uint nodeGlobalRank1, uint nodeGlobalRank2,
		uint cellNodesThreshold) {
	if (nodeGlobalRank1 < cellNodesThreshold
			&& nodeGlobalRank2 < cellNodesThreshold) {
		return true;
	} else {
		return false;
	}
}

__device__ bool isSameCell(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if (nodeGlobalRank1 < cellNodeBeginPos
			|| nodeGlobalRank2 < cellNodeBeginPos) {
		return false;
	}
	if ((nodeGlobalRank1 - cellNodeBeginPos) / nodeCountPerCell
			== (nodeGlobalRank2 - cellNodeBeginPos) / nodeCountPerCell) {
		return true;
	} else {
		return false;
	}
}

__device__ bool isSameECM(uint nodeGlobalRank1, uint nodeGlobalRank2) {
	if ((nodeGlobalRank1 - ECMbeginPos) / nodeCountPerECM
			== (nodeGlobalRank2 - ECMbeginPos) / nodeCountPerECM) {
		return true;
	} else {
		return false;
	}
}

__device__ bool isNeighborECMNodes(uint nodeGlobalRank1, uint nodeGlobalRank2) {
// this means that two nodes are from the same ECM
	if ((nodeGlobalRank1 - ECMbeginPos) / nodeCountPerECM
			== (nodeGlobalRank2 - ECMbeginPos) / nodeCountPerECM) {
		// this means that two nodes are actually close to each other
		// seems to be strange because of unsigned int.
		if ((nodeGlobalRank1 > nodeGlobalRank2
				&& nodeGlobalRank1 - nodeGlobalRank2 == 1)
				|| (nodeGlobalRank2 > nodeGlobalRank1
						&& nodeGlobalRank2 - nodeGlobalRank1 == 1)) {
			return true;
		}
	}
	return false;
}

__device__ bool isNeighborProfileNodes(uint nodeGlobalRank1,
		uint nodeGlobalRank2) {
	if ((nodeGlobalRank1 > nodeGlobalRank2
			&& nodeGlobalRank1 - nodeGlobalRank2 == 1)
			|| (nodeGlobalRank2 > nodeGlobalRank1
					&& nodeGlobalRank2 - nodeGlobalRank1 == 1)) {
		return true;
	}
	return false;
}

__device__ bool ofSameType(uint cellType1, uint cellType2) {
	if (cellType1 == cellType2) {
		return true;
	} else {
		return false;
	}
}

__device__ bool bothCellNodes(SceNodeType &type1, SceNodeType &type2) {
	if ((type1 == MX || type1 == FNM) && (type2 == MX || type2 == FNM)) {
		return true;
	} else {
		return false;
	}
}

__device__
void calculateForceBetweenLinkNodes(double &xLoc, double &yLoc, double &zLoc,
		double &xLocLeft, double &yLocLeft, double &zLocLeft, double &xLocRight,
		double &yLocRight, double &zLocRight, double &xVel, double &yVel,
		double &zVel) {
	double linkLengthLeft = computeDist(xLoc, yLoc, zLoc, xLocLeft, yLocLeft,
			zLocLeft);
	double forceValueLeft = sceProfilePara[5]
			* (linkLengthLeft - sceProfilePara[6]);
	xVel = xVel + forceValueLeft * (xLocLeft - xLoc) / linkLengthLeft;
	yVel = yVel + forceValueLeft * (yLocLeft - yLoc) / linkLengthLeft;
	zVel = zVel + forceValueLeft * (zLocLeft - zLoc) / linkLengthLeft;

	double linkLengthRight = computeDist(xLoc, yLoc, zLoc, xLocRight, yLocRight,
			zLocRight);
	double forceValueRight = sceProfilePara[5]
			* (linkLengthRight - sceProfilePara[6]);
	xVel = xVel + forceValueRight * (xLocRight - xLoc) / linkLengthRight;
	yVel = yVel + forceValueRight * (yLocRight - yLoc) / linkLengthRight;
	zVel = zVel + forceValueRight * (zLocRight - zLoc) / linkLengthRight;

}
__device__
void handleForceBetweenNodes(uint &nodeRank1, SceNodeType &type1,
		uint &nodeRank2, SceNodeType &type2, double &xPos, double &yPos,
		double &zPos, double &xPos2, double &yPos2, double &zPos2, double &xRes,
		double &yRes, double &zRes, double &maxForce, double* _nodeLocXAddress,
		double* _nodeLocYAddress, double* _nodeLocZAddress,
		uint beginPosOfCells) {
// this means that both nodes come from cells
	if (bothCellNodes(type1, type2)) {
		// this means that nodes come from different type of cell, apply differential adhesion
		if (type1 != type2) {
			// TODO: apply differential adhesion here.
			// It should be a different type of inter force.
			calculateAndAddDiffInterCellForce(xPos, yPos, zPos,
					_nodeLocXAddress[nodeRank2], _nodeLocYAddress[nodeRank2],
					_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
		} else {
			if (isSameCell(nodeRank1, nodeRank2)) {
				calculateAndAddIntraForce(xPos, yPos, zPos,
						_nodeLocXAddress[nodeRank2],
						_nodeLocYAddress[nodeRank2],
						_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
			} else {
				double xPre = xRes;
				double yPre = yRes;
				double zPre = zRes;
				calculateAndAddInterForce(xPos, yPos, zPos,
						_nodeLocXAddress[nodeRank2],
						_nodeLocYAddress[nodeRank2],
						_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
				double xDiff = xRes - xPre;
				double yDiff = yRes - yPre;
				double zDiff = zRes - zPre;
				double force = sqrt(
						xDiff * xDiff + yDiff * yDiff + zDiff * zDiff);
				if (force > maxForce) {
					maxForce = force;
				}
			}
		}
	}

// this means that both nodes come from ECM and from same ECM
	else if (type1 == ECM && type2 == ECM && isSameECM(nodeRank1, nodeRank2)) {
		if (isNeighborECMNodes(nodeRank1, nodeRank2)) {
			calculateAndAddECMForce(xPos, yPos, zPos,
					_nodeLocXAddress[nodeRank2], _nodeLocYAddress[nodeRank2],
					_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
		}
		// if both nodes belong to same ECM but are not neighbors they shouldn't interact.
	}
// this means that both nodes come from profile ( Epithilum layer).
	else if (type1 == Profile && type2 == Profile) {
		if (isNeighborProfileNodes(nodeRank1, nodeRank2)) {
			// TODO: need a set of parameters for calculating linking force between profile nodes
			//calculateAndAddProfileForce(xPos, yPos, zPos,
			//		_nodeLocXAddress[nodeRank2], _nodeLocYAddress[nodeRank2],
			//		_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
		}
		// if both nodes belong to Profile but are not neighbors they shouldn't interact.

	} else {
		// for now, we assume that interaction between other nodes are the same as inter-cell force.
		calculateAndAddInterForce(xPos, yPos, zPos, _nodeLocXAddress[nodeRank2],
				_nodeLocYAddress[nodeRank2], _nodeLocZAddress[nodeRank2], xRes,
				yRes, zRes);
	}

}

void SceNodes::extendBuckets2D() {
	static const uint extensionFactor2D = 9;
	uint valuesCount = auxVecs.bucketValues.size();
	auxVecs.bucketKeysExpanded.resize(valuesCount * extensionFactor2D);
	auxVecs.bucketValuesIncludingNeighbor.resize(
			valuesCount * extensionFactor2D);

	/**
	 * beginning of constant iterator
	 */
	thrust::constant_iterator<uint> first(extensionFactor2D);
	/**
	 * end of constant iterator.
	 * the plus sign only indicate movement of position, not value.
	 * e.g. movement is 5 and first iterator is initialized as 9
	 * result array is [9,9,9,9,9];
	 */
	thrust::constant_iterator<uint> last = first + valuesCount;

	expand(first, last,
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeys.begin(),
							auxVecs.bucketValues.begin())),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							auxVecs.bucketValuesIncludingNeighbor.begin())));

	thrust::counting_iterator<uint> countingBegin(0);
	thrust::counting_iterator<uint> countingEnd = countingBegin
			+ valuesCount * extensionFactor2D;

//std::cout << "number of values for array holding extended value= "
//		<< valuesCount * extensionFactor2D << std::endl;
//thrust::for_each(
//		thrust::make_zip_iterator(
//				make_tuple(bucketKeysExpanded.begin(), countingBegin)),
//		thrust::make_zip_iterator(
//				make_tuple(bucketKeysExpanded.end(), countingEnd)),
//		NeighborFunctor2D(numOfBucketsInXDim, numOfBucketsInYDim));

	thrust::transform(
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.end(), countingEnd)),
			make_zip_iterator(
					make_tuple(auxVecs.bucketKeysExpanded.begin(),
							countingBegin)),
			NeighborFunctor2D(domainPara.numOfBucketsInXDim,
					domainPara.numOfBucketsInYDim));

	int numberOfOutOfRange = thrust::count(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(), UINT_MAX);
//std::cout << "number out of range = " << numberOfOutOfRange << std::endl;
	int sizeBeforeShrink = auxVecs.bucketKeysExpanded.size();
	int numberInsideRange = sizeBeforeShrink - numberOfOutOfRange;
	thrust::sort_by_key(auxVecs.bucketKeysExpanded.begin(),
			auxVecs.bucketKeysExpanded.end(),
			auxVecs.bucketValuesIncludingNeighbor.begin());
	auxVecs.bucketKeysExpanded.erase(
			auxVecs.bucketKeysExpanded.begin() + numberInsideRange,
			auxVecs.bucketKeysExpanded.end());
	auxVecs.bucketValuesIncludingNeighbor.erase(
			auxVecs.bucketValuesIncludingNeighbor.begin() + numberInsideRange,
			auxVecs.bucketValuesIncludingNeighbor.end());
}
void SceNodes::applySceForces() {

// reason for casting these pointers every time is for flexibility.
// Because these are thrust vectors, the begin address could possibly change.
// Therefore, we re-define the begin address of these pointers everytime.
	uint* valueAddress = thrust::raw_pointer_cast(
			&auxVecs.bucketValuesIncludingNeighbor[0]);
	double* nodeLocXAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocX[0]);
	double* nodeLocYAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocY[0]);
	double* nodeLocZAddress = thrust::raw_pointer_cast(&infoVecs.nodeLocZ[0]);
	uint* nodeRankAddress = thrust::raw_pointer_cast(&infoVecs.nodeCellRank[0]);
	SceNodeType* nodeTypeAddress = thrust::raw_pointer_cast(
			&infoVecs.nodeCellType[0]);

//std::cout << "begin transformation" << std::endl;

	thrust::transform(
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.begin()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.begin()),
							auxVecs.bucketValues.begin(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeLocZ.begin(),
									auxVecs.bucketValues.begin()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(auxVecs.keyBegin.begin(),
									auxVecs.bucketKeys.end()),
							make_permutation_iterator(auxVecs.keyEnd.begin(),
									auxVecs.bucketKeys.end()),
							auxVecs.bucketValues.end(),
							make_permutation_iterator(infoVecs.nodeLocX.begin(),
									auxVecs.bucketValues.end()),
							make_permutation_iterator(infoVecs.nodeLocY.begin(),
									auxVecs.bucketValues.end()),
							make_permutation_iterator(infoVecs.nodeLocZ.begin(),
									auxVecs.bucketValues.end()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(infoVecs.nodeVelX.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelY.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(infoVecs.nodeVelZ.begin(),
									auxVecs.bucketValues.begin()),
							make_permutation_iterator(
									infoVecs.nodeMaxForce.begin(),
									auxVecs.bucketValues.begin()))),
			AddSceForce(valueAddress, nodeLocXAddress, nodeLocYAddress,
					nodeLocZAddress, nodeRankAddress, nodeTypeAddress,
					allocPara.maxTotalCellNodeCount, allocPara.startPosCells,
					allocPara.maxNodeOfOneCell, allocPara.maxNodePerECM));

//std::cout << "after transformation" << std::endl;
}

void SceNodes::calculateAndApplySceForces() {
	prepareSceForceComputation();
	applySceForces();
	applyProfileForces();
}

const SceDomainPara& SceNodes::getDomainPara() const {
	return domainPara;
}

void SceNodes::setDomainPara(const SceDomainPara& domainPara) {
	this->domainPara = domainPara;
}

const NodeAllocPara& SceNodes::getAllocPara() const {
	return allocPara;
}

void SceNodes::setAllocPara(const NodeAllocPara& allocPara) {
	this->allocPara = allocPara;
}

const NodeAuxVecs& SceNodes::getAuxVecs() const {
	return auxVecs;
}

void SceNodes::setAuxVecs(const NodeAuxVecs& auxVecs) {
	this->auxVecs = auxVecs;
}

NodeInfoVecs& SceNodes::getInfoVecs() {
	return infoVecs;
}

void SceNodes::setInfoVecs(const NodeInfoVecs& infoVecs) {
	this->infoVecs = infoVecs;
}

void SceNodes::allocSpaceForNodes(uint maxTotalNodeCount) {
	infoVecs.nodeLocX.resize(maxTotalNodeCount);
	infoVecs.nodeLocY.resize(maxTotalNodeCount);
	infoVecs.nodeLocZ.resize(maxTotalNodeCount);
	infoVecs.nodeVelX.resize(maxTotalNodeCount);
	infoVecs.nodeVelY.resize(maxTotalNodeCount);
	infoVecs.nodeVelZ.resize(maxTotalNodeCount);
	infoVecs.nodeMaxForce.resize(maxTotalNodeCount);
	infoVecs.nodeCellType.resize(maxTotalNodeCount);
	infoVecs.nodeCellRank.resize(maxTotalNodeCount);
	infoVecs.nodeIsActive.resize(maxTotalNodeCount);
}

void SceNodes::initNodeAllocPara(uint totalBdryNodeCount,
		uint maxProfileNodeCount, uint maxTotalECMCount, uint maxNodeInECM,
		uint maxTotalCellCount, uint maxNodeInCell) {
	allocPara.maxCellCount = maxTotalCellCount;
	allocPara.maxNodeOfOneCell = maxNodeInCell;
	allocPara.maxNodePerECM = maxNodeInECM;
	allocPara.maxECMCount = maxTotalECMCount;
	allocPara.maxProfileNodeCount = maxProfileNodeCount;
	//std::cout << "break point 1" << std::endl;
	allocPara.currentActiveProfileNodeCount = 0;
	allocPara.BdryNodeCount = totalBdryNodeCount;
	allocPara.currentActiveCellCount = 0;
	allocPara.maxTotalECMNodeCount = allocPara.maxECMCount
			* allocPara.maxNodePerECM;
	allocPara.currentActiveECM = 0;
	//std::cout << "break point 2" << std::endl;
	// will need to change this value after we have more detail about ECM
	allocPara.maxTotalCellNodeCount = maxTotalCellCount
			* allocPara.maxNodeOfOneCell;

	allocPara.startPosProfile = totalBdryNodeCount;
	allocPara.startPosECM = allocPara.startPosProfile
			+ allocPara.maxProfileNodeCount;
	allocPara.startPosCells = allocPara.startPosECM
			+ allocPara.maxTotalECMNodeCount;
}

void SceNodes::initNodeAllocPara_v2(uint totalBdryNodeCount,
		uint maxProfileNodeCount, uint maxCartNodeCount, uint maxTotalECMCount,
		uint maxNodeInECM, uint maxTotalCellCount, uint maxNodeInCell) {
	allocPara.maxCellCount = maxTotalCellCount;
	allocPara.maxNodeOfOneCell = maxNodeInCell;
	allocPara.maxNodePerECM = maxNodeInECM;
	allocPara.maxECMCount = maxTotalECMCount;
	allocPara.maxProfileNodeCount = maxProfileNodeCount;
	allocPara.maxCartNodeCount = maxCartNodeCount;
	//std::cout << "break point 1" << std::endl;
	allocPara.currentActiveProfileNodeCount = 0;
	allocPara.currentActiveCartNodeCount = 0;
	allocPara.BdryNodeCount = totalBdryNodeCount;
	allocPara.currentActiveCellCount = 0;
	allocPara.maxTotalECMNodeCount = allocPara.maxECMCount
			* allocPara.maxNodePerECM;
	allocPara.currentActiveECM = 0;
	//std::cout << "break point 2" << std::endl;
	// will need to change this value after we have more detail about ECM
	allocPara.maxTotalCellNodeCount = maxTotalCellCount
			* allocPara.maxNodeOfOneCell;

	allocPara.startPosProfile = totalBdryNodeCount;
	allocPara.startPosCart = allocPara.startPosProfile
			+ allocPara.maxProfileNodeCount;
	allocPara.startPosECM = allocPara.startPosCart + allocPara.maxCartNodeCount;
	allocPara.startPosCells = allocPara.startPosECM
			+ allocPara.maxTotalECMNodeCount;
}

