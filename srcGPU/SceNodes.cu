#include "SceNodes.h"

__constant__ double sceInterPara[5];
__constant__ double sceIntraPara[4];
__constant__ double sceInterDiffPara[5];
__constant__ double sceProfilePara[7];
__constant__ double sceECMPara[5];
__constant__ double sceDiffPara[5];

double sceInterParaCPU[5];
double sceIntraParaCPU[4];
double sceInterDiffParaCPU[5];
double sceProfileParaCPU[7];
double sceECMParaCPU[5];
double sceDiffParaCPU[5];

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

SceNodes::SceNodes(uint totalBdryNodeCount, uint maxProfileNodeCount,
		uint maxTotalECMCount, uint maxNodeInECM, uint maxTotalCellCount,
		uint maxNodeInCell) {
	std::cout << "start creating SceNodes object" << std::endl;
	maxCellCount = maxTotalCellCount;
	maxNodeOfOneCell = maxNodeInCell;
	maxNodePerECM = maxNodeInECM;
	maxECMCount = maxTotalECMCount;
	this->maxProfileNodeCount = maxProfileNodeCount;
	currentActiveProfileNodeCount = 0;
	BdryNodeCount = totalBdryNodeCount;
	currentActiveCellCount = 0;
	maxTotalECMNodeCount = maxECMCount * maxNodePerECM;
	currentActiveECM = 0;

	// will need to change this value after we have more detail about ECM
	maxTotalCellNodeCount = maxTotalCellCount * maxNodeOfOneCell;

	//cellRanks.resize(maxTotalNodeCount);
	//nodeRanks.resize(maxTotalNodeCount);
	//std::cout << "before resizing vectors" << std::endl;
	uint maxTotalNodeCount = totalBdryNodeCount + maxProfileNodeCount
			+ maxTotalECMNodeCount + maxTotalCellNodeCount;
	//std::cout << "maxTotalNodeCount = " << maxTotalNodeCount << std::endl;
	//thrust::host_vector<bool> nodeIsActiveHost

	nodeLocX.resize(maxTotalNodeCount);
	nodeLocY.resize(maxTotalNodeCount);
	nodeLocZ.resize(maxTotalNodeCount);
	nodeVelX.resize(maxTotalNodeCount);
	nodeVelY.resize(maxTotalNodeCount);
	nodeVelZ.resize(maxTotalNodeCount);
	nodeCellType.resize(maxTotalNodeCount);
	nodeCellRank.resize(maxTotalNodeCount);
	nodeIsActive.resize(maxTotalNodeCount);

	startPosProfile = totalBdryNodeCount;
	startPosECM = startPosProfile + maxProfileNodeCount;
	startPosCells = startPosECM + maxTotalECMNodeCount;

	std::cout << "start pos Profile = " << startPosProfile << ", startPosECM = "
			<< startPosECM << ", startPosCells = " << startPosCells
			<< std::endl;
	//int jj;
	//std::cin >> jj;

	thrust::host_vector<CellType> hostTmpVector(maxTotalNodeCount);
	thrust::host_vector<bool> hostTmpVector2(maxTotalNodeCount);
	thrust::host_vector<int> hostTmpVector3(maxTotalNodeCount);
	for (int i = 0; i < maxTotalNodeCount; i++) {
		if (i < startPosProfile) {
			hostTmpVector[i] = Boundary;
			hostTmpVector3[i] = 0;
		} else if (i < startPosECM) {
			hostTmpVector[i] = Profile;
			hostTmpVector3[i] = 0;
		} else if (i < startPosCells) {
			hostTmpVector[i] = ECM;
			hostTmpVector3[i] = (i - startPosECM) / maxNodeInECM;
		} else {
			// all initialized as FNM
			hostTmpVector[i] = FNM;
			hostTmpVector3[i] = (i - startPosCells) / maxNodeOfOneCell;
		}
		nodeIsActive[i] = false;
	}
	nodeCellType = hostTmpVector;
	nodeIsActive = hostTmpVector2;
	nodeCellRank = hostTmpVector3;
	copyParaToGPUConstMem();
}

void SceNodes::copyParaToGPUConstMem() {
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

	sceInterParaCPU[0] = U0;
	sceInterParaCPU[1] = V0;
	sceInterParaCPU[2] = k1;
	sceInterParaCPU[3] = k2;
	sceInterParaCPU[4] = interLinkEffectiveRange;

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
	sceIntraParaCPU[0] = U0_Intra;
	sceIntraParaCPU[1] = V0_Intra;
	sceIntraParaCPU[2] = k1_Intra;
	sceIntraParaCPU[3] = k2_Intra;

	//std::cout << "in SceNodes, before cuda memory copy to symbol:" << std::endl;
	cudaMemcpyToSymbol(sceInterPara, sceInterParaCPU, 5 * sizeof(double));
	cudaMemcpyToSymbol(sceIntraPara, sceIntraParaCPU, 4 * sizeof(double));
	cudaMemcpyToSymbol(ProfilebeginPos, &startPosProfile, sizeof(uint));
	cudaMemcpyToSymbol(ECMbeginPos, &startPosECM, sizeof(uint));
	cudaMemcpyToSymbol(cellNodeBeginPos, &startPosCells, sizeof(uint));
	cudaMemcpyToSymbol(nodeCountPerECM, &maxNodePerECM, sizeof(uint));
	cudaMemcpyToSymbol(nodeCountPerCell, &maxNodeOfOneCell, sizeof(uint));

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
	sceInterDiffParaCPU[0] = U0_Diff;
	sceInterDiffParaCPU[1] = V0_Diff;
	sceInterDiffParaCPU[2] = k1_Diff;
	sceInterDiffParaCPU[3] = k2_Diff;
	sceInterDiffParaCPU[4] = interLinkEffectiveRange;

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

	sceProfileParaCPU[0] = U0_Bdry;
	sceProfileParaCPU[1] = V0_Bdry;
	sceProfileParaCPU[2] = k1_Bdry;
	sceProfileParaCPU[3] = k2_Bdry;
	sceProfileParaCPU[4] = interLinkEffectiveRange;
	sceProfileParaCPU[5] = linearParameter;
	sceProfileParaCPU[6] = neutralLength;

	std::cout << "linear parameter = " << linearParameter << std::endl;
	std::cout << "neutralLength  =" << neutralLength << std::endl;
	//int jj;
	//std::cin >> jj;

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
	sceECMParaCPU[0] = U0_ECM;
	sceECMParaCPU[1] = V0_ECM;
	sceECMParaCPU[2] = k1_ECM;
	sceECMParaCPU[3] = k2_ECM;
	sceECMParaCPU[4] = interLinkEffectiveRange;

	cudaMemcpyToSymbol(sceProfilePara, sceProfileParaCPU, 7 * sizeof(double));

	cudaMemcpyToSymbol(sceInterDiffPara, sceInterDiffParaCPU,
			5 * sizeof(double));

	cudaMemcpyToSymbol(sceECMPara, sceECMParaCPU, 5 * sizeof(double));
	//std::cout << "finished SceNodes:" << std::endl;
}

void SceNodes::initDimension(double domainMinX, double domainMaxX,
		double domainMinY, double domainMaxY, double domainBucketSize) {
	minX = domainMinX;
	maxX = domainMaxX;
	minY = domainMinY;
	maxY = domainMaxY;
	bucketSize = domainBucketSize;
	numOfBucketsInXDim = (maxX - minX) / bucketSize + 1;
	numOfBucketsInYDim = (maxY - minY) / bucketSize + 1;
	totalBucketCount = numOfBucketsInXDim * numOfBucketsInYDim;

	keyBegin.resize(totalBucketCount);
	keyEnd.resize(totalBucketCount);

	/*
	 std::cout << "after initialization, values:" << std::endl;
	 std::cout << "minX = " << minX << ", maxX = " << maxX << std::endl;
	 std::cout << "minX = " << minX << ", maxX = " << maxX << std::endl;
	 std::cout << "numOfBucketsInXDim = " << numOfBucketsInXDim
	 << ", numOfBucketsInYDim = " << numOfBucketsInYDim << std::endl;
	 std::cout << "totalBucketCount= " << totalBucketCount << std::endl;
	 */

	//int jj;
	//std::cin >> jj;
}

std::vector<std::pair<uint, uint> > SceNodes::obtainNeighborPairs() {
	std::vector<std::pair<uint, uint> > result;
	thrust::host_vector<uint> keyBeginCPU = keyBegin;
	thrust::host_vector<uint> keyEndCPU = keyEnd;
	thrust::host_vector<uint> bucketKeysCPU = bucketKeys;
	thrust::host_vector<uint> bucketValuesCPU = bucketValues;
	thrust::host_vector<uint> bucketValuesExtendedCPU =
			bucketValuesIncludingNeighbor;

	int size = bucketKeysCPU.size();
	for (int i = 0; i < size; i++) {
		for (int j = keyBeginCPU[bucketKeysCPU[i]];
				j < keyEndCPU[bucketKeysCPU[i]]; j++) {
			//std::cout << "pair node 1: " << bucketValues[i] << ",pair node2: "
			//		<< bucketValuesIncludingNeighbor[j] << std::endl;
			result.push_back(
					std::make_pair<uint, uint>(bucketValues[i],
							bucketValuesIncludingNeighbor[j]));
		}
	}

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

	uint beginAddressOfProfile = startPosProfile;
	// find the begining position of ECM.
	uint beginAddressOfECM = startPosECM;
	// find the begining position of FNM cells.
	uint beginAddressOfFNM = startPosCells;
	// find the begining position of MX cells.
	uint beginAddressOfMX = beginAddressOfFNM + FNMNodeCountX;

	//std::cerr << "before copying arrays" << endl;

	thrust::copy(initBdryCellNodePosX.begin(), initBdryCellNodePosX.end(),
			nodeLocX.begin());
	thrust::copy(initBdryCellNodePosY.begin(), initBdryCellNodePosY.end(),
			nodeLocY.begin());

	//std::cerr << "copy 1" << endl;

	// copy x and y position of nodes of Profile to actual node position.
	thrust::copy(initProfileNodePosX.begin(), initProfileNodePosX.end(),
			nodeLocX.begin() + beginAddressOfProfile);
	thrust::copy(initProfileNodePosY.begin(), initProfileNodePosY.end(),
			nodeLocY.begin() + beginAddressOfProfile);

	//std::cerr << "copy 2" << endl;

	// copy x and y position of nodes of ECM to actual node position.
	thrust::copy(initECMNodePosX.begin(), initECMNodePosX.end(),
			nodeLocX.begin() + beginAddressOfECM);
	thrust::copy(initECMNodePosY.begin(), initECMNodePosY.end(),
			nodeLocY.begin() + beginAddressOfECM);

	// debug
	for (int i = 0; i < initECMNodePosX.size(); i++) {
		std::cout << "i + beginAddressOfECM = " << (i + beginAddressOfECM)
				<< "nodeLocX =" << nodeLocX[i + beginAddressOfECM] << std::endl;
		assert(nodeLocX[i + beginAddressOfECM] == initECMNodePosX[i]);
		assert(!isnan(initECMNodePosX[i]));
	}

	// std::cerr << "copy 3" << endl;

	// copy x and y position of nodes of FNM cells to actual node position.
	thrust::copy(initFNMCellNodePosX.begin(), initFNMCellNodePosX.end(),
			nodeLocX.begin() + beginAddressOfFNM);
	thrust::copy(initFNMCellNodePosY.begin(), initFNMCellNodePosY.end(),
			nodeLocY.begin() + beginAddressOfFNM);

	// std::cerr << "copy 4" << endl;

	thrust::fill(nodeCellType.begin() + beginAddressOfFNM,
			nodeCellType.begin() + beginAddressOfMX, FNM);

	// copy x and y position of nodes of MX cells to actual node position.
	thrust::copy(initMXCellNodePosX.begin(), initMXCellNodePosX.end(),
			nodeLocX.begin() + beginAddressOfMX);
	thrust::copy(initMXCellNodePosY.begin(), initMXCellNodePosY.end(),
			nodeLocY.begin() + beginAddressOfMX);

	//std::cerr << "after copying arrays" << endl;

	thrust::fill(nodeCellType.begin() + beginAddressOfMX,
			nodeCellType.begin() + beginAddressOfMX + MXNodeCountX, MX);

	//std::cout << "initial MX cell numbers: " << mxQuotient << std::endl;
}

void SceNodes::applyProfileForces() {
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(
			currentActiveProfileNodeCount);

	double* nodeLocXAddressEpiBegin = thrust::raw_pointer_cast(
			&nodeLocX[startPosProfile]);
	double* nodeLocYAddressEpiBegin = thrust::raw_pointer_cast(
			&nodeLocY[startPosProfile]);
	double* nodeLocZAddressEpiBegin = thrust::raw_pointer_cast(
			&nodeLocZ[startPosProfile]);

	double* nodeVelXAddressEpiBegin = thrust::raw_pointer_cast(
			&nodeVelX[startPosProfile]);
	double* nodeVelYAddressEpiBegin = thrust::raw_pointer_cast(
			&nodeVelY[startPosProfile]);
	double* nodeVelZAddressEpiBegin = thrust::raw_pointer_cast(
			&nodeVelZ[startPosProfile]);

	thrust::transform(countingIterBegin, countingIterEnd,
			thrust::make_zip_iterator(
					thrust::make_tuple(nodeVelX.begin(), nodeVelY.begin(),
							nodeVelZ.begin())) + startPosProfile,
			AddLinkForces(nodeLocXAddressEpiBegin, nodeLocYAddressEpiBegin,
					nodeLocZAddressEpiBegin, nodeVelXAddressEpiBegin,
					nodeVelYAddressEpiBegin, nodeVelZAddressEpiBegin,
					currentActiveProfileNodeCount));
}

void SceNodes::addNewlyDividedCells(
		thrust::device_vector<double> &nodeLocXNewCell,
		thrust::device_vector<double> &nodeLocYNewCell,
		thrust::device_vector<double> &nodeLocZNewCell,
		thrust::device_vector<bool> &nodeIsActiveNewCell,
		thrust::device_vector<CellType> &nodeCellTypeNewCell) {

	// data validation
	uint nodesSize = nodeLocXNewCell.size();
	assert(nodesSize % maxNodeOfOneCell == 0);
	uint addCellCount = nodesSize / maxNodeOfOneCell;

	// position that we will add newly divided cells.
	uint shiftStartPosNewCell = startPosCells
			+ currentActiveCellCount * maxNodeOfOneCell;

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
					thrust::make_tuple(nodeLocX.begin(), nodeLocY.begin(),
							nodeLocZ.begin(), nodeIsActive.begin(),
							nodeCellType.begin())) + shiftStartPosNewCell);

	// total number of cells has increased.
	currentActiveCellCount = currentActiveCellCount + addCellCount;
}

void SceNodes::buildBuckets2D() {
	int totalActiveNodes = startPosCells
			+ currentActiveCellCount * maxNodeOfOneCell;

	bucketKeys.resize(totalActiveNodes);
	bucketValues.resize(totalActiveNodes);
	thrust::counting_iterator<uint> countingIterBegin(0);
	thrust::counting_iterator<uint> countingIterEnd(totalActiveNodes);

	// takes counting iterator and coordinates
	// return tuple of keys and values

	// transform the points to their bucket indices
	thrust::transform(
			make_zip_iterator(
					make_tuple(nodeLocX.begin(), nodeLocY.begin(),
							nodeLocZ.begin(), nodeIsActive.begin(),
							countingIterBegin)),
			make_zip_iterator(
					make_tuple(nodeLocX.begin(), nodeLocY.begin(),
							nodeLocZ.begin(), nodeIsActive.begin(),
							countingIterBegin)) + totalActiveNodes,
			make_zip_iterator(
					make_tuple(bucketKeys.begin(), bucketValues.begin())),
			pointToBucketIndex2D(minX, maxX, minY, maxY, bucketSize));

	// sort the points by their bucket index
	thrust::sort_by_key(bucketKeys.begin(), bucketKeys.end(),
			bucketValues.begin());
	// for those nodes that are inactive, key value of UINT_MAX will be returned.
	// we need to removed those keys along with their values.
	int numberOfOutOfRange = thrust::count(bucketKeys.begin(), bucketKeys.end(),
			UINT_MAX);
	bucketKeys.erase(bucketKeys.end() - numberOfOutOfRange, bucketKeys.end());
	bucketValues.erase(bucketValues.end() - numberOfOutOfRange,
			bucketValues.end());
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
	/*
	 if (linkLength > sceProfilePara[4]) {
	 forceValue = 0;
	 } else {
	 forceValue = -sceProfilePara[0] / sceProfilePara[2]
	 * exp(-linkLength / sceProfilePara[2])
	 + sceProfilePara[1] / sceProfilePara[3]
	 * exp(-linkLength / sceProfilePara[3]);
	 // positive value means force is attraction
	 if (linkLength > sceProfilePara[6]) {
	 forceValue = sceProfilePara[5] * (linkLength - sceProfilePara[6]);
	 //if (forceValue < 0) {
	 //	forceValue = 0;
	 //}
	 }
	 }
	 */

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
		forceValue = 0;
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

__device__ bool bothCellNodes(CellType &type1, CellType &type2) {
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

	/*
	 if (linkLength > sceProfilePara[4]) {
	 forceValue = 0;
	 } else {
	 forceValue = -sceProfilePara[0] / sceProfilePara[2]
	 * exp(-linkLength / sceProfilePara[2])
	 + sceProfilePara[1] / sceProfilePara[3]
	 * exp(-linkLength / sceProfilePara[3]);
	 // positive value means force is attraction
	 if (linkLength > sceProfilePara[6]) {
	 forceValue = sceProfilePara[5] * (linkLength - sceProfilePara[6]);
	 //if (forceValue < 0) {
	 //	forceValue = 0;
	 //}
	 }
	 }
	 */
}

__device__
void handleForceBetweenNodes(uint &nodeRank1, CellType &type1, uint &nodeRank2,
		CellType &type2, double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes, double* _nodeLocXAddress, double* _nodeLocYAddress,
		double* _nodeLocZAddress, uint beginPosOfCells) {
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
			// TODO: this function needs to be modified.
			// (1) nodeCountPerCell need to be stored in constant memory.
			// (2) begin address of cell nodes need to be stored in constant memory.
			if (isSameCell(nodeRank1, nodeRank2)) {
				calculateAndAddIntraForce(xPos, yPos, zPos,
						_nodeLocXAddress[nodeRank2],
						_nodeLocYAddress[nodeRank2],
						_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
			} else {
				calculateAndAddInterForce(xPos, yPos, zPos,
						_nodeLocXAddress[nodeRank2],
						_nodeLocYAddress[nodeRank2],
						_nodeLocZAddress[nodeRank2], xRes, yRes, zRes);
			}
		}
	}
	// this means that both nodes come from ECM and from same ECM
	else if (type1 == ECM && type2 == ECM && isSameECM(nodeRank1, nodeRank2)) {
		if (isNeighborECMNodes(nodeRank1, nodeRank2)) {
			// TODO: need to create another two vectors that holds the neighbor information for ECM.
			// TODO: alternatively, try to store ECM begin address and number of node per ECM in constant memory.
			// TODO: implement this function.
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
	uint valuesCount = bucketValues.size();
	bucketKeysExpanded.resize(valuesCount * extensionFactor2D);
	bucketValuesIncludingNeighbor.resize(valuesCount * extensionFactor2D);

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
					make_tuple(bucketKeys.begin(), bucketValues.begin())),
			make_zip_iterator(
					make_tuple(bucketKeysExpanded.begin(),
							bucketValuesIncludingNeighbor.begin())));

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
					make_tuple(bucketKeysExpanded.begin(), countingBegin)),
			make_zip_iterator(
					make_tuple(bucketKeysExpanded.end(), countingEnd)),
			make_zip_iterator(
					make_tuple(bucketKeysExpanded.begin(), countingBegin)),
			NeighborFunctor2D(numOfBucketsInXDim, numOfBucketsInYDim));

	int numberOfOutOfRange = thrust::count(bucketKeysExpanded.begin(),
			bucketKeysExpanded.end(), UINT_MAX);
//std::cout << "number out of range = " << numberOfOutOfRange << std::endl;
	int sizeBeforeShrink = bucketKeysExpanded.size();
	int numberInsideRange = sizeBeforeShrink - numberOfOutOfRange;
	thrust::sort_by_key(bucketKeysExpanded.begin(), bucketKeysExpanded.end(),
			bucketValuesIncludingNeighbor.begin());
	bucketKeysExpanded.erase(bucketKeysExpanded.begin() + numberInsideRange,
			bucketKeysExpanded.end());
	bucketValuesIncludingNeighbor.erase(
			bucketValuesIncludingNeighbor.begin() + numberInsideRange,
			bucketValuesIncludingNeighbor.end());
}
void SceNodes::applySceForces() {
	std::cout << "begin apply sce forces" << std::endl;
	std::cout << "size of lower = " << keyBegin.size() << std::endl;
	thrust::counting_iterator<unsigned int> search_begin(0);
	thrust::lower_bound(bucketKeysExpanded.begin(), bucketKeysExpanded.end(),
			search_begin, search_begin + totalBucketCount, keyBegin.begin());
	thrust::upper_bound(bucketKeysExpanded.begin(), bucketKeysExpanded.end(),
			search_begin, search_begin + totalBucketCount, keyEnd.begin());

	thrust::host_vector<uint> lowerCPU = keyBegin;

	std::cout << "finished finding bounds" << std::endl;

	//int test1 = lowerCPU[0];
	//int test2 = lowerCPU[0];

	//std::cout << "test 1 =" << test1 << ", test 2 = " << test2 << std::endl;
	//std::cout.flush();

	//int test3 = keyBegin[totalBucketCount - 1];
	//int test4 = keyEnd[totalBucketCount - 1];

	//std::cout << "test 3 =" << test3 << ", test 4 = " << test4 << std::endl;
	uint* valueAddress = thrust::raw_pointer_cast(
			&bucketValuesIncludingNeighbor[0]);

	std::cout << "begin pointer casting" << std::endl;

	double* nodeLocXAddress = thrust::raw_pointer_cast(&nodeLocX[0]);
	double* nodeLocYAddress = thrust::raw_pointer_cast(&nodeLocY[0]);
	double* nodeLocZAddress = thrust::raw_pointer_cast(&nodeLocZ[0]);
	uint* nodeRankAddress = thrust::raw_pointer_cast(&nodeCellRank[0]);
	CellType* nodeTypeAddress = thrust::raw_pointer_cast(&nodeCellType[0]);

	std::cout << "begin transformation" << std::endl;

	thrust::transform(
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(keyBegin.begin(),
									bucketKeys.begin()),
							make_permutation_iterator(keyEnd.begin(),
									bucketKeys.begin()), bucketValues.begin(),
							make_permutation_iterator(nodeLocX.begin(),
									bucketValues.begin()),
							make_permutation_iterator(nodeLocY.begin(),
									bucketValues.begin()),
							make_permutation_iterator(nodeLocZ.begin(),
									bucketValues.begin()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(keyBegin.begin(),
									bucketKeys.end()),
							make_permutation_iterator(keyEnd.begin(),
									bucketKeys.end()), bucketValues.end(),
							make_permutation_iterator(nodeLocX.begin(),
									bucketValues.end()),
							make_permutation_iterator(nodeLocY.begin(),
									bucketValues.end()),
							make_permutation_iterator(nodeLocZ.begin(),
									bucketValues.end()))),
			make_zip_iterator(
					make_tuple(
							make_permutation_iterator(nodeVelX.begin(),
									bucketValues.begin()),
							make_permutation_iterator(nodeVelY.begin(),
									bucketValues.begin()),
							make_permutation_iterator(nodeVelZ.begin(),
									bucketValues.begin()))),
			AddSceForce(valueAddress, nodeLocXAddress, nodeLocYAddress,
					nodeLocZAddress, nodeRankAddress, nodeTypeAddress,
					maxTotalCellNodeCount, startPosCells, maxNodeOfOneCell,
					maxNodePerECM));

	std::cout << "after transformation" << std::endl;
}

void SceNodes::calculateAndApplySceForces() {
	//const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1;
	//const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1;
	std::cout << "in SceNodes, before build buckets 2D:" << std::endl;
	buildBuckets2D();
	std::cout << "in SceNodes, before extend buckets 2D:" << std::endl;
	extendBuckets2D();
	std::cout << "in SceNodes, before apply sce forces:" << std::endl;
	applySceForces();
	std::cout << "in SceNodes, finished apply sce forces:" << std::endl;
	applyProfileForces();
	std::cout << "in SceNodes, finished apply sce forces:" << std::endl;
}

