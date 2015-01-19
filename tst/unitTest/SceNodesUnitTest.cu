#include <iostream>
//#include "gtest/gtest.h"
#include "gtest/gtest.h"
#include "SceNodes.h"
#include <vector>
#include <cuda_runtime.h>
#include <algorithm>
#include <tr1/unordered_set>
#include "commonData.h"
using namespace std;

const double errTol = 1.0e-12;

GlobalConfigVars globalConfigVars;

double computeDistInTest(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2)
					+ (zPos - zPos2) * (zPos - zPos2));
}

void calculateAndAddIntraForceInTest(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes, vector<double> &sceIntraPara) {
	double linkLength = computeDistInTest(xPos, yPos, zPos, xPos2, yPos2,
			zPos2);
	double forceValue = -sceIntraPara[0] / sceIntraPara[2]
			* exp(-linkLength / sceIntraPara[2])
			+ sceIntraPara[1] / sceIntraPara[3]
					* exp(-linkLength / sceIntraPara[3]);
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}
void calculateAndAddInterForceInTest(double &xPos, double &yPos, double &zPos,
		double &xPos2, double &yPos2, double &zPos2, double &xRes, double &yRes,
		double &zRes, vector<double> &sceInterPara) {
	double linkLength = computeDistInTest(xPos, yPos, zPos, xPos2, yPos2,
			zPos2);
	double forceValue = 0;
	if (linkLength > sceInterPara[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterPara[0] / sceInterPara[2]
				* exp(-linkLength / sceInterPara[2])
				+ sceInterPara[1] / sceInterPara[3]
						* exp(-linkLength / sceInterPara[3]);
	}
	//if (forceValue > 0) {
	//	forceValue = 0;
	//}
	xRes = xRes + forceValue * (xPos2 - xPos) / linkLength;
	yRes = yRes + forceValue * (yPos2 - yPos) / linkLength;
	zRes = zRes + forceValue * (zPos2 - zPos) / linkLength;
}

void computeResultFromCPUAllIntra2D(vector<double> &xPoss,
		vector<double> &yPoss, vector<double> &zPoss, vector<double> &xVels,
		vector<double> &yVels, vector<double> &zVels, vector<bool> &isActive,
		vector<double> &paraSet, double bucketSize, double minX, double maxX,
		double minY, double maxY) {
	unsigned int i, j;
	for (i = 0; i < xVels.size(); i++) {
		xVels[i] = 0.0;
		yVels[i] = 0.0;
		zVels[i] = 0.0;
	}
	for (i = 0; i < xPoss.size(); i++) {
		if (isActive[i] == true) {
			int xBucketPos = (xPoss[i] - minX) / bucketSize;
			int yBucketPos = (yPoss[i] - minY) / bucketSize;
			for (j = 0; j < xPoss.size(); j++) {
				if (j != i) {
					if (isActive[j] == true) {
						int xBucketPos2 = (xPoss[j] - minX) / bucketSize;
						int yBucketPos2 = (yPoss[j] - minY) / bucketSize;
						if (abs(xBucketPos - xBucketPos2) <= 1
								&& abs(yBucketPos - yBucketPos2) <= 1) {
							calculateAndAddIntraForceInTest(xPoss[i], yPoss[i],
									zPoss[i], xPoss[j], yPoss[j], zPoss[i],
									xVels[i], yVels[i], zVels[i], paraSet);
						}
					}
				}
			}
		}
	}
}

void computeResultFromCPUAllIntraAndInter2D(vector<double> &xPoss,
		vector<double> &yPoss, vector<double> &zPoss, vector<double> &xVels,
		vector<double> &yVels, vector<double> &zVels, vector<bool> &isActive,
		vector<double> &paraSet1, vector<double> &paraSet2, double bucketSize,
		uint activeCellCount, uint maxNodesPerCell, double minX, double maxX,
		double minY, double maxY) {
	unsigned int i, j;
	uint numberOfActiveNodes = activeCellCount * maxNodesPerCell;

	for (i = 0; i < numberOfActiveNodes; i++) {
		xVels[i] = 0.0;
		yVels[i] = 0.0;
		zVels[i] = 0.0;
	}
	for (i = 0; i < numberOfActiveNodes; i++) {
		if (isActive[i] == true) {
			int xBucketPos = (xPoss[i] - minX) / bucketSize;
			int yBucketPos = (yPoss[i] - minY) / bucketSize;
			for (j = 0; j < numberOfActiveNodes; j++) {
				if (j != i) {
					if (isActive[j] == true) {
						int xBucketPos2 = (xPoss[j] - minX) / bucketSize;
						int yBucketPos2 = (yPoss[j] - minY) / bucketSize;
						if (abs(xBucketPos - xBucketPos2) <= 1
								&& abs(yBucketPos - yBucketPos2) <= 1) {
							if (i / maxNodesPerCell == j / maxNodesPerCell) {
								calculateAndAddIntraForceInTest(xPoss[i],
										yPoss[i], zPoss[i], xPoss[j], yPoss[j],
										zPoss[j], xVels[i], yVels[i], zVels[i],
										paraSet1);

							} else {
								calculateAndAddInterForceInTest(xPoss[i],
										yPoss[i], zPoss[i], xPoss[j], yPoss[j],
										zPoss[j], xVels[i], yVels[i], zVels[i],
										paraSet2);
								//std::cout << "inter logic:" << std::endl;
								//cin >> j;
							}
						}
					}
				} else {
					continue;
				}
			}
		} else {
			continue;
		}
	}
}

int calculateBucketKey(double minX, double maxX, double minY, double maxY,
		double bucketSize, double xcoord, double ycoord) {
	uint width = (maxX - minX) / bucketSize + 1;
	unsigned int x = (unsigned int) ((xcoord - minX) / bucketSize);
	unsigned int y = (unsigned int) ((ycoord - minY) / bucketSize);
	return (y * width + x);
}

class SceNodeTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		ConfigParser parser;
		std::string configFileName = "./resources/unitTest.cfg";
		globalConfigVars = parser.parseConfigFile(configFileName);
		cudaSetDevice(
				globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	}
};

uint Test_totalBdryNodeCount = 10;
uint Test_maxProfileNodeCount = 13;
uint Test_maxTotalECMCount = 5;
uint Test_maxNodeInECM = 3;
uint Test_maxTotalCellCount = 4;
uint Test_maxNodeInCell = 20;
double Test_minX = -1.0;
double Test_maxX = 34.4;
double Test_minY = -1.2;
double Test_maxY = 25.9;
double Test_bucketSize = 0.97;

TEST(DummyTest, SanityTest) {
	cudaSetDevice(0);
	EXPECT_EQ(32, 32);
	int size = 256;
	thrust::device_vector<unsigned int> dv(size);
}

TEST_F(SceNodeTest, SceNodeInitTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
}

TEST_F(SceNodeTest, ParameterInitTest) {
	//cout << " point 1 , before everything starts" << endl;
	int deviceID = globalConfigVars.getConfigValue("GPUDeviceNumber").toInt();
	int totalDeviceCount;
	cudaGetDeviceCount(&totalDeviceCount);
	EXPECT_TRUE(deviceID >= 0 && deviceID < totalDeviceCount);
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	//cout << " point 1.1 " << endl;

	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);

	//cout << " point 1.2 " << endl;
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);

	// following part tests parameters of the neighbor grid.
	int Expected_BucketXDim = (Test_maxX - Test_minX) / Test_bucketSize + 1;
	int Expected_BucketYDim = (Test_maxY - Test_minY) / Test_bucketSize + 1;
	int Expected_TotalBucket = Expected_BucketXDim * Expected_BucketYDim;
	EXPECT_EQ(Expected_BucketXDim, nodes.getDomainPara().numOfBucketsInXDim);
	EXPECT_EQ(Expected_BucketYDim, nodes.getDomainPara().numOfBucketsInYDim);
	EXPECT_EQ(Expected_TotalBucket, nodes.getDomainPara().totalBucketCount);
	//cout << " point 2, middle of no where" << endl;

	// following part tests initial parameters for node location info
	int startPosOfProfile = Test_totalBdryNodeCount;
	int startPosOfECM = Test_maxProfileNodeCount + startPosOfProfile;
	int startPosOfCells = Test_maxTotalECMCount * Test_maxNodeInECM
			+ startPosOfECM;
	EXPECT_EQ(startPosOfProfile, nodes.getAllocPara().startPosProfile);
	EXPECT_EQ(startPosOfECM, nodes.getAllocPara().startPosECM);
	EXPECT_EQ(startPosOfCells, nodes.getAllocPara().startPosCells);
	EXPECT_EQ(Test_totalBdryNodeCount, nodes.getAllocPara().BdryNodeCount);
	EXPECT_EQ(Test_maxTotalECMCount, nodes.getAllocPara().maxECMCount);
	EXPECT_EQ(Test_maxTotalCellCount, nodes.getAllocPara().maxCellCount);
	EXPECT_EQ(Test_maxProfileNodeCount,
			nodes.getAllocPara().maxProfileNodeCount);
	EXPECT_EQ(Test_maxNodeInCell, nodes.getAllocPara().maxNodeOfOneCell);
	EXPECT_EQ(Test_maxNodeInECM, nodes.getAllocPara().maxNodePerECM);

	//cout << " point 3, everything has finished" << endl;
}

TEST_F(SceNodeTest, MemSizeTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);
	int totalNodeSize = Test_totalBdryNodeCount + Test_maxProfileNodeCount
			+ Test_maxTotalECMCount * Test_maxNodeInECM
			+ Test_maxTotalCellCount * Test_maxNodeInCell;

	// size of buckets should be same with keyBeing and KeyEnd as they are
	// recording information for buckets
	EXPECT_EQ(nodes.getDomainPara().totalBucketCount,
			nodes.getAuxVecs().keyBegin.size());
	EXPECT_EQ(nodes.getDomainPara().totalBucketCount,
			nodes.getAuxVecs().keyEnd.size());

	// size of these node information should be the same with max node number.
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeCellRank.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeCellType.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeIsActive.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeLocX.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeLocY.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeLocZ.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeVelX.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeVelY.size());
	EXPECT_EQ(totalNodeSize, nodes.getInfoVecs().nodeVelZ.size());
}

TEST_F(SceNodeTest, InitialValueTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);
	//nodes.initValues();
}

TEST_F(SceNodeTest,BuildBucketFixedTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);

	const uint testTotalNodeCount = 4;

	thrust::host_vector<double> nodeLocXHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocYHost(testTotalNodeCount);
	thrust::host_vector<bool> nodeIsActiveHost(testTotalNodeCount);
	thrust::host_vector<uint> nodeExpectedBucket(testTotalNodeCount);
	nodeLocXHost[0] = 0.0;
	nodeLocYHost[0] = 0.0;
	nodeIsActiveHost[0] = 1;
	nodeExpectedBucket[0] = calculateBucketKey(Test_minX, Test_maxX, Test_minY,
			Test_maxY, Test_bucketSize, nodeLocXHost[0], nodeLocYHost[0]);

	nodeLocXHost[1] = 0.199;
	nodeLocYHost[1] = 0.4;
	nodeIsActiveHost[1] = 1;
	nodeExpectedBucket[1] = calculateBucketKey(Test_minX, Test_maxX, Test_minY,
			Test_maxY, Test_bucketSize, nodeLocXHost[1], nodeLocYHost[1]);

	nodeLocXHost[2] = 0.2;
	nodeLocYHost[2] = 3.4;
	nodeIsActiveHost[2] = 1;
	nodeExpectedBucket[2] = calculateBucketKey(Test_minX, Test_maxX, Test_minY,
			Test_maxY, Test_bucketSize, nodeLocXHost[2], nodeLocYHost[2]);

	nodeLocXHost[3] = 9.5;
	nodeLocYHost[3] = 5.212;
	nodeIsActiveHost[3] = 0;
	nodeExpectedBucket[3] = calculateBucketKey(Test_minX, Test_maxX, Test_minY,
			Test_maxY, Test_bucketSize, nodeLocXHost[3], nodeLocYHost[3]);

	nodes.getInfoVecs().nodeLocX = nodeLocXHost;
	nodes.getInfoVecs().nodeLocY = nodeLocYHost;
	nodes.getInfoVecs().nodeIsActive = nodeIsActiveHost;

	NodeAllocPara para = nodes.getAllocPara();
	para.startPosCells = 0;
	para.currentActiveCellCount = 2;
	para.maxNodeOfOneCell = 2;
	nodes.setAllocPara(para);
	nodes.prepareSceForceComputation();
	thrust::host_vector<uint> keysFromGPU = nodes.getAuxVecs().bucketKeys;
	thrust::host_vector<uint> valuesFromGPU = nodes.getAuxVecs().bucketValues;
	uint activeNodeCount = 0;
	for (uint i = 0; i < nodeIsActiveHost.size(); i++) {
		if (nodeIsActiveHost[i] == true) {
			activeNodeCount++;
		}
	}
	EXPECT_EQ(keysFromGPU.size(), activeNodeCount);
	EXPECT_EQ(keysFromGPU.size(), valuesFromGPU.size());
	for (uint i = 0; i < keysFromGPU.size(); i++) {
		uint nodeRank = valuesFromGPU[i];
		uint expectedResultFromCPU = nodeExpectedBucket[nodeRank];
		EXPECT_EQ(expectedResultFromCPU, keysFromGPU[i]);
	}
}

TEST_F(SceNodeTest,BuildBucketRandomTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
	const double minX = 0.5;
	const double maxX = 1.5;
	const double minY = 0.1;
	const double maxY = 2.2;
	const double bucketSize = 0.1;
	nodes.initDimension(minX, maxX, minY, maxY, bucketSize);
	const uint maxNodeCount = 40000;
	thrust::host_vector<double> nodeLocXHost(maxNodeCount);
	thrust::host_vector<double> nodeLocYHost(maxNodeCount);
	thrust::host_vector<double> nodeLocZHost(maxNodeCount);
	thrust::host_vector<bool> nodeIsActiveHost(maxNodeCount);

	const int width = (maxX - minX) / bucketSize + 1;
	thrust::counting_iterator<unsigned int> index_sequence_begin(0);
	thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
			nodeLocXHost.begin(), Prg(minX, maxX));
	thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
			nodeLocYHost.begin(), Prg(minY, maxY));

	for (uint i = 0; i < maxNodeCount; i++) {
		if (i % 2 == 0) {
			nodeIsActiveHost[i] = true;
		} else {
			nodeIsActiveHost[i] = false;
		}
	}
	nodes.getInfoVecs().nodeLocX = nodeLocXHost;
	nodes.getInfoVecs().nodeLocY = nodeLocYHost;
	nodes.getInfoVecs().nodeLocZ = nodeLocZHost;
	nodes.getInfoVecs().nodeIsActive = nodeIsActiveHost;

	NodeAllocPara para = nodes.getAllocPara();
	para.startPosCells = maxNodeCount;
	para.currentActiveCellCount = 0;
	para.maxNodeOfOneCell = 0;
	nodes.setAllocPara(para);

	nodes.prepareSceForceComputation();
	thrust::host_vector<uint> keysFromGPU = nodes.getAuxVecs().bucketKeys;
	thrust::host_vector<uint> valuesFromGPU = nodes.getAuxVecs().bucketValues;
	uint activeNodeCount = 0;
	for (uint i = 0; i < nodeIsActiveHost.size(); i++) {
		if (nodeIsActiveHost[i] == true) {
			activeNodeCount++;
		}
	}
	EXPECT_EQ(keysFromGPU.size(), activeNodeCount);
	EXPECT_EQ(valuesFromGPU.size(), activeNodeCount);
	for (uint i = 0; i < activeNodeCount; i++) {
		uint nodeRank = valuesFromGPU[i];
		uint resultFromCPU =
				(int) ((nodeLocXHost[nodeRank] - minX) / bucketSize)
						+ (int) ((nodeLocYHost[nodeRank] - minY) / bucketSize)
								* width;
		EXPECT_EQ(resultFromCPU, keysFromGPU[i]);
	}
}

TEST_F(SceNodeTest, ExtendBucketfixedTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
	const uint testCellCount = 2;
	const uint testNodePerCell = 2;
	const uint testTotalNodeCount = testCellCount * testNodePerCell;

	NodeAllocPara para = nodes.getAllocPara();
	para.maxNodeOfOneCell = testCellCount;
	para.currentActiveCellCount = testCellCount;
	nodes.setAllocPara(para);

	thrust::host_vector<double> nodeLocXHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocYHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocZHost(testTotalNodeCount);
	thrust::host_vector<bool> nodeIsActiveHost(testTotalNodeCount);
	thrust::host_vector<uint> nodeExpectedBucket(testTotalNodeCount);
	const double minX = 0.0;
	const double maxX = 0.99999;
	const double minY = 0.0;
	const double maxY = 0.99999;
	//const double minZ = 0.0;
	//const double maxZ = 0.0;
	const double bucketSize = 0.1;
	nodes.initDimension(minX, maxX, minY, maxY, bucketSize);
	nodeLocXHost[0] = 0.0;
	nodeLocYHost[0] = 0.0;
	nodeIsActiveHost[0] = true;
	// 0
	nodeLocXHost[1] = 0.51;
	nodeLocYHost[1] = 0.212;
	nodeIsActiveHost[1] = true;
	// 25
	nodeLocXHost[2] = 0.52;
	nodeLocYHost[2] = 0.211;
	nodeIsActiveHost[2] = true;
	// 25
	nodeLocXHost[3] = 0.63;
	nodeLocYHost[3] = 0.207;
	nodeIsActiveHost[3] = false;
	// 26
	nodes.getInfoVecs().nodeLocX = nodeLocXHost;
	nodes.getInfoVecs().nodeLocY = nodeLocYHost;
	nodes.getInfoVecs().nodeLocZ = nodeLocZHost;
	nodes.getInfoVecs().nodeIsActive = nodeIsActiveHost;
	nodes.prepareSceForceComputation();

	thrust::host_vector<uint> extendedKeysFromGPU =
			nodes.getAuxVecs().bucketKeysExpanded;
	thrust::host_vector<uint> extendValuesFromGPU =
			nodes.getAuxVecs().bucketValuesIncludingNeighbor;
	EXPECT_EQ(extendedKeysFromGPU.size(), (uint )22);
	EXPECT_EQ(extendValuesFromGPU.size(), (uint )22);
	int expectedKeys[] = { 0, 1, 10, 11, 14, 14, 15, 15, 16, 16, 24, 24, 25, 25,
			26, 26, 34, 34, 35, 35, 36, 36 };
	std::vector<uint> expectedResultsKeys(expectedKeys, expectedKeys + 22);
	for (int i = 0; i < 22; i++) {
		EXPECT_EQ(expectedResultsKeys[i], extendedKeysFromGPU[i]);
	}
}

/*
 * Expected size of the extended buckets equals to computed
 * all results fits the requirement
 * no duplicate
 */ //
TEST(SceExtendBucket2D, extendBucketRandomTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			0, Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell, false);
	uint maxNodeCount = 5253;

	thrust::host_vector<double> nodeLocXHost(maxNodeCount);
	thrust::host_vector<double> nodeLocYHost(maxNodeCount);
	thrust::host_vector<double> nodeLocZHost(maxNodeCount);
	thrust::host_vector<bool> nodeIsActiveHost(maxNodeCount);
	const double minX = 0.8;
	const double maxX = 1.9;
	const double minY = 0.6;
	const double maxY = 3.2;
	const double minZ = 0.0;
	const double maxZ = 0.0;
	const double bucketSize = 0.12;

	nodes.initDimension(minX, maxX, minY, maxY, bucketSize);

	NodeAllocPara para = nodes.getAllocPara();
	para.startPosCells = maxNodeCount;
	para.currentActiveCellCount = 0;
	para.maxNodeOfOneCell = 0;
	nodes.setAllocPara(para);

	thrust::counting_iterator<unsigned int> index_sequence_begin(0);
	thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
			nodeLocXHost.begin(), Prg(minX, maxX));
	thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
			nodeLocYHost.begin(), Prg(minY, maxY));
	thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
			nodeLocZHost.begin(), Prg(minZ, maxZ));
	for (uint i = 0; i < maxNodeCount; i++) {
		if (i % 3 == 0) {
			nodeIsActiveHost[i] = false;
		} else {
			nodeIsActiveHost[i] = true;
		}
	}
	nodes.getInfoVecs().nodeLocX = nodeLocXHost;
	nodes.getInfoVecs().nodeLocY = nodeLocYHost;
	nodes.getInfoVecs().nodeLocZ = nodeLocZHost;
	nodes.getInfoVecs().nodeIsActive = nodeIsActiveHost;
	const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1;
	const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1;
	nodes.prepareSceForceComputation();

	thrust::host_vector<uint> extendedKeysFromGPU =
			nodes.getAuxVecs().bucketKeysExpanded;
	thrust::host_vector<uint> extendValuesFromGPU =
			nodes.getAuxVecs().bucketValuesIncludingNeighbor;
	uint expectedResultKeyCount = 0;
	for (uint i = 0; i < maxNodeCount; i++) {
		if (nodeIsActiveHost[i] == true) {
			int xPos = (int) ((nodeLocXHost[i] - minX) / bucketSize);
			int yPos = (int) ((nodeLocYHost[i] - minY) / bucketSize);
			int xQuota = 3;
			int yQuota = 3;
			if (xPos == 0) {
				xQuota--;
			}
			if (xPos == numberOfBucketsInXDim - 1) {
				xQuota--;
			}
			if (yPos == 0) {
				yQuota--;
			}
			if (yPos == numberOfBucketsInYDim - 1) {
				yQuota--;
			}
			expectedResultKeyCount += xQuota * yQuota;
		}
	}
	// verify if size is correct
	EXPECT_EQ(expectedResultKeyCount, extendedKeysFromGPU.size());
	EXPECT_EQ(expectedResultKeyCount, extendValuesFromGPU.size());
	// verify if all key- value pairs fits our the requirement
	for (uint i = 0; i < extendValuesFromGPU.size(); i++) {
		int nodeRank = extendValuesFromGPU[i];
		int xPos = (int) ((nodeLocXHost[nodeRank] - minX) / bucketSize);
		int yPos = (int) ((nodeLocYHost[nodeRank] - minY) / bucketSize);
		int bucketXPos = extendedKeysFromGPU[i] % numberOfBucketsInXDim;
		int bucketYPos = extendedKeysFromGPU[i] / numberOfBucketsInXDim;
		EXPECT_TRUE(abs(xPos - bucketXPos) <= 1);
		EXPECT_TRUE(abs(yPos - bucketYPos) <= 1);
	}
	//verify for each key, there is no duplicate for its values
	std::vector<uint> previousValues;
	bool startNewFlag = 1;
	for (uint i = 0; i < extendedKeysFromGPU.size(); i++) {
		if (startNewFlag == 1) {
			previousValues.clear();
			previousValues.push_back(extendValuesFromGPU[i]);
			startNewFlag = 0;
		} else {
			for (uint j = 0; j < previousValues.size(); j++) {
				EXPECT_TRUE(extendValuesFromGPU[i] != previousValues[j]);
			}
			previousValues.push_back(extendValuesFromGPU[i]);
		}

		if (i < extendedKeysFromGPU.size() - 1) {
			if (extendedKeysFromGPU[i] != extendedKeysFromGPU[i + 1]) {
				startNewFlag = 1;
			}
		}
	}

}

/**
 * This test function was setup to verify if all pair of nodes that are relatively close
 * can be found out using the find pair method.
 */ //
TEST_F(SceNodeTest, FindingPossiblePairsTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	// means that we have four nodes.
	SceNodes nodes = SceNodes(1, 1, 0, 1, 1, 1, 1, false);
	const uint testCellCount = 1;
	const uint testNodePerCell = 1;
	const uint testTotalNodeCount = 4;

	NodeAllocPara para = nodes.getAllocPara();
	para.maxNodeOfOneCell = testNodePerCell;
	para.currentActiveCellCount = testCellCount;
	nodes.setAllocPara(para);

	thrust::host_vector<double> nodeLocXHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocYHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocZHost(testTotalNodeCount);
	thrust::host_vector<bool> nodeIsActiveHost(testTotalNodeCount);
	thrust::host_vector<uint> nodeExpectedBucket(testTotalNodeCount);
	const double minX = -1.34;
	const double maxX = 2;
	const double minY = -0.5;
	const double maxY = 2.6;
	//const double minZ = 0.0;
	//const double maxZ = 0.0;
	const double bucketSize = 0.3;
	nodes.initDimension(minX, maxX, minY, maxY, bucketSize);
	nodeLocXHost[0] = 0.09;      // (0.09 - (-1.34)) / 0.3 = 4
	nodeLocYHost[0] = -0.09;     // (-0.09 - (-0.5)) / 0.3 = 1
	nodeIsActiveHost[0] = true;
	// 1 * 12 + 4 = 16
	nodeLocXHost[1] = 0.21;       // (0.21 - (-1.34)) / 0.3 = 5
	nodeLocYHost[1] = 0.2;        // (0.2 - (-0.5)) / 0.3 = 2
	nodeIsActiveHost[1] = true;
	// 2 * 12 + 5 = 29
	nodeLocXHost[2] = -0.1;      // (-0.1 - (-1.34)) / 0.3 = 4
	nodeLocYHost[2] = 0.5;      // (0.5 - (-0.5)) / 0.3 = 3
	nodeIsActiveHost[2] = true;
	// 3 * 12 + 4 = 40
	nodeLocXHost[3] = 0.63;      // (0.63 - (-1.34)) / 0.3 = 6
	nodeLocYHost[3] = 0.207;     // (0.207 - (-0.5)) / 0.3 = 2
	nodeIsActiveHost[3] = true;
	// 2 * 12 + 6 = 30

	// 2D bucket:
	// 3   4   5
	// 15  16  17  18  19
	// 27  28  29  30  31
	// 39  40  41  42  43
	// 51  52  53

	nodes.getInfoVecs().nodeLocX = nodeLocXHost;
	nodes.getInfoVecs().nodeLocY = nodeLocYHost;
	nodes.getInfoVecs().nodeLocZ = nodeLocZHost;
	nodes.getInfoVecs().nodeIsActive = nodeIsActiveHost;
	nodes.prepareSceForceComputation();
	//nodes.buildBuckets2D();
	//const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1; // (2 - (-1.34)) / 0.3 + 1 = 12
	//const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1; // (2.6 - (-0.5)) / 0.3 + 1 = 11
	//nodes.extendBuckets2D();

	thrust::host_vector<uint> extendedKeysFromGPU =
			nodes.getAuxVecs().bucketKeysExpanded;
	thrust::host_vector<uint> extendValuesFromGPU =
			nodes.getAuxVecs().bucketValuesIncludingNeighbor;
	EXPECT_EQ(extendedKeysFromGPU.size(), (uint )36);
	EXPECT_EQ(extendValuesFromGPU.size(), (uint )36);
	int expectedKeys[] = { 3, 4, 5, 15, 16, 16, 17, 17, 17, 18, 18, 19, 27, 27,
			28, 28, 28, 29, 29, 29, 29, 30, 30, 31, 39, 40, 40, 41, 41, 41, 42,
			42, 43, 51, 52, 53 };
	std::vector<uint> expectedResultsKeys(expectedKeys, expectedKeys + 36);
	for (int i = 0; i < 36; i++) {
		EXPECT_EQ(expectedResultsKeys[i], extendedKeysFromGPU[i]);
	}

	// keyBegin and keyEnd will be built in this step
	//nodes.applySceForces();

	// make sure that size matches our expectation
	std::vector<std::pair<uint, uint> > possiblePairs =
			nodes.obtainPossibleNeighborPairs();
	EXPECT_EQ((uint )3, possiblePairs.size());

	sort(possiblePairs.begin(), possiblePairs.end());
	std::vector<std::pair<uint, uint> > expectedPairs;
	expectedPairs.push_back(std::make_pair<uint, uint>(0, 1));
	expectedPairs.push_back(std::make_pair<uint, uint>(1, 2));
	expectedPairs.push_back(std::make_pair<uint, uint>(1, 3));

	// make sure that the possible pairs data matches our expectation.
	for (uint i = 0; i < possiblePairs.size(); i++) {
		EXPECT_EQ(possiblePairs[i], expectedPairs[i]);
	}
}

/**
 * This test function was setup to verify if output animation method was properly setup.
 */ //
TEST_F(SceNodeTest, outputAnimationLinks) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	// means that we have four nodes.
	SceNodes nodes = SceNodes(2, 0, 0, 0, 0, 2, 2, false);
	const uint testCellCount = 2;
	const uint testNodePerCell = 2;
	const uint testTotalNodeCount = 6;
	//nodes.setCurrentActiveCellCount(testCellCount);
	NodeAllocPara para = nodes.getAllocPara();
	para.maxNodeOfOneCell = testNodePerCell;
	para.currentActiveCellCount = testCellCount;
	nodes.setAllocPara(para);
	thrust::host_vector<double> nodeLocXHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocYHost(testTotalNodeCount);
	thrust::host_vector<double> nodeLocZHost(testTotalNodeCount);
	thrust::host_vector<bool> nodeIsActiveHost(testTotalNodeCount);
	thrust::host_vector<uint> nodeExpectedBucket(testTotalNodeCount);
	const double minX = -1.34;
	const double maxX = 2;
	const double minY = -0.5;
	const double maxY = 2.6;
	//const double minZ = 0.0;
	//const double maxZ = 0.0;
	const double bucketSize = 0.3;
	nodes.initDimension(minX, maxX, minY, maxY, bucketSize);
	nodeLocXHost[0] = -0.3;      // (-0.3 - (-1.34)) / 0.3 = 3
	nodeLocYHost[0] = -0.39;     // (-0.39 - (-0.5)) / 0.3 = 0
	nodeIsActiveHost[0] = true;
	// 0 * 12 + 3 = 3
	nodeExpectedBucket[0] = calculateBucketKey(minX, maxX, minY, maxY,
			bucketSize, nodeLocXHost[0], nodeLocYHost[0]);

	nodeLocXHost[1] = 0.0;       // (0.0 - (-1.34)) / 0.3 = 4
	nodeLocYHost[1] = -0.3;      // (-0.3 - (-0.5)) / 0.3 = 0
	nodeIsActiveHost[1] = true;
	// 0 * 12 + 4 = 4
	nodeExpectedBucket[1] = calculateBucketKey(minX, maxX, minY, maxY,
			bucketSize, nodeLocXHost[1], nodeLocYHost[1]);

	nodeLocXHost[2] = 0.3;       // (0.3 - (-1.34)) / 0.3 = 5
	nodeLocYHost[2] = -0.09;     // (-0.09 - (-0.5)) / 0.3 = 1
	nodeIsActiveHost[2] = true;
	// 1 * 12 + 5 = 17
	nodeExpectedBucket[2] = calculateBucketKey(minX, maxX, minY, maxY,
			bucketSize, nodeLocXHost[2], nodeLocYHost[2]);

	nodeLocXHost[3] = 0.21;      // (0.21 - (-1.34)) / 0.3 = 5
	nodeLocYHost[3] = -0.1;      // (-0.1 - (-0.5)) / 0.3 = 1
	nodeIsActiveHost[3] = true;
	// 1 * 12 + 5 = 17
	nodeExpectedBucket[3] = calculateBucketKey(minX, maxX, minY, maxY,
			bucketSize, nodeLocXHost[3], nodeLocYHost[3]);

	nodeLocXHost[4] = 0.82;       // (0.82 - (-1.34)) / 0.3 = 7
	nodeLocYHost[4] = 0.51;      // (0.51 - (-0.5)) / 0.3 = 3
	nodeIsActiveHost[4] = false; // intentionally mark it inactive
	// 3 * 12 + 7 = 43
	nodeExpectedBucket[4] = calculateBucketKey(minX, maxX, minY, maxY,
			bucketSize, nodeLocXHost[4], nodeLocYHost[4]);

	nodeLocXHost[5] = 0.5;       // (0.5 - (-1.34)) / 0.3 = 6
	nodeLocYHost[5] = 0.25;      // (0.25 - (-0.5)) / 0.3 = 2
	nodeIsActiveHost[5] = true;
	// 2 * 12 + 6 = 30
	nodeExpectedBucket[5] = calculateBucketKey(minX, maxX, minY, maxY,
			bucketSize, nodeLocXHost[5], nodeLocYHost[5]);

	//const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1; // (2 - (-1.34)) / 0.3 + 1 = 12
	//const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1; // (2.6 - (-0.5)) / 0.3 + 1 = 11

	nodes.getInfoVecs().nodeLocX = nodeLocXHost;
	nodes.getInfoVecs().nodeLocY = nodeLocYHost;
	nodes.getInfoVecs().nodeLocZ = nodeLocZHost;
	nodes.getInfoVecs().nodeIsActive = nodeIsActiveHost;
	//nodes.buildBuckets2D();
	//nodes.extendBuckets2D();
	// This step is necessary for
	//nodes.applySceForces();
	nodes.prepareSceForceComputation();
	thrust::host_vector<uint> keysFromGPU = nodes.getAuxVecs().bucketKeys;
	thrust::host_vector<uint> valuesFromGPU = nodes.getAuxVecs().bucketValues;
	uint activeNodeCount = 0;
	for (uint i = 0; i < nodeIsActiveHost.size(); i++) {
		if (nodeIsActiveHost[i] == true) {
			activeNodeCount++;
		}
	}
	EXPECT_EQ(keysFromGPU.size(), activeNodeCount);
	EXPECT_EQ(keysFromGPU.size(), valuesFromGPU.size());
	for (uint i = 0; i < keysFromGPU.size(); i++) {
		uint nodeRank = valuesFromGPU[i];
		uint expectedResultFromCPU = nodeExpectedBucket[nodeRank];
		EXPECT_EQ(expectedResultFromCPU, keysFromGPU[i]);
	}

	// make sure that size matches our expectation
	std::vector<std::pair<uint, uint> > possiblePairs =
			nodes.obtainPossibleNeighborPairs();
	EXPECT_EQ((uint )6, possiblePairs.size());

	sort(possiblePairs.begin(), possiblePairs.end());
	std::vector<std::pair<uint, uint> > expectedPairs;
	expectedPairs.push_back(std::make_pair<uint, uint>(0, 1));
	expectedPairs.push_back(std::make_pair<uint, uint>(1, 2));
	expectedPairs.push_back(std::make_pair<uint, uint>(1, 3));
	expectedPairs.push_back(std::make_pair<uint, uint>(2, 3));
	expectedPairs.push_back(std::make_pair<uint, uint>(2, 5));
	expectedPairs.push_back(std::make_pair<uint, uint>(3, 5));

	// make sure that the possible pairs data matches our expectation.
	for (uint i = 0; i < possiblePairs.size(); i++) {
		EXPECT_EQ(possiblePairs[i], expectedPairs[i]);
	}

	AnimationCriteria aniCri;
	aniCri.defaultEffectiveDistance = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();
	aniCri.isStressMap = false;
	VtkAnimationData vtkData = nodes.obtainAnimationData(aniCri);

	EXPECT_EQ((uint )4, vtkData.pointsAniData.size());
	EXPECT_EQ((uint )2, vtkData.linksAniData.size());

	std::vector<PointAniData> expectedAniPoints;
	PointAniData ptData;
	ptData.pos = CVector(nodeLocXHost[0], nodeLocYHost[0], 0.0);
	ptData.colorScale = nodeTypeToScale(Boundary);
	expectedAniPoints.push_back(ptData);
	ptData.pos = CVector(nodeLocXHost[1], nodeLocYHost[1], 0.0);
	ptData.colorScale = nodeTypeToScale(Boundary);
	expectedAniPoints.push_back(ptData);
	ptData.pos = CVector(nodeLocXHost[2], nodeLocYHost[2], 0.0);
	ptData.colorScale = nodeTypeToScale(FNM);
	expectedAniPoints.push_back(ptData);
	ptData.pos = CVector(nodeLocXHost[3], nodeLocYHost[3], 0.0);
	ptData.colorScale = nodeTypeToScale(FNM);
	expectedAniPoints.push_back(ptData);

	std::vector<LinkAniData> expectedAniPairs;
	LinkAniData linkData;
	linkData.node1Index = 0;
	linkData.node2Index = 1;
	expectedAniPairs.push_back(linkData);
	linkData.node1Index = 2;
	linkData.node2Index = 3;
	expectedAniPairs.push_back(linkData);

	for (uint i = 0; i < expectedAniPairs.size(); i++) {
		EXPECT_EQ(expectedAniPairs[i].node1Index,
				vtkData.linksAniData[i].node1Index);
		EXPECT_EQ(expectedAniPairs[i].node2Index,
				vtkData.linksAniData[i].node2Index);
	}

	for (uint i = 0; i < expectedAniPoints.size(); i++) {
		EXPECT_NEAR(expectedAniPoints[i].pos.x, vtkData.pointsAniData[i].pos.x,
				errTol);
		EXPECT_NEAR(expectedAniPoints[i].pos.y, vtkData.pointsAniData[i].pos.y,
				errTol);
		EXPECT_NEAR(expectedAniPoints[i].pos.z, vtkData.pointsAniData[i].pos.z,
				errTol);
		EXPECT_NEAR(expectedAniPoints[i].colorScale,
				vtkData.pointsAniData[i].colorScale, errTol);
	}
}

