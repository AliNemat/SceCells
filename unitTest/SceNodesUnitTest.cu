#include <iostream>
#include "gtest/gtest.h"
#include "SceNodes.h"
#include <vector>
#include <cuda_runtime.h>
#include <algorithm>
using namespace std;

/**
 * The CPU parameters are declared as extern while GPU parameters are locally effective.
 * Reason is that CUDA const memory cannot be shared between files, at least before CUDA 5.0
 */

__constant__ double sceInterPara[5];
__constant__ double sceIntraPara[4];
__constant__ double sceInterDiffPara[5];
__constant__ double sceProfilePara[7];
__constant__ double sceECMPara[5];
__constant__ double sceDiffPara[5];

extern double sceInterParaCPU[5];
extern double sceIntraParaCPU[4];
extern double sceInterDiffParaCPU[5];
extern double sceProfileParaCPU[7];
extern double sceECMParaCPU[5];
extern double sceDiffParaCPU[5];

__constant__ uint ProfilebeginPos;
__constant__ uint ECMbeginPos;
__constant__ uint cellNodeBeginPos;
__constant__ uint nodeCountPerECM;
__constant__ uint nodeCountPerCell;

const int myDeviceId = 0;
const uint maxCellCount = 10;
const uint maxNodePerCell = 100;
const uint maxNodeCount = maxCellCount * maxNodePerCell;
const double dt = 0.1;
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
	//double xMax = *std::max_element(xPoss.begin(), xPoss.end());
	//double xMin = *std::min_element(xPoss.begin(), xPoss.end());
	//double yMax = *std::max_element(yPoss.begin(), yPoss.end());
	//double yMin = *std::min_element(yPoss.begin(), yPoss.end());
	//uint xBucketCount = (xMax - xMin) / bucketSize + 1;
	//uint yBuckeyCount = (yMax - yMin) / bucketSize + 1;
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

class SceNodeTest: public ::testing::Test {
protected:
	double sceInterParaCPU[5];
	double sceIntraParaCPU[4];
	virtual void SetUp() {
		ConfigParser parser;
		std::string configFileName = "../sceCell.config";
		globalConfigVars = parser.parseConfigFile(configFileName);
		cudaSetDevice(
				globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
		static const double U0 =
				globalConfigVars.getConfigValue("InterCell_U0_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"InterCell_U0_DivFactor").toDouble();
		static const double V0 =
				globalConfigVars.getConfigValue("InterCell_V0_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"InterCell_V0_DivFactor").toDouble();
		static const double k1 =
				globalConfigVars.getConfigValue("InterCell_k1_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"InterCell_k1_DivFactor").toDouble();
		static const double k2 =
				globalConfigVars.getConfigValue("InterCell_k2_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"InterCell_k2_DivFactor").toDouble();
		static const double interLinkEffectiveRange =
				globalConfigVars.getConfigValue("InterCellLinkBreakRange").toDouble();

		sceInterParaCPU[0] = U0;
		sceInterParaCPU[1] = V0;
		sceInterParaCPU[2] = k1;
		sceInterParaCPU[3] = k2;
		sceInterParaCPU[4] = interLinkEffectiveRange;

		static const double U0_Intra =
				globalConfigVars.getConfigValue("IntraCell_U0_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_U0_DivFactor").toDouble();
		static const double V0_Intra =
				globalConfigVars.getConfigValue("IntraCell_V0_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_V0_DivFactor").toDouble();
		static const double k1_Intra =
				globalConfigVars.getConfigValue("IntraCell_k1_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_k1_DivFactor").toDouble();
		static const double k2_Intra =
				globalConfigVars.getConfigValue("IntraCell_k2_Original").toDouble()
						/ globalConfigVars.getConfigValue(
								"IntraCell_k2_DivFactor").toDouble();
		sceIntraParaCPU[0] = U0_Intra;
		sceIntraParaCPU[1] = V0_Intra;
		sceIntraParaCPU[2] = k1_Intra;
		sceIntraParaCPU[3] = k2_Intra;

		//std::cout << "in SceNodes, before cuda memory copy to symbol:" << std::endl;
		cudaMemcpyToSymbol(sceInterPara, sceInterParaCPU, 5 * sizeof(double));
		cudaMemcpyToSymbol(sceIntraPara, sceIntraParaCPU, 4 * sizeof(double));

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
						/ k2_Bdry
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

		cudaMemcpyToSymbol(sceProfilePara, sceProfileParaCPU,
				7 * sizeof(double));

		cudaMemcpyToSymbol(sceInterDiffPara, sceInterDiffParaCPU,
				5 * sizeof(double));

		cudaMemcpyToSymbol(sceECMPara, sceECMParaCPU, 5 * sizeof(double));
		//std::cout << "finished SceNodes:" << std::endl;
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

TEST_F(SceNodeTest, ParameterInitTest) {
        cout<<" point 1 , before everything starts"<<endl;
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
        cout<<" point 1.1 "<<endl;
        
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
        cout<<" point 1.2 "<<endl;
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);

	// following part tests parameters of the neighbor grid.
	int Expected_BucketXDim = (Test_maxX - Test_minX) / Test_bucketSize + 1;
	int Expected_BucketYDim = (Test_maxY - Test_minY) / Test_bucketSize + 1;
	int Expected_TotalBucket = Expected_BucketXDim * Expected_BucketYDim;
	EXPECT_EQ(Expected_BucketXDim, nodes.numOfBucketsInXDim);
	EXPECT_EQ(Expected_BucketYDim, nodes.numOfBucketsInYDim);
	EXPECT_EQ(Expected_TotalBucket, nodes.totalBucketCount);
        cout<<" point 2, middle of no where"<<endl;

	// following part tests initial parameters for node location info
	int startPosOfProfile = Test_totalBdryNodeCount;
	int startPosOfECM = Test_maxProfileNodeCount + startPosOfProfile;
	int startPosOfCells = Test_maxTotalECMCount * Test_maxNodeInECM
			+ startPosOfECM;
	EXPECT_EQ(startPosOfProfile, nodes.startPosProfile);
	EXPECT_EQ(startPosOfECM, nodes.startPosECM);
	EXPECT_EQ(startPosOfCells, nodes.startPosCells);
	EXPECT_EQ(Test_totalBdryNodeCount, nodes.BdryNodeCount);
	EXPECT_EQ(Test_maxTotalECMCount, nodes.maxECMCount);
	EXPECT_EQ(Test_maxTotalCellCount, nodes.maxCellCount);
	EXPECT_EQ(Test_maxProfileNodeCount, nodes.maxProfileNodeCount);
	EXPECT_EQ(Test_maxNodeInCell, nodes.maxNodeOfOneCell);
	EXPECT_EQ(Test_maxNodeInECM, nodes.maxNodePerECM);

        cout<<" point 3, everything has finished"<<endl;
}

TEST_F(SceNodeTest, GPUConstMemTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
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

	double sceInterParaFromGPU[5];
	//std::cout << "in SceNodes, before cuda memory copy to symbol:" << std::endl;
	cudaMemcpyFromSymbol(sceInterParaFromGPU, sceInterPara, 5 * sizeof(double));
	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(sceInterParaFromGPU[i], sceInterParaCPU[i], errTol);
	}
	double sceIntraParaFromGPU[4];
	cudaMemcpyFromSymbol(sceIntraParaFromGPU, sceIntraPara, 4 * sizeof(double));
	for (int i = 0; i < 4; i++) {
		EXPECT_NEAR(sceIntraParaFromGPU[i], sceIntraParaCPU[i], errTol);
	}

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
	static const double neutralLength = globalConfigVars.getConfigValue(
			"Bdry_base_neutral_dist").toDouble() / k2_Bdry
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

	double sceProfileParaFromGPU[7];
	cudaMemcpyFromSymbol(sceProfileParaFromGPU, sceProfilePara,
			7 * sizeof(double));
	for (int i = 0; i < 7; i++) {
		EXPECT_NEAR(sceProfileParaFromGPU[i], sceProfileParaCPU[i], errTol);
	}

	double sceInterDiffParaFromGPU[5];
	cudaMemcpyFromSymbol(sceInterDiffParaFromGPU, sceInterDiffPara,
			5 * sizeof(double));
	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(sceInterDiffParaFromGPU[i], sceInterDiffParaFromGPU[i],
				errTol);
	}

	double sceECMParaFromGPU[5];
	cudaMemcpyFromSymbol(sceECMParaFromGPU, sceECMPara, 5 * sizeof(double));
	for (int i = 0; i < 5; i++) {
		EXPECT_NEAR(sceECMParaCPU[i], sceECMParaFromGPU[i], errTol);
	}
}

TEST_F(SceNodeTest, MemSizeTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);
	int totalNodeSize = Test_totalBdryNodeCount + Test_maxProfileNodeCount
			+ Test_maxTotalECMCount * Test_maxNodeInECM
			+ Test_maxTotalCellCount * Test_maxNodeInCell;

	// size of buckets should be same with keyBeing and KeyEnd as they are
	// recording information for buckets
	EXPECT_EQ(nodes.totalBucketCount, nodes.keyBegin.size());
	EXPECT_EQ(nodes.totalBucketCount, nodes.keyEnd.size());

	// size of these node information should be the same with max node number.
	EXPECT_EQ(totalNodeSize, nodes.nodeCellRank.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeCellType.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeIsActive.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeLocX.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeLocY.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeLocZ.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeVelX.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeVelY.size());
	EXPECT_EQ(totalNodeSize, nodes.nodeVelZ.size());
}

void generateInitVectors() {

}

TEST_F(SceNodeTest, InitialValueTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
	nodes.initDimension(Test_minX, Test_maxX, Test_minY, Test_maxY,
			Test_bucketSize);
	//nodes.initValues();
}

int calculateBucketKey(double minX, double maxX, double minY, double maxY,
		double bucketSize, double xcoord, double ycoord) {
	uint width = (maxX - minX) / bucketSize + 1;
	unsigned int x = (unsigned int) ((xcoord - minX) / bucketSize);
	unsigned int y = (unsigned int) ((ycoord - minY) / bucketSize);
	return (y * width + x);
}

TEST_F(SceNodeTest,BuildBucketFixedTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
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

	nodes.nodeLocX = nodeLocXHost;
	nodes.nodeLocY = nodeLocYHost;

	nodes.nodeIsActive = nodeIsActiveHost;
	nodes.startPosCells = 0;
	nodes.currentActiveCellCount = 2;
	nodes.maxNodeOfOneCell = 2;
	nodes.buildBuckets2D();
	thrust::host_vector<uint> keysFromGPU = nodes.bucketKeys;
	thrust::host_vector<uint> valuesFromGPU = nodes.bucketValues;
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
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
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
	nodes.nodeLocX = nodeLocXHost;
	nodes.nodeLocY = nodeLocYHost;
	nodes.nodeLocZ = nodeLocZHost;
	nodes.nodeIsActive = nodeIsActiveHost;
	nodes.startPosCells = maxNodeCount;
	nodes.currentActiveCellCount = 0;
	nodes.maxNodeOfOneCell = 0;
	nodes.buildBuckets2D();
	thrust::host_vector<uint> keysFromGPU = nodes.bucketKeys;
	thrust::host_vector<uint> valuesFromGPU = nodes.bucketValues;
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
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
	const uint testCellCount = 2;
	const uint testNodePerCell = 2;
	const uint testTotalNodeCount = testCellCount * testNodePerCell;
	nodes.setCurrentActiveCellCount(testCellCount);
	nodes.maxNodeOfOneCell = testCellCount;
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
	nodes.nodeLocX = nodeLocXHost;
	nodes.nodeLocY = nodeLocYHost;
	nodes.nodeLocZ = nodeLocZHost;
	nodes.nodeIsActive = nodeIsActiveHost;
	nodes.buildBuckets2D();
	const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1;
	const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1;
	nodes.extendBuckets2D();

	thrust::host_vector<uint> extendedKeysFromGPU = nodes.bucketKeysExpanded;
	thrust::host_vector<uint> extendValuesFromGPU =
			nodes.bucketValuesIncludingNeighbor;
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
 * no duplicate */

TEST(SceExtendBucket2D, extendBucketRandomTest) {
	cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
	SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
			Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
			Test_maxNodeInCell);
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
	nodes.startPosCells = maxNodeCount;
	nodes.currentActiveCellCount = 0;
	nodes.maxNodeOfOneCell = 0;
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
	nodes.nodeLocX = nodeLocXHost;
	nodes.nodeLocY = nodeLocYHost;
	nodes.nodeLocZ = nodeLocZHost;
	nodes.nodeIsActive = nodeIsActiveHost;
	const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1;
	const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1;
	nodes.buildBuckets2D();
	nodes.extendBuckets2D();

	thrust::host_vector<uint> extendedKeysFromGPU = nodes.bucketKeysExpanded;
	thrust::host_vector<uint> extendValuesFromGPU =
			nodes.bucketValuesIncludingNeighbor;
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
/*
 TEST_F(SceNodeTest, addForceFixedNeighborTest) {
 cudaSetDevice(globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
 SceNodes nodes = SceNodes(Test_totalBdryNodeCount, Test_maxProfileNodeCount,
 Test_maxTotalECMCount, Test_maxNodeInECM, Test_maxTotalCellCount,
 Test_maxNodeInCell);
 const uint testCellCount = 1;
 const uint testNodePerCell = 4;
 const uint testTotalNodeCount = testCellCount * testNodePerCell;

 thrust::host_vector<double> nodeLocXHost(testTotalNodeCount);
 thrust::host_vector<double> nodeLocYHost(testTotalNodeCount);
 thrust::host_vector<double> nodeLocZHost(testTotalNodeCount);
 thrust::host_vector<bool> nodeIsActiveHost(testTotalNodeCount);
 const double minX = 0.0;
 const double maxX = 3.0 - 1.0e-10;
 const double minY = 0.0;
 const double maxY = 2.0 - 1.0e-10;
 //const double minZ = 0.0;
 //const double maxZ = 0.0;
 const double bucketSize = 1.0;
 nodes.initDimension(minX, maxX, minY, maxY, bucketSize);

 nodeLocXHost[0] = 0.2;
 nodeLocYHost[0] = 0.5;
 nodeLocZHost[0] = 0.0;
 nodeIsActiveHost[0] = 1;
 // 0
 nodeLocXHost[1] = 1.2;
 nodeLocYHost[1] = 0.2;
 nodeLocZHost[1] = 0.0;
 nodeIsActiveHost[1] = 1;
 // 1
 nodeLocXHost[2] = 1.3;
 nodeLocYHost[2] = 0.5;
 nodeLocZHost[2] = 0.0;
 nodeIsActiveHost[2] = 1;
 // 1
 nodeLocXHost[3] = 2.7;
 nodeLocYHost[3] = 1.1;
 nodeLocZHost[3] = 0.0;
 nodeIsActiveHost[3] = 1;
 // 5
 nodes.nodeLocX = nodeLocXHost;
 nodes.nodeLocY = nodeLocYHost;
 nodes.nodeLocZ = nodeLocZHost;
 nodes.nodeIsActive = nodeIsActiveHost;
 nodes.startPosCells = 0;
 nodes.currentActiveCellCount = testCellCount;
 nodes.maxNodeOfOneCell = testNodePerCell;
 //nodes.buildBuckets2D(minX, maxX, minY, maxY, bucketSize);
 nodes.calculateAndApplySceForces();
 //const int numberOfBucketsInXDim = (maxX - minX) / bucketSize + 1;
 //const int numberOfBucketsInYDim = (maxY - minY) / bucketSize + 1;
 //nodes.extendBuckets2D(numberOfBucketsInXDim, numberOfBucketsInYDim);
 //std::cout << "before applying forces:" << std::endl;
 //thrust::host_vector<double> nodeVelXFromGPU_init = nodes.nodeVelX;
 //thrust::host_vector<double> nodeVelYFromGPU_init = nodes.nodeVelY;
 //thrust::host_vector<double> nodeVelZFromGPU_init = nodes.nodeVelZ;
 //for (uint i = 0; i < nodeVelXFromGPU_init.size(); i++) {
 //	std::cout << nodeVelXFromGPU_init[i] << ", " << nodeVelYFromGPU_init[i]
 //			<< ", " << nodeVelZFromGPU_init[i] << std::endl;
 //}
 //std::cout << std::endl;
 //nodes.applySceForces(numberOfBucketsInXDim, numberOfBucketsInYDim);
 //thrust::host_vector<uint> bucketsKeysFromGPU = nodes.bucketKeys;
 //thrust::host_vector<uint> bucketsValuesFromGPU = nodes.bucketValues;
 //std::cout << "printing key-value pairs:" << std::endl;
 //for (uint i = 0; i < bucketsKeysFromGPU.size(); i++) {
 //	std::cout << "Key :" << bucketsKeysFromGPU[i] << ", value: "
 //			<< bucketsValuesFromGPU[i] << std::endl;
 //}
 thrust::host_vector<double> nodeVelXFromGPU = nodes.nodeVelX;
 thrust::host_vector<double> nodeVelYFromGPU = nodes.nodeVelY;
 thrust::host_vector<double> nodeVelZFromGPU = nodes.nodeVelZ;
 for (uint i = 0; i < nodeVelXFromGPU.size(); i++) {
 //std::cout << nodeVelXFromGPU[i] << ", " << nodeVelYFromGPU[i] << ", "
 //		<< nodeVelZFromGPU[i] << std::endl;
 }

 vector<double> xPoss(testTotalNodeCount, 0.0);
 vector<double> yPoss(testTotalNodeCount, 0.0);
 vector<double> zPoss(testTotalNodeCount, 0.0);
 vector<double> xVels(testTotalNodeCount, 0.0);
 vector<double> yVels(testTotalNodeCount, 0.0);
 vector<double> zVels(testTotalNodeCount, 0.0);
 vector<bool> isActive(testTotalNodeCount, 0);
 for (uint i = 0; i < testTotalNodeCount; i++) {
 xPoss[i] = nodeLocXHost[i];
 yPoss[i] = nodeLocYHost[i];
 zPoss[i] = nodeLocZHost[i];
 isActive[i] = nodeIsActiveHost[i];
 }
 vector<double> paraSetIntra(4, 0.0);
 for (uint i = 0; i < 4; i++) {
 paraSetIntra[i] = sceIntraParaCPU[i];
 }
 computeResultFromCPUAllIntra2D(xPoss, yPoss, zPoss, xVels, yVels, zVels,
 isActive, paraSetIntra, bucketSize, minX, maxX, minY, maxY);
 for (uint i = 0; i < nodeVelXFromGPU.size(); i++) {
 EXPECT_NEAR(xVels[i], nodeVelXFromGPU[i], errTol);
 EXPECT_NEAR(yVels[i], nodeVelYFromGPU[i], errTol);
 EXPECT_NEAR(zVels[i], nodeVelZFromGPU[i], errTol);
 }

 //std::cout << std::endl;
 }
 */
/*
 TEST_F(SceNodeTest, addForceRandomTest) {
 cudaSetDevice(myDeviceId);
 SceNodes nodes(maxCellCount, maxNodePerCell);
 int currentActiveCellCount = maxCellCount - 4;
 nodes.setCurrentActiveCellCount(currentActiveCellCount);
 thrust::host_vector<double> nodeLocXHost(maxNodeCount);
 thrust::host_vector<double> nodeLocYHost(maxNodeCount);
 thrust::host_vector<double> nodeLocZHost(maxNodeCount);
 thrust::host_vector<bool> nodeIsActiveHost(maxNodeCount);
 const double minX = 0.9;
 const double maxX = 2.4;
 const double minY = 0.8;
 const double maxY = 3.7;
 const double minZ = 0.0;
 const double maxZ = 0.0;
 const double bucketSize = 0.15;
 thrust::counting_iterator<unsigned int> index_sequence_begin(0);
 thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
 nodeLocXHost.begin(), Prg(minX, maxX));
 thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
 nodeLocYHost.begin(), Prg(minY, maxY));
 thrust::transform(index_sequence_begin, index_sequence_begin + maxNodeCount,
 nodeLocZHost.begin(), Prg(minZ, maxZ));
 for (uint i = 0; i < maxNodeCount; i++) {
 if (i % maxNodePerCell < maxNodePerCell / 2) {
 nodeIsActiveHost[i] = true;
 } else {
 nodeIsActiveHost[i] = false;
 }
 }
 nodes.nodeLocX = nodeLocXHost;
 nodes.nodeLocY = nodeLocYHost;
 nodes.nodeLocZ = nodeLocZHost;
 nodes.nodeIsActive = nodeIsActiveHost;
 nodes.calculateAndApplySceForces(minX, maxX, minY, maxY, bucketSize);

 thrust::host_vector<double> nodeVelXFromGPU = nodes.nodeVelX;
 thrust::host_vector<double> nodeVelYFromGPU = nodes.nodeVelY;
 thrust::host_vector<double> nodeVelZFromGPU = nodes.nodeVelZ;

 vector<double> xPoss(maxNodeCount, 0.0);
 vector<double> yPoss(maxNodeCount, 0.0);
 vector<double> zPoss(maxNodeCount, 0.0);
 vector<double> xVels(maxNodeCount, 0.0);
 vector<double> yVels(maxNodeCount, 0.0);
 vector<double> zVels(maxNodeCount, 0.0);
 vector<bool> isActive(maxNodeCount, 0.0);

 for (uint i = 0; i < maxNodeCount; i++) {
 xPoss[i] = nodeLocXHost[i];
 yPoss[i] = nodeLocYHost[i];
 zPoss[i] = nodeLocZHost[i];
 isActive[i] = nodeIsActiveHost[i];
 }

 vector<double> paraSetIntra(4, 0.0);
 for (uint i = 0; i < 4; i++) {
 paraSetIntra[i] = sceIntraParaCPU[i];
 }
 vector<double> paraSetInter(5, 0.0);
 for (uint i = 0; i < 5; i++) {
 paraSetInter[i] = sceInterParaCPU[i];
 }
 computeResultFromCPUAllIntraAndInter2D(xPoss, yPoss, zPoss, xVels, yVels,
 zVels, isActive, paraSetIntra, paraSetInter, bucketSize,
 currentActiveCellCount, maxNodePerCell, minX, maxX, minY, maxY);

 for (uint i = 0; i < currentActiveCellCount * maxNodePerCell; i++) {
 //std::cout << "xVel expected = " << xVels[i] << " xVel from GPU = "
 //		<< nodeVelXFromGPU[i] << std::endl;
 EXPECT_NEAR(xVels[i], nodeVelXFromGPU[i], errTol);
 EXPECT_NEAR(yVels[i], nodeVelYFromGPU[i], errTol);
 EXPECT_NEAR(zVels[i], nodeVelZFromGPU[i], errTol);
 }
 }

 */

/*
 TEST(AddNewCellTest, addCellFixedTest) {
 cudaSetDevice(myDeviceId);
 const uint maxCellCount = 4;
 const uint initCellCount = 2;
 const uint maxNodePerCell = 2;
 const uint maxECMCount = 2;
 const uint maxNodeInECM = 1;
 const uint addNodeSize = 4;
 const uint maxTotalNodeCount = maxCellCount * maxNodePerCell
 + maxECMCount * maxNodeInECM;
 SceNodes nodes(maxCellCount, maxNodePerCell, maxECMCount, maxNodeInECM);
 nodes.setCurrentActiveCellCount(initCellCount);
 nodes.setCurrentActiveEcm(maxNodeInECM);

 double nodeXInput[] = { 1.2, 3, 2, 1.5, 0.3, 1.1, 9.9, 0.0, 0.0, 0.0 };
 double nodeYInput[] = { 2.3, 1, 2, 5.6, 0.9, 8.6, 2.3, 0.0, 0.0, 0.0 };
 double nodeZInput[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 bool nodeIsActInput[] = { true, true, true, false, true, true, true, true,
 true, false };

 thrust::host_vector<double> nodeLocXHost(nodeXInput,
 nodeXInput + maxTotalNodeCount);
 thrust::host_vector<double> nodeLocYHost(nodeYInput,
 nodeYInput + maxTotalNodeCount);
 thrust::host_vector<double> nodeLocZHost(nodeZInput,
 nodeZInput + maxTotalNodeCount);
 thrust::host_vector<bool> nodeIsActiveHost(nodeIsActInput,
 nodeIsActInput + maxTotalNodeCount);
 nodes.nodeLocX = nodeLocXHost;
 nodes.nodeLocY = nodeLocYHost;
 nodes.nodeLocZ = nodeLocZHost;
 nodes.nodeIsActive = nodeIsActiveHost;

 double nodeXInputToAdd[] = { 8.2, 7.4, 9.5, 6.7 };
 double nodeYInputToAdd[] = { 3.2, 5.6, 8.8, 7.7 };
 double nodeZInputToAdd[] = { 0, 0, 0, 0 };
 bool nodeIsActInputToAdd[] = { true, false, true, false };

 thrust::host_vector<double> locXToBeAddedHost(nodeXInputToAdd,
 nodeXInputToAdd + addNodeSize);
 thrust::host_vector<double> locYToBeAddedHost(nodeYInputToAdd,
 nodeYInputToAdd + addNodeSize);
 thrust::host_vector<double> locZToBeAddedHost(nodeZInputToAdd,
 nodeZInputToAdd + addNodeSize);
 thrust::host_vector<bool> locIsActiveToBeAddedHost(nodeIsActInputToAdd,
 nodeIsActInputToAdd + addNodeSize);

 thrust::device_vector<double> locXToBeAdded = locXToBeAddedHost;
 thrust::device_vector<double> locYToBeAdded = locYToBeAddedHost;
 thrust::device_vector<double> locZToBeAdded = locZToBeAddedHost;
 thrust::device_vector<bool> locIsActiveToBeAdded = locIsActiveToBeAddedHost;

 double nodeXExpecteRes[] = { 1.2, 3, 2, 1.5, 8.2, 7.4, 9.5, 6.7, 0.3, 1.1 };
 double nodeYExpecteRes[] = { 2.3, 1, 2, 5.6, 3.2, 5.6, 8.8, 7.7, 0.9, 8.6 };
 double nodeZExpecteRes[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
 bool nodeIsActExpecteRes[] = { true, true, true, false, true, false, true,
 false, true, true };

 nodes.addNewlyDividedCells(locXToBeAdded, locYToBeAdded, locZToBeAdded,
 locIsActiveToBeAdded);
 thrust::host_vector<double> locXAfterExpandFromGPU = nodes.nodeLocX;
 thrust::host_vector<double> locYAfterExpandFromGPU = nodes.nodeLocY;
 thrust::host_vector<double> locZAfterExpandFromGPU = nodes.nodeLocZ;
 thrust::host_vector<double> isActiveAfterExpandFromGPU = nodes.nodeIsActive;
 EXPECT_EQ(locXAfterExpandFromGPU.size(), maxTotalNodeCount);
 EXPECT_EQ(locYAfterExpandFromGPU.size(), maxTotalNodeCount);
 EXPECT_EQ(locZAfterExpandFromGPU.size(), maxTotalNodeCount);
 EXPECT_EQ(isActiveAfterExpandFromGPU.size(), maxTotalNodeCount);
 for (uint i = 0; i < maxTotalNodeCount; i++) {
 EXPECT_EQ(locXAfterExpandFromGPU[i], nodeXExpecteRes[i]);
 EXPECT_EQ(locYAfterExpandFromGPU[i], nodeYExpecteRes[i]);
 EXPECT_EQ(locZAfterExpandFromGPU[i], nodeZExpecteRes[i]);
 EXPECT_EQ(isActiveAfterExpandFromGPU[i], nodeIsActExpecteRes[i]);
 }
 }
 */

