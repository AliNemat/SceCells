#include <iostream>
#include "gtest/gtest.h"
#include "SceNodes.h"
#include "SceCells.h"
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

//const int myDeviceId = 0;
//const uint maxCellCount = 10;
//const uint maxNodePerCell = 100;
//const uint maxNodeCount = maxCellCount * maxNodePerCell;
//const double dt = 0.1;
//const double errTol = 1.0e-12;

extern GlobalConfigVars globalConfigVars;

class SingleCellTest: public ::testing::Test {
protected:
	virtual void SetUp() {
	}
};

/*
 TEST(SceCellsInitTest, sizeTest) {
 cudaSetDevice(myDeviceId);
 SceNodes nodes(maxCellCount, maxNodePerCell);
 SceCells cells(&nodes);

 EXPECT_EQ(cells.growthProgress.size(), maxCellCount);
 EXPECT_EQ(cells.activeNodeCountOfThisCell.size(), maxCellCount);
 EXPECT_EQ(cells.lastCheckPoint.size(), maxCellCount);
 EXPECT_EQ(cells.isDivided.size(), maxCellCount);
 EXPECT_EQ(cells.centerCoordX.size(), maxCellCount);
 EXPECT_EQ(cells.centerCoordY.size(), maxCellCount);
 EXPECT_EQ(cells.centerCoordZ.size(), maxCellCount);

 //EXPECT_EQ(cells.xCoordTmp.size(), maxNodeCount);
 //EXPECT_EQ(cells.yCoordTmp.size(), maxNodeCount);
 //EXPECT_EQ(cells.zCoordTmp.size(), maxNodeCount);
 }

 TEST(SceCellsDistriIsActiveInfoTest,fixedTest) {
 cudaSetDevice(myDeviceId);
 const uint maxCellCount = 2;
 const uint initCellCount = 2;
 const uint maxNodePerCell = 4;
 const uint maxECMCount = 2;
 const uint maxNodeInECM = 1;
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

 SceCells cells(&nodes);
 thrust::host_vector<uint> activeNodeCount(2);
 activeNodeCount[0] = 4;
 activeNodeCount[1] = 3;
 cells.activeNodeCountOfThisCell = activeNodeCount;
 cells.distributeIsActiveInfo();
 bool expectedNodeIsActiveOutput[] = { true, true, true, true, true, true,
 true, false, true, false };
 thrust::host_vector<bool> nodeIsActiveOutputFromGPU = nodes.nodeIsActive;
 for (uint i = 0; i < nodes.getCurrentActiveCellCount() * maxNodePerCell;
 i++) {
 EXPECT_EQ(expectedNodeIsActiveOutput[i], nodeIsActiveOutputFromGPU[i]);
 }
 }

 TEST(SceCellsCompCelLCenterTest,fixedTest) {
 cudaSetDevice(myDeviceId);
 const uint maxCellCount = 2;
 const uint initCellCount = 2;
 const uint maxNodePerCell = 4;
 const uint maxECMCount = 2;
 const uint maxNodeInECM = 1;
 const uint maxTotalNodeCount = maxCellCount * maxNodePerCell
 + maxECMCount * maxNodeInECM;
 SceNodes nodes(maxCellCount, maxNodePerCell, maxECMCount, maxNodeInECM);
 nodes.setCurrentActiveCellCount(initCellCount);
 nodes.setCurrentActiveEcm(maxNodeInECM);

 double nodeXInput[] = { 1.2, 3, 2, 1.5, 0.3, 1.1, 9.9, 4.2, 0.0, 0.0 };
 double nodeYInput[] = { 2.3, 1, 2, 5.6, 0.9, 8.6, 2.3, 5.9, 0.0, 0.0 };
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

 SceCells cells(&nodes);
 thrust::host_vector<uint> activeNodeCount(2);
 activeNodeCount[0] = 4;
 activeNodeCount[1] = 3;
 cells.activeNodeCountOfThisCell = activeNodeCount;
 cells.distributeIsActiveInfo();
 cells.computeCenterPos();
 thrust::host_vector<double> centerXFromGPU = cells.centerCoordX;
 thrust::host_vector<double> centerYFromGPU = cells.centerCoordY;
 double expectedCenterX[] = { 7.7 / 4, 11.3 / 3 };
 double expectedCenterY[] = { 10.9 / 4, 11.8 / 3 };
 EXPECT_NEAR(centerXFromGPU[0], expectedCenterX[0], errTol);
 EXPECT_NEAR(centerXFromGPU[1], expectedCenterX[1], errTol);
 EXPECT_NEAR(centerYFromGPU[0], expectedCenterY[0], errTol);
 EXPECT_NEAR(centerYFromGPU[1], expectedCenterY[1], errTol);
 }

 TEST(GrowthTest, fixedTest) {
 cudaSetDevice(myDeviceId);
 const uint maxCellCount = 2;
 const uint initCellCount = 2;
 const uint maxNodePerCell = 4;
 const uint maxECMCount = 2;
 const uint maxNodeInECM = 1;
 const uint maxTotalNodeCount = maxCellCount * maxNodePerCell
 + maxECMCount * maxNodeInECM;
 SceNodes nodes(maxCellCount, maxNodePerCell, maxECMCount, maxNodeInECM);
 nodes.setCurrentActiveCellCount(initCellCount);
 nodes.setCurrentActiveEcm(maxNodeInECM);

 double nodeXInput[] = { 1.2, 3, 2, 1.5, 7.3, 6.1, 9.9, 4.2, 0.0, 0.0 };
 double nodeYInput[] = { 4.3, 9, 8, 9.6, 0.9, 3.6, 2.3, 5.9, 0.0, 0.0 };
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

 uint growthVectorSize = 4;
 SceCells cells(&nodes);
 thrust::host_vector<uint> activeNodeCount(2);
 activeNodeCount[0] = 3;
 activeNodeCount[1] = 4;
 cells.activeNodeCountOfThisCell = activeNodeCount;
 cells.distributeIsActiveInfo();
 cells.computeCenterPos();
 double expectedCenterXCoord[2];
 double expectedCenterYCoord[2];
 expectedCenterXCoord[0] = (nodeXInput[0] + nodeXInput[1] + nodeXInput[2])
 / 3.0;
 expectedCenterXCoord[1] = (nodeXInput[4] + nodeXInput[5] + nodeXInput[6]
 + nodeXInput[7]) / 4.0;
 expectedCenterYCoord[0] = (nodeYInput[0] + nodeYInput[1] + nodeYInput[2])
 / 3.0;
 expectedCenterYCoord[1] = (nodeYInput[4] + nodeYInput[5] + nodeYInput[6]
 + nodeYInput[7]) / 4.0;
 thrust::host_vector<double> centerXFromGPU = cells.centerCoordX;
 thrust::host_vector<double> centerYFromGPU = cells.centerCoordY;
 EXPECT_NEAR(expectedCenterXCoord[0], centerXFromGPU[0], errTol);
 EXPECT_NEAR(expectedCenterXCoord[1], centerXFromGPU[1], errTol);
 EXPECT_NEAR(expectedCenterYCoord[0], centerYFromGPU[0], errTol);
 EXPECT_NEAR(expectedCenterYCoord[1], centerYFromGPU[1], errTol);
 thrust::device_vector<double> growthMag, growthDirX, growthDirY;
 double growMagArr[] = { 1.2, 3.4, 5.6, 7.8 };
 double growDirXArr[] = { 1.0, 0.0, 0.0, -1.0 };
 double growDirYArr[] = { 0.0, -1.0, 1.0, 0.0 };
 thrust::host_vector<double> growthMagHost(growMagArr,
 growMagArr + growthVectorSize);
 thrust::host_vector<double> growthDirXHost(growDirXArr,
 growDirXArr + growthVectorSize);
 thrust::host_vector<double> growthDirYHost(growDirYArr,
 growDirYArr + growthVectorSize);
 growthMag = growthMagHost;
 growthDirX = growthDirXHost;
 growthDirY = growthDirYHost;
 uint gridDimX = 2, gridDimY = 2;
 double gridSpacing = 5.0, dt = 0.1;
 cells.grow2DSimplified(dt, growthMag, growthDirX, growthDirY, gridDimX,
 gridDimY, gridSpacing);
 double expectedGrowthMag[2];
 double expectedGrowthXDir[2];
 double expectedGrowthYDir[2];
 double expectedGrowthProgress[2];
 bool   expectedIsScheduleToGrow[2];
 // cell 1 : X Pos = 0 Y Pos = 1
 // cell 2 : X Pos = 1 Y Pos = 0
 expectedGrowthMag[0] = growMagArr[2];
 expectedGrowthMag[1] = growMagArr[1];
 expectedGrowthXDir[0] = growDirXArr[2];
 expectedGrowthXDir[1] = growDirXArr[1];
 expectedGrowthYDir[0] = growDirYArr[2];
 expectedGrowthYDir[1] = growDirYArr[1];
 expectedGrowthProgress[0] = expectedGrowthMag[0] * dt;
 expectedGrowthProgress[1] = expectedGrowthMag[1] * dt;
 expectedIsScheduleToGrow[0] = true;
 expectedIsScheduleToGrow[1] = false;
 thrust::host_vector<double> growthSpeedFromGPU = cells.growthSpeed;
 thrust::host_vector<double> growthXDirFromGPU = cells.growthXDir;
 thrust::host_vector<double> growthYDirFromGPU = cells.growthYDir;
 thrust::host_vector<double> growthProgressFromGPU = cells.growthProgress;
 thrust::host_vector<bool> isGoingToAddFromGPU = cells.isScheduledToGrow;
 for (uint i = 0; i < 2; i++) {
 EXPECT_NEAR(expectedGrowthMag[i], growthSpeedFromGPU[i], errTol);
 EXPECT_NEAR(expectedGrowthXDir[i], growthXDirFromGPU[i], errTol);
 EXPECT_NEAR(expectedGrowthYDir[i], growthYDirFromGPU[i], errTol);
 EXPECT_NEAR(expectedGrowthProgress[i], growthProgressFromGPU[i],
 errTol);
 //EXPECT_EQ(expectedIsScheduleToGrow[i],isGoingToAddFromGPU[i]);
 }
 }
 */
