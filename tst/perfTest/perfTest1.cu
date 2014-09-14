#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/generate.h>
#include <cuda_runtime.h>
#include "gtest/gtest.h"
#include "SceNodes.h"
#include "commonData.h"
#include "TimedObject.h"

double errTol = 1.0e-9;

// size = 2^24
uint vectorSize = 1 << 10;
GlobalConfigVars globalConfigVars;

TEST(DummyTest, dummyTest) {
	cudaEvent_t start, stop;
	float elapsedTime;
	thrust::host_vector<double> h_vec(vectorSize);
	thrust::generate(h_vec.begin(), h_vec.end(), rand);
	std::vector<double> stdVector(h_vec.size());
	for (uint i = 0; i < vectorSize; i++) {
		stdVector[i] = h_vec[i];
	}
	thrust::device_vector<double> d_vec = h_vec;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	thrust::sort(d_vec.begin(), d_vec.end());
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	cudaEventRecord(start, 0);
	sort(stdVector.begin(), stdVector.end());
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);

	thrust::host_vector<double> d_res = d_vec;
	for (uint i = 0; i < vectorSize; i++) {
		EXPECT_NEAR(d_res[i], stdVector[i], errTol);
	}

	cout << "Purpose of this dummy test is to make sure that " << endl;
	cout << "(1) sorting result is correct and " << endl;
	cout << "(2) warm up GPU to reduce overhead" << endl;
}

class PerfTest: public ::testing::Test {
protected:

	std::vector<CVector> initNodesPos;
	std::vector<CVector> initCentersPos;

	virtual void SetUp() {
		ConfigParser parser;
		std::string configFileName = "./resources/perfTest.cfg";
		PerfTestUtils utils;
		//CellInitHelper helper;
		globalConfigVars = parser.parseConfigFile(configFileName);
		cudaSetDevice(
				globalConfigVars.getConfigValue("GPUDeviceNumber").toInt());
		double cellRadius =
				globalConfigVars.getConfigValue("TestCellRadius").toDouble();
		double cellNodeInterval = globalConfigVars.getConfigValue(
				"TestCellNodeInterval").toDouble();
		initNodesPos = utils.obtainCellInitPointForTesting(cellRadius,
				cellNodeInterval);
		initNodesPos.push_back(CVector(0.2, 0.0, 0.0));
		initCentersPos = utils.obtainCellInitCentersForTesting(80, 0.85,
				CVector(25.0, 25.0, 0));
	}
};

TEST(DummyTest, sortTest) {
	uint sortVectorSize = 1 << 23;
	double sortRandomRange = 2.0;
	RunConfig cpuConfig, gpuConfig;
	cpuConfig.size = sortVectorSize;
	cpuConfig.randomRange = sortRandomRange;
	cpuConfig.type = SortCPU;
	gpuConfig.size = sortVectorSize;
	gpuConfig.randomRange = sortRandomRange;
	gpuConfig.type = SortGPU;
	TimedObject sortObj;
	double sortTimeCPU = sortObj.countExecutionTime(cpuConfig);
	double sortTimeGPU = sortObj.countExecutionTime(gpuConfig);
	cout << "gpu sorting time = " << sortTimeGPU << endl;
	cout << "cpu sorting time = " << sortTimeCPU << endl;
	cout << "GPU speed up for sorting: " << (sortTimeCPU / sortTimeGPU)
			<< " times" << endl;
}

TEST_F(PerfTest, sceNodePerfTest128Cells) {
	uint size1 = 128;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 128 cells = " << sceExecTime1 << endl;
}

TEST_F(PerfTest, sceNodePerfTest250Cells) {
	uint size1 = 250;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 250 cells = " << sceExecTime1 << endl;
}

TEST_F(PerfTest, sceNodePerfTest500Cells) {
	uint size1 = 500;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 500 cells = " << sceExecTime1 << endl;
}

TEST_F(PerfTest, sceNodePerfTest1000Cells) {
	uint size1 = 1000;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 1000 cells = " << sceExecTime1 << endl;
}

TEST_F(PerfTest, sceNodePerfTest5000Cells) {
	uint size1 = 5000;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 5000 cells = " << sceExecTime1 << endl;
}

TEST_F(PerfTest, sceNodePerfTest10000Cells) {
	uint size1 = 10000;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 10000 cells = " << sceExecTime1 << endl;
}

TEST_F(PerfTest, sceNodePerfTest20000Cells) {
	uint size1 = 20000;
	RunConfig gpuConfig1;
	gpuConfig1.size = size1;
	gpuConfig1.type = SceNodeGPU;
	TimedObject sceNodeTimerObj;
	sceNodeTimerObj.loadInitSceData(initNodesPos, initCentersPos);
	double sceExecTime1 = sceNodeTimerObj.countExecutionTime(gpuConfig1);
	cout << "gpu sceNode time for 20000 cells = " << sceExecTime1 << endl;
}

