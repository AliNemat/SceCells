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
	double gpuSortingTime = elapsedTime;
	cudaEventRecord(start, 0);
	sort(stdVector.begin(), stdVector.end());
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	double cpuSortingTime = elapsedTime;

	thrust::host_vector<double> d_res = d_vec;
	for (uint i = 0; i < vectorSize; i++) {
		EXPECT_NEAR(d_res[i], stdVector[i], errTol);
	}

	cout << "Purpose of this dummy test is to make sure that " << endl;
	cout << "(1) sorting result is correct and " << endl;
	cout << "(2) warm up GPU to reduce overhead" << endl;
}

TEST(DummyTest, sortTest) {
	uint sortVectorSize = 1 << 25;
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
	cout << " GPU speed up: " << (sortTimeCPU / sortTimeGPU) << " times"
			<< endl;
}
