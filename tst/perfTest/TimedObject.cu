/*
 * TimedObject.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: wsun2
 */

#include "TimedObject.h"

TimedObject::TimedObject() {
	isSortingDataInitialized = 0;
}

double TimedObject::countExecutionTime(RunConfig cfg) {
	double result;
	switch (cfg.type) {
	case SortGPU:
		return countGPUSorting(cfg.size, cfg.randomRange);
	case SortCPU:
		return countCPUSorting(cfg.size, cfg.randomRange);
	}
}

double TimedObject::countCPUSorting(uint size, double randomRange) {
	if ((isSortingDataInitialized && initialVectorForSorting.size() != size)
			|| !isSortingDataInitialized) {
		initialVectorForSorting = utils.obtainRandomVector(size, randomRange);
		isSortingDataInitialized = true;
	}

	std::clock_t start;
	double duration;
	start = std::clock();
	sort(initialVectorForSorting.begin(), initialVectorForSorting.end());
	duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;

	return duration;
}

double TimedObject::countGPUSorting(uint size, double randomRange) {
	if ((isSortingDataInitialized && initialVectorForSorting.size() != size)
			|| !isSortingDataInitialized) {
		initialVectorForSorting = utils.obtainRandomVector(size, randomRange);
		isSortingDataInitialized = true;
	}
	cudaEvent_t start, stop;
	float elapsedTime;
	thrust::host_vector<double> h_vec(size);
	for (uint i = 0; i < size; i++) {
		h_vec[i] = initialVectorForSorting[i];
	}

	thrust::device_vector<double> d_vec = h_vec;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	thrust::sort(d_vec.begin(), d_vec.end());
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	// elapsedTime is in unit of ms
	double gpuSortingTime = elapsedTime / 1000.0;

	return gpuSortingTime;

}

TimedObject::~TimedObject() {
// TODO Auto-generated destructor stub
}

