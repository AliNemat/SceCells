/*
 * TimedObject.h
 *
 *  Created on: Sep 12, 2014
 *      Author: wsun2
 */
#include <vector>
#include "PerfTestUtils.h"
#include <ctime>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/fill.h>
#include <thrust/generate.h>
#include <cuda_runtime.h>
#include <algorithm>
#include "SceNodes.h"

#ifndef TIMEDOBJECT_H_
#define TIMEDOBJECT_H_

enum PerfTestType {
	SortGPU, SortCPU, SceNodeGPU
};

/**
 * config for timedObject
 */
struct RunConfig {
	PerfTestType type;
	uint size;
	double randomRange;
};

class TimedObject {

	double countCPUSorting(uint size, double randomRange);
	double countGPUSorting(uint size, double randomRange);
	double countSceMove(uint size);
	bool isSortingDataInitialized;
	std::vector<double> initialVectorForSorting;
	std::vector<CVector> nodeInitPos;
	std::vector<CVector> centerInitPos;
	PerfTestUtils utils;

public:
	TimedObject();
	void loadInitSceData(std::vector<CVector> nodeInitPosData,
			std::vector<CVector> centerInitPosData);
	double countExecutionTime(RunConfig cfg);
	virtual ~TimedObject();
};

#endif /* TIMEDOBJECT_H_ */
