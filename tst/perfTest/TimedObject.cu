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
	double result = -1.0;
	switch (cfg.type) {
	case SortGPU:
		return countGPUSorting(cfg.size, cfg.randomRange);
	case SortCPU:
		return countCPUSorting(cfg.size, cfg.randomRange);
	case SceNodeGPU:
		return countSceMove(cfg.size);
	}
	return result;
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

void TimedObject::loadInitSceData(std::vector<CVector> nodeInitPosData,
		std::vector<CVector> centerInitPosData) {
	nodeInitPos = nodeInitPosData;
	centerInitPos = centerInitPosData;
}

double TimedObject::countSceMove(uint size) {

	double minX = globalConfigVars.getConfigValue("DOMAIN_XMIN").toDouble();
	double maxX = globalConfigVars.getConfigValue("DOMAIN_XMAX").toDouble();
	double minY = globalConfigVars.getConfigValue("DOMAIN_YMIN").toDouble();
	double maxY = globalConfigVars.getConfigValue("DOMAIN_YMAX").toDouble();

	double gridSpacing =
			globalConfigVars.getConfigValue("DOMAIN_GRID_SPACING").toDouble();

	std::vector<double> emptyVector;
	std::vector<double> nodeXVector(size * 20);
	std::vector<double> nodeYVector(size * 20);

	utils.transformVals(nodeXVector, nodeYVector, nodeInitPos, centerInitPos);

	SceNodes nodes = SceNodes(0, 0, 0, 0, 0, size, 20, false);
	nodes.initDimension(minX, maxX, minY, maxY, gridSpacing);

	NodeAllocPara para = nodes.getAllocPara();
	para.currentActiveCellCount = size;
	nodes.setAllocPara(para);

	nodes.initValues(emptyVector, emptyVector, emptyVector, emptyVector,
			emptyVector, emptyVector, emptyVector, emptyVector, nodeXVector,
			nodeYVector);
	thrust::fill(nodes.getInfoVecs().nodeIsActive.begin(),
			nodes.getInfoVecs().nodeIsActive.end(), true);

	cudaEvent_t start, stop;
	float elapsedTime;

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	for (uint i = 0; i < 3000; i++) {
		nodes.calculateAndApplySceForces();
	}

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	// elapsedTime is in unit of ms

	double gpuSceMoveTime = elapsedTime / 1000.0;
	return gpuSceMoveTime;

}

TimedObject::~TimedObject() {
// TODO Auto-generated destructor stub
}

