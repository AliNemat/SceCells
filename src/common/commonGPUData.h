/*
 * commonGPUData.h
 *
 *  Created on: Sep 18, 2014
 *      Author: wsun2
 */

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cuda_runtime.h>
#include <vector>
#include "commonData.h"

#ifndef COMMONGPUDATA_H_
#define COMMONGPUDATA_H_

/**
 * This function converts a Thrust vector (host) to c++ std vector
 */
template<class T>
inline void thrustVecToStd(thrust::host_vector<T> &thrustVec,
		std::vector<T> &stdVec) {
	uint vecSize = thrustVec.size();
	stdVec.resize(vecSize);
	for (uint i = 0; i < vecSize; i++) {
		stdVec[i] = thrustVec[i];
	}
}

/**
 * This function converts a Thrust vector (device) to c++ std vector
 */
template<class T>
inline void thrustVecToStd(thrust::device_vector<T> &thrustVec,
		std::vector<T> &stdVec) {
	uint vecSize = thrustVec.size();
	thrust::host_vector<T> vecCPU = thrustVec;
	stdVec.resize(vecSize);
	for (uint i = 0; i < vecSize; i++) {
		stdVec[i] = vecCPU[i];
	}
}

#endif /* COMMONGPUDATA_H_ */
