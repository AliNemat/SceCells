/*
 * PerfTestUtils.cpp
 *
 *  Created on: Sep 12, 2014
 *      Author: wsun2
 */

#include "PerfTestUtils.h"

PerfTestUtils::PerfTestUtils() {
	// TODO Auto-generated constructor stub

}

std::vector<double> PerfTestUtils::obtainRandomVector(uint size, double range) {
	std::vector<double> result(size);
	std::generate(result.begin(), result.end(), std::rand);
	for (uint i = 0; i < size; i++) {
		result[i] = result[i] * range / RAND_MAX;
	}
	return result;
}

PerfTestUtils::~PerfTestUtils() {
	// TODO Auto-generated destructor stub
}

