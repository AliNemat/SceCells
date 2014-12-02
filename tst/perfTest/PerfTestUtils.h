/*
 * PerfTestUtils.h
 *
 *  Created on: Sep 12, 2014
 *      Author: wsun2
 */

#include <vector>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include "GeoVector.h"
#include <assert.h>

#ifndef PERFTESTUTILS_H_
#define PERFTESTUTILS_H_

class PerfTestUtils {
public:
	PerfTestUtils();
	std::vector<double> obtainRandomVector(uint size, double range);
	std::vector<CVector> obtainCellInitPointForTesting(double diskRadius,
			double circleRadius);
	std::vector<CVector> obtainCellInitCentersForTesting(uint cellCount,
			double dist, CVector initCenter);
	void transformVals(std::vector<double> &nodeXVector,
			std::vector<double> &nodeYVector, std::vector<CVector> &nodeInitPos,
			std::vector<CVector> &centerInitPos);
	void transformVals(std::vector<CVector> &nodePosVector,
			std::vector<CVector> &nodeInitPos,
			std::vector<CVector> &centerInitPos);
	virtual ~PerfTestUtils();
};

#endif /* PERFTESTUTILS_H_ */
