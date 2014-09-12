/*
 * PerfTestUtils.h
 *
 *  Created on: Sep 12, 2014
 *      Author: wsun2
 */

#include <vector>
#include <algorithm>
#include <cstdlib>

#ifndef PERFTESTUTILS_H_
#define PERFTESTUTILS_H_

class PerfTestUtils {
public:
	PerfTestUtils();
	std::vector<double> obtainRandomVector(uint size,double range);
	virtual ~PerfTestUtils();
};

#endif /* PERFTESTUTILS_H_ */
