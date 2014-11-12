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

std::vector<CVector> PerfTestUtils::obtainCellInitPointForTesting(
		double diskRadius, double circleRadius) {
	const double PI = acos(-1.0);
	int sideCount = (int) (diskRadius / circleRadius);
	std::vector<CVector> result;
	CVector startPoint(-circleRadius * sideCount, 0, 0);
	int maxNum = 2 * sideCount + 1;
	CVector topIncrease = CVector(circleRadius * sin(PI / 6.0),
			circleRadius * cos(PI / 6.0), 0.0);
	CVector bottomIncrease = CVector(circleRadius * sin(PI / 6.0),
			-circleRadius * cos(PI / 6.0), 0.0);
	CVector rightIncrease = CVector(circleRadius, 0, 0);
	CVector topIter = startPoint;
	int topCount = 0;
	while (topCount <= sideCount) {
		int num = maxNum - topCount;
		for (int i = 0; i < num; i++) {
			CVector rightIter = topIter + i * rightIncrease;
			result.push_back(rightIter);
		}
		topCount++;
		topIter = topIter + topIncrease;
	}
	int botCount = 1;
	CVector botIter = startPoint + bottomIncrease;
	while (botCount <= sideCount) {
		int num = maxNum - botCount;
		for (int i = 0; i < num; i++) {
			CVector rightIter = botIter + i * rightIncrease;
			result.push_back(rightIter);
		}
		botCount++;
		botIter = botIter + bottomIncrease;
	}
	return result;
}

std::vector<CVector> PerfTestUtils::obtainCellInitCentersForTesting(
		uint cellPerSizeCount, double dist, CVector initCenter) {
	const double PI = acos(-1.0);
	int sideCount = cellPerSizeCount;
	double circleRadius = dist;
	std::vector<CVector> result;
	CVector startPoint = CVector(-circleRadius * sideCount, 0, 0) + initCenter;
	int maxNum = 2 * sideCount + 1;
	CVector topIncrease = CVector(circleRadius * sin(PI / 6.0),
			circleRadius * cos(PI / 6.0), 0.0);
	CVector bottomIncrease = CVector(circleRadius * sin(PI / 6.0),
			-circleRadius * cos(PI / 6.0), 0.0);
	CVector rightIncrease = CVector(circleRadius, 0, 0);
	CVector topIter = startPoint;
	int topCount = 0;
	while (topCount <= sideCount) {
		int num = maxNum - topCount;
		for (int i = 0; i < num; i++) {
			CVector rightIter = topIter + i * rightIncrease;
			result.push_back(rightIter);
		}
		topCount++;
		topIter = topIter + topIncrease;
	}
	int botCount = 1;
	CVector botIter = startPoint + bottomIncrease;
	while (botCount <= sideCount) {
		int num = maxNum - botCount;
		for (int i = 0; i < num; i++) {
			CVector rightIter = botIter + i * rightIncrease;
			result.push_back(rightIter);
		}
		botCount++;
		botIter = botIter + bottomIncrease;
	}
	cout << "when generating, result size = " << result.size() << endl;
	return result;
}

void PerfTestUtils::transformVals(std::vector<double> &nodeXVector,
		std::vector<double> &nodeYVector, std::vector<CVector>& nodeInitPos,
		std::vector<CVector>& centerInitPos) {
	//assert(nodeXVector.size() == nodeYVector.size());
	cout << "node total size = " << nodeXVector.size() << endl;
	cout << "single cell init pos size =" << nodeInitPos.size() << endl;
	cout << "all cell centers pos size= " << centerInitPos.size() << endl;
	uint maxIndex = nodeXVector.size();
	uint cellCount = centerInitPos.size();
	uint nodeInCellCount = nodeInitPos.size();
	uint index = 0;
	for (uint i = 0; i < cellCount; i++) {
		for (uint j = 0; j < nodeInCellCount; j++) {
			index = i * nodeInCellCount + j;
			nodeXVector[index] = centerInitPos[i].x + nodeInitPos[j].x;
			nodeYVector[index] = centerInitPos[i].y + nodeInitPos[j].y;
			if (index >= maxIndex) {
				break;
			}
		}
		if (index >= maxIndex) {
			break;
		}
	}
}

PerfTestUtils::~PerfTestUtils() {
// TODO Auto-generated destructor stub
}

