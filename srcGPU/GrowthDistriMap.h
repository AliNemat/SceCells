#ifndef GrowthDistriMap_H_
#define GrowthDistriMap_H_

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>

/**
 * Provides information for cell growth influenced by chemical concentration gradient.
 */
class GrowthDistriMap {
public:
	thrust::device_vector<double> growthFactorMag;
	thrust::device_vector<double> growthFactorDirXComp;
	thrust::device_vector<double> growthFactorDirYComp;
	uint gridDimensionX;
	uint gridDimensionY;
	double gridSpacing;
	double lowerLeftPtX;
	double lowerLeftPtY;

	GrowthDistriMap() {
	}
	GrowthDistriMap(uint xDim, uint yDim, double spacing);
	void initialize(double lowerLeftPointX, double lowerLeftPointY,
			double centerX, double centerY, double high, double low,
			double slope);
};

#endif
