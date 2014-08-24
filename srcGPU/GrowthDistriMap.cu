#include "GrowthDistriMap.h"

GrowthDistriMap::GrowthDistriMap(uint xDim, uint yDim, double spacing) {
	gridDimensionX = xDim;
	gridDimensionY = yDim;
	gridSpacing = spacing;
	uint vectorSize = gridDimensionX * gridDimensionY;
	growthFactorMag.resize(vectorSize);
	growthFactorDirXComp.resize(vectorSize);
	growthFactorDirYComp.resize(vectorSize);
}

void GrowthDistriMap::initialize(double lowerLeftPointX, double lowerLeftPointY,
		double centerX, double centerY, double high, double low, double slope) {
	uint vectorSize = gridDimensionX * gridDimensionY;
	thrust::host_vector<double> tmpMag(vectorSize);
	thrust::host_vector<double> tmpXDir(vectorSize);
	thrust::host_vector<double> tmpYDir(vectorSize);
	uint i, j, posInVector;
	double tmpXPos, tmpYPos;
	for (i = 0; i < gridDimensionX; i++) {
		for (j = 0; j < gridDimensionY; j++) {
			tmpXPos = lowerLeftPointX + i * gridSpacing + 0.5 * gridSpacing;
			tmpYPos = lowerLeftPointY + j * gridSpacing + 0.5 * gridSpacing;
			posInVector = i + j * gridDimensionX;
			double distanceToCenter = sqrt(
					(centerX - tmpXPos) * (centerX - tmpXPos)
							+ (centerY - tmpYPos) * (centerY - tmpYPos));
			double dirX = (centerX - tmpXPos) / distanceToCenter;
			double dirY = (centerY - tmpYPos) / distanceToCenter;
			double mag = high - distanceToCenter * slope;
			if (mag < 0) {
				mag = 0;
			}
			tmpMag[posInVector] = mag;
			tmpXDir[posInVector] = dirX;
			tmpYDir[posInVector] = dirY;
		}
	}
	growthFactorMag = tmpMag;
	growthFactorDirXComp = tmpXDir;
	growthFactorDirYComp = tmpYDir;
}
