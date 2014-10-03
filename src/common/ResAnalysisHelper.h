/*
 * ResAnalysisHelper.h
 *
 *  Created on: Oct 2, 2014
 *      Author: wenzhao
 */

#include "commonData.h"

#ifndef RESANALYSISHELPER_H_
#define RESANALYSISHELPER_H_

struct Index2D {
	uint indexX;
	uint indexY;
};

struct PixelizePara {
	uint pixelXDim;
	uint pixelYDim;
	double xMin, xMax;
	double yMin, yMax;
	double effectiveRange;
};

struct NodeWithLabel {
	CVector position;
	uint cellRank;
};

struct LabelWithDist {
	uint label;
	double dist;
};

class ResAnalysisHelper {
public:
	ResAnalysisHelper();
	static CVector obtainCenterLoc(Index2D index2D);
	static std::vector<Index2D> obtainNeighborPixels(CVector &pos,
			PixelizePara& paras);
	static void updateRawMatrix(
			std::vector<std::vector<std::vector<LabelWithDist> > > &rawMatrix,
			std::vector<Index2D> &indicies2D, CVector &pos);
	static void updateLabelMatrix(std::vector<std::vector<uint> > &resultMatrix,
			std::vector<std::vector<std::vector<LabelWithDist> > > &rawMatrix);
	static std::vector<std::vector<uint> > outputLabelMatrix(
			std::vector<NodeWithLabel> &nodeLabels, PixelizePara &paras);
	virtual ~ResAnalysisHelper();
};

#endif /* RESANALYSISHELPER_H_ */
