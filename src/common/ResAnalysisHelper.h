/*
 * ResAnalysisHelper.h
 *
 *  Created on: Oct 2, 2014
 *      Author: wenzhao
 */

#include "commonData.h"
#include "ConfigParser.h"

extern GlobalConfigVars globalConfigVars;

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
	double effectiveRange_single;
	double allowedAbsoluteError;
	void initFromConfigFile();
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
	PixelizePara _pixelPara;
	double _pixelSpacing;
	/**
	 * This variable is used for
	 */
	uint _integerRadius;
	double computeDist(NodeWithLabel& nodeLabel, Index2D &index2D);
	CVector obtainCenterLoc(Index2D &index2D);
	Index2D obtainIndex2D(CVector &pos);
	std::vector<Index2D> obtainNeighborPixels(NodeWithLabel &nodeLabel);
	void updateRawMatrix(
			std::vector<std::vector<std::vector<LabelWithDist> > > &rawMatrix,
			NodeWithLabel &nodeLabel);
	void updateLabelMatrix(std::vector<std::vector<int> > &resultMatrix,
			std::vector<std::vector<std::vector<LabelWithDist> > > &rawMatrix);
public:
	ResAnalysisHelper();
	void setPixelPara(PixelizePara &pixelPara);
	std::vector<std::vector<int> > outputLabelMatrix(
			std::vector<NodeWithLabel> &nodeLabels);
	virtual ~ResAnalysisHelper();
};

#endif /* RESANALYSISHELPER_H_ */
