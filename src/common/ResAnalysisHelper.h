/*
 * ResAnalysisHelper.h
 *
 *  Created on: Oct 2, 2014
 *      Author: wenzhao
 */

#include "commonData.h"
#include "ConfigParser.h"
#include "string.h"
#include <fstream>
#include <tr1/unordered_set>
#include <tr1/unordered_map>

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
	void generateRGBMatrix(std::vector<std::vector<int> > &labelMatrix,
			std::vector<std::vector<double> > &red,
			std::vector<std::vector<double> > &green,
			std::vector<std::vector<double> > &blue);
	void transformToRGB(int &labelValue, int &maxLabelValue, double &rValue,
			double &gValue, double &bValue);
public:
	ResAnalysisHelper();
	void setPixelPara(PixelizePara &pixelPara);
	std::vector<std::vector<int> > outputLabelMatrix(
			std::vector<NodeWithLabel> &nodeLabels);
	/**
	 * convert a label matrix to a bmp image.
	 * labels with value -1 are treated as empty.
	 */
	void outputImg_formatBMP(std::string fileName,
			std::vector<std::vector<int> > &labelMatrix);

	/**
	 * output stat file for polygon counting.
	 *
	 * @param growthProVec vector with information of cell's growth progress. optional.
	 *    supplying growth progress vector indicates statistics will include mitiotic shift effect.
	 *    by supplying such vector, the method will print out statistics for dividing and non-dividing cells.
	 */
	void outputStat_PolygonCounting(std::string fileName, uint step,
			std::vector<std::vector<int> > &labelMatrix,
			std::vector<double> growthProVec = std::vector<double>());

	virtual ~ResAnalysisHelper();
};

#endif /* RESANALYSISHELPER_H_ */
