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

/**
 * Used for analyze result from the simulation.
 */
class ResAnalysisHelper {
	PixelizePara _pixelPara;
	double _pixelSpacing;
	/**
	 * This integer indicates how many pixels should be extended.
	 */
	uint _integerRadius;
	/**
	 * compute distance from node to a pixel.
	 */
	double computeDist(NodeWithLabel& nodeLabel, Index2D &index2D);
	/**
	 * obtain the center location of a pixel.
	 */
	CVector obtainCenterLoc(Index2D &index2D);

	/**
	 * given current exact position,
	 */
	Index2D obtainIndex2D(CVector &pos);

	/**
	 * given a node label, return all of the neighbor pixels.
	 */
	std::vector<Index2D> obtainNeighborPixels(NodeWithLabel &nodeLabel);
	/**
	 * updates raw matrix by inserting a data point which contains information of node label and
	 * distance of the label from pixel.
	 */
	void updateRawMatrix(
			std::vector<std::vector<std::vector<LabelWithDist> > > &rawMatrix,
			NodeWithLabel &nodeLabel);
	/**
	 * Generate a result label matrix given raw matrix.
	 * Raw matrix is a matrix which all entries are list of possible labels with their shortest distance.
	 */
	void updateLabelMatrix(std::vector<std::vector<int> > &resultMatrix,
			std::vector<std::vector<std::vector<LabelWithDist> > > &rawMatrix);
	/**
	 * Convert the matrix which is labeled by cell rank to three matrices,
	 * Red, Green, Blue, respectively.
	 */
	void generateRGBMatrix(std::vector<std::vector<int> > &labelMatrix,
			std::vector<std::vector<double> > &red,
			std::vector<std::vector<double> > &green,
			std::vector<std::vector<double> > &blue);
	/**
	 * Given a label value and maximum possible label value, generate three values,
	 * representing red, green, blue weight respectively.
	 */
	void transformToRGB(int &labelValue, int &maxLabelValue, double &rValue,
			double &gValue, double &bValue);
public:
	ResAnalysisHelper();
	void setPixelPara(PixelizePara &pixelPara);
	/**
	 * outputs label matrix given vector of node labels.
	 */
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
	 * @param fileName name of the polygon counting statistical data file.
	 * @param step current step.
	 * @param labelMatrix the label matrix for differentiate cells.
	 */
	void outputStat_PolygonCounting(std::string fileName, uint step,
			std::vector<std::vector<int> > &labelMatrix,
			std::vector<double> growthProVec = std::vector<double>());

	virtual ~ResAnalysisHelper();
};

#endif /* RESANALYSISHELPER_H_ */
