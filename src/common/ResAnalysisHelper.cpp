/*
 * ResAnalysisHelper.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: wenzhao
 */

#include "ResAnalysisHelper.h"

std::vector<std::vector<int> > ResAnalysisHelper::outputLabelMatrix(
		std::vector<NodeWithLabel>& nodeLabels) {
	std::vector<std::vector<int> > result;
	std::vector<std::vector<std::vector<LabelWithDist> > > matrixRawData;
	result.resize(_pixelPara.pixelYDim);
	matrixRawData.resize(_pixelPara.pixelYDim);
	for (uint i = 0; i < _pixelPara.pixelYDim; i++) {
		result[i].resize(_pixelPara.pixelXDim, -1);
		matrixRawData[i].resize(_pixelPara.pixelXDim);
	}

	//cout << "bb in outputLabelMatrix" << endl;
	//cout.flush();

	for (uint i = 0; i < nodeLabels.size(); i++) {
		//std::vector<Index2D> neighborPixels = obtainNeighborPixels(
		//		nodeLabels[i]);
		//cout << "bb " << i << "in outputLabelMatrix" << endl;
		//cout.flush();
		updateRawMatrix(matrixRawData, nodeLabels[i]);
	}
	updateLabelMatrix(result, matrixRawData);
	return result;
}

CVector ResAnalysisHelper::obtainCenterLoc(Index2D &index2D) {
	double xMultiplier = index2D.indexX + 0.5;
	double yMultiplier = index2D.indexY + 0.5;
	CVector result(xMultiplier * _pixelSpacing + _pixelPara.xMin,
			yMultiplier * _pixelSpacing + _pixelPara.yMin, 0.0);
	return result;
}

std::vector<Index2D> ResAnalysisHelper::obtainNeighborPixels(
		NodeWithLabel& nodeLabel) {
	std::vector<Index2D> result;
	Index2D index2D = obtainIndex2D(nodeLabel.position);
	Index2D tmpIndex;
	uint minX = std::max((uint) 0, index2D.indexX - _integerRadius);
	uint maxX = std::min(_pixelPara.pixelXDim - 1,
			index2D.indexX + _integerRadius);
	uint minY = std::max((uint) 0, index2D.indexY - _integerRadius);
	uint maxY = std::min(_pixelPara.pixelYDim - 1,
			index2D.indexY + _integerRadius);
	// note the for loop is inclusive
	for (uint i = minX; i <= maxX; i++) {
		// note the for loop is inclusive
		for (uint j = minY; j <= maxY; j++) {
			tmpIndex.indexX = i;
			tmpIndex.indexY = j;
			CVector tmpPos = obtainCenterLoc(tmpIndex);
			CVector tmpDistVec = nodeLabel.position - tmpPos;
			//cout << "node pos: ";
			//nodeLabel.position.Print();
			//cout << "tmp node pos: ";
			//tmpPos.Print();
			double dist = tmpDistVec.getModul();
			//cout << "distance =  " << dist << endl;
			if (dist
					< _pixelPara.effectiveRange
							+ _pixelPara.allowedAbsoluteError) {
				result.push_back(tmpIndex);
				//cout << "pixel X: " << tmpIndex.indexX << "pixel Y:"
				//		<< tmpIndex.indexY << endl;
			}
		}
	}
	return result;
}

void ResAnalysisHelper::updateRawMatrix(
		std::vector<std::vector<std::vector<LabelWithDist> > >& rawMatrix,
		NodeWithLabel& nodeLabel) {
	//cout << "bb 0 in updateRawMatrix" << endl;
	//cout.flush();
	std::vector<Index2D> indicies2D = obtainNeighborPixels(nodeLabel);
	//cout << "bb 1 in updateRawMatrix" << endl;
	//cout.flush();
	for (uint i = 0; i < indicies2D.size(); i++) {
		//cout << "bbb" << i << " in updateRawMatrix" << endl;
		//cout.flush();
		LabelWithDist labelDist;
		labelDist.dist = computeDist(nodeLabel, indicies2D[i]);
		labelDist.label = nodeLabel.cellRank;
		rawMatrix[indicies2D[i].indexY][indicies2D[i].indexX].push_back(
				labelDist);
	}
}

void ResAnalysisHelper::updateLabelMatrix(
		std::vector<std::vector<int> >& resultMatrix,
		std::vector<std::vector<std::vector<LabelWithDist> > >& rawMatrix) {
	for (uint i = 0; i < rawMatrix.size(); i++) {
		for (uint j = 0; j < rawMatrix[i].size(); j++) {
			if (rawMatrix[i][j].size() == 0) {
				resultMatrix[i][j] = -1;
			} else {
				double minDist = rawMatrix[i][j][0].dist;
				uint label = rawMatrix[i][j][0].label;

				uint originalLabel = label;
				bool allSameLabel = true;

				for (uint k = 1; k < rawMatrix[i][j].size(); k++) {

					if (rawMatrix[i][j][k].label != originalLabel) {
						allSameLabel = false;
					}

					if (rawMatrix[i][j][k].dist < minDist) {
						minDist = rawMatrix[i][j][k].dist;
						label = rawMatrix[i][j][k].label;

					}
				}
				if (!allSameLabel) {
					resultMatrix[i][j] = label;
				} else {
					if (minDist < _pixelPara.effectiveRange_single) {
						resultMatrix[i][j] = label;
					} else {
						resultMatrix[i][j] = -1;
					}
				}
			}
		}
	}
}

Index2D ResAnalysisHelper::obtainIndex2D(CVector &pos) {
	Index2D result;
	double xFromMin = pos.GetX() - _pixelPara.xMin;
	double yFromMin = pos.GetY() - _pixelPara.yMin;
	uint positionX = xFromMin / _pixelSpacing;
	uint positionY = yFromMin / _pixelSpacing;
	result.indexX = positionX;
	result.indexY = positionY;
	return result;
}

double ResAnalysisHelper::computeDist(NodeWithLabel& nodeLabel,
		Index2D& index2D) {
	CVector posForIndex = obtainCenterLoc(index2D);
	CVector tmpDist = nodeLabel.position - posForIndex;
	return tmpDist.getModul();
}

ResAnalysisHelper::ResAnalysisHelper() {
	_pixelPara.initFromConfigFile();
	double pixelSpacingX = (_pixelPara.xMax - _pixelPara.xMin)
			/ _pixelPara.pixelXDim;
	double pixelSpacingY = (_pixelPara.yMax - _pixelPara.yMin)
			/ _pixelPara.pixelYDim;
	if (fabs(pixelSpacingX - pixelSpacingY) > _pixelPara.allowedAbsoluteError) {
		throw SceException("domain size and pixel size does not fit",
				OutputAnalysisDataException);
	}
	_pixelSpacing = pixelSpacingX;
	_integerRadius = _pixelPara.effectiveRange / _pixelSpacing + 1;
}

void ResAnalysisHelper::setPixelPara(PixelizePara& pixelPara) {
	_pixelPara = pixelPara;
	double pixelSpacingX = (_pixelPara.xMax - _pixelPara.xMin)
			/ _pixelPara.pixelXDim;
	double pixelSpacingY = (_pixelPara.yMax - _pixelPara.yMin)
			/ _pixelPara.pixelYDim;
	if (fabs(pixelSpacingX - pixelSpacingY) > _pixelPara.allowedAbsoluteError) {
		throw SceException("domain size and pixel size does not fit",
				OutputAnalysisDataException);
	}
	_pixelSpacing = pixelSpacingX;
	_integerRadius = _pixelPara.effectiveRange / _pixelSpacing + 1;
}

ResAnalysisHelper::~ResAnalysisHelper() {
// TODO Auto-generated destructor stub
}

void PixelizePara::initFromConfigFile() {
	pixelXDim = globalConfigVars.getConfigValue("Pixel_Para_X_DIM").toInt();
	pixelYDim = globalConfigVars.getConfigValue("Pixel_Para_Y_DIM").toInt();

	xMin = globalConfigVars.getConfigValue("Pixel_Para_X_MIN").toDouble();
	xMax = globalConfigVars.getConfigValue("Pixel_Para_X_MAX").toDouble();

	yMin = globalConfigVars.getConfigValue("Pixel_Para_Y_MIN").toDouble();
	yMax = globalConfigVars.getConfigValue("Pixel_Para_Y_MAX").toDouble();

	effectiveRange = globalConfigVars.getConfigValue(
			"Pixel_Para_Effective_Range").toDouble();

	effectiveRange_single = globalConfigVars.getConfigValue(
			"Pixel_Para_Effective_Range_Single").toDouble();

	allowedAbsoluteError = globalConfigVars.getConfigValue(
			"Pixel_Para_Allowed_Error").toDouble();
}
