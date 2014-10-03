/*
 * ResAnalysisHelper.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: wenzhao
 */

#include "ResAnalysisHelper.h"

std::vector<std::vector<uint> > ResAnalysisHelper::outputLabelMatrix(
		std::vector<NodeWithLabel>& nodeLabels) {
	std::vector<std::vector<uint> > result;
	std::vector<std::vector<std::vector<LabelWithDist> > > matrixRawData;
	result.resize(_pixelPara.pixelYDim);
	matrixRawData.resize(_pixelPara.pixelYDim);
	for (uint i = 0; i < _pixelPara.pixelYDim; i++) {
		result[i].resize(_pixelPara.pixelXDim, -1);
		matrixRawData[i].resize(_pixelPara.pixelXDim);
	}
	for (uint i = 0; i < nodeLabels.size(); i++) {
		std::vector<Index2D> neighborPixels = obtainNeighborPixels(
				nodeLabels[i]);
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
			double dist = tmpDistVec.getModul();
			if (dist < _pixelPara.effectiveRange) {
				result.push_back(tmpIndex);
			}
		}
	}
	return result;
}

void ResAnalysisHelper::updateRawMatrix(
		std::vector<std::vector<std::vector<LabelWithDist> > >& rawMatrix,
		NodeWithLabel& nodeLabel) {
	std::vector<Index2D> indicies2D = obtainNeighborPixels(nodeLabel);
	for (uint i = 0; i < indicies2D.size(); i++) {
		LabelWithDist labelDist;
		labelDist.dist = computeDist(nodeLabel, indicies2D[i]);
		labelDist.label = nodeLabel.cellRank;
		rawMatrix[indicies2D[i].indexX][indicies2D[i].indexY].push_back(
				labelDist);
	}
}

void ResAnalysisHelper::updateLabelMatrix(
		std::vector<std::vector<uint> >& resultMatrix,
		std::vector<std::vector<std::vector<LabelWithDist> > >& rawMatrix) {
	for (uint i = 0; i < rawMatrix.size(); i++) {
		for (uint j = 0; j < rawMatrix[i].size(); j++) {
			if (rawMatrix[i][j].size() == 0) {
				resultMatrix[i][j] = -1;
			} else {
				double minDist = rawMatrix[i][j][0].dist;
				uint label = rawMatrix[i][j][0].label;
				for (uint k = 1; k < rawMatrix[i][j].size(); k++) {
					if (rawMatrix[i][j][k].dist < minDist) {
						minDist = rawMatrix[i][j][k].dist;
						label = rawMatrix[i][j][k].label;
					}
				}
				resultMatrix[i][j] = label;
			}
		}
	}
}

ResAnalysisHelper::ResAnalysisHelper(PixelizePara& pixelPara) {
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

ResAnalysisHelper::~ResAnalysisHelper() {
// TODO Auto-generated destructor stub
}

