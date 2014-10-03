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
		updateRawMatrix(matrixRawData, neighborPixels, nodeLabels[i].position);
	}
	updateLabelMatrix(result, matrixRawData);
	return result;
}

CVector ResAnalysisHelper::obtainCenterLoc(Index2D index2D) {
	double xMultiplier = index2D.indexX + 0.5;
	double yMultiplier = index2D.indexY + 0.5;
	CVector result(xMultiplier * _pixelSpacing, yMultiplier * _pixelSpacing,
			0.0);
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
	for (uint i = minX; i <= maxX; i++) {
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
		std::vector<Index2D>& indicies2D, CVector& pos) {
}

void ResAnalysisHelper::updateLabelMatrix(
		std::vector<std::vector<uint> >& resultMatrix,
		std::vector<std::vector<std::vector<LabelWithDist> > >& rawMatrix) {
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

Index2D ResAnalysisHelper::obtainIndex2D(CVector pos) {
	Index2D result;
	double xFromMin = pos.GetX() - _pixelPara.xMin;
	return result;

}

ResAnalysisHelper::~ResAnalysisHelper() {
// TODO Auto-generated destructor stub
}

