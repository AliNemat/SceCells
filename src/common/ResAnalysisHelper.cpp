/*
 * ResAnalysisHelper.cpp
 *
 *  Created on: Oct 2, 2014
 *      Author: wenzhao
 */

#include "ResAnalysisHelper.h"

ResAnalysisHelper::ResAnalysisHelper() {
	// TODO Auto-generated constructor stub

}

std::vector<std::vector<uint> > ResAnalysisHelper::outputLabelMatrix(
		std::vector<NodeWithLabel>& nodeLabels, PixelizePara& paras) {
	std::vector<std::vector<uint> > result;
	std::vector<std::vector<std::vector<LabelWithDist> > > matrixRawData;
	result.resize(paras.pixelYDim);
	matrixRawData.resize(paras.pixelYDim);
	for (uint i = 0; i < paras.pixelYDim; i++) {
		result[i].resize(paras.pixelXDim, -1);
		matrixRawData[i].resize(paras.pixelXDim);
	}
	for (uint i = 0; i < nodeLabels.size(); i++) {
		std::vector<Index2D> neighborPixels = obtainNeighborPixels(
				nodeLabels[i].position, paras);
		updateRawMatrix(matrixRawData, neighborPixels, nodeLabels[i].position);
	}
	updateLabelMatrix(result, matrixRawData);
	return result;
}

CVector ResAnalysisHelper::obtainCenterLoc(Index2D index2D) {
	CVector result;
	return result;
}

std::vector<Index2D> ResAnalysisHelper::obtainNeighborPixels(CVector& pos,
		PixelizePara& paras) {
	std::vector<Index2D> result;
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

ResAnalysisHelper::~ResAnalysisHelper() {
// TODO Auto-generated destructor stub
}

