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

	for (uint i = 0; i < nodeLabels.size(); i++) {

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
			if (dist
					< _pixelPara.effectiveRange
							+ _pixelPara.allowedAbsoluteError) {
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
				/*
				 if (!allSameLabel) {
				 resultMatrix[i][j] = label;
				 } else {
				 if (minDist < _pixelPara.effectiveRange_single) {
				 resultMatrix[i][j] = label;
				 } else {
				 resultMatrix[i][j] = -1;
				 }
				 }
				 */
				// changed logic
				resultMatrix[i][j] = label;
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

void ResAnalysisHelper::outputImg_formatBMP(std::string fileName,
		std::vector<std::vector<int> >& labelMatrix) {
	if (labelMatrix.size() == 0) {
		return;
	}

	int imgWidth = labelMatrix[0].size();
	int imgHeight = labelMatrix.size();
	FILE *f;
	int imgIndexX, imgIndexY;
	int r, g, b;
	unsigned char *img = NULL;
	int filesize = 54 + 3 * imgWidth * imgHeight;
	img = (unsigned char *) malloc(3 * imgWidth * imgHeight);
	memset(img, 0, sizeof(3 * imgWidth * imgHeight));

	vector<vector<double> > red, green, blue;
	generateRGBMatrix(labelMatrix, red, green, blue);

	for (int i = 0; i < imgWidth; i++) {
		for (int j = 0; j < imgHeight; j++) {
			imgIndexX = i;
			imgIndexY = (imgHeight - 1) - j;
			r = red[i][j] * 255;
			g = green[i][j] * 255;
			b = blue[i][j] * 255;
			if (r > 255)
				r = 255;
			if (g > 255)
				g = 255;
			if (b > 255)
				b = 255;
			img[(imgIndexX + imgIndexY * imgWidth) * 3 + 2] =
					(unsigned char) (r);
			img[(imgIndexX + imgIndexY * imgWidth) * 3 + 1] =
					(unsigned char) (g);
			img[(imgIndexX + imgIndexY * imgWidth) * 3 + 0] =
					(unsigned char) (b);
		}
	}

	unsigned char bmpfileheader[14] = { 'B', 'M', 0, 0, 0, 0, 0, 0, 0, 0, 54, 0,
			0, 0 };
	unsigned char bmpinfoheader[40] = { 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
			0, 24, 0 };
	unsigned char bmppad[3] = { 0, 0, 0 };

	bmpfileheader[2] = (unsigned char) (filesize);
	bmpfileheader[3] = (unsigned char) (filesize >> 8);
	bmpfileheader[4] = (unsigned char) (filesize >> 16);
	bmpfileheader[5] = (unsigned char) (filesize >> 24);

	bmpinfoheader[4] = (unsigned char) (imgWidth);
	bmpinfoheader[5] = (unsigned char) (imgWidth >> 8);
	bmpinfoheader[6] = (unsigned char) (imgWidth >> 16);
	bmpinfoheader[7] = (unsigned char) (imgWidth >> 24);
	bmpinfoheader[8] = (unsigned char) (imgHeight);
	bmpinfoheader[9] = (unsigned char) (imgHeight >> 8);
	bmpinfoheader[10] = (unsigned char) (imgHeight >> 16);
	bmpinfoheader[11] = (unsigned char) (imgHeight >> 24);

	f = fopen(fileName.c_str(), "wb");
	fwrite(bmpfileheader, 1, 14, f);
	fwrite(bmpinfoheader, 1, 40, f);
	for (int i = 0; i < imgHeight; i++) {
		fwrite(img + (imgWidth * (imgHeight - i - 1) * 3), 3, imgWidth, f);
		fwrite(bmppad, 1, (4 - (imgWidth * 3) % 4) % 4, f);
	}
	fclose(f);
	free(img);
}

void ResAnalysisHelper::generateRGBMatrix(
		std::vector<std::vector<int> >& labelMatrix,
		std::vector<std::vector<double> >& red,
		std::vector<std::vector<double> >& green,
		std::vector<std::vector<double> >& blue) {
	// check if input size is valid.
	if (labelMatrix.size() == 0) {
		throw SceException("input label matrix must not be empty",
				OutputAnalysisDataException);
	}

	uint dimX = labelMatrix.size();
	uint dimY = labelMatrix[0].size();

	red.resize(dimY);
	green.resize(dimY);
	blue.resize(dimY);
	for (uint i = 0; i < dimY; i++) {
		red[i].resize(dimX);
		green[i].resize(dimX);
		blue[i].resize(dimX);
	}

	std::tr1::unordered_set<int> cellLabelSet;
	for (uint i = 0; i < dimX; i++) {
		for (uint j = 0; j < dimY; j++) {
			cellLabelSet.insert(labelMatrix[i][j]);
		}
	}

	int labelSetSize = cellLabelSet.size();

	for (uint i = 0; i < dimX; i++) {
		for (uint j = 0; j < dimY; j++) {
			transformToRGB(labelMatrix[i][j], labelSetSize, red[j][i],
					green[j][i], blue[j][i]);
		}
	}
}

void ResAnalysisHelper::transformToRGB(int& labelValue, int& maxLabelValue,
		double& rValue, double& gValue, double& bValue) {
	double relPos = (double) labelValue / (double) maxLabelValue;
	// set as between Red and Green
	if (relPos >= 0 && relPos < 1.0 / 3) {
		rValue = 1.0 - relPos * 3;
		gValue = relPos * 3;
		bValue = 0;
	}
	// set as between Green and Blue
	else if (relPos >= 1.0 / 3 && relPos < 2.0 / 3) {
		rValue = 0.0;
		gValue = (relPos - 1.0 / 3) * 3 * 3;
		bValue = 1.0 - (relPos - 1.0 / 3) * 3;
	}
	// set as between Blue and Red
	// slightly bigger than 1.0 to prevent possible numerical error.
	else if (relPos >= 2.0 / 3 && relPos <= 1.000000001) {
		rValue = 1.0 - (relPos - 2.0 / 3) * 3;
		gValue = 0.0;
		bValue = (relPos - 2.0 / 3) * 3;
	}
	// set as white
	else {
		rValue = 1.0;
		gValue = 1.0;
		bValue = 1.0;
	}
	if (rValue < 0) {
		rValue = 0;
	}
	if (gValue < 0) {
		gValue = 0;
	}
	if (bValue < 0) {
		bValue = 0;
	}
}

void ResAnalysisHelper::outputStat_PolygonCounting(std::string fileName,
		uint step, std::vector<std::vector<int> >& labelMatrix,
		std::vector<double> growthProVec) {
	// check if input size is valid.
	if (labelMatrix.size() == 0) {
		throw SceException("input label matrix must not be empty",
				OutputAnalysisDataException);
	}

	uint dimX = labelMatrix.size();
	uint dimY = labelMatrix[0].size();

	std::tr1::unordered_set<int> cellLabelSet;
	int maxLabel = -1;
	for (uint i = 0; i < dimX; i++) {
		for (uint j = 0; j < dimY; j++) {
			cellLabelSet.insert(labelMatrix[i][j]);
			if (labelMatrix[i][j] > maxLabel) {
				maxLabel = labelMatrix[i][j];
			}
		}
	}
	uint labelSetSize = cellLabelSet.size() - 1;
	std::cout << "Updating polygon stat file " << fileName << ", max label = "
			<< maxLabel << ", labelSetsize = " << labelSetSize << std::endl;
	assert(maxLabel == (int )labelSetSize - 1);

	std::tr1::unordered_set<int> neighborCellSet[labelSetSize];
	std::tr1::unordered_set<int> boundaryCellsSet;
	for (uint i = 0; i < dimX; i++) {
		for (uint j = 0; j < dimY; j++) {
			int label = labelMatrix[i][j];
			if (label == -1) {
				continue;
			}
			if (i != 0) {
				if (labelMatrix[i - 1][j] == -1) {
					boundaryCellsSet.insert(label);
				}
				neighborCellSet[label].insert(labelMatrix[i - 1][j]);
			}
			if (i != dimX - 1) {
				if (labelMatrix[i + 1][j] == -1) {
					boundaryCellsSet.insert(label);
				}
				neighborCellSet[label].insert(labelMatrix[i + 1][j]);
			}
			if (j != 0) {
				if (labelMatrix[i][j - 1] == -1) {
					boundaryCellsSet.insert(label);
				}
				neighborCellSet[label].insert(labelMatrix[i][j - 1]);
			}
			if (j != dimY - 1) {
				if (labelMatrix[i][j + 1] == -1) {
					boundaryCellsSet.insert(label);
				}
				neighborCellSet[label].insert(labelMatrix[i][j + 1]);
			}
		}
	}

	uint sideCount[labelSetSize];
	for (uint i = 0; i < labelSetSize; i++) {
		// the cell is not a boundary cell
		if (boundaryCellsSet.find(i) == boundaryCellsSet.end()) {
			// minus one because I need to remove the label of the cell itself.
			sideCount[i] = neighborCellSet[i].size() - 1;
		} else {
			// for boundary cells, assign value -1.
			sideCount[i] = -1;
		}
	}

	// Size =0 means the vector is not supplied and we do not need to consider mitiotic difference
	if (growthProVec.size() == 0) {
		std::tr1::unordered_map<int, int> polygonCountingRes;
		std::tr1::unordered_map<int, int>::iterator it;
		for (uint i = 0; i < labelSetSize; i++) {
			it = polygonCountingRes.find(sideCount[i]);
			if (it == polygonCountingRes.end()) {
				polygonCountingRes.insert(std::pair<int, int>(sideCount[i], 1));
			} else {
				it->second = it->second + 1;
			}
		}

		// -1 comes from boundary cells and should be removed.
		it = polygonCountingRes.find(-1);
		if (it != polygonCountingRes.end()) {
			polygonCountingRes.erase(it);
		}

		ofstream ofs;
		ofs.open(fileName.c_str(), ios::app);
		ofs << step << " ";
		for (it = polygonCountingRes.begin(); it != polygonCountingRes.end();
				++it) {
			ofs << it->first << "," << it->second << " ";
		}
		ofs << std::endl;
		ofs.close();
	} else {
		assert(growthProVec.size() == labelSetSize);
		double mitiThreshold = globalConfigVars.getConfigValue(
				"MitioticThreshold").toDouble();
		std::tr1::unordered_map<int, int> polyCountResMiti, polyCountResNonMiti;
		std::tr1::unordered_map<int, int>::iterator it;
		for (uint i = 0; i < labelSetSize; i++) {
			if (growthProVec[i] > mitiThreshold) {
				it = polyCountResMiti.find(sideCount[i]);
				if (it == polyCountResMiti.end()) {
					polyCountResMiti.insert(
							std::pair<int, int>(sideCount[i], 1));
				} else {
					it->second = it->second + 1;
				}
			} else {
				it = polyCountResNonMiti.find(sideCount[i]);
				if (it == polyCountResNonMiti.end()) {
					polyCountResNonMiti.insert(
							std::pair<int, int>(sideCount[i], 1));
				} else {
					it->second = it->second + 1;
				}
			}
		}

		// -1 comes from boundary cells and should be removed.
		it = polyCountResMiti.find(-1);
		if (it != polyCountResMiti.end()) {
			polyCountResMiti.erase(it);
		}
		it = polyCountResNonMiti.find(-1);
		if (it != polyCountResNonMiti.end()) {
			polyCountResNonMiti.erase(it);
		}

		ofstream ofs;
		ofs.open(fileName.c_str(), ios::app);
		ofs << step << " ";
		for (it = polyCountResMiti.begin(); it != polyCountResMiti.end();
				++it) {
			ofs << it->first << "," << it->second << " ";
		}
		ofs << "# ";
		for (it = polyCountResNonMiti.begin(); it != polyCountResNonMiti.end();
				++it) {
			ofs << it->first << "," << it->second << " ";
		}
		ofs << std::endl;
		ofs.close();
	}

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

	//effectiveRange_single = globalConfigVars.getConfigValue(
	//		"Pixel_Para_Effective_Range_Single").toDouble();

	allowedAbsoluteError = globalConfigVars.getConfigValue(
			"Pixel_Para_Allowed_Error").toDouble();
}

ResAnalysisHelper::~ResAnalysisHelper() {
}
