#include "commonData.h"

const std::string SceInputPoint::delimiter = " ";

std::string toString(SceExceptionType type) {
	std::string result("Undefined type");
	switch (type) {
	case BaseException:
		result = "Base exception (usage discouraged)";
		break;
	case InputInitException:
		result = "Initialization Exception";
		break;
	case ConfigFileNotFound:
		result = "Configuration file not found Exception";
		break;
	case ConfigValueException:
		result = "Configuration value Exception";
		break;
	case OutputAnalysisDataException:
		result = "Exception while outputting analysis result";
		break;
	case FileIOException:
		result = "Exception while file IO";
		break;
	case MemoryInvalidAccess:
		result = "Exception while executing memory operation";
		break;
	case InvalidInput:
		result = "Exception processing method input";
		break;
	case AlgorithmBug:
		result = "Possible bug in algorithm used";
		break;
	}
	return result;
}

std::string toString(SceNodeType type) {
	std::string result("Undefined type");
	switch (type) {
	case Boundary:
		result = "Boundary";
		break;
	case Profile:
		result = "Profile";
		break;
	case ECM:
		result = "ECM";
		break;
	case FNM:
		result = "FNM";
		break;
	case MX:
		result = "MX";
		break;
	case Cart:
		result = "Cartilage";
		break;
	case Base:
		result = "Base";
		break;
	}
	return result;
}

double nodeTypeToScale(SceNodeType type) {
	double result = 0.0;
	switch (type) {
	case Boundary:
		result = 1;
		break;
	case Profile:
		result = 7;
		break;
	case ECM:
		result = 3;
		break;
	case FNM:
		result = 4;
		break;
	case MX:
		result = 5;
		break;
	case Cart:
		result = 2;
		break;
	case Base:
		result = 6;
		break;
	}
	return result;
}

void inputInitialData::addNewPoints(std::vector<SceInputPoint>& newPoints) {
	for (unsigned int i = 0; i < newPoints.size(); i++) {
		inputPoints.push_back(newPoints[i]);
	}
}

bool AnimationCriteria::isPairQualify(uint seq1, uint seq2, double x1,
		double y1, double z1, SceNodeType t1, uint r1, double x2, double y2,
		double z2, SceNodeType t2, uint r2) {
	bool condi1 = false, condi2 = false;
	if (t1 == t2 && r1 == r2) {
		if (t1 == Boundary || t1 == ECM) {
			if (abs(seq1 - seq2) == 1) {
				condi1 = true;
			}
		} else if (t1 == MX || t1 == FNM) {
			condi1 = true;
		}
		if (condi1) {
			double dist = compuDistHost(x1, y1, z1, x2, y2, z2);
			if (dist < defaultEffectiveDistance) {
				condi2 = true;
			}
		}
	}
	return condi1 && condi2;
}

double compuDistHost(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2)
					+ (zPos - zPos2) * (zPos - zPos2));
}

void VtkAnimationData::outputVtkAni(std::string scriptNameBase, int rank) {
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string scriptNameRank = ss.str();
	std::string vtkFileName = scriptNameBase + scriptNameRank + ".vtk";
	std::cout << "start to create vtk file" << vtkFileName << std::endl;
	std::ofstream fs;
	fs.open(vtkFileName.c_str());
	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << "Lines and points representing subcelluar element cells "
			<< std::endl;
	fs << "ASCII" << std::endl;
	fs << std::endl;
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << pointsAniData.size() << " float" << std::endl;
	for (uint i = 0; i < pointsAniData.size(); i++) {
		fs << pointsAniData[i].pos.x << " " << pointsAniData[i].pos.y << " "
				<< pointsAniData[i].pos.z << std::endl;
	}

	fs << std::endl;
	fs << "CELLS " << linksAniData.size() << " " << 3 * linksAniData.size()
			<< std::endl;
	for (uint i = 0; i < linksAniData.size(); i++) {
		fs << 2 << " " << linksAniData[i].node1Index << " "
				<< linksAniData[i].node2Index << std::endl;
	}
	fs << "CELL_TYPES " << linksAniData.size() << endl;
	for (uint i = 0; i < linksAniData.size(); i++) {
		fs << "3" << endl;
	}
	fs << "POINT_DATA " << pointsAniData.size() << endl;
	fs << "SCALARS relative_compression float" << endl;
	fs << "LOOKUP_TABLE default" << endl;

	for (uint i = 0; i < pointsAniData.size(); i++) {
		fs << pointsAniData[i].colorScale << endl;
	}

	fs.close();
}

std::vector<double> getArrayXComp(std::vector<CVector>& nodePosVec) {
	std::vector<double> result;
	for (uint i = 0; i < nodePosVec.size(); i++) {
		result.push_back(nodePosVec[i].GetX());
	}
	return result;
}

std::vector<double> getArrayYComp(std::vector<CVector>& nodePosVec) {
	std::vector<double> result;
	for (uint i = 0; i < nodePosVec.size(); i++) {
		result.push_back(nodePosVec[i].GetY());
	}
	return result;
}

std::vector<double> getArrayZComp(std::vector<CVector>& nodePosVec) {
	std::vector<double> result;
	for (uint i = 0; i < nodePosVec.size(); i++) {
		result.push_back(nodePosVec[i].GetZ());
	}
	return result;
}

SimulationType parseTypeFromConfig(int configValue) {
	if (configValue == 0) {
		return Beak;
	} else if (configValue == 1) {
		return Disc;
	} else {
		throw SceException("Simulation type in config file is not defined",
				ConfigValueException);
	}
}

bool valueToType(int value) {
	if (value != 0) {
		return true;
	} else {
		return false;
	}
}

uint findClosestArrIndexGivenPos(std::vector<CVector>& vecArr, CVector& pos) {
	if (vecArr.size() == 0) {
		std::string errorMsg =
				"while finding closest array index, the input array has zero elements";
		throw SceException(errorMsg, InvalidInput);
	} else {
		uint index = 0;
		CVector tmpVec = vecArr[0] - pos;
		double minDis = tmpVec.getModul();
		for (uint i = 1; i < vecArr.size(); i++) {
			tmpVec = vecArr[i] - pos;
			double dist = tmpVec.getModul();
			if (dist < minDis) {
				index = i;
				minDis = dist;
			}
		}
		return index;
	}
}
