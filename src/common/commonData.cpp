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
		result = 2;
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
	case Base:
		result = 6;
		break;
	}
	return result;
}

void inputInitialData::initFromFile(std::string fileName) {
	std::ifstream infile(fileName.c_str());
	if (!infile.is_open()) {
		throw SceException("Fatal error: input file not found!",
				InputInitException);
	}
	std::string line;
	while (std::getline(infile, line)) {
		SceInputPoint pt(line);
		inputPoints.push_back(pt);
	}
}

void SceInputPoint::initFromString(std::string inputLine) {
	std::string s = inputLine;
	std::vector<std::string> strs;
	size_t pos = 0;
	std::string token;
	while ((pos = s.find(delimiter)) != std::string::npos) {
		token = s.substr(0, pos);
		strs.push_back(token);
		s.erase(0, pos + delimiter.length());
	}
	strs.push_back(s);
	cellRank = atoi(strs[0].c_str());
	cellType = atoi(strs[1].c_str());
	xCoord = atof(strs[2].c_str());
	xCoord = atof(strs[3].c_str());
	xCoord = atof(strs[4].c_str());
}

SceInputPoint::SceInputPoint(std::string inputLine) {
	initFromString(inputLine);
}

void SceInputPoint::outputToString(std::string &outputLine) {
	std::stringstream s;
	s << cellRank << delimiter;
	s << cellType << delimiter;
	s << xCoord << delimiter;
	s << yCoord << delimiter;
	s << zCoord << delimiter;
	s << std::endl;
	outputLine = s.str();
}

void inputInitialData::addNewPoints(std::vector<SceInputPoint>& newPoints) {
	for (unsigned int i = 0; i < newPoints.size(); i++) {
		inputPoints.push_back(newPoints[i]);
	}
}

void inputInitialData::outputToFile(std::string fileName) {
	std::ofstream outFile(fileName.c_str());
	std::string outLine;
	for (unsigned int i = 0; i < inputPoints.size(); i++) {
		inputPoints[i].outputToString(outLine);
		outFile << outLine;
	}
	outFile.close();
}

bool AnimationCriteria::isPairQualify(uint seq1, uint seq2, double x1,
		double y1, double z1, SceNodeType t1, uint r1, double x2, double y2,
		double z2, SceNodeType t2, uint r2) {
	bool condi1 = false, condi2 = false;
	if (t1 == t2 && r1 == r2) {
		if (t1 == Boundary || t1 == Profile || t1 == ECM) {
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
