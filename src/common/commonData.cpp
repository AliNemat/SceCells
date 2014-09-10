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
