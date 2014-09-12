#include <vector>
#include "GeoVector.h"
#include <fstream>
#include <exception>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sstream>
#include <map>

#ifndef COMMONDATA_H_
#define COMMONDATA_H_

typedef unsigned int uint;
typedef std::map<uint, uint> IndexMap;

enum SceExceptionType {
	BaseException, InputInitException, ConfigFileNotFound, ConfigValueException
};

std::string toString(SceExceptionType type);
double compuDistHost(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2);

enum SceNodeType {
	Boundary, Profile, ECM, FNM, MX, Base
};

double nodeTypeToScale(SceNodeType type);

class SceException: public std::exception {
private:
	std::string _message;
	SceExceptionType _exceptionType;
public:
	SceException(const std::string& message) :
			_message(message), _exceptionType(BaseException) {
	}
	SceException(const std::string& message, const SceExceptionType type) :
			_message(message), _exceptionType(type) {
	}
	~SceException() throw () {
	}
	virtual const char* what() const throw () {
		return std::string(
				_message + ", Exception type: " + toString(_exceptionType)).c_str();
	}
};

struct RawDataInput {
	std::vector<CVector> bdryNodes;
	std::vector<CVector> profileNodes;
	std::vector<CVector> FNMCellCenters;
	std::vector<CVector> MXCellCenters;
	std::vector<CVector> ECMCenters;
	std::vector<double> ECMAngles;
	std::vector<CVector> initCellNodePoss;
	std::vector<CVector> initECMNodePoss;
};

/**
 * an data structure that was specifically designed for Beak project.
 */
struct SimulationInitData {
	std::vector<SceNodeType> cellTypes;
	std::vector<uint> numOfInitActiveNodesOfCells;
	std::vector<double> initBdryCellNodePosX;
	std::vector<double> initBdryCellNodePosY;
	std::vector<double> initProfileNodePosX;
	std::vector<double> initProfileNodePosY;
	std::vector<double> initECMNodePosX;
	std::vector<double> initECMNodePosY;
	std::vector<double> initFNMCellNodePosX;
	std::vector<double> initFNMCellNodePosY;
	std::vector<double> initMXCellNodePosX;
	std::vector<double> initMXCellNodePosY;
};

struct SceInputPoint {
	static const std::string delimiter;
	uint cellRank;
	SceNodeType cellType;
	double xCoord;
	double yCoord;
	double zCoord;
	SceInputPoint(std::string inputLine);
	void initFromString(std::string inputLine);
	void outputToString(std::string& outputLine);
};

struct SceMechanicalParameters {

};

struct SceMemoryParameters {

};

struct inputInitialData {
	SceMemoryParameters memParas;
	SceMechanicalParameters mechParas;
	std::vector<SceInputPoint> inputPoints;
	void initFromFile(std::string fileName);
	void addNewPoints(std::vector<SceInputPoint> &newPoints);
	void outputToFile(std::string fileName);
};

/**
 * Data structure that controls the animation criteria.
 */
struct AnimationCriteria {
	// If this varible is set to be true, output stress map;
	// otherwise, output normal animation.
	bool isStressMap;
	// We will only animate links that are close enough.
	double defaultEffectiveDistance;
	// determines if a potential pair is qualified for animation.
	bool isPairQualify(uint seq1, uint seq2, double x1, double y1, double z1,
			SceNodeType t1, uint r1, double x2, double y2, double z2,
			SceNodeType t2, uint r2);

};

/**
 * Necessary information to animate a point in VTK.
 */
struct PointAniData {
	// position of the node.
	CVector pos;
	// In VTK animation software, color scale represents the relative distance to red and blue;
	// bigger value means close to red. smaller value means close to blue.
	double colorScale;
};

/**
 * Necessary information to animate a link in VTK.
 */
struct LinkAniData {
	uint node1Index, node2Index;
};

/**
 * contains all information need to be
 */
struct VtkAnimationData {
	std::vector<PointAniData> pointsAniData;
	std::vector<LinkAniData> linksAniData;
};

#endif /* COMMONDATA_H_ */
