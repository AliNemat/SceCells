/*
 * CellInitHelper.h
 *
 *  Created on: Sep 22, 2013
 *      Author: wsun2
 */

#ifndef CELLINITHELPER_H_
#define CELLINITHELPER_H_

#include <vector>
#include "GeoVector.h"
#include <cmath>
#include "commonData.h"
#include "ConfigParser.h"
#include <assert.h>
#include <time.h>
#include <stdlib.h>
#include <string>
#include "MeshGen.h"

const double numericalErrorEps = 1.0e-10;

using namespace std;

/*
 * This class helps the simulation domain to determine
 * center position of a cell and its cell type.
 */
class CellPlacementInfo {
public:
	CVector centerLocation;
	SceNodeType cellNodeType;
};

class CellInitHelper2 {
public:
	std::vector<CellPlacementInfo> obtainCentersInCircle(double radius,
			int precision, Criteria criteria);
};

struct CellInitHelperException: public std::exception {
	std::string errorMessage;
	CellInitHelperException(std::string errMsg) :
			errorMessage(errMsg) {
	}
	~CellInitHelperException() throw () {
	}
	const char* what() const throw () {
		return errorMessage.c_str();
	}
};

class CellInitHelper {
	CVector getPointGivenAngle(double currentAngle, double r,
			CVector centerPos);

	void generateRandomAngles(vector<double> &randomAngles,
			int initProfileNodeSize);
	void generateCellInitNodeInfo(vector<CVector> &initPos,
			std::string meshInput);
	void generateCellInitNodeInfo_v2(vector<CVector> &initPos);
	void generateECMInitNodeInfo(vector<CVector> &initECMNodePoss,
			int initNodeCountPerECM);
	void generateECMCenters(vector<CVector> &ECMCenters,
			vector<CVector> &CellCenters, vector<CVector> &bdryNodes);

	bool anyECMCenterTooClose(vector<CVector> &ecmCenters, CVector position);
	bool anyCellCenterTooClose(vector<CVector> &cellCenters, CVector position);
	bool anyBoundaryNodeTooClose(vector<CVector> &bdryNodes, CVector position);
	bool isInitNodesInitializedFlag;

	double getRandomNum(double min, double max);
	vector<CVector> generateInitCellNodes();
	vector<CVector> attemptGeenerateInitCellNodes();
	bool isPositionQualify(vector<CVector> &poss);

	void initInternalBdry();
	vector<CVector> internalBdryPts;

	void transformRawCartData(CartilageRawData &cartRawData, CartPara &cartPara,
			std::vector<CVector> &initNodePos);

	bool isMXType_v2(CVector position);

public:

	CellInitHelper();

	virtual ~CellInitHelper();

	vector<CellPlacementInfo> obtainPreciseCellInfoArray(double interval,
			double deformRatio);

	vector<CVector> rotate2D(vector<CVector> &initECMNodePoss, double angle);

	void generateBoundaryCellNodesArray(vector<CVector> &bdryNodes,
			double interval);
	void generateThreeInputCellInfoArrays(vector<CVector> &bdryNodes,
			vector<CVector> &FNMCellCenters, vector<CVector> &MXCellCenters,
			double cellCenterInterval, double bdryNodeInterval);
	void initInputsFromCellInfoArray(vector<SceNodeType> &cellTypes,
			vector<uint> &numOfInitNodesOfCells,
			vector<double> &initBdryCellNodePosX,
			vector<double> &initBdryCellNodePosY,
			vector<double> &initFNMCellNodePosX,
			vector<double> &initFNMCellNodePosY,
			vector<double> &initMXCellNodePosX,
			vector<double> &initMXCellNodePosY, vector<CVector> &bdryNodes,
			vector<CVector> &FNMCellCenters, vector<CVector> &MXCellCenters,
			vector<CVector> &initCellNodePoss);

	//RawDataInput generateRawInput(std::string meshInput);
	//SimulationInitData generateInput(std::string meshInput);

	vector<CVector> generateCircleCentersInDisk(double diskRadius,
			double circleRadius);
	vector<CVector> generateCircleCentersInDisk2(double diskRadius,
			double circleRadius);
	RawDataInput generateDiskRawInput(std::string meshInput);
	RawDataInput generateRawInput_stab();
	RawDataInput generateRawInputWithProfile(
			std::vector<CVector> &cellCenterPoss, bool isInnerBdryIncluded =
					true);
	RawDataInput generateRawInput_V2(std::vector<CVector> &cellCenterPoss);
	SimulationInitData generateDiskInput(std::string meshInput);

	SimulationInitData initInputsV2(RawDataInput &rawData);
	SimulationInitData_V2 initInputsV3(RawDataInput &rawData);
};

#endif /* CELLINITHELPER_H_ */
