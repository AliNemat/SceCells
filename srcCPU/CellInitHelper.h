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
#include "SceCells.h"
#include <time.h>
#include <stdlib.h>

const double numericalErrorEps = 1.0e-10;

using namespace std;

struct RawDataInput {
	vector<CVector> bdryNodes;
	vector<CVector> profileNodes;
	vector<CVector> FNMCellCenters;
	vector<CVector> MXCellCenters;
	vector<CVector> ECMCenters;
	vector<double> ECMAngles;
	vector<CVector> initCellNodePoss;
	vector<CVector> initECMNodePoss;
};

struct SimulationInitData {
	std::vector<CellType> cellTypes;
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

/*
 * This class helps the simulation domain to determine
 * center position of a cell and its cell type.
 */
class CellPlacementInfo {
public:
	CVector centerLocation;
	CellType cellType;
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

	enum ConditionType {
		Outside, Left, Right, Up, Down
	};
	class BoundaryLine {
	public:
		ConditionType condiType;
		virtual bool isFitCondition(CVector position) {
			return false;
		}
		virtual void printWhoAmI() {
			cout
					<< "I am an uninitialized boundaryLine and you shouldn't use me"
					<< endl;
		}
	};
	// dealing with general arc is unnecessary.
	// Here arc is initialized in a special way, other methods are not acceptable.
	class Arc: public BoundaryLine {
	public:
		void printWhoAmI() {
			cout << "I am an arc my radius is " << this->r
					<< ", and my center position is: ";
			this->centerPos.Print();
		}
		CVector centerPos;
		double r;
		CVector p1;
		CVector p2;
		Arc() {
			centerPos = CVector(0, 0, 0);
			r = 0;
		}
		Arc(CVector point1, CVector point2, CVector point3) {
			if (point1.x != point3.x) {
				throw CellInitHelperException(
						"point 1 and point 3 must have same x coordinate!");
			}
			if (fabs((point1.y + point3.y) / 2 - point2.y) > 1.0e-9) {
				throw CellInitHelperException(
						"point 2 must locate in the middle of point 1 and point 3");
			}
			double tmpDiff1 = fabs(point2.x - point1.x);
			double tmpDiff2 = fabs(point2.y - point1.y);
			r = (tmpDiff1 * tmpDiff1 + tmpDiff2 * tmpDiff2) / 2.0 / tmpDiff1;
			centerPos = point2 + CVector(-r, 0, 0);
			p1 = point1;
			p2 = point3;
		}
		bool isWithInRange(CVector position) {
			if (this->condiType == Outside) {
				double minY = p1.y;
				double maxY = p2.y;
				if (p2.y < p1.y) {
					minY = p2.y;
					maxY = p1.y;
				}
				if (position.y <= maxY && position.y >= minY) {
					return true;
				} else {
					return false;
				}
			} else {
				throw CellInitHelperException("Arc cannot deal with direction "
						"other than outside");
			}
		}
		Arc getOutside(double distance) {
			if (distance < 0) {
				throw CellInitHelperException(
						"distance must be larger than 0!");
			}
			double newR = this->r + distance;
			Arc result;
			result.centerPos = this->centerPos;
			result.r = newR;
			result.condiType = this->condiType;
			return result;
		}
		bool isVeryCloseToLine(CVector position) {
			CVector vecToCenter = position - centerPos;
			double dist = Modul(vecToCenter);
			if (fabs(dist - this->r) < numericalErrorEps) {
				return true;
			} else {
				return false;
			}
		}
		bool isInside(CVector position) {
			CVector vecToCenter = position - centerPos;
			double dist = Modul(vecToCenter);
			if (dist > this->r) {
				//cout << "not inside arc range" << endl;
				return false;
			} else {
				double difference = dist - this->r;
				//cout << " inside arc range" << "distance = " << difference
				//		<< endl;
				return true;
			}
		}
		bool isFitCondition(CVector position) {

			//cout
			//		<< "determine if position fits the condition of outside of an arc"
			//		<< endl;
			//cout << "arc radius =" << this->r << "and arc center: ";
			if (!isWithInRange(position)) {
				//cout << "Not in range, count as pass!" << endl;
				return true;
			}
			//this->centerPos.Print();
			if (isVeryCloseToLine(position)) {
				//cout << "arc is very close to line, count as fit condition!"
				//		<< endl;
				return true;
			}
			return !isInside(position);
		}
	};

	/* General purpose class for parabola is difficult and not worthwhile.
	 * Use two straight lines instead
	 */
	/*
	 // y = a*(x-b)^2+c  , a must be negative
	 class DownwardsParabola: public BoundaryLine {
	 double a, b, c;
	 bool isInsideOf(CVector position) {
	 if (position.y < c) {
	 return false;
	 } else {
	 double temp = sqrt((position.y - c) / a);
	 double xLeft = b - temp;
	 double xRight = b + temp;
	 if (position.x >= xLeft && position.x <= xRight) {
	 return true;
	 } else {
	 return false;
	 }
	 }
	 }
	 };
	 */

	class StraightLineEquationVertical: public BoundaryLine {
	public:
		void printWhoAmI() {
			cout << "I am a vertical line, my x position is " << this->xPos
					<< endl;
		}
		double xPos;
		StraightLineEquationVertical(double xP) {
			xPos = xP;
		}
		bool isVeryCloseToLine(CVector position) {
			if (fabs(position.x - this->xPos) < numericalErrorEps) {
				return true;
			} else {
				return false;
			}
		}
		bool isRightOfLine(CVector position) {
			if (isVeryCloseToLine(position)) {
				return true;
			}
			return position.x > xPos;
		}
		StraightLineEquationVertical getRightOf(double distance) {
			StraightLineEquationVertical result = *this;
			result.xPos = this->xPos + distance;
			result.condiType = this->condiType;
			return result;
		}
		bool isFitCondition(CVector position) {
			if (isVeryCloseToLine(position)) {
				return true;
			}
			if (condiType == Right) {
				return isRightOfLine(position);
			} else if (condiType == Left) {

				return !isRightOfLine(position);
			} else {
				string errorMessage1 = " condition type is" + condiType;
				string errorMessage2 =
						" Vertical straight line cannot handle condition other than Right and Left";
				throw CellInitHelperException(errorMessage1);
			}
		}
	};
	class StraightLineEquationNoneVertical: public BoundaryLine {

	public:
		double k, b;
		CVector p1, p2;
		// equation is of form y-kx-b = 0
		// because we don't plan to include corner cases of
		// k = infinite or 0 this is sufficient
		void printWhoAmI() {
			//cout << "I am a none vertical line, my k is " << this->k
			//		<< ", and my b is: " << this->b << endl;
		}

		StraightLineEquationNoneVertical(CVector n1, CVector n2) {
			k = (n2.y - n1.y) / (n2.x - n1.x);
			b = n1.y - k * n1.x;
			p1 = n1;
			p2 = n2;
		}
		StraightLineEquationNoneVertical() {
			k = 0;
			b = 0;
		}
		bool isWithInRange(CVector position) {
			if (this->condiType == Left || this->condiType == Right) {
				double minY = p1.y;
				double maxY = p2.y;
				if (p2.y < p1.y) {
					minY = p2.y;
					maxY = p1.y;
				}
				if (position.y <= maxY && position.y >= minY) {
					return true;
				} else {
					return false;
				}
			} else if (this->condiType == Up || this->condiType == Down) {
				double minX = p1.x;
				double maxX = p2.x;
				if (p2.x < p1.x) {
					minX = p2.x;
					maxX = p1.x;
				}
				if (position.x <= maxX && position.x >= minX) {
					return true;
				} else {
					return false;
				}
			} else {
				throw CellInitHelperException(
						"none vertical line cannot deal with direction"
								" other than left, right, down and up");
			}

		}
		bool isVeryCloseToLine(CVector position) {
			double yOfSameX = k * position.x + b;
			if (fabs(yOfSameX - position.y) < numericalErrorEps) {
				return true;
			} else {
				return false;
			}
		}
		bool isLeftOfLine(CVector position) {
			double xOfSameY = (position.y - b) / k;
			//cout << "determining if point is left of line:";
			//position.Print();
			//cout << "line k = " << k << " b = " << b << endl;
			if (xOfSameY > position.x) {
				//cout << "x of same y = " << xOfSameY << endl;
				//cout << "result: true" << endl;
				return true;
			} else {
				//cout << "x of same y = " << xOfSameY << endl;
				//cout << "result: false" << endl;
				return false;
			}
		}

		bool isRightOfLine(CVector position) {
			return !isLeftOfLine(position);
		}
		bool isDownOfLine(CVector position) {
			double yOfSameX = k * position.x + b;
			//cout << "determining if a point is downside of a line: k = " << k
			//		<< " b = " << b << endl;
			//position.Print();
			if (yOfSameX > position.y) {
				return true;
			} else {
				return false;
			}
		}
		bool isUpOfLine(CVector position) {
			return !isDownOfLine(position);
		}
		double getDistance(CVector position) {
			double result = fabs(position.y - k * position.x - b)
					/ sqrt(1.0 * 1.0 + k * k);
			return result;
		}
		StraightLineEquationNoneVertical getUpOf(double distance) {
			StraightLineEquationNoneVertical result =
			//StraightLineEquationNoneVertical();
					*this;
			result.k = this->k;
			result.b = this->b + (sqrt(1 + this->k * this->k) * distance);
			result.condiType = this->condiType;
			return result;
		}
		StraightLineEquationNoneVertical getDownOf(double distance) {
			StraightLineEquationNoneVertical result = *this;
			result.k = this->k;
			result.b = this->b - (sqrt(1 + this->k * this->k) * distance);
			result.condiType = this->condiType;
			return result;
		}
		bool isFitCondition(CVector position) {
			if (!isWithInRange(position)) {
				//cout << "Not in range, count as pass!" << endl;
				return true;
			}
			if (isVeryCloseToLine(position)) {
				//cout << "Point Very Close to Line, count as pass!" << endl;
				return true;
			}
			if (condiType == Right) {
				return isRightOfLine(position);
			} else if (condiType == Left) {
				return isLeftOfLine(position);
			} else if (condiType == Up) {
				return isUpOfLine(position);
			} else if (condiType == Down) {
				return isDownOfLine(position);
			} else {
				string errorMessage1 = " condition type is" + condiType;
				throw CellInitHelperException(
						errorMessage1
								+ ", straight line cannot handle condition other than Up,Down,Right and Left");
			}
		}
	};

	StraightLineEquationNoneVertical linEqn1, linEqn2, linEqn3, linEqn4;
	bool isInsideFNMRegion(CVector position);
	bool isInsideMXRegion(CVector position);
	double getMinDistanceFromFNMBorder(CVector position);
	double getMinDistanceFromMXBorder(CVector position);

	vector<CVector> initPoints;
	vector<CVector> transformedPoints;
	vector<BoundaryLine *> boundaryLines;
	vector<BoundaryLine *> boundariesForCellCenter;
	vector<BoundaryLine *> internalBoundaryLines;
	void initPrecisionBoundaryPoints();
	void transformBoundaryPoints();
	// init both boundary lines and internal boundary lines
	// given interval, we also need to generate a set of boundaries for cell center
	void initBoundaryLines(double interval);
	double getStartingXGivenY(double yPos);
	vector<CVector> getCellCentersInside(double interval);

	//bool isInsidePreciseFNMRegion(CVector position);
	//bool isInsidePreciseMXRegion(CVector position);

	bool isCellCenterInsidePreciseRegion(CVector position);
	bool isMXType(CVector position);

	CVector getPointGivenAngle(double currentAngle, double r,
			CVector centerPos);

	void generateProfileNodesArray(vector<CVector> &initProfileNodes,
			double profileNodeInterval);

	void generateRandomAngles(vector<double> &randomAngles,
			int initProfileNodeSize);
	void generateCellInitNodeInfo(vector<CVector> &initPos,
			std::string meshInput);
	void generateECMInitNodeInfo(vector<CVector> &initECMNodePoss,
			int initNodeCountPerECM);
	void generateECMCenters(vector<CVector> &ECMCenters,
			vector<CVector> &CellCenters, vector<CVector> &bdryNodes);

	bool anyECMCenterTooClose(vector<CVector> &ecmCenters, CVector position);
	bool anyCellCenterTooClose(vector<CVector> &cellCenters, CVector position);
	bool anyBoundaryNodeTooClose(vector<CVector> &bdryNodes, CVector position);
	bool isInitNodesInitializedFlag;

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
	void initInputsFromCellInfoArray(vector<CellType> &cellTypes,
			vector<uint> &numOfInitNodesOfCells,
			vector<double> &initBdryCellNodePosX,
			vector<double> &initBdryCellNodePosY,
			vector<double> &initFNMCellNodePosX,
			vector<double> &initFNMCellNodePosY,
			vector<double> &initMXCellNodePosX,
			vector<double> &initMXCellNodePosY, vector<CVector> &bdryNodes,
			vector<CVector> &FNMCellCenters, vector<CVector> &MXCellCenters,
			vector<CVector> &initCellNodePoss);

	RawDataInput generateRawInput(std::string meshInput);
	SimulationInitData initInputsV2(RawDataInput &rawData);
	SimulationInitData generateInput(std::string meshInput);
};

#endif /* CELLINITHELPER_H_ */
