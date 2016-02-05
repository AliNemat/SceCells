/*
 * Point2D.h
 *
 *  Created on: Aug 20, 2014
 *      Author: wsun2
 */

#ifndef POINT2D_H_
#define POINT2D_H_

namespace GEOMETRY {

class Point2D {
	bool isOnBdry;
	double x;
	double y;
public:
	Point2D();
	Point2D(double inputX, double inputY);

	virtual ~Point2D();

	double getX() const {
		return x;
	}

	void setX(double x) {
		this->x = x;
	}

	double getY() const {
		return y;
	}

	void setY(double y) {
		this->y = y;
	}

	bool isIsOnBdry() const {
		return isOnBdry;
	}

	void setIsOnBdry(bool isOnBdry) {
		this->isOnBdry = isOnBdry;
	}

        void Assign_M2(double inputx1, double inputy1) ; 
                 

};

} /* namespace GEOMETRY */
#endif /* POINT2D_H_ */
