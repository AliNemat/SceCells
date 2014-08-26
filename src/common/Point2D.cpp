/*
 * Point2D.cpp
 *
 *  Created on: Aug 20, 2014
 *      Author: wsun2
 */

#include "Point2D.h"

namespace GEOMETRY {

Point2D::Point2D() {
	x = 0.0;
	y = 0.0;
}

Point2D::~Point2D() {

}

} /* namespace GEOMETRY */

GEOMETRY::Point2D::Point2D(double inputX, double inputY) {
	x = inputX;
	y = inputY;
}
