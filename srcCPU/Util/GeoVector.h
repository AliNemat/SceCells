/*
 * GeoVector.h
 *
 *  Created on: Sep 13, 2013
 *      Author: wsun2
 *      different from vector in C++ std, this CVector class
 *      is used for geometry
 */

#ifndef GEOVECTOR_H_
#define GEOVECTOR_H_

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <exception>

using namespace std;

struct UnitVectorCalculationException: public std::exception {
	std::string errorMessage;
	UnitVectorCalculationException(std::string errMsg) :
			errorMessage(errMsg) {
	}
	~UnitVectorCalculationException() throw () {
	}
	const char* what() const throw () {
		return errorMessage.c_str();
	}
};

class CVector {
public:
	double x, y, z;
	CVector();
	CVector(double a);
	CVector(double a, double b, double c);
	double GetX() {
		return x;
	}
	double GetY() {
		return y;
	}
	double GetZ() {
		return z;
	}
	double getModul() {
		return sqrt(x * x + y * y + z * z);
	}
	void Printf(char *name);
	void FilePrint(const char *);
	CVector getUnitVector(double tolerance = 1.0e-14);

	friend CVector operator+(const CVector& a, const CVector& b);
	friend CVector operator-(const CVector& a, const CVector& b);
	friend double operator*(const CVector& a, const CVector& b);
	friend CVector operator*(const double& a, const CVector& b);
	friend CVector operator*(const CVector& a, const double& b);
	friend CVector operator/(const double& a, const CVector& b);
	friend CVector operator/(const CVector& a, const double& b);
	friend CVector Cross(const CVector& a, const CVector& b);
	friend double Modul(const CVector& a);
	friend ostream& operator<<(ostream& os, const CVector& dt);
	void Print() const;
};

class Cvector {
public:
	Cvector(void);
	~Cvector(void);
};

#endif /* GEOVECTOR_H_ */
