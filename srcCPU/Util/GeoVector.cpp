/*
 * GeoVector.cpp
 *
 *  Created on: Sep 13, 2013
 *      Author: wsun2
 */

#include "GeoVector.h"

CVector::CVector() {
	x = y = z = 0.0;
}

CVector::CVector(double a) {
	x = y = z = a;
}

CVector::CVector(double a, double b, double c) {
	x = a;
	y = b;
	z = c;
}

double Modul(const CVector& a) {
	return sqrt(a * a);
}

CVector CVector::getUnitVector(double tolerance) {
	CVector result;
	if (fabs(this->GetX()) < tolerance && fabs(this->GetY()) < tolerance
			&& fabs(this->GetZ()) < tolerance) {
		result = CVector(1, 0, 0);
		throw UnitVectorCalculationException(
				"Warning ! two points are too close when calculating unit "
						"vector of a vector");
	} else {
		result = *this/ Modul(*this);
	}
	return result;
}

void CVector::Print() const {
	std::cout << "(" << x << "," << y << "," << z << ")\n";
}
void CVector::FilePrint(const char *str) {
	ofstream myfile;
	myfile.open(str, ios_base::out | ios_base::app);
	myfile << "(" << x << "," << y << "," << z << ")\n";
	myfile.close();
}

void CVector::Printf(char *s) {
	std::cout << s << "(" << x << "," << y << "," << z << ")\n";
}

CVector operator+(const CVector& a, const CVector& b) {
	return CVector(a.x + b.x, a.y + b.y, a.z + b.z);
}

CVector operator-(const CVector& a, const CVector& b) {
	return CVector(a.x - b.x, a.y - b.y, a.z - b.z);
}

double operator*(const CVector& a, const CVector& b) {
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}

CVector operator*(const double& a, const CVector& b) {
	return CVector(a * b.x, a * b.y, a * b.z);
}

CVector operator*(const CVector& a, const double& b) {
	return CVector(a.x * b, a.y * b, a.z * b);
}

CVector operator/(const double& a, const CVector& b) {
	return CVector(a / b.x, a / b.y, a / b.z);
}

CVector operator/(const CVector& a, const double& b) {
	return CVector(a.x / b, a.y / b, a.z / b);
}

CVector Cross(const CVector& a, const CVector& b) {
	return CVector(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
			a.x * b.y - a.y * b.x);
}

ostream& operator<<(ostream& os, const CVector& vec){
	os<<"("<<vec.x<<","<<vec.y<<","<<vec.z<<")";
	return os;
}

Cvector::Cvector(void) {
}

Cvector::~Cvector(void) {
}

