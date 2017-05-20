#include <vector>
#include <iostream>
#include <cmath>
#include "GeoVector.h"
using namespace std ; 

vector<double> updateSignal(vector<double> & dppLevels, const vector<CVector> & cellCentersHost, int cellMax, double MinX, double MaxX, double MinY, double MaxY, double dt, double InitTimeStage, double curTime, int & plotSignal) ; 
