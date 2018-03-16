#include <vector>
#include <iostream>
#include <cmath>
#include <string>
#include <fstream> 
#include "GeoVector.h"
#include "commonData.h" 
using namespace std ; 

double NormalCDFInverse(double p);
double RationalApproximation(double t);

class Signal {

	public :
	uint maxAllNodePerCell ;
	uint maxMembrNodePerCell ; 
	uint maxTotalNodes ; 
	uint maxCellCount ; 
	int periodCount ; 

	std::vector<bool> nodeIsActiveHost ; 
	std::vector<double> nodeLocXHost, nodeLocYHost ; 

	vector<double> updateSignal(const vector<CVector> & cellCentersHost, double MinX, double MaxX, double MinY, double MaxY, double curTime, int maxTotalNumActiveNodes) ;

	void Initialize(uint maxAllNodePerCell, uint maxMembrNodePerCellECM, uint maxTotalNodes, uint maxCellCount) ; 

} ; 

