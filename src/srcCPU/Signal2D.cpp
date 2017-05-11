#include <fstream>
#include "Signal2D.h"


std::vector<double> updateSignal(vector<CVector> CellCentersHost, int cellMax) {

cout << "I am in update signal" << std::endl ; 
vector <double> DPPLevel ;
 
for (int i=0 ; i<cellMax ; i++) {
         DPPLevel.push_back(0.0) ; 
        }

   for (int i=0 ; i<CellCentersHost.size() ; i++) {
         DPPLevel[i]=5.0 ; 
        }

 
return  DPPLevel ; 
}
