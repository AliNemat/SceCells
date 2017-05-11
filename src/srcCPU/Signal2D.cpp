#include <fstream>
#include "Signal2D.h"

void printAli (int a, int b)
{
int c ;  
c=a+b ; 
std::cout << a+b << std::endl ; 
}


std::vector<double> updateSignal(vector<CVector> CellCentersHost) {

cout << "I am in update signal" << std::endl ; 
vector <double> DPPLevel ; 

   for (int i=0 ; i<CellCentersHost.size() ; i++) {
         DPPLevel.push_back(5.0) ; 
        }

 
return  DPPLevel ; 
}
