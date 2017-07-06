#include "SceECM.h"

/*
struct InitECMX 
{
  const double  a ; 
 
  InitECMX (double _a) : a(_a) {}

  __host__ __device__ 

   double operator() (const double & x , const double & y) const {
   return (x*a-50) ; 
  }
} ; 
*/
void SceECM::applyECMForce_M() {
double Ali=5 ; 
} 



void SceECM::Initialize() {

double eCMMinX= -50;
double eCMMaxX= 50;
// eCMY=24.0;
double eCMMinDist=0.04;
eCMY=23.7 ; 
int  numNodesECM= (eCMMaxX-eCMMinX)/eCMMinDist ; 
thrust:: device_vector<double> tmp ; 
cout<<" I am inside ECM initialization0  " << endl ; 
tmp.resize(numNodesECM,0.0) ;
nodeECMLocX.resize(numNodesECM,0.0) ;
nodeECMLocY.resize(numNodesECM,eCMY) ;

thrust::sequence (tmp.begin(),tmp.begin()+numNodesECM);

cout<<" I am inside ECM initialization1 " << endl ; 
cout<<" I am inside ECM initialization2 " << endl ; 
 
cout<<" I am inside ECM initialization3 " << endl ; 
thrust::fill (nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,0.0); 
thrust::fill (nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,eCMY); 
thrust::transform(tmp.begin(),tmp.begin()+numNodesECM,nodeECMLocY.begin(),nodeECMLocX.begin(),InitECMLoc(eCMMinX,eCMMinDist));

for (int i=0;  i<nodeECMLocX.size() ;  i++) {
cout<< nodeECMLocX[i]<<endl; 
} 
cout<<" I am inside ECM initialization4  " << endl ; 
}

void SceECM:: ApplyECMConstrain(int totalNodeCountForActiveCellsECM){   

thrust::counting_iterator<int> iBegin(0) ; 
nodeDeviceTmpLocX.resize(totalNodeCountForActiveCellsECM,0.0) ;
nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
 
thrust::copy(nodeDeviceLocX.begin(),nodeDeviceLocX.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocX.begin()) ; 
thrust::copy(nodeDeviceLocY.begin(),nodeDeviceLocY.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocY.begin()) ; 


double yBound=10 ;
 eCMY=23.7 ;  
int maxAllNodePerCell=360 ;
int maxMembrNodePerCell=280 ;  

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					nodeDeviceTmpLocX.begin(),
					nodeDeviceTmpLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					 nodeDeviceTmpLocX.begin(),
                                         nodeDeviceTmpLocY.begin()))+totalNodeCountForActiveCellsECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeDeviceLocX.begin(),
					nodeDeviceLocY.begin())),
				MyFunctor2(eCMY,yBound,maxMembrNodePerCell)); 


}

