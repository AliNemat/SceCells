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
int  numNodesECM= (eCMMaxX-eCMMinX)/eCMMinDist ; 
thrust:: device_vector<double> tmp ; 
thrust::sequence (tmp.begin(),tmp.begin()+numNodesECM); 
thrust::fill (nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,eCMY); 
thrust::fill (nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,0.0); 
thrust::transform(tmp.begin(),tmp.begin()+numNodesECM,nodeECMLocX.begin(),nodeECMLocX.begin(),InitECMLocX(eCMMinX,eCMMinDist)); 
}

void SceECM:: ApplyECMConstrain(int totalNodeCountForActiveCellsECM){   

//int totalNodeCountForActiveCells=allocPara_M.currentActiveCellCount*allocPara_M.maxAllNodePerCell ;   
//int totalnodecountforactivecells=   
  //                               cells.getAllocParaM().currentActiveCellCount+
    //                             cells.getAllocParaM().maxAllNodePerCell ; 
            
//int totalNodeCountForActiveCellsECM=1680;   

nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
 
thrust::copy(nodeDeviceLocY.begin(),nodeDeviceLocY.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocY.begin()) ; 

cout << "the temporary values are "<< endl ; 
//thrust::copy (nodeDeviceTmpLocY.begin(),nodeDeviceTmpLocY.begin()+160,std::ostream_iterator<double>(std::cout,"\n")) ; 
cout << "the permanent values are "<< endl ; 
//thrust::copy (nodeDeviceLocY.begin(),nodeDeviceLocY.begin()+160,std::ostream_iterator<double>(std::cout,"\n")) ; 
/*	thrust::transform(nodeDeviceLocX.begin(),
			  nodeDeviceLocX.begin()+totalNodeCountForActiveCells,
			  nodeDeviceTmpLocY.begin(),
			  nodeDeviceLocY.begin(),
			  AddECMForce(eCMY));

*/

double yBound=10 ;
 eCMY=23.5 ;  
thrust:: transform(nodeDeviceTmpLocY.begin(),
                    nodeDeviceTmpLocY.begin()+totalNodeCountForActiveCellsECM,
                    nodeDeviceLocY.begin(),
                    MyFunction(eCMY,yBound));   

cout<< "I am after transform function  tmp value is "<< endl ; 
//thrust::copy (nodeDeviceTmpLocY.begin(),nodeDeviceTmpLocY.begin()+160,std::ostream_iterator<double>(std::cout,"\n")) ; 
cout<< "I am after transform function  permanent value is "<< endl ; 
//thrust::copy (nodeDeviceLocY.begin(),nodeDeviceLocY.begin()+160,std::ostream_iterator<double>(std::cout,"\n")) ; 
//nodeECMLocX.resize(1680,2.0) ; 
//nodeECMLocY.resize(1680,1.0) ; 
//thrust::copy(nodes->getInfoVecs().nodeLocY.begin(), 
  //           nodes->getInfoVecs().nodeLocY.begin()+totalNodeCountForActiveCells,
    //         nodeECMLocX.begin()) ; 



/*	thrust::transform(nodeECMLocY.begin(),
			  nodeECMLocY.begin()+1680,
                          nodeECMLocX.begin(),
			  nodeECMLocY.begin(),
			  AddECMForce(eCMY));
*/
cout << "I am in SceECM.cu"<< endl ; 
}

