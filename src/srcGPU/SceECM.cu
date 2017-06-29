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
//thrust::fill (nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,0.0); 
//thrust::fill (nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,eCMY); 
thrust::transform(tmp.begin(),tmp.begin()+numNodesECM,nodeECMLocX.begin(),nodeECMLocX.begin(),InitECMLoc(eCMMinX,eCMMinDist)); 
cout<<" I am inside ECM initialization4  " << endl ; 
}

void SceECM:: ApplyECMConstrain(int totalNodeCountForActiveCellsECM){   

//int totalNodeCountForActiveCells=allocPara_M.currentActiveCellCount*allocPara_M.maxAllNodePerCell ;   
//int totalnodecountforactivecells=   
  //                               cells.getAllocParaM().currentActiveCellCount+
    //                             cells.getAllocParaM().maxAllNodePerCell ; 
            
//int totalNodeCountForActiveCellsECM=1680;   
thrust::counting_iterator<int> iBegin(0) ; 
nodeDeviceTmpLocX.resize(totalNodeCountForActiveCellsECM,0.0) ;
nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
 
thrust::copy(nodeDeviceLocX.begin(),nodeDeviceLocX.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocX.begin()) ; 
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
 eCMY=23.7 ;  
/*thrust:: transform(nodeDeviceTmpLocY.begin(),
                    nodeDeviceTmpLocY.begin()+totalNodeCountForActiveCellsECM,
                    nodeDeviceLocY.begin(),
                    MyFunction(eCMY,yBound));   

*/
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
int maxAllNodePerCell=240 ;
int maxMembrNodePerCell=200 ;  

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

