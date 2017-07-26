#include "SceECM.h"


namespace patch
{
template <typename  T> std::string to_string (const T& n) 
{
std:: ostringstream stm ; 
stm << n ; 
return stm.str() ; 
}
}



void SceECM::applyECMForce_M() {
double Ali=5 ; 
} 



void SceECM::Initialize() {

double eCMMinX= -50;
double eCMMaxX= 50;
double eCMMinDist=0.04;

lastPrintECM=0 ; 
outputFrameECM=0 ; 
restLenECMSpring=0.03 ; //eCMMinDist ;  
eCMY=23.7 ; 
numNodesECM= (eCMMaxX-eCMMinX)/eCMMinDist ; 
//thrust:: device_vector<double> tmp ; 
cout<<" I am inside ECM initialization0  " << endl ; 
nodeECMLocX.resize(numNodesECM,0.0) ;
nodeECMLocY.resize(numNodesECM,eCMY) ;

indexECM.resize(numNodesECM,0) ;
//tmp.resize(numNodesECM,0.0) ;

//thrust::sequence (tmp.begin(),tmp.begin()+numNodesECM);
thrust::sequence (indexECM.begin(),indexECM.begin()+numNodesECM);



linSpringForceECMX.resize(numNodesECM,0.0); 
linSpringForceECMY.resize(numNodesECM,0.0); 


bendSpringForceECMX.resize(numNodesECM,0.0); 
bendSpringForceECMY.resize(numNodesECM,0.0); 
 
memMorseForceECMX.resize(numNodesECM,0.0); 
memMorseForceECMY.resize(numNodesECM,0.0);
 
totalForceECMX.resize(numNodesECM,0.0); 
totalForceECMY.resize(numNodesECM,0.0);


 
thrust::fill (nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,0.0); 
thrust::fill (nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,eCMY); 
thrust::transform(indexECM.begin(),indexECM.begin()+numNodesECM,nodeECMLocY.begin(),nodeECMLocX.begin(),InitECMLoc(eCMMinX,eCMMinDist));

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


double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 


double eCMLinSpringStiff=200 ;   

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
 

double* nodeCellLocXAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocX[0]) ; 
double* nodeCellLocYAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocY[0]) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                                         nodeECMLocY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin())),
				LinSpringForceECM(numNodesECM,restLenECMSpring,eCMLinSpringStiff,nodeECMLocXAddr,nodeECMLocYAddr));


//double* nodeCellLocXAddr= thrust::raw_pointer_cast (
//			&nodeDeviceTmpLocX[0]) ; 
//double* nodeCellLocYAddr= thrust::raw_pointer_cast (
//			&nodeDeviceTmpLocY[0]) ;

bool* nodeIsActive_CellAddr= thrust::raw_pointer_cast (
			&nodeIsActive_Cell[0]) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					nodeECMLocX.begin(),
					nodeECMLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 nodeECMLocX.begin(),
                                         nodeECMLocY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin())),
				MorseForceECM(totalNodeCountForActiveCellsECM,maxAllNodePerCell,maxMembrNodePerCell,nodeCellLocXAddr,nodeCellLocYAddr,nodeIsActive_CellAddr));


double dummy=0.0 ;

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					linSpringForceECMX.begin(),
					linSpringForceECMY.begin(),
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin(),
					memMorseForceECMX.begin(),
					memMorseForceECMY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					totalForceECMX.begin(),
					totalForceECMY.begin())),
				TotalECMForceCompute(dummy));


double dt=0.05; //0.005  
double dampECM=36.0 ; 

nodeECMTmpLocX.resize(numNodesECM,0.0) ;
nodeECMTmpLocY.resize(numNodesECM,0.0) ;
 
thrust::copy(nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,nodeECMTmpLocX.begin()) ; 
thrust::copy(nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,nodeECMTmpLocY.begin()) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMTmpLocX.begin(),
					nodeECMTmpLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMTmpLocX.begin(),
					nodeECMTmpLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())),
				MoveNodeECM(dt,dampECM));




lastPrintECM=lastPrintECM+1 ; 
               if (lastPrintECM==10000) { 
			outputFrameECM++ ; 
			lastPrintECM=0 ; 
			std::string vtkFileName = "ECM_" + patch::to_string(outputFrameECM) + ".vtk";
			ofstream ECMOut;
			ECMOut.open(vtkFileName.c_str());
			ECMOut<< "# vtk DataFile Version 3.0" << endl;
			ECMOut<< "Result for paraview 2d code" << endl;
			ECMOut << "ASCII" << endl;
			ECMOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
			ECMOut << "POINTS " << nodeECMLocX.size() << " float" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
				<< 0.0 << std::endl;
			}
			ECMOut<< std::endl;
			 ECMOut<< "CELLS " << nodeECMLocX.size()-1<< " " << 3 *(nodeECMLocX.size()-1)<< std::endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << 2 << " " << indexECM[i] << " "
				<< indexECM[i+1] << std::endl;
			}

			ECMOut << "CELL_TYPES " << nodeECMLocX.size()-1<< endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << "3" << endl;
			}

			ECMOut.close(); 
			}

}





