#include "SceECM.h"
// task: frequency of plotting the ECM should be imported. Right now is given explicitly
__constant__ double sceInterCell_ECM[5]; 
__constant__ double wLCPara_ECM[4]; 
__constant__ double restLenECMAdhSpringGPU  ;  
__constant__ double maxLenECMAdhSpringGPU ;
__constant__ double kAdhECMGPU ;

namespace patch{
	template <typename  T> std::string to_string (const T& n) 
	{
	std:: ostringstream stm ; 
	stm << n ; 
	return stm.str() ; 
	}
}




__device__
double calMorse_ECM(const double& linkLength ) {

	double forceValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		forceValue = 0;
	} else {
		forceValue = -sceInterCell_ECM[0] / sceInterCell_ECM[2]
				* exp(-linkLength / sceInterCell_ECM[2])
				+ sceInterCell_ECM[1] / sceInterCell_ECM[3]
						* exp(-linkLength / sceInterCell_ECM[3]);
//		if (forceValue > 0) {
//			forceValue = 0;
//		}
	}

	return (forceValue) ; 
}

__device__
double calWLC_ECM(const double& linkLength ) {

	double x=linkLength/wLCPara_ECM[0] ;
	return (wLCPara_ECM[1]*( 6*x+ ( x*x*(3.0-2*x))/( (1-x)*(1-x) ) )
	       -wLCPara_ECM[2]/pow(linkLength,wLCPara_ECM[3]) ) ; 	
}

__device__
bool IsValidAdhPair(const double& dist ) {
		if (dist > restLenECMAdhSpringGPU  && dist < maxLenECMAdhSpringGPU){ 
			return true ;
		}
		else {
			return false ;
			}
	}
__device__
bool IsValidAdhPairForNotInitPhase(const double& dist ) {
		if (dist > restLenECMAdhSpringGPU){ 
			return true ;
		}
		else {
			return false ;
			}
	}


__device__
double  CalAdhECM(const double& dist ) {
	return (kAdhECMGPU*(dist-restLenECMAdhSpringGPU)); 
	// in the function IsValid pair, distance already checked to be greater than neutral length
			}


void SceECM::Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM) {

maxAllNodePerCell=maxAllNodePerCellECM ; 
maxMembrNodePerCell= maxMembrNodePerCellECM ; 
maxTotalNodes=maxTotalNodesECM ; //Ali 


std::fstream readCoord_ECM ;
std::fstream readInput_ECM ;  
int numberNodes_ECM ; 
double tmpPosX_ECM,tmpPosY_ECM ; 
vector<double> posXIni_ECM,posYIni_ECM ;
 
readCoord_ECM.open("./resources/coordinate_ECM8.txt") ;
if (readCoord_ECM.is_open()) {
	cout << "ECM coordinates file opened successfully" <<endl ; 
}
else {
	cout << "ECM coordinates file is not opened successfully" << endl ; 
}




readCoord_ECM>>numberNodes_ECM ;
for (int i=0 ; i<numberNodes_ECM ; i++){
	readCoord_ECM>>tmpPosX_ECM>>tmpPosY_ECM ; 
	posXIni_ECM.push_back(tmpPosX_ECM) ; 
	posYIni_ECM.push_back(tmpPosY_ECM) ; 
}



readInput_ECM.open("./resources/input_ECM.txt") ;
if (readInput_ECM.is_open()) {
	cout << "ECM Mech input opened successfully" <<endl ; 
}
else {
	cout << "ECM Mech input is not opened successfully" << endl ; 
}


 for (int i=0 ; i<5; i++) {
 	readInput_ECM>> mechPara_ECM.sceInterCellCPU_ECM[i] ; //=39.0 ; 
 }
 
 readInput_ECM>>restLenECMSpring ;
 readInput_ECM>>eCMLinSpringStiff ;    
 readInput_ECM>>restLenECMAdhSpring  ;  
 readInput_ECM>>maxLenECMAdhSpring ;
 readInput_ECM>>kAdhECM ;
 for ( int i=0 ; i<4 ; i++) {
	readInput_ECM>>mechPara_ECM.wLCParaCPU_ECM[i] ;
 }    




cout<< "numer of ECM nodes is"<< numberNodes_ECM <<endl ; 
for (int i=0 ; i<posXIni_ECM.size() ; i++){
	cout << "ECM nodes read in cpu"<<endl; 
	cout << posXIni_ECM[i] <<posYIni_ECM[i]<<endl;  ; 
}




 for (int i=0 ; i<5; i++) {
	cout <<"Morse parameter number"<<i<<" is " <<mechPara_ECM.sceInterCellCPU_ECM[i]<<endl ; 

} 

 cout <<"rest length of ECM spring is "<<restLenECMSpring<<endl ;   

 cout <<"ECM spring stiffness is "<<eCMLinSpringStiff<<endl ; 

 cout <<"ECM Membrane neutral adhesion length is "<<restLenECMAdhSpring<<endl ;  
 cout <<"ECM Membrane max adhesion length is "<<maxLenECMAdhSpring<<endl ;
 cout <<"ECM Membrane adhesion stiffness is "<<kAdhECM<<endl ;
 cout << "ECM only applies adhesvie force" << endl ; 

for ( int i=0 ; i<4 ; i++) {
	cout<<"wLC parameter "<< i << " is "<<mechPara_ECM.wLCParaCPU_ECM[i]<<endl ;  ;
}    

cudaMemcpyToSymbol(sceInterCell_ECM,mechPara_ECM.sceInterCellCPU_ECM
			,5*sizeof(double)); 
cudaMemcpyToSymbol(wLCPara_ECM,mechPara_ECM.wLCParaCPU_ECM
			,4*sizeof(double)); 

cudaMemcpyToSymbol(restLenECMAdhSpringGPU, &restLenECMAdhSpring,sizeof(double));

cudaMemcpyToSymbol(maxLenECMAdhSpringGPU, &maxLenECMAdhSpring,sizeof(double));
cudaMemcpyToSymbol(kAdhECMGPU, &kAdhECM,sizeof(double));

lastPrintECM=1000000 ; // large number 
outputFrameECM=0 ; 
numNodesECM= numberNodes_ECM ; //(eCMMaxX-eCMMinX)/eCMMinDist ; 




indexECM.resize(numNodesECM,0) ;

nodeECMLocX.resize(numNodesECM,0.0) ;
nodeECMLocY.resize(numNodesECM,0.0) ;

linSpringForceECMX.resize(numNodesECM,0.0); 
linSpringForceECMY.resize(numNodesECM,0.0); 
linSpringAvgTension.resize(numNodesECM,0.0); 


bendSpringForceECMX.resize(numNodesECM,0.0); 
bendSpringForceECMY.resize(numNodesECM,0.0); 
 
memMorseForceECMX.resize(numNodesECM,0.0); 
memMorseForceECMY.resize(numNodesECM,0.0);
 
fBendCenterX.resize(numNodesECM,0.0); 
fBendCenterY.resize(numNodesECM,0.0); 
fBendLeftX.resize(numNodesECM,0.0); 
fBendLeftY.resize(numNodesECM,0.0); 
fBendRightX.resize(numNodesECM,0.0); 
fBendRightY.resize(numNodesECM,0.0); 
 
totalForceECMX.resize(numNodesECM,0.0); 
totalForceECMY.resize(numNodesECM,0.0);

memNodeType.resize(maxTotalNodes,notAssigned1) ; 

thrust::sequence (indexECM.begin(),indexECM.begin()+numNodesECM);
 
thrust::copy(posXIni_ECM.begin(),posXIni_ECM.end(),nodeECMLocX.begin()) ; 
thrust::copy(posYIni_ECM.begin(),posYIni_ECM.end(),nodeECMLocY.begin()) ; 
for (int i=0;  i<nodeECMLocX.size() ;  i++) {
	cout<< nodeECMLocX[i]<<endl; 
} 


} //initilaization function finished






void SceECM:: ApplyECMConstrain(int totalNodeCountForActiveCellsECM, double curTime, double dt, double Damp_Coef, bool cellPolar, bool subCellPolar, bool isInitPhase){   

thrust::counting_iterator<int> iBegin(0) ; 
nodeDeviceTmpLocX.resize(totalNodeCountForActiveCellsECM,0.0) ;
nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
//isBasalMemNode.resize(totalNodeCountForActiveCellsECM,false) ;
adhPairECM_Cell.resize(totalNodeCountForActiveCellsECM,-1) ;
 
thrust::copy(nodeDeviceLocX.begin(),nodeDeviceLocX.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocX.begin()) ; 
thrust::copy(nodeDeviceLocY.begin(),nodeDeviceLocY.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocY.begin()) ; 
cout << " max all node per cell in ECM module is " << maxAllNodePerCell << endl ; 
cout<< "max membrane node per cell in ECM module is " << maxMembrNodePerCell<< endl ; 

double eCMBendStiff=0.0 ; // need to be an input

//if (cellPolar) {eCMLinSpringStiff=100 ; }
//if (subCellPolar) {eCMLinSpringStiff=100 ; }


double* nodeECMLocXAddr= thrust::raw_pointer_cast (
			&nodeECMLocX[0]) ; 
double* nodeECMLocYAddr= thrust::raw_pointer_cast (
			&nodeECMLocY[0]) ; 



// move the nodes of epithelial cells 

 thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					nodeDeviceTmpLocX.begin(),
					nodeDeviceTmpLocY.begin(), 
					nodeIsActive_Cell.begin(),
					memNodeType.begin()
					)), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					 nodeDeviceTmpLocX.begin(),
                     nodeDeviceTmpLocY.begin(),
					 nodeIsActive_Cell.begin(),
					 memNodeType.begin()
					 ))+totalNodeCountForActiveCellsECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeDeviceLocX.begin(),
					nodeDeviceLocY.begin(),
					memNodeType.begin(),
					adhPairECM_Cell.begin())),
				MoveNodes2_Cell(nodeECMLocXAddr,nodeECMLocYAddr,maxMembrNodePerCell,numNodesECM,dt,Damp_Coef,isInitPhase));


double* nodeCellLocXAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocX[0]) ; 
double* nodeCellLocYAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocY[0]) ;
 
bool* nodeIsActive_CellAddr= thrust::raw_pointer_cast (
			&nodeIsActive_Cell[0]) ; 

int * adhPairECM_CellAddr= thrust::raw_pointer_cast (
			&adhPairECM_Cell[0]) ; 

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
					linSpringForceECMY.begin(),
					linSpringAvgTension.begin())),
				LinSpringForceECM(numNodesECM,restLenECMSpring,eCMLinSpringStiff,nodeECMLocXAddr,nodeECMLocYAddr));

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
					fBendCenterX.begin(),
					fBendCenterY.begin(),
					fBendLeftX.begin(),
					fBendLeftY.begin(),
					fBendRightX.begin(),
					fBendRightY.begin())),
				CalBendECM(nodeECMLocXAddr,nodeECMLocYAddr,numNodesECM,eCMBendStiff));

double* fBendLeftXAddr= thrust::raw_pointer_cast (
			&fBendLeftX[0]) ; 
double* fBendLeftYAddr= thrust::raw_pointer_cast (
			&fBendLeftY[0]) ; 
double* fBendRightXAddr= thrust::raw_pointer_cast (
			&fBendRightX[0]) ; 
double* fBendRightYAddr= thrust::raw_pointer_cast (
			&fBendRightY[0]) ; 


thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					indexECM.begin(),
					fBendCenterX.begin(),
					fBendCenterY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 indexECM.begin(),
					 fBendCenterX.begin(),
                                         fBendCenterY.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					bendSpringForceECMX.begin(),
					bendSpringForceECMY.begin())),
				SumBendForce(fBendLeftXAddr,fBendLeftYAddr,fBendRightXAddr,fBendRightYAddr,numNodesECM));


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
				MorseAndAdhForceECM(totalNodeCountForActiveCellsECM,maxAllNodePerCell,maxMembrNodePerCell,nodeCellLocXAddr,nodeCellLocYAddr,nodeIsActive_CellAddr,adhPairECM_CellAddr));


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



nodeECMTmpLocX.resize(numNodesECM,0.0) ;
nodeECMTmpLocY.resize(numNodesECM,0.0) ;
 
thrust::copy(nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,nodeECMTmpLocX.begin()) ; 
thrust::copy(nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,nodeECMTmpLocY.begin()) ; 

double Damp_CoefECM=Damp_Coef ; 
// move the nodes of ECM 
thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMTmpLocX.begin(),
					nodeECMTmpLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					indexECM.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMTmpLocX.begin(),
					nodeECMTmpLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					indexECM.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())),
				MoveNodeECM(dt,Damp_CoefECM,numNodesECM));




lastPrintECM=lastPrintECM+1 ; 
               if (lastPrintECM>=5000) { 
			outputFrameECM++ ; 
			lastPrintECM=0 ; 
			std::string vtkFileName = "ECM_" + patch::to_string(outputFrameECM-1) + ".vtk";
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
			 ECMOut<< "CELLS " << nodeECMLocX.size()<< " " << 3 *nodeECMLocX.size()<< std::endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << 2 << " " << indexECM[i] << " "
				<< indexECM[i+1] << std::endl;
			}
			ECMOut << 2 << " " << indexECM[nodeECMLocX.size()-1] << " "<< indexECM[0] << std::endl; //last point to the first point


			ECMOut << "CELL_TYPES " << nodeECMLocX.size()<< endl;
			for (uint i = 0; i < nodeECMLocX.size() ; i++) {
				ECMOut << "3" << endl;
			}
			ECMOut << "POINT_DATA "<<nodeECMLocX.size() <<endl ; 
			ECMOut << "SCALARS Avg_Tension " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<linSpringAvgTension[i] <<endl ; 
			}
			ECMOut.close(); 
			}

}





