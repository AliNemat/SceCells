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
double calMorseEnergy_ECM(const double& linkLength ) {

	double energyValue=0.0 ; 
	if (linkLength > sceInterCell_ECM[4]) {
		energyValue = 0;
	} else {
		energyValue = sceInterCell_ECM[0]* exp(-linkLength / sceInterCell_ECM[2])
				    - sceInterCell_ECM[1]* exp(-linkLength / sceInterCell_ECM[3]);
	}

	return (energyValue) ; 
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

__device__
double  CalAdhEnergy(const double& dist ) {
	return (0.5*kAdhECMGPU*(dist-restLenECMAdhSpringGPU)*(dist-restLenECMAdhSpringGPU)); 
	// in the function IsValid pair, distance already checked to be greater than neutral length
			}


EType SceECM:: ConvertStringToEType(string eNodeRead) {
	if (eNodeRead=="perip")  {
		return perip ; 
	}
	else if (eNodeRead=="bc") {
		return bc2 ; 
	}
	else if (eNodeRead=="excm") {
		return excm ; 
	}
	else {
		cout << "Error in defining type of external nodes" << endl ; 
		return excm ;// To just return something to avoid compiler complain 
	}
} 


void SceECM::Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM, int freqPlotData) {

maxAllNodePerCell=maxAllNodePerCellECM ; 
maxMembrNodePerCell= maxMembrNodePerCellECM ; 
maxTotalNodes=maxTotalNodesECM ; //Ali 
this->freqPlotData=freqPlotData ; 


std::fstream readCoord_ECM ;
std::fstream readInput_ECM ;  
int numberNodes_ECM ; 
double tmpPosX_ECM,tmpPosY_ECM ; 
vector<double> posXIni_ECM,posYIni_ECM ;
vector <EType> eNodeVec ;  
readCoord_ECM.open("./resources/coordinate_ECM12.txt") ;
if (readCoord_ECM.is_open()) {
	cout << "ECM coordinates file opened successfully" <<endl ; 
}
else {
	cout << "ECM coordinates file is not opened successfully" << endl ; 
}



string eNodeRead ; 
readCoord_ECM>>numberNodes_ECM ;
for (int i=0 ; i<numberNodes_ECM ; i++){
	readCoord_ECM>>tmpPosX_ECM>>tmpPosY_ECM>>eNodeRead  ; 
	posXIni_ECM.push_back(tmpPosX_ECM) ; 
	posYIni_ECM.push_back(tmpPosY_ECM) ; 
	EType eNode=ConvertStringToEType(eNodeRead) ; 
	eNodeVec.push_back(eNode) ; 
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




cout<< "number of ECM nodes is"<< numberNodes_ECM <<endl ; 
for (int i=0 ; i<posXIni_ECM.size() ; i++){
	cout << "ECM nodes read in cpu"<<endl; 
	cout << posXIni_ECM[i] <<", "<<posYIni_ECM[i]<<", " << eNodeVec[i] <<endl;  ; 
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
peripORexcm.resize(numNodesECM,perip) ;

nodeECMLocX.resize(numNodesECM,0.0) ;
nodeECMLocY.resize(numNodesECM,0.0) ;

stiffLevel.resize(numNodesECM,eCMLinSpringStiff) ;
sponLen.resize(numNodesECM,restLenECMSpring) ;

linSpringForceECMX.resize(numNodesECM,0.0); 
linSpringForceECMY.resize(numNodesECM,0.0); 
linSpringAvgTension.resize(numNodesECM,0.0); 
linSpringEnergy.resize(numNodesECM,0.0); 
morseEnergy.resize(numNodesECM,0.0); 
adhEnergy.resize(numNodesECM,0.0); 


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
//thrust::fill (peripORexcm.begin(),peripORexcm.begin()+604,excm );
//thrust::fill (peripORexcm.begin()+1166,peripORexcm.end(),excm );
//thrust::fill (peripORexcm.begin(),peripORexcm.begin()+int(numNodesECM/2),excm );
//thrust::fill (peripORexcm.begin(),peripORexcm.begin()+int(numNodesECM/2),excm );
 
thrust::copy(posXIni_ECM.begin(),posXIni_ECM.end(),nodeECMLocX.begin()) ; 
thrust::copy(posYIni_ECM.begin(),posYIni_ECM.end(),nodeECMLocY.begin()) ; 
thrust::copy(eNodeVec.begin(),eNodeVec.end(),peripORexcm.begin()) ;

//cout << "GPU level initial coordinates and type of external nodes are: " << endl ; 
//for (int i=0;  i<nodeECMLocX.size() ;  i++) {
//	cout<< nodeECMLocX[i]<<", "<<nodeECMLocY[i]<<", "<<peripORexcm[i] << endl; 
//}

PrintECM(0.0) ; 
for (int i=0 ; i<maxTotalNodes ; i++) {
	int nodeRankPerCell=i%maxAllNodePerCell ;
	if (nodeRankPerCell<7) {
		memNodeType[i]=lateral1 ;
	}
	else if (nodeRankPerCell<21) {
		memNodeType[i]=apical1 ;
	}
	else if (nodeRankPerCell<35) {
		memNodeType[i]=lateral1 ;
	}
	else if (nodeRankPerCell<49) {
		memNodeType[i]=basal1 ;
	}
	else if (nodeRankPerCell<56) {
		memNodeType[i]=lateral1 ;
	}
	else {
		memNodeType[i]=notAssigned1;
	}
}

std::string cSVFileName = "EnergyExport.CSV";
			ofstream EnergyExport ;
			EnergyExport.open(cSVFileName.c_str());

			EnergyExport <<"Time,"<<"TotalMorseEnergyCell," << "TotalAdhEnergyCell,"<<  "TotalMorseEnergy,"<<"TotalAdhEnergy,"<< "TotalLinSpringEnergy" <<" TotalEnergy"<< std::endl;





} //initilaization function finished






void SceECM:: ApplyECMConstrain(int totalNodeCountForActiveCellsECM, double curTime, double dt, double Damp_Coef, bool cellPolar, bool subCellPolar, bool isInitPhase){   

thrust::counting_iterator<int> iBegin(0) ; 
nodeDeviceTmpLocX.resize(totalNodeCountForActiveCellsECM,0.0) ;
nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
//isBasalMemNode.resize(totalNodeCountForActiveCellsECM,false) ;
adhPairECM_Cell.resize(totalNodeCountForActiveCellsECM,-1) ;
 
morseEnergyCell.resize(totalNodeCountForActiveCellsECM,0.0); 
adhEnergyCell.resize(totalNodeCountForActiveCellsECM,0.0); 
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

EType* peripORexcmAddr= thrust::raw_pointer_cast (
			&peripORexcm[0]) ; 

memNodeType.resize(maxTotalNodes,notAssigned1) ; 

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
					adhPairECM_Cell.begin(),
					morseEnergyCell.begin(),
					adhEnergyCell.begin())),
				MoveNodes2_Cell(nodeECMLocXAddr,nodeECMLocYAddr,maxMembrNodePerCell,numNodesECM,dt,Damp_Coef,isInitPhase,peripORexcmAddr,curTime));

totalMorseEnergyCell = thrust::reduce( morseEnergyCell.begin(),morseEnergyCell.begin()+totalNodeCountForActiveCellsECM,(double) 0.0, thrust::plus<double>() ); 
totalAdhEnergyCell   = thrust::reduce( adhEnergyCell.begin()  ,adhEnergyCell.begin()  +totalNodeCountForActiveCellsECM,(double) 0.0, thrust::plus<double>() );

double* nodeCellLocXAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocX[0]) ; 
double* nodeCellLocYAddr= thrust::raw_pointer_cast (
			&nodeDeviceTmpLocY[0]) ;
 
bool* nodeIsActive_CellAddr= thrust::raw_pointer_cast (
			&nodeIsActive_Cell[0]) ; 

int * adhPairECM_CellAddr= thrust::raw_pointer_cast (
			&adhPairECM_Cell[0]) ; 


thrust:: transform (peripORexcm.begin(), peripORexcm.begin()+numNodesECM,thrust::make_zip_iterator (thrust::make_tuple (
											stiffLevel.begin(),sponLen.begin())),MechProp(isInitPhase,eCMLinSpringStiff,restLenECMSpring));


double* stiffLevelAddr=thrust::raw_pointer_cast (
			&stiffLevel[0]) ; 

double*  sponLenAddr =thrust::raw_pointer_cast (
			&sponLen[0]) ; 


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
					linSpringAvgTension.begin(),
					linSpringEnergy.begin())),
				LinSpringForceECM(numNodesECM,restLenECMSpring,eCMLinSpringStiff,nodeECMLocXAddr,nodeECMLocYAddr,stiffLevelAddr,sponLenAddr));


totalLinSpringEnergy = thrust::reduce( linSpringEnergy.begin(),linSpringEnergy.begin()+numNodesECM,(double) 0.0, thrust::plus<double>() ); 
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
					memMorseForceECMY.begin(),
					morseEnergy.begin(),
					adhEnergy.begin())),
				MorseAndAdhForceECM(totalNodeCountForActiveCellsECM,maxAllNodePerCell,maxMembrNodePerCell,nodeCellLocXAddr,nodeCellLocYAddr,nodeIsActive_CellAddr,adhPairECM_CellAddr));

totalMorseEnergy = thrust::reduce( morseEnergy.begin(),morseEnergy.begin()+numNodesECM,(double) 0.0, thrust::plus<double>() ); 
totalAdhEnergy   = thrust::reduce( adhEnergy.begin()  ,adhEnergy.begin()  +numNodesECM,(double) 0.0, thrust::plus<double>() );


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
double Damp_CoefPerip=Damp_Coef ; 
// move the nodes of ECM 
thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMTmpLocX.begin(),
					nodeECMTmpLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					indexECM.begin(),
					peripORexcm.begin())),
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					nodeECMTmpLocX.begin(),
					nodeECMTmpLocY.begin(),
					totalForceECMX.begin(),
					totalForceECMY.begin(),
					indexECM.begin(),
					peripORexcm.begin()))+numNodesECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeECMLocX.begin(),
					nodeECMLocY.begin())),
				MoveNodeECM(dt,Damp_CoefECM,Damp_CoefPerip,numNodesECM,curTime));

//cout << "GPU level first updated coordinates and type of external nodes are: " << endl ; 
//for (int i=0;  i<nodeECMLocX.size() ;  i++) {
//	cout<< nodeECMLocX[i]<<", "<<nodeECMLocY[i]<<", "<<peripORexcm[i] << endl; 
//}

PrintECM(curTime); 

}

void  SceECM:: PrintECM(double curTime) {
		lastPrintECM=lastPrintECM+1 ; 
               if (lastPrintECM>=freqPlotData) { 
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
//			ECMOut << "POINT_DATA "<<nodeECMLocX.size() <<endl ; 
///			ECMOut << "SCALARS Avg_Tension2 " << "float"<< endl;
//			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
//			for (uint i = 0; i < nodeECMLocX.size(); i++) {
//				ECMOut<<linSpringAvgTension[i] <<endl ; 
//			}

//			ECMOut << "POINT_DATA "<<nodeECMLocX.size() <<endl ; 
			ECMOut << "SCALARS Node_Type " << "float"<< endl;
			ECMOut << "LOOKUP_TABLE " << "default"<< endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut<<peripORexcm[i] <<endl ; 
			}

			ECMOut.close();
			// second output file for curvature estimation //
			std::string txtFileName = "ECMExport_" + patch::to_string(outputFrameECM-1) + ".txt";
			ofstream ECMExport ;
			ECMExport.open(txtFileName.c_str());
			ECMExport << "ECM pouch coordinates" << std::endl;

			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				if (peripORexcm[i]==excm) {
					ECMExport<< nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
					<< 0.0 << std::endl;
				}
			}

			ECMExport << "ECM lumen side coordinates" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				if (peripORexcm[i]==perip) {
					ECMExport << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
					<< 0.0 << std::endl;
				}
			}

			ECMExport.close();
			
			double totalEnergyECM=0.5*(totalMorseEnergyCell+totalMorseEnergy)+ 0.5*(totalAdhEnergyCell+totalAdhEnergy)+ 0.5*totalLinSpringEnergy ;


			std::string cSVFileName = "EnergyExport.CSV";
			ofstream EnergyExport ;
			EnergyExport.open(cSVFileName.c_str(),ofstream::app);
			
			//EnergyExport <<"totalMorseEnergyCell " << "totalAdhEnergyCell "<<  "totalMorseEnergy "<<"totalAdhEnergy "<< "totalLinSpringEnergy " << std::endl;
			EnergyExport <<curTime<<","<<totalMorseEnergyCell << "," << totalAdhEnergyCell<< "," << totalMorseEnergy<< "," << totalAdhEnergy<< "," << totalLinSpringEnergy <<"," << totalEnergyECM << std::endl;


			}

}




