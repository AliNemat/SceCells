#include <fstream>
#include <sstream>
#include "Signal2D.h"
#include <unistd.h>
#include <time.h>  
#include <stdlib.h>
#include <stdio.h>

namespace patch
{
template  < typename T> std::string to_string (const T & n)
{
std::ostringstream stm ; 
stm << n ; 
return stm.str(); 
}



}


void Signal::Initialize (uint maxAllNodePerCell, uint maxMembrNodePerCell, uint maxTotalNodes, uint maxCellCount) {
	
	this->maxAllNodePerCell=maxAllNodePerCell ; 
	this->maxMembrNodePerCell=maxMembrNodePerCell ;
	this->maxCellCount=maxCellCount ; 
	periodCount=0 ; 
	nodeLocXHost.resize(maxTotalNodes, 0.0) ; 
	nodeLocYHost.resize(maxTotalNodes, 0.0) ;
	nodeIsActiveHost.resize(maxTotalNodes,false) ; 
	cellCenterX.resize(maxCellCount,0.0) ; 
	cellCenterY.resize(maxCellCount,0.0) ; 
	dppLevel.resize(maxCellCount,0.0) ; 

	minResol=0.1 ;// max distance of the first imported coordinate of DPP from the tissue center to accept it for that cell    
	resol=501 ; // the number of imported DPP values
	cout << "I am at the end of signal initialization function" << endl ; 
	cout << "size of node is active in signal module initialization is " << nodeIsActiveHost.size() << endl ; 
	cout << "max of all nodes per cell in signal module initialization is " << maxAllNodePerCell << endl ; 
}
void Signal::updateSignal(double minX, double maxX, double minY, double maxY, double curTime, int maxTotalNumActiveNodes, int numActiveCells)  {
	this->maxX=maxX;
	this->maxY=maxY;
	this->minX=minX;
	this->minY=minY;
	this->curTime=curTime ; 
	this->maxTotalNumActiveNodes=maxTotalNumActiveNodes ; 
	this->numActiveCells=numActiveCells ; 
	cout << "I am in update signal started" << std::endl ; 
	exportGeometryInfo()	; 
	importSignalInfoCellLevel()	; 
	processSignalInfoCellLevel()	; 

}

void Signal::exportGeometryInfo() {
 	double Center_X=minX+0.5*(maxX-minX); 
 	int cellRank ; 


	srand(time(NULL));

 	double max_Rx=max (maxX-Center_X,Center_X-minX) ; 
	cout << "The information is read from "<< periodCount<<" data" <<endl ; 
		
		 

	cout << "size of node is active in signal module is " << nodeIsActiveHost.size() << endl ; 
	cout << "max total number of active nodes in signal module is " << maxTotalNumActiveNodes  << endl ; 
	cout << "max of all nodes per cell in signal module is " << maxAllNodePerCell << endl ; 
	//std :: string  txtFileName="ExportTisuProp_" + patch::to_string(periodCount)+".txt" ; 
	std :: string  txtFileName="ExportCellProp_" + patch::to_string(periodCount)+".txt" ; 
	ofstream ExportOut ; 
	ExportOut.open(txtFileName.c_str()); 
	//ExportOut << "Time (s), Tissue_CenterX(micro meter),Max_Length_X(micro meter)"<<endl; 
	//ExportOut<<curTime<<","<<Center_X<<","<<max_Rx<<endl   ;
	//ExportOut << "CellRank,CellCenterX,CellCenterY"<<endl;
	for (int k=0; k<numActiveCells; k++){
		ExportOut<<k<<","<<cellCenterX[k]<<","<<cellCenterY[k]<<endl   ;
	}
 	//ExportOut << "CellRank,MembraneNodeX,MembraneNodeY"<<endl;

	for ( uint i=0 ; i< maxTotalNumActiveNodes ; i++) {

		cellRank= i/maxAllNodePerCell ; 
		if (nodeIsActiveHost[i] && (i%maxAllNodePerCell)<maxMembrNodePerCell) {
			ExportOut<<cellRank<<","<<nodeLocXHost[i]<<","<<nodeLocYHost[i]<<endl   ;
		}
	}
	int lastLineIndicator=123456789 ; 
	ExportOut<<lastLineIndicator<<","<<lastLineIndicator<<","<<lastLineIndicator<<endl   ;
	ExportOut.flush() ;
	cout << "I exported  the data for signaling model"<< endl ; 
	ExportOut.close() ;  
		
}


void Signal::importSignalInfoCellLevel() {


		float dppDist,dppLevelTmp ; 
		dppLevelV.clear(); 
		dppDistV.clear();

		std:: string importDppFileName= "Dpp_cell_T" + patch::to_string(periodCount) + ".txt" ;


		std:: ifstream inputDpp ;
		bool fileIsOpen=false ;
			
		cout << "The file name I am looking for is " << importDppFileName <<endl ;
		
		sleep(300) ; 
		while (fileIsOpen==false) {
			inputDpp.open(importDppFileName.c_str()) ;
			if (inputDpp.is_open()){
				cout<< "dpp file opened successfully"<<endl; 
				cout << "the opened file namei is " << importDppFileName <<endl ;
				fileIsOpen=true ; 
				 clock_t t;
				t = clock();
				cout << "I start to sleep. Time is"  <<endl ;
				sleep(30) ; 
				 t = clock() - t;
				cout << "Sleep takes "<< ((float)t)/CLOCKS_PER_SEC  <<endl ;
				cout << "the opened file name is " << importDppFileName <<endl ;
			}
			else {
		//		cout << "failed openining dpp file"<<endl ;
			}	
		}
		if (inputDpp.good()) {
		cout << " I passed opening the file in the while loop"<< endl ;
		}

 		periodCount+= 1 ;// abs(floor((curTime-InitTimeStage)/exchPeriod)) ;
		for (int i=0; i<numActiveCells ; i++) {
			inputDpp >> dppLevelTmp ;
			cout<<"zeroth dpp is"<<dppLevelTmp<< endl ; 
			dppLevelV.push_back(dppLevelTmp) ;  
		}	
		cout <<"first dpp value is"<< dppLevelV.at(0)<< endl ; 	
}





void Signal::importSignalInfoTissueLevel() {


		float dppDist,dppLevelTmp ; 
		dppLevelV.clear(); 
		dppDistV.clear();

		std:: string importDppFileName= "DppImport_" + patch::to_string(periodCount) + ".txt" ;


		std:: ifstream inputDpp ;
		bool fileIsOpen=false ;
			
		cout << "the file name I am looking for is " << importDppFileName <<endl ;
		
		sleep(200) ; 
		while (fileIsOpen==false) {
			inputDpp.open(importDppFileName.c_str()) ;
			if (inputDpp.is_open()){
				cout<< "dpp file opened successfully"<<endl; 
				cout << "the opened file nameis " << importDppFileName <<endl ;
				fileIsOpen=true ; 
				 clock_t t;
				t = clock();
				cout << "I start to sleep time is"  <<endl ;
				sleep(30) ; 
				 t = clock() - t;
				cout << "Sleep takes"<< ((float)t)/CLOCKS_PER_SEC  <<endl ;
				cout << "the opened file name is " << importDppFileName <<endl ;
			}
			else {
				//cout << "failed openining dpp file"<<endl ;
			}	
		}
		if (inputDpp.good()) {
		cout << " I passed opening the file in the while loop"<< endl ;
		}
 		periodCount+= 1 ;// abs(floor((curTime-InitTimeStage)/exchPeriod)) ;
		for (int i=0; i<resol ; i++) {
			inputDpp >> dppDist >> dppLevelTmp ;
			//cout<<"zeroth dpp is"<<dppDist<<dppLevel<< endl ; 
			dppDistV.push_back(dppDist) ; 
			dppLevelV.push_back(dppLevelTmp) ;  
		}	
		cout <<"first dpp value is"<< dppLevelV.at(0)<< endl ; 	
}


void Signal::processSignalInfoTissueLevel() {


	vector<double> dppLevels_Cell ;
	dppLevels_Cell.clear() ;


 	double Center_X=minX+0.5*(maxX-minX);
 	double Center_Y=minY+0.5*(maxY-minY);

  		for (int k=0; k<numActiveCells; k++){
     			double DistXCell=abs(cellCenterX[k]-Center_X); 
     			int StoredI=-1 ;  //We will figure out that there is an error if stays -1  
  
        	for (int i=0 ; i<resol ; i++){
        		if (DistXCell<dppDistV[i] || DistXCell<minResol) {

        			StoredI=i ;
           			break;  
        		}
				if (i==resol-1) {
        			StoredI=i ;
				}
       		}  
       		dppLevels_Cell.push_back(dppLevelV[StoredI]);
       	}

       	for (int k=numActiveCells; k<maxCellCount ; k++){
       		dppLevels_Cell.push_back(0.0) ;   //these cells are not active
       	}


		for (int k=0 ;  k<numActiveCells; k++) {
			double distYAbs=abs (cellCenterY[k]-Center_Y); 
          		
			double dummy = (static_cast<double>(rand()) / RAND_MAX);
			double ranNum = NormalCDFInverse(dummy);	

			dppLevel[k]=dppLevels_Cell[k]; //+
					  //dppLevels_Cell[k]*(0.1*sin(0.2*3.141592*distYAbs)+0.12*ranNum); 
		}	
}

void Signal::processSignalInfoCellLevel() {


	vector<double> dppLevels_Cell ;
	dppLevels_Cell.clear() ;

 	double Center_Y=minY+0.5*(maxY-minY);


  		for (int k=0; k<numActiveCells; k++){
        	dppLevels_Cell.push_back(dppLevelV[k]);
       	}

       	for (int k=numActiveCells; k<maxCellCount ; k++){
       		dppLevels_Cell.push_back(0.0) ;   //these cells are not active
       	}


		for (int k=0 ;  k<numActiveCells; k++) {
			double distYAbs=abs (cellCenterY[k]-Center_Y); 
          		
			double dummy = (static_cast<double>(rand()) / RAND_MAX);
			double ranNum = NormalCDFInverse(dummy);	

			dppLevel[k]=dppLevels_Cell[k] ;  //+ dppLevels_Cell[k]*(0.1*sin(0.2*3.141592*distYAbs)+0.12*ranNum); 
		}	
}





double NormalCDFInverse(double p) {
	
	if (p < 0.5) {
		return -RationalApproximation( sqrt(-2.0*log(p)));
	}
	else {
		return RationalApproximation( sqrt(-2.0*log(1-p)));
	}
}

double RationalApproximation(double t) {

	double c[] = {2.515517, 0.802853, 0.010328};
	double d[] = {1.432788, 0.189269, 0.001308};
	return (t - ((c[2]*t + c[1])*t + c[0]) / (((d[2]*t + d[1])*t + d[0])*t + 1.0));

}


