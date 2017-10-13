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
vector<double> updateSignal(vector<double> & dppLevelsV,const vector<CVector> & CellCentersHost, int cellMax, double MinX, double MaxX, double MinY, double MaxY, double dt, double InitTimeStage, double curTime, int & plotSignal, double & lastTimeExchang, int & periodCount)  {

cout << "I am in update signal started" << std::endl ; 

 vector<double> dppLevels_Cell ;
 double R_x=0.5*(MaxX-MinX); 
 double R_y=0.5*(MaxY-MinY); 
 double Center_X=MinX+0.5*(MaxX-MinX); 
 double Center_Y=MinY+0.5*(MaxY-MinY);
 bool importData=true ; 

srand(time(NULL));

 double max_Rx=max (MaxX-Center_X,Center_X-MinX) ; 
float dppDist,dppLevel ; 
vector<double> dppDistV,dppLevelV ;
//double exchPeriod=2*60*60 ;  
double minResol=0.1 ; 
int resol=501 ; 
if (importData) {
	dppLevels_Cell.clear(); 
	cout << "last time data exchanged is" << lastTimeExchang <<endl ;  
	//lastTimeExchang=lastTimeExchang+dt ;
	cout << "The iformation is read from "<< periodCount<<" data" <<endl ; 
		
		lastTimeExchang=0 ; 
		 

		//exporting data//
		std :: string  txtFileName="ExportTisuProp_" + patch::to_string(periodCount)+".txt" ; 
		ofstream ExportOut ; 
		ExportOut.open(txtFileName.c_str()); 
		ExportOut << "Time (s), Tissue_CenterX(micro meter),Max_Length_X(micro meter)"<<endl; 
		ExportOut<<curTime<<","<<Center_X<<","<<max_Rx<<endl   ;
		ExportOut.flush() ;
		cout << "I export the data"<< endl ; 
		ExportOut.close() ;  
		
		//Importing data //
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
				cout << "the opened file nameis " << importDppFileName <<endl ;
			}
			else {
				cout << "failed openining dpp file"<<endl ;
			}	
		}
		if (inputDpp.good()) {
		cout << " I passed opening the file in the while loop"<< endl ;
		}
 		periodCount+= 1 ;// abs(floor((curTime-InitTimeStage)/exchPeriod)) ;
		for (int i=0; i<resol ; i++) {
			inputDpp >> dppDist >> dppLevel ;
			//cout<<"zeroth dpp is"<<dppDist<<dppLevel<< endl ; 
			dppDistV.push_back(dppDist) ; 
			dppLevelV.push_back(dppLevel) ;  
		}	
		cout <<"first dpp value is"<< dppLevelV.at(0)<< endl ; 	

  		for (int k=0; k<CellCentersHost.size(); k++){
     			double DistXCell=abs(CellCentersHost[k].x-Center_X); 
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

       		for (int k=CellCentersHost.size(); k<cellMax ; k++){
       			dppLevels_Cell.push_back(0.0) ;   //these cells are not active
       		}


		for (int k=0 ;  k<CellCentersHost.size() ; k++) {
			double distYAbs=abs (CellCentersHost[k].y-Center_Y); 
          		
			double dummy = (static_cast<double>(rand()) / RAND_MAX);
			double ranNum = NormalCDFInverse(dummy);	

			dppLevels_Cell[k]=dppLevels_Cell[k]+
					  dppLevels_Cell[k]*(0.1*sin(0.2*3.141592*distYAbs)+0.12*ranNum); 
		}	
//		cout <<"second dpp value is"<< dppLevels_Cell.at(0)<< endl ; 	
		


     // }
    //  else {
//		for (int k=0; k<cellMax ; k++){
  //     			dppLevels_Cell.push_back(0.0) ;   //these cells are not active
    //  		}
    //  }

  }
	//std :: string  normFileName="normDistr_" + patch::to_string(periodCount)+".txt" ; 
	//ofstream ofs(normFileName.c_str());

	//double r, norm;

	//for (int i = 0; i < 1000; i++) {
	//	r = (static_cast<double>(rand()) / RAND_MAX);
	//	norm = NormalCDFInverse(r);	
	//	ofs << norm << endl;
//	}

return  dppLevels_Cell ; 
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


