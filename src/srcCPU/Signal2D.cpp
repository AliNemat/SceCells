#include <fstream>
#include <sstream>
#include "Signal2D.h"
#include <unistd.h>
#include <time.h>  
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

 if (importData==false) {
	const double MinFinalX=-40; 
 	const double MaxFinalX=100; 
 	const double MinFinalY=-40; 
 	const double MaxFinalY=100; 
 
 	const int nx=201; 
 	const int ny=201; 
 	double dx=(MaxFinalX-MinFinalX)/(nx-1); 
 	double dy=(MaxFinalY-MinFinalY)/(ny-1);
 	double x[nx+1],y[ny+1] ; 
 	double dppLevelsOld[nx+1][ny+1]; 
 	double dppLevels[nx+1][ny+1]; 
 	double C_Tisu[nx+1][ny+1]; 
 	double Prd[nx+1][ny+1]; 

 	double dppS_Width=5 ; 
 	double miu=1.0 ; 

	cout << "I am in update signal 0 " << std::endl ; 
 	for (int i=0 ; i<=nx ; i++) {
		x[i]=MinFinalX+0.5*dx+(i-1)*dx; 
        }
 	for (int j=0 ; j<=ny ; j++){
		y[j]=MinFinalY+0.5*dy+(j-1)*dy; 
     	}

        if (curTime==(InitTimeStage+dt)){

        	for (int i=0 ; i<=nx; i++){
        		for (int j=0 ; j<=ny ; j++){
        			dppLevelsOld[i][j]=0.0; 
        		}
       		}

       }
       else { 
       		int kk=0 ; 
     		for (int i=0 ; i<=nx; i++){
        		for (int j=0 ; j<=ny ; j++){
        			dppLevelsOld[i][j]=dppLevelsV[kk]; 
        			kk=kk+1 ; 
        		}
      		}
      }

      double Anydist; 
      for (int i=0; i<=nx; i++){
     		for (int j=0; j<=ny ; j++){
        		Anydist=pow((x[i]-Center_X),2)/pow(R_x,2)+pow((y[j]-Center_Y),2)/pow(R_y,2) ; 
       			if (Anydist<=1.0){
         			C_Tisu[i][j]=1.0 ; 
         		}
       			else {
         			C_Tisu[i][j]=0.0; 
         		}

       		}
    	}


    	for (int i=0; i<=nx; i++){
    		for (int j=0; j<=ny ; j++){
        		if (C_Tisu[i][j]==1){
         			Anydist=abs(x[i]-Center_X) ; 
         			if (Anydist<dppS_Width){
         				Prd[i][j]=100.0 ;     
         			}
         			else
           			{		  
           				Prd[i][j]=0.0 ;
           			}     
         		}
         
         		else
         		{
	 			Prd[i][j]=0.0; 
         		}

     		}
   	}

  	for ( int i=1 ; i<=nx-1 ; i++){
    		for (int j=1 ; j<ny-1 ; j++){
    
    			dppLevels[i][j]=dppLevelsOld[i][j] +dt*Prd[i][j] 
                                       +dt*miu*(
                                               (dppLevelsOld[i+1][j]-2*dppLevelsOld[i][j]+dppLevelsOld[i-1][j])/(dx*dx)
                                              +(dppLevelsOld[i][j+1]-2*dppLevelsOld[i][j]+dppLevelsOld[i][j-1])/(dy*dy)
                                               ); 
    		}
  	}

 // treating boundary
  	for (int i=1 ; i<=nx-1; i++){
    		dppLevels[i][0]=dppLevels[i][1]; 
    		dppLevels[i][ny]=dppLevels[i][ny-1]; 
   	}
  	for (int j=0 ; j <=ny ; j++){
    		dppLevels[0][j]=dppLevels[1][j]; 
    		dppLevels[nx][j]=dppLevels[nx-1][j]; 
   	}


  	double AnyX; 
  	double AnyY; 
 
  //  if (curTime==(InitTimeStage+dt)){
	dppLevels_Cell.clear(); 
  	for (int k=0; k<CellCentersHost.size(); k++)
   	{
     		AnyX=CellCentersHost[k].x; 
     		AnyY=CellCentersHost[k].y; 
     		int StoredI=-1 ;  //We will figure out that there is an error if stays -1  
     		int StoredJ=-1 ; 
  
     for (i=0 ; i<nx ; i++)
     {
       if (x[i]<=AnyX && x[i+1]>=AnyX)
        {  StoredI=i ;
           break;  
        }
     }  
     for (j=0 ; j<ny ; j++)
     {
       if (y[j]<=AnyY && y[j+1]>=AnyY)
        {  StoredJ=j ;
           break;  
        }
     }
cout<< "stored I and J is"<<StoredI<<StoredJ<<endl; 
cout<< "dpp level for stored I and J is"<<dppLevels[StoredI][StoredJ]<<endl; 
cout << "I am in first time loop for dppLevels_cell vector pushing " << std::endl ; 
     dppLevels_Cell.push_back(dppLevels[StoredI][StoredJ]);
   }

for (int k=CellCentersHost.size(); k<cellMax ; k++){
     dppLevels_Cell.push_back(0.0) ;   //these cells are not active
   }



cout << "I cellmax is " <<cellMax <<std::endl ; 
cout << "dppLevels_cell size is" << dppLevels_Cell.size()<<std::endl ; 

cout << "I am in update signal end-2 " << std::endl ; 
    if (curTime==(InitTimeStage+dt)){
        plotSignal=0 ; 
      for (int i=0 ; i<=nx; i++)
      {
        for (int j=0 ; j<=ny ; j++)
        {
        dppLevelsV.push_back(dppLevels[i][j]); 
        }
       }

     }

cout << "I am in update signal end-1 " << std::endl ; 
     int kk=0 ; 
     for (int i=0 ; i<=nx; i++)
      {
      for (int j=0 ; j<=ny ; j++)
        {
        dppLevelsV[kk]=dppLevels[i][j]; 
        kk=kk+1 ; 
        }
      }


cout << "I am in update signal end " << std::endl ;
 double z[2]; 
 z[1]=1 ; //for output purpose
 int nz=2 ; 
            plotSignal++ ; 
	if (plotSignal == 7000) {
                int curTimePlot=curTime*100 ; 
		std::string vtkFileName = "DPP_" + patch::to_string(curTimePlot) + ".vtk";
//		std::string vtkFileName = "DPP.vtk";
		ofstream SignalOut;
		SignalOut.open(vtkFileName.c_str());
		SignalOut << "# vtk DataFile Version 2.0" << endl;
		SignalOut << "Result for paraview 2d code" << endl;
		SignalOut << "ASCII" << endl;
		SignalOut << "DATASET RECTILINEAR_GRID" << endl;
		SignalOut << "DIMENSIONS" << " " << nx - 1 << " " << " " << ny - 1 << " " << nz - 1 << endl;

			SignalOut << "X_COORDINATES " << nx - 1 << " float" << endl;
			//write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
 			for (int i = 1; i <= nx - 1; i++) {
 				SignalOut << x[i] << endl;
 			}

 			SignalOut << "Y_COORDINATES " << ny - 1 << " float" << endl;
 			//write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
 			for (int j = 1; j <= ny - 1; j++) {
 				SignalOut << y[j] << endl;
 			}
 
 			SignalOut << "Z_COORDINATES " << nz - 1 << " float" << endl;
 			//write(tp + 10000, 106) 'X_COORDINATES ', Nx - 1, ' float'
 			for (int k = 1; k <= nz - 1; k++) {
 				SignalOut << z[k] << endl;
 			}
 
 			SignalOut << "POINT_DATA " << (nx - 1)*(ny - 1)*(nz - 1) << endl;
 			SignalOut << "SCALARS DPP float 1" << endl;
 			SignalOut << "LOOKUP_TABLE default" << endl;
 
 			for (int k = 1; k <= nz - 1; k++) {
 				for (int j = 1; j <= ny - 1; j++) {
 					for (int i = 1; i <= nx - 1; i++) {
 						SignalOut << dppLevels[i][j] << endl;
 
 					}
 				}
 			}
 
 			SignalOut << "SCALARS Wing float 1" << endl;
 			SignalOut << "LOOKUP_TABLE default" << endl;
 
 			for (int k = 1; k <= nz - 1; k++) {
 				for (int j = 1; j <= ny - 1; j++) {
 					for (int i = 1; i <= nx - 1; i++) {
 						SignalOut << C_Tisu[i][j] << endl;
 
 					}
 				}
 			}
 
 			plotSignal = 0; 
 		}
 
}


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
		
		sleep(20) ; 
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
          
				
		cout <<"second dpp value is"<< dppLevels_Cell.at(0)<< endl ; 	
		


     // }
    //  else {
//		for (int k=0; k<cellMax ; k++){
  //     			dppLevels_Cell.push_back(0.0) ;   //these cells are not active
    //  		}
    //  }

  }



 
return  dppLevels_Cell ; 
}
