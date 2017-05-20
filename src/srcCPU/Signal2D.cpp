#include <fstream>
#include "Signal2D.h"


vector<double> updateSignal(vector<double> & dppLevelsV,const vector<CVector> & CellCentersHost, int cellMax, double MinX, double MaxX, double MinY, double MaxY, double dt, double InitTimeStage, double curTime)  {

cout << "I am in update signal started" << std::endl ; 

 vector<double> dppLevels_Cell ; 
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

cout << "I am in update signal started 2" << std::endl ; 
 double R_x=0.5*(MaxX-MinX); 
 double R_y=0.5*(MaxY-MinY); 
 double Center_X=MinX+0.5*(MaxX-MinX); 
 double Center_Y=MinY+0.5*(MaxY-MinY); 
 double dppS_Width=5 ; 
 double miu=1.0 ; 

cout << "I am in update signal 0 " << std::endl ; 
 for (int i=0 ; i<=nx ; i++) 
     {
	x[i]=0.5*dx+(i-1)*dx; 
     }
 for (int j=0 ; j<=ny ; j++)
     {
	y[j]=0.5*dx+(j-1)*dy; 
     }

cout << "I am in update signal 0.5  " << std::endl ; 
cout << "current time is " <<curTime <<std::endl ; 
cout << "Init time stage is"<<InitTimeStage << std::endl ; 
    if (curTime==(InitTimeStage+dt)){

      for (int i=0 ; i<=nx; i++)
      {
        for (int j=0 ; j<=ny ; j++)
        {
        dppLevelsOld[i][j]=0.0; 
        }
       }

cout << "I am in update signal 1 " << std::endl ; 
     }
    else { 
     int kk=0 ; 
     for (int i=0 ; i<=nx; i++)
      {
      for (int j=0 ; j<=ny ; j++)
        {
        dppLevelsOld[i][j]=dppLevelsV[kk]; 
        kk=kk+1 ; 
        }
      }
     }

  double Anydist; 
  for (int i=0; i<=nx; i++)
    {
    for (int j=0; j<=ny ; j++)
       {
       Anydist=pow((x[i]-Center_X),2)/pow(R_x,2)+pow((y[j]-Center_Y),2)/pow(R_y,2) ; 
       if (Anydist<=1.0)
         {
         C_Tisu[i][j]=1.0 ; 
         }
       else 
         {
         C_Tisu[i][j]=0.0; 
         }

       }
    }


  for (int i=0; i<=nx; i++)
     {
     for (int j=0; j<=ny ; j++)
       {
       if (C_Tisu[i][j]==1)
         {
         Anydist=abs(x[i]-Center_X) ; 
         if (Anydist<dppS_Width)
           {
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

  for ( int i=1 ; i<=nx-1 ; i++)
  {
    for (int j=1 ; j<ny-1 ; j++)
    {
    
    dppLevels[i][j]=dppLevelsOld[i][j] +dt*Prd[i][j] 
                                       +dt*miu*(
                                               (dppLevelsOld[i+1][j]-2*dppLevelsOld[i][j]+dppLevelsOld[i-1][j])/(dx*dx)
                                              +(dppLevelsOld[i][j+1]-2*dppLevelsOld[i][j]+dppLevelsOld[i][j-1])/(dy*dy)
                                               ); 
    }
  }

cout << "I am in update signal 2 " << std::endl ; 
 // treating boundary
  for (int i=1 ; i<=nx-1; i++)
   {
    dppLevels[i][0]=dppLevels[i][1]; 
    dppLevels[i][ny]=dppLevels[i][ny-1]; 
   }
  for (int j=0 ; j <=ny ; j++)
   {
    dppLevels[0][j]=dppLevels[1][j]; 
    dppLevels[nx][j]=dppLevels[nx-1][j]; 
   }

cout << "I am in update signal 3 " << std::endl ; 

  double AnyX; 
  double AnyY; 
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

cout << "I am in loop for dppLevels_cell vector pushing " << std::endl ; 
     dppLevels_Cell.push_back(dppLevels[StoredI][StoredJ]);
   }


   for (int k=CellCentersHost.size(); k<cellMax ; k++){
     dppLevels_Cell.push_back(0.0);  //these cells are not active
   }
cout << "I am in update signal end-2 " << std::endl ; 
    if (curTime==(InitTimeStage+dt)){

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
return  dppLevels_Cell ; 
}
