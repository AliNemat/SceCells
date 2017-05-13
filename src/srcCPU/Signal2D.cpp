#include <fstream>
#include "Signal2D.h"


std::vector<double> updateSignal(const vector<CVector> & CellCentersHost, int cellMax, double MinX, double MaxX, double MinY, double MaxY, double dt, double InitTimeStage, double curTime)  {

cout << "I am in update signal" << std::endl ; 
vector <double> DPPLevel ;

 const double MinFinalX=-40; 
 const double MaxFinalX=100; 
 const double MinFinalY=-40; 
 const double MaxFinalY=100; 
 
 const int nx=201
 const int ny=201

 double dx=(MaxFinalX-MinFinalX)/(nx-1); 
 double dy=(MaxFinalY-MinFinalY)/(ny-1);

 for (int i=0 ; i<=nx+1 ; i++){
	x[i]=0.5*dx+(i-1)*dx
     }
 for (int j=0 ; i<=ny+1 ; j++){
	y[j]=0.5*dx+(j-1)*dy
     }
    if (curTime==InitTimeStage){
}

  for (int i=0 ; i<=nx; i++){
    for (int j=0 ; j<=ny ; j++){
     DPPLevelOld[i][j]=DPPLevel[i][j]
    }
  }

  double R_x=0.5*(MaxX-MinX); 
  double R_y=0.5*(MaxY-MinY); 
  double Center_X=MinX+0.5*(MaxX-MinX); 
  double Center_Y=MinY+0.5*(MaxY-MinY); 

  double Anydist; 
  for (int i=0; i<nx; i++){
     for (int j=0; j<ny ; j++){
       Anydist=pow((x[i]-Center_X),2)/pow(R_x,2)+pow((y[j]-Center_Y),2)/pow(R_y,2) ; 
       if (disx<=1) {
         c_Tisu[i]=1.0 ; 
       }
       else {
	 c_Tisu[i][j]=0.0; 
       }

}
}


  for (int i=0; i<=nx; i++){
     for (int j=0; j<=ny ; j++)
     {
       if (c_Tisu[i][j]==1)
       {
         AnyDist=abs(x[i]-CenterX) ; 
         if (AnyDist<5.0)
          {
          Prd[i][j]=1.0 ;     
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
    
    DppLevel[i][j]=DppLevelOld[i][j] +dt*Prd[i][j] 
                                     +dt*miu*(
                                              (DppLevel[i+1][j]-2*DppLevel[i][j]+DppLevel[i-1][j])/(dx*dx)
                                              (DppLevel[i][j+1]-2*DppLevel[i][j]+DppLevel[i][j-1])/(dx*dx)
                                             ); 
    }
 }

 // treating boundary
  for (int i=1 ; i<=nx-1; i++)
   {
    DppLevel[i][0]=DppLevel[i][1]; 
    DppLevel[i][ny]=DppLevel[i][ny-1]; 
   }
  for (int j=0 ; j <=ny ; j++)
   {
    DppLevel[0][j]=DppLevel[1][j]; 
    DppLevel[nx][j]=DppLevel[nx-1][j]; 
   }


  double AnyX
  double AnyY
  for (int k=0; i<cellCentersHost.Size(); k++)
   {
     AnyX=cellCentersHost[k].x
     AnyY=cellCentersHost[k].y
     StoredI=-1 ;  //We will figure out that there is an error if stays -1  
     StoredJ=-1 ; 
  
    for (i=0 ; i<nx-1 ; i++)
     {
       if (x[i]<=AnyX && x[i+1]>=AnyX)
        {  double StoredI=i ;
           break;  
        }
    }  
    for (j=0 ; j<ny-1 ; j++)
     {
       if (y[j]<=AnyY && y[j+1]>=AnyY)
        {  double StoredJ=i ;
           break;  
        }
    }

   DPPLevel.push_back(DppLevel[StoredI][StoredJ]; 
  }
  
    

return  DPPLevel ; 
}
