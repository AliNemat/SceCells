#include "commonData.h" 
#include "SceNodes.h" 
#include <string>
#include <sstream> 
//#include "SimulationDomainGPU.h"


typedef thrust ::tuple<int,double,double> IDD ; 
typedef thrust ::tuple<double,double> DD ; 

class SceECM {
//	SceNodes* nodes;

public:
        void applyECMForce_M() ; 
        void ApplyECMConstrain(int totalNodeCountForActiveCellsECM) ; 
 double eCMY ; 
double restLenECMSpring ; 
int outputFrameECM ;  
int lastPrintECM ;
int numNodesECM ; 
thrust::device_vector<int> indexECM ;
thrust::device_vector<double> nodeECMLocX ; 
thrust::device_vector<double> nodeECMLocY ; 
thrust::device_vector<double> nodeECMVelX ; 
thrust::device_vector<double> nodeECMVelY ;


thrust::device_vector<double> nodeDeviceLocX ; 
thrust::device_vector<double> nodeDeviceLocY ; 
thrust::device_vector<double> nodeDeviceTmpLocX ; 
thrust::device_vector<double> nodeDeviceTmpLocY ; 

thrust::device_vector<double> linSpringForceECMX ; 
thrust::device_vector<double> linSpringForceECMY ; 
 
 
        void Initialize(); 
}; 



struct InitECMLoc
{
     double  _MinLocX;
     double  _MinDist;

    InitECMLoc (double MinLocX, double MinDist) : _MinLocX(MinLocX), _MinDist(MinDist) {}

   __host__  __device__

        double operator()(const int & x, const double & y)  {
        return (_MinLocX+x*_MinDist) ; 

  }
};


struct AddECMForce
{
    const double  eCMY;


    AddECMForce(double _eCMY) : eCMY(_eCMY) {}

   __host__  __device__

        double operator()( const double & x, const double & y) const {
        if (y<eCMY) {
        return (eCMY) ; 
        }
        else {
        return (y) ; 
         }
  }
};

struct MyFunction: public thrust::unary_function <double,double>{
       double _k, _bound ; 
      

       __host__ __device__ MyFunction(double k , double bound ) :
                             _k(k),_bound(bound) {
	}
       __host__ __device__
       double operator()(const double &y) {
	if (y<=_k && y>_bound) {
		return _k ; 
	}
        else {
                return y; 
	}
	}
} ; 

struct ModuloFunctor2: public thrust::unary_function <int,int>{
       int _dividend ;  
      

       __host__ __device__ ModuloFunctor2(int dividend) :
                             _dividend(dividend) {
	}
       __host__ __device__
       int operator()(const int &num) {
                return num %_dividend;
	} 
	} ; 


struct MyFunctor2: public thrust::unary_function<IDD,DD> {

	 double  _eCMY, _activeBound ;
         int _maxMembrNodePerCell ; 

	__host__ __device__ MyFunctor2 (double eCMY , double activeBound, int maxMembrNodePerCell) :
				_eCMY(eCMY),_activeBound(activeBound),_maxMembrNodePerCell(maxMembrNodePerCell) {
	}
	__host__ __device__ DD operator()(const IDD & iDD) const {
	
	int nodeRank= thrust::get<0>(iDD) ; 
	double  LocX=     thrust::get<1>(iDD) ; 
	double  LocY=     thrust::get<2>(iDD) ; 

	if (nodeRank> _maxMembrNodePerCell) {
		return thrust::make_tuple (LocX,LocY) ; 
	}

	else if (LocY<=_eCMY && LocY>_activeBound) {
			return thrust::make_tuple (LocX,_eCMY) ; 
	}	
		else if (LocY>(_eCMY+0.5)){
		
		return thrust::make_tuple (LocX,LocY) ; 
		}
        	else {
                	return thrust::make_tuple (LocX,(LocY-10.25*(LocY-_eCMY)*0.005/36.0))  ; 
		}

	}

	
} ;

struct LinSpringForceECM: public thrust::unary_function<IDD,DD> {
         double  _numNodes ; 	
	 double  _restLen ; 
	 double  _linSpringStiff ;
         double  *_locXAddr; 
         double  *_locYAddr; 
	  

	__host__ __device__ LinSpringForceECM (double numNodes, double restLen, double linSpringStiff , double * locXAddr, double * locYAddr) :
	 _numNodes(numNodes),_restLen(restLen),_linSpringStiff(linSpringStiff),_locXAddr(locXAddr),_locYAddr(locYAddr) {
	}
	__host__ __device__ DD operator()(const IDD & iDD) const {
	
	int     index=    thrust::get<0>(iDD) ; 
	double  locX=     thrust::get<1>(iDD) ; 
	double  locY=     thrust::get<2>(iDD) ; 

	double locXLeft  ; 
	double locYLeft  ;

	double distLeft ;
	double forceLeft  ; 
	double forceLeftX ; 
	double forceLeftY ; 

	double locXRight ; 
	double locYRight ;
 
	double distRight ; 
	double forceRight ; 
	double forceRightX ; 
	double forceRightY ; 


        if (index != 0) {

		locXLeft=_locXAddr[index-1] ; 
		locYLeft=_locYAddr[index-1] ;
		distLeft=sqrt( ( locX-locXLeft )*(locX-locXLeft) +( locY-locYLeft )*(locY-locYLeft) ) ;
		forceLeft=_linSpringStiff*(distLeft-_restLen) ; 
		forceLeftX=forceLeft*(locXLeft-locX)/distLeft ; 
		forceLeftY=forceLeft*(locYLeft-locY)/distLeft ; 
	}
	else {
		forceLeftX=0.0 ; 
		forceLeftY=0.0 ; 
	}

	if (index != _numNodes-1) {

		
		locXRight=_locXAddr[index+1] ; 
		locYRight=_locYAddr[index+1] ;
		distRight=sqrt( ( locX-locXRight )*(locX-locXRight) +( locY-locYRight )*(locY-locYRight) ) ; 
		forceRight=_linSpringStiff*(distRight-_restLen) ; 
		forceRightX=forceRight*(locXRight-locX) /distRight ; 
		forceRightY=forceRight*(locYRight-locY) /distRight ; 
	}
	else {

		 forceRightX=0.0 ; 
		 forceRightY=0.0 ; 
		
	}

	return thrust::make_tuple(forceLeftX+forceRightX,forceLeftY+forceRightY) ; 
        

		}

	
} ; 
 





