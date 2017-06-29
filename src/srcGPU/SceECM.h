#include "commonData.h" 
#include "SceNodes.h" 
//#include "SimulationDomainGPU.h"


typedef thrust ::tuple<int,double,double> IDD ; 
typedef thrust ::tuple<double,double> DD ; 

class SceECM {
//	SceNodes* nodes;

public:
        void applyECMForce_M() ; 
        void ApplyECMConstrain(int totalNodeCountForActiveCellsECM) ; 
 double eCMY ; 
thrust::device_vector<double> nodeECMLocX ; 
thrust::device_vector<double> nodeECMLocY ; 
thrust::device_vector<double> nodeECMVelX ; 
thrust::device_vector<double> nodeECMVelY ;


thrust::device_vector<double> nodeDeviceLocX ; 
thrust::device_vector<double> nodeDeviceLocY ; 
thrust::device_vector<double> nodeDeviceTmpLocX ; 
thrust::device_vector<double> nodeDeviceTmpLocY ; 
 
        void Initialize(); 
}; 



struct InitECMLoc
{
    const double  _MinLocX;
    const double  _MinDist;

    InitECMLoc (double MinLocX, double MinDist) : _MinLocX(MinLocX), _MinDist(MinDist) {}

   __host__  __device__

        double operator()(const double & x, const double & y) const {
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
		else if (LocY>(_eCMY+0.4)){
		
		return thrust::make_tuple (LocX,LocY) ; 
		}
        	else {
                	return thrust::make_tuple (LocX,(LocY-2.5*(LocY-_eCMY)*0.003/36.0))  ; 
		}

	}

	
} ; 





