#include "commonData.h" 
#include "SceNodes.h" 
//#include "SimulationDomainGPU.h"

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
thrust::device_vector<double> nodeDeviceTmpLocY ; 
 
        void Initialize(); 
}; 



struct InitECMLocX
{
    const double  MinLocX;
    const double  MinDist;

    InitECMLocX(double _MinLocX, double _MinDist) : MinLocX(_MinLocX), MinDist(_MinDist) {}

   __host__  __device__

        double operator()(const double & x, const double & y) const {
        return (MinLocX+x*MinDist) ; 

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


