#include "commonData.h" 
#include "SceNodes.h" 
#include <string>
#include <sstream>
#include <fstream>  
//#include "SimulationDomainGPU.h"


typedef thrust ::tuple<int,double,double> IDD ; 
typedef thrust ::tuple<int,double,double,bool> IDDB ; 
typedef thrust ::tuple<double,double> DD ; 
typedef thrust ::tuple<double,double,double,double,double,double> DDDDDD ;
typedef thrust ::tuple<double,double,double,double> DDDD ; 


struct MechPara_ECM {
	double sceInterCellCPU_ECM[5] ;
	double linSpringCPU_ECM ; 
	double linSpringRestLenCPU_ECM ; 
}; 


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

MechPara_ECM mechPara_ECM ; 
 
thrust::device_vector<int> indexECM ;
thrust::device_vector<double> nodeECMLocX ; 
thrust::device_vector<double> nodeECMLocY ; 
thrust::device_vector<double> nodeECMVelX ; 
thrust::device_vector<double> nodeECMVelY ;

thrust::device_vector<double> nodeECMTmpLocX ; 
thrust::device_vector<double> nodeECMTmpLocY ; 

thrust::device_vector<double> nodeDeviceLocX ; 
thrust::device_vector<double> nodeDeviceLocY ; 
thrust::device_vector<double> nodeDeviceTmpLocX ; 
thrust::device_vector<double> nodeDeviceTmpLocY ;
 
thrust::device_vector<bool> nodeIsActive_Cell ; 

thrust::device_vector<double> linSpringForceECMX ; 
thrust::device_vector<double> linSpringForceECMY ; 
 
 
thrust::device_vector<double> bendSpringForceECMX ; 
thrust::device_vector<double> bendSpringForceECMY ;

 
thrust::device_vector<double> memMorseForceECMX ; 
thrust::device_vector<double> memMorseForceECMY ;
 
thrust::device_vector<double> totalForceECMX ; 
thrust::device_vector<double> totalForceECMY ;

        void Initialize(); 
};
 
__device__
double calMorse_ECM (const double & linkLength); 




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


struct MoveNodes_Cell: public thrust::unary_function<IDDB,DD> {

	 double  _eCMY, _activeBound ;
         int _maxMembrNodePerCell ; 

	__host__ __device__ MoveNodes_Cell (double eCMY , double activeBound, int maxMembrNodePerCell) :
				_eCMY(eCMY),_activeBound(activeBound),_maxMembrNodePerCell(maxMembrNodePerCell) {
	}
	__host__ __device__ DD operator()(const IDDB & iDDB) const {
	
	int nodeRank=     thrust::get<0>(iDDB) ; 
	double  LocX=     thrust::get<1>(iDDB) ; 
	double  LocY=     thrust::get<2>(iDDB) ; 
	bool    nodeIsActive= thrust::get<3>(iDDB) ; 

//	if (nodeRank> _maxMembrNodePerCell) {
//		return thrust::make_tuple (LocX,LocY) ; 
//	}

	 if (LocY<=_eCMY && nodeIsActive) {
			return thrust::make_tuple (LocX,_eCMY+0.1) ; 
	}	
		else if (LocY>(_eCMY+0.5)){
		
		return thrust::make_tuple (LocX,LocY) ; 
		}
        	else {
                	return thrust::make_tuple (LocX,(LocY-10.25*(LocY-_eCMY)*0.005/36.0))  ; 
		}

	}

	
} ;

struct MoveNodes2_Cell: public thrust::unary_function<IDDB,DD> {
	 double  *_locXAddr_ECM; 
         double  *_locYAddr_ECM; 
         int _maxMembrNodePerCell ; 

	__host__ __device__ MoveNodes2_Cell (double * locXAddr_ECM, double * locYAddr_ECM, int maxMembrNodePerCell) :
				_locXAddr_ECM(locXAddr_ECM),_locYAddr_ECM(locYAddr_ECM),_maxMembrNodePerCell(maxMembrNodePerCell) {
	}
	__device__ DD operator()(const IDDB & iDDB) const {
	
	int nodeRankInOneCell=     thrust::get<0>(iDDB) ; 
	double  locX=     thrust::get<1>(iDDB) ; 
	double  locY=     thrust::get<2>(iDDB) ; 
	bool    nodeIsActive= thrust::get<3>(iDDB) ; 
	
	double locX_ECM, locY_ECM ; 
	double dist ;
	double fMorse ;  
	double fTotalMorseX=0.0 ; 
	double fTotalMorseY=0.0 ;
	double fTotalMorse=0.0 ;
	int _numNodes_ECM=2500 ;  




	 if ( nodeIsActive && nodeRankInOneCell<_maxMembrNodePerCell) {

		for (int i=0 ; i<_numNodes_ECM ; i++) {
			locX_ECM=_locXAddr_ECM[i]; 
			locY_ECM=_locYAddr_ECM[i];
			dist=sqrt((locX-locX_ECM)*(locX-locX_ECM)+(locY-locY_ECM)*(locY-locY_ECM)) ;
			fMorse=calMorse_ECM(dist);  
			fTotalMorseX=fTotalMorseX+fMorse*(locX_ECM-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_ECM-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ; 
		}
	
                return thrust::make_tuple ((locX+fTotalMorseX*0.005/36.0),(locY+fTotalMorseY*0.005/36.0))  ; 
	}
		
	else {
		return thrust::make_tuple (locX,locY) ; 
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
	//return thrust::make_tuple(0.0,0.0) ; 
        

		}

	
} ;

 
 struct MorseForceECM: public thrust::unary_function<IDD,DD> {
         int  _numActiveNodes_Cell ; 	
         int  _maxNodePerCell ; 	
         int  _maxMembrNodePerCell ; 	
         double  *_locXAddr_Cell; 
         double  *_locYAddr_Cell; 
	 bool    *_nodeIsActive_Cell ;  

	__host__ __device__ MorseForceECM (int numActiveNodes_Cell, int maxNodePerCell, double maxMembrNodePerCell, double * locXAddr_Cell, double * locYAddr_Cell, bool * nodeIsActive_Cell) :
	_numActiveNodes_Cell(numActiveNodes_Cell), _maxNodePerCell(maxNodePerCell), _maxMembrNodePerCell(maxMembrNodePerCell),_locXAddr_Cell(locXAddr_Cell),_locYAddr_Cell(locYAddr_Cell),_nodeIsActive_Cell(nodeIsActive_Cell) {
	}
	 __device__ 
	DD operator()(const IDD & iDD) const {
	
	int     index=    thrust::get<0>(iDD) ; 
	double  locX=     thrust::get<1>(iDD) ; 
	double  locY=     thrust::get<2>(iDD) ; 

	double fMorse ; 
	double locX_C, locY_C ; 
	double dist ; 
	double fTotalMorseX=0.0 ; 
	double fTotalMorseY=0.0 ;
	double fTotalMorse=0.0 ;
	double kStiff=1.0 ; 
	double restLen=0.03 ; 
	// we are already in active cells. Two more conditions: 1-it is membrane 2-it is active node 
        for (int i=0 ; i<_numActiveNodes_Cell ; i++) {
		if (_nodeIsActive_Cell[i] && (i%_maxNodePerCell)<_maxMembrNodePerCell){
		//if (_nodeIsActive_Cell[i]){
				
		 
			locX_C=_locXAddr_Cell[i]; 
			locY_C=_locYAddr_Cell[i];
		//	if (locY_C>10 && locY_C<30) {
			dist=sqrt((locX-locX_C)*(locX-locX_C)+(locY-locY_C)*(locY-locY_C)) ;
		//	if (0.0001<dist<2.0) { 
		//		fMorse=kStiff*(dist-restLen) ;
		//	}
		///	else {
		//		fMorse=0.0 ; 
		//	}
			fMorse=calMorse_ECM(dist);  
			fTotalMorseX=fTotalMorseX+fMorse*(locX_C-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_C-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ; 
		}
	}
	return thrust::make_tuple(fTotalMorseX,fTotalMorseY) ; 
	//return thrust::make_tuple(5,5) ; 

	}

}; 




struct TotalECMForceCompute: public thrust::unary_function<DDDDDD,DD> {

	double _dummy ; 

	__host__ __device__ TotalECMForceCompute(double dummy):_dummy(dummy) {
	}

	__host__ __device__ DD operator() (const DDDDDD & dDDDDD) const {

	double fLinSpringX=  thrust:: get<0>(dDDDDD); 
	double fLinSpringY=  thrust:: get<1>(dDDDDD); 
	double fBendSpringX= thrust:: get<2>(dDDDDD); 
	double fBendSpringY= thrust:: get<3>(dDDDDD); 
	double fMembX       = thrust:: get<4>(dDDDDD); 
	double fMembY       = thrust:: get<5>(dDDDDD); 


	return thrust::make_tuple(fLinSpringX+fBendSpringX+fMembX,fLinSpringY+fBendSpringY+fMembY); 
	}
}; 



struct MoveNodeECM: public thrust::unary_function<DDDD,DD> {

	double _dt ; 
	double _damp ; 
	__host__ __device__ MoveNodeECM (double dt, double damp): _dt(dt), _damp(damp) {
	}

	__host__ __device__ DD operator() (const DDDD & dDDD) const  {


	double locXOld= thrust:: get <0> (dDDD) ;
	double locYOld= thrust:: get <1> (dDDD) ;
	double fx= 	thrust:: get <2> (dDDD) ;
	double fy= 	thrust:: get <3> (dDDD) ;
 
	return thrust::make_tuple (locXOld+fx*_dt/_damp, locYOld+fy*_dt/_damp) ;
 
	}
}; 
