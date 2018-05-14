#include "commonData.h"
#include "SceNodes.h" 
#include <string>
#include <sstream>
#include <fstream>  
//#include "SimulationDomainGPU.h"


typedef thrust ::tuple<int,double,double> IDD ; 
typedef thrust ::tuple<int,double,double,bool,MembraneType1> IDDBT ; 
typedef thrust ::tuple<double,double> DD ; 
typedef thrust ::tuple<double,double,double,double> DDDD ; 
typedef thrust ::tuple<double,double,bool> DDB ; 
typedef thrust ::tuple<double,double,MembraneType1,int,double,double> DDTIDD ; 
typedef thrust ::tuple<double,double,double,double,double,double> DDDDDD ;
typedef thrust ::tuple<double,double,double,double,int,EType> DDDDIT ; 


struct MechPara_ECM {
	double sceInterCellCPU_ECM[5] ;
	double wLCParaCPU_ECM[4] ;
	double linSpringCPU_ECM ; 
	double linSpringRestLenCPU_ECM ; 
}; 


class SceECM {
//	SceNodes* nodes;

public:
        void ApplyECMConstrain(int totalNodeCountForActiveCellsECM, double curTime, double dt, double Damp_Coef, bool cellPolar, bool subCellPolar, bool isInitPhase) ; 
		void Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM, int freqPlotData); 
		EType ConvertStringToEType (string eNodeRead) ;
		void PrintECM(double curTime); 
double restLenECMSpring ;
double eCMLinSpringStiff ; 
double restLenECMAdhSpring ; 
double maxLenECMAdhSpring ; 
double kAdhECM ;
double totalLinSpringEnergy,totalMorseEnergy, totalAdhEnergy, totalMorseEnergyCell, totalAdhEnergyCell ; 

int outputFrameECM ;  
int lastPrintECM ;
int numNodesECM ;
int freqPlotData ; 
uint maxAllNodePerCell ; 
uint maxMembrNodePerCell ;
uint maxTotalNodes ; 
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
thrust::device_vector<MembraneType1> memNodeType ;
thrust::device_vector<int>   adhPairECM_Cell ;
 
thrust::device_vector<bool> nodeIsActive_Cell ; 

thrust::device_vector<double> linSpringForceECMX ; 
thrust::device_vector<double> linSpringForceECMY ; 
thrust::device_vector<double> linSpringAvgTension  ; 
thrust::device_vector<double> linSpringEnergy  ; 
thrust::device_vector<double> morseEnergy  ; 
thrust::device_vector<double> adhEnergy  ;

thrust::device_vector<double> morseEnergyCell ; //it should be cell size 
thrust::device_vector<double> adhEnergyCell  ; // it should be cell size
 
 
thrust::device_vector<double> bendSpringForceECMX ; 
thrust::device_vector<double> bendSpringForceECMY ;

 
thrust::device_vector<double> memMorseForceECMX ; 
thrust::device_vector<double> memMorseForceECMY ;
 
thrust::device_vector<double> fBendCenterX ;
thrust::device_vector<double> fBendCenterY ;
thrust::device_vector<double> fBendLeftX ;
thrust::device_vector<double> fBendLeftY ;
thrust::device_vector<double> fBendRightX ;
thrust::device_vector<double> fBendRightY ;

thrust::device_vector<double> totalForceECMX ; 
thrust::device_vector<double> totalForceECMY ;
thrust::device_vector<EType>  peripORexcm ;

thrust::device_vector<double> stiffLevel ;
thrust::device_vector<double> sponLen ;
};
 
__device__
double calMorse_ECM (const double & linkLength); 

__device__
double calMorseEnergy_ECM (const double & linkLength); 
__device__
double calWLC_ECM (const double & linkLength); 


__device__
bool IsValidAdhPair (const double & dist); 

__device__
bool IsValidAdhPairForNotInitPhase (const double & dist); 

__device__
double CalAdhECM (const double & dist);

__device__
double CalAdhEnergy (const double & dist);

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


struct MoveNodes2_Cell: public thrust::unary_function<IDDBT,DDTIDD> {
	 double  *_locXAddr_ECM; 
         double  *_locYAddr_ECM; 
        uint _maxMembrNodePerCell ; 
	 int _numNodes_ECM ;
	 double _dt ; 
	 double _Damp_Coef ;
	 bool _isInitPhase ;
	 EType*  _peripORexcmAddr ;
	 double _curTime ; 
	__host__ __device__ MoveNodes2_Cell (double * locXAddr_ECM, double * locYAddr_ECM, uint maxMembrNodePerCell, int numNodes_ECM, double dt, double Damp_Coef, bool isInitPhase, EType * peripORexcmAddr, double curTime) :
				_locXAddr_ECM(locXAddr_ECM),_locYAddr_ECM(locYAddr_ECM),_maxMembrNodePerCell(maxMembrNodePerCell),_numNodes_ECM(numNodes_ECM),_dt(dt),
			    _Damp_Coef(Damp_Coef), _isInitPhase (isInitPhase), _peripORexcmAddr(peripORexcmAddr),_curTime (curTime)	{
	}
	__device__ DDTIDD  operator()(const IDDBT & iDDBT) const {
	
	int nodeRankInOneCell=          thrust::get<0>(iDDBT) ; 
	double            locX=         thrust::get<1>(iDDBT) ; 
	double            locY=         thrust::get<2>(iDDBT) ; 
	bool              nodeIsActive= thrust::get<3>(iDDBT) ; 
	MembraneType1     nodeType=     thrust::get<4>(iDDBT) ; 
	
	double locX_ECM, locY_ECM ; 
	double dist ;
	double fMorse ;  
	double fTotalMorseX=0.0 ; 
	double fTotalMorseY=0.0 ;
	double fTotalMorse=0.0 ;
	double eMorseCell=0 ; 
	double eAdhCell= 0 ; 
	double distMin=10000 ; // large number
	double distMinX, distMinY ; 
	//double kStifMemECM=3.0 ; // need to take out 
	//double distMaxAdh=0.78125; // need to take out
	//double distSponAdh=0.0625 ;  // need to take out
	double fAdhMemECM=0.0   ; 
	double fAdhMemECMX=0.0 ; 
	double fAdhMemECMY=0.0  ; 
	int    adhPairECM=-1 ; //no adhere Pair
	int   iPair=-1 ;
	double smallNumber=0.000001;

		if ( nodeIsActive && nodeRankInOneCell<_maxMembrNodePerCell ) {
			for (int i=0 ; i<_numNodes_ECM ; i++) {
				locX_ECM=_locXAddr_ECM[i]; 
				locY_ECM=_locYAddr_ECM[i];
				dist=sqrt((locX-locX_ECM)*(locX-locX_ECM)+(locY-locY_ECM)*(locY-locY_ECM)) ;
				fMorse=calMorse_ECM(dist);
				eMorseCell=eMorseCell + calMorseEnergy_ECM(dist);  
				fTotalMorseX=fTotalMorseX+fMorse*(locX_ECM-locX)/dist ; 
				fTotalMorseY=fTotalMorseY+fMorse*(locY_ECM-locY)/dist ; 
				fTotalMorse=fTotalMorse+fMorse ; 
				
				if ( dist < distMin && (nodeType==basal1 || nodeType==apical1) ) {  //adhesion only for basal and apical nodes
					distMin=dist ; 
					distMinX=(locX_ECM-locX) ;
					distMinY=(locY_ECM-locY) ; 
					iPair=i ; 
				}

			}

			if (IsValidAdhPairForNotInitPhase(distMin)&& iPair!=-1) {

        		fAdhMemECM=CalAdhECM(distMin) ; 
				eAdhCell=CalAdhEnergy(distMin) ; 

				fAdhMemECMX=fAdhMemECM*distMinX/distMin ;  
				fAdhMemECMY=fAdhMemECM*distMinY/distMin ; 
				adhPairECM=iPair ; 
			}
	 }

	 return thrust::make_tuple ((locX+(fTotalMorseX+fAdhMemECMX)*_dt/_Damp_Coef),
	 							(locY+(fTotalMorseY+fAdhMemECMY)*_dt/_Damp_Coef),
								 nodeType,adhPairECM,eMorseCell,eAdhCell )  ; 
		
}
	

	
} ;

        



struct LinSpringForceECM: public thrust::unary_function<IDD,DDDD> {
         int   _numNodes ; 	
	 double  _restLen ; 
	 double  _linSpringStiff ;
         double  *_locXAddr; 
         double  *_locYAddr;
		 double *_stiffAddr; 
		 double * _sponLenAddr ; 


	__host__ __device__ LinSpringForceECM (double numNodes, double restLen, double linSpringStiff , double * locXAddr, double * locYAddr, 
										   double * stiffAddr, double * sponLenAddr ) :
	 _numNodes(numNodes),_restLen(restLen),_linSpringStiff(linSpringStiff),_locXAddr(locXAddr),_locYAddr(locYAddr),_stiffAddr (stiffAddr), _sponLenAddr (sponLenAddr) {
	}
	 __device__ DDDD operator()(const IDD & iDD) const {
	
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

	double energyLeft, energyRight ; 
	int indexLeft, indexRight ; 

        if (index != 0) {

		locXLeft=_locXAddr[index-1] ; 
		locYLeft=_locYAddr[index-1] ;
		indexLeft=index-1 ; 
		}
	else {
		locXLeft=_locXAddr[_numNodes-1] ;
		locYLeft=_locYAddr[_numNodes-1] ;
		indexLeft=_numNodes-1 ; 
	}

	distLeft=sqrt( ( locX-locXLeft )*(locX-locXLeft) +( locY-locYLeft )*(locY-locYLeft) ) ;
	double kStiffLeft =0.5*(  _stiffAddr[indexLeft]  +  _stiffAddr[index]  ) ; 
	double sponLenLeft=0.5*(_sponLenAddr[indexLeft]  +_sponLenAddr[index]  ) ; 
	//	forceLeft=calWLC_ECM(distLeft) ; 
		//forceLeft=_linSpringStiff*(distLeft-_restLen) ; 
		forceLeft=kStiffLeft*(distLeft-sponLenLeft) ; 
		energyLeft=0.5*kStiffLeft*(distLeft-sponLenLeft)*(distLeft-sponLenLeft) ; 
		forceLeftX=forceLeft*(locXLeft-locX)/distLeft ; 
		forceLeftY=forceLeft*(locYLeft-locY)/distLeft ; 

	if (index != _numNodes-1) {

		
		locXRight=_locXAddr[index+1] ; 
		locYRight=_locYAddr[index+1] ;
		indexRight=index+1 ; 
		}
	else {

		locXRight=_locXAddr[0] ; 
		locYRight=_locYAddr[0] ;
		indexRight=0  ; 
	}
	distRight=sqrt( ( locX-locXRight )*(locX-locXRight) +( locY-locYRight )*(locY-locYRight) ) ; 
	double kStiffRight =0.5*(_stiffAddr  [indexRight]  +  _stiffAddr[index]  ) ; 
	double sponLenRight=0.5*(_sponLenAddr[indexRight]  +_sponLenAddr[index]  ) ; 
   	    	//forceRight=_linSpringStiff*(distRight-_restLen) ; 
   	    	forceRight=kStiffRight*(distRight-sponLenRight) ;
			energyRight=0.5*kStiffRight*(distRight-sponLenRight)*(distRight-sponLenRight)  ;

        //  	forceRight=calWLC_ECM(distRight) ; 
		forceRightX=forceRight*(locXRight-locX)/distRight ; 
		forceRightY=forceRight*(locYRight-locY)/distRight ; 
		//for open ECM.
//	if (index == 0 || index==int(_numNodes/2) ) {
//		return thrust::make_tuple(forceRightX,forceRightY,forceRight) ;
  // }
 //  else if (index ==_numNodes-1 || index==(int(_numNodes/2)-1) ) {
//		return thrust::make_tuple(forceLeftX,forceLeftY,forceLeft) ;
//	}
//	else {
		return thrust::make_tuple(forceLeftX+forceRightX,forceLeftY+forceRightY,0.5*(forceLeft+forceRight), energyLeft+energyRight) ;
//	}
        

}

	
} ;

 
 struct MorseAndAdhForceECM: public thrust::unary_function<IDD,DDDD> {
         int  _numActiveNodes_Cell ; 	
         uint  _maxNodePerCell ; 	
         uint  _maxMembrNodePerCell ; 	
         double  *_locXAddr_Cell; 
         double  *_locYAddr_Cell; 
	 bool    *_nodeIsActive_Cell ;  
	 int     *_adhPairECM_Cell ; 
	__host__ __device__ MorseAndAdhForceECM (int numActiveNodes_Cell, uint maxNodePerCell, uint maxMembrNodePerCell, double * locXAddr_Cell, double * locYAddr_Cell, bool * nodeIsActive_Cell, int * adhPairECM_Cell) :
	_numActiveNodes_Cell(numActiveNodes_Cell), _maxNodePerCell(maxNodePerCell), _maxMembrNodePerCell(maxMembrNodePerCell),_locXAddr_Cell(locXAddr_Cell),_locYAddr_Cell(locYAddr_Cell),_nodeIsActive_Cell(nodeIsActive_Cell),_adhPairECM_Cell(adhPairECM_Cell) {
	}
	 __device__ 
	DDDD operator()(const IDD & iDD) const {
	
	int     index=    thrust::get<0>(iDD) ; 
	double  locX=     thrust::get<1>(iDD) ; 
	double  locY=     thrust::get<2>(iDD) ; 

	double fMorse ; 
	double locX_C, locY_C ; 
	double dist ;
	double distMin=100000 ; //large number
	double distMinX,distMinY ; 
	double fAdhX=0 ; 
	double fAdhY=0 ; 
	double eAdh=0 ; 
	double eMorse=0 ; 
	double fTotalMorseX=0.0 ; 
	double fTotalMorseY=0.0 ;
	double fTotalMorse=0.0 ;
	//double kAdhECM=3 ;  //need to take out 
	//double distAdhSpon=0.0625 ; // need to take out 
	//double distAdhMax=0.78125 ; // need to take out
	double fAdh ; 
	// we are already in active cells. Two more conditions: 1-it is membrane 2-it is active node 
        for (int i=0 ; i<_numActiveNodes_Cell ; i++) {
		if (_nodeIsActive_Cell[i] && (i%_maxNodePerCell)<_maxMembrNodePerCell){
		//if (_nodeIsActive_Cell[i]){
		
		 
			locX_C=_locXAddr_Cell[i]; 
			locY_C=_locYAddr_Cell[i];
			dist=sqrt((locX-locX_C)*(locX-locX_C)+(locY-locY_C)*(locY-locY_C)) ;
			fMorse=calMorse_ECM(dist);  
			eMorse=eMorse + calMorseEnergy_ECM(dist);  
			fTotalMorseX=fTotalMorseX+fMorse*(locX_C-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_C-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ;
			if (_adhPairECM_Cell[i]==index) {
				fAdh=CalAdhECM(dist) ; 
				eAdh=eAdh+CalAdhEnergy(dist) ; 
				fAdhX=fAdhX+fAdh*(locX_C-locX)/dist ; 
				fAdhY=fAdhY+fAdh*(locY_C-locY)/dist ; 

			}			 
		}
	}
	
	return thrust::make_tuple(fTotalMorseX+fAdhX,fTotalMorseY+fAdhY,eMorse,eAdh) ;

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
//	return thrust::make_tuple(fLinSpringX+fMembX,fLinSpringY+fMembY); 
	}
}; 

struct MechProp: public thrust::unary_function<EType,DD> {

	bool _isInitPhase ;
	double _stiffness ; 
	double _sponLen ; 

	__host__ __device__ MechProp(bool isInitPhase, double stiffness, double sponLen): _isInitPhase(isInitPhase), _stiffness(stiffness),_sponLen (sponLen) {
	}

	__host__ __device__ DD operator() (const EType  & nodeType) const {

	double stiffness= _stiffness  ;
	double sponLen=0 ;   ; 
	if (_isInitPhase == false ) {
		if (nodeType==excm) {
			stiffness=1.0* _stiffness ; 
			sponLen=0.08  ; 
		}
		if (nodeType==perip) {
			stiffness=_stiffness ; 
			sponLen=0.08 ; // 0.1 ; 
		}

		if (nodeType==bc2) {
			stiffness=0.1*_stiffness ; 
			sponLen=0.08 ;// _sponLen ; 
		}
	}

	return thrust::make_tuple(stiffness,sponLen); 
	}
}; 



struct MoveNodeECM: public thrust::unary_function<DDDDIT,DD> {

	double _dt ; 
	double _dampECM ; 
	double _dampPeri ; 
	int _numNodes ;
	double _curTime ; 
	__host__ __device__ MoveNodeECM (double dt, double dampECM, double dampPeri, int numNodes, double curTime): _dt(dt), _dampECM(dampECM),_dampPeri(dampPeri), _numNodes(numNodes),_curTime(curTime) {
	}
	__host__ __device__ DD operator() (const DDDDIT & dDDDIT) const  {
	double 			  locXOld= thrust:: get <0> (dDDDIT) ;
	double 			  locYOld= thrust:: get <1> (dDDDIT) ;
	double 			  fx= 	   thrust:: get <2> (dDDDIT) ;
	double 			  fy= 	   thrust:: get <3> (dDDDIT) ;
	int    			  index=   thrust:: get <4> (dDDDIT) ; 
	EType             nodeType=thrust:: get <5> (dDDDIT) ; 
	//if (index == 0 || index==_numNodes-1 || index==( int (_numNodes/2)-1) || index == int (_numNodes/2) ) {
	//if (( index == 0  || index==( int (_numNodes/2)-1) ) && (_curTime<=200 ) ) {
		
//	if (( index == 0  || index==560 ) && (_curTime<=200 ) ) {
//		return thrust::make_tuple (locXOld+fx*_dt/_dampECM, locYOld) ;
//	}
//	else {
		if (nodeType==excm) {		
			return thrust::make_tuple (locXOld+fx*_dt/_dampECM, locYOld+fy*_dt/_dampECM) ;
		}
		else {
//			if (_curTime<=200) {
//				return thrust::make_tuple (locXOld, locYOld+0.36*_dt/_dampPeri) ;
//			} 
//			else {
				return thrust::make_tuple (locXOld+fx*_dt/_dampPeri, locYOld+fy*_dt/_dampPeri) ;
		//	}
		}
//	}
}
}; 


struct CalBendECM: public thrust::unary_function<IDD, DDDDDD> {
	double* _locXAddr;
	double* _locYAddr;
	int  _numECMNodes ;
	double _eCMBendStiff ;  

	__host__ __device__ CalBendECM (double* locXAddr, double* locYAddr, int numECMNodes, double eCMBendStiff) :
				_locXAddr(locXAddr), _locYAddr(locYAddr),_numECMNodes(numECMNodes),_eCMBendStiff(eCMBendStiff){
	}
	
	__host__ __device__ DDDDDD operator()(const IDD &  iDD) const {

		int   nodeRank = thrust::get<0>(iDD);
		double locX = thrust::get<1>(iDD);
		double locY = thrust::get<2>(iDD);

		double kb=_eCMBendStiff ; 
		double leftPosX,leftPosY;
		double lenLeft;

		double rightPosX,rightPosY;
		double lenRight;
			
			int index_left = nodeRank - 1;
			if (index_left == -1) {
				index_left = _numECMNodes - 1;
			}
		
				leftPosX = _locXAddr[index_left];
				leftPosY = _locYAddr[index_left];
				lenLeft = sqrt((leftPosX - locX) * (leftPosX - locX) + (leftPosY - locY) * (leftPosY - locY) );


			int index_right = nodeRank + 1;
			if (index_right ==  _numECMNodes) {
				index_right = 0;
			}
				rightPosX = _locXAddr[index_right];
				rightPosY = _locYAddr[index_right];
				lenRight = sqrt((rightPosX - locX) * (rightPosX - locX) + (rightPosY - locY) * (rightPosY - locY) );

			double cosTheta=( (leftPosX-locX)*(rightPosX-locX)+(leftPosY-locY)*(rightPosY-locY) )/(lenRight*lenLeft) ; 

			
			double 	bendLeftX= -kb*(rightPosX-locX)/(lenRight*lenLeft)+
						            kb*cosTheta/(lenLeft*lenLeft)*(leftPosX-locX) ; 

			double 	bendRightX=-kb*(leftPosX-locX)/(lenRight*lenLeft)+
						            kb*cosTheta/(lenRight*lenRight)*(rightPosX-locX) ; 

			double  bendLeftY =-kb*(rightPosY-locY)/(lenRight*lenLeft)+
						            kb* cosTheta/(lenLeft*lenLeft)*(leftPosY-locY) ; 

		    double 	bendRightY=-kb*(leftPosY-locY)/(lenRight*lenLeft)+
						            kb* cosTheta/(lenRight*lenRight)*(rightPosY-locY) ;

			double bendCenterX=-(bendLeftX+bendRightX) ; 
			double bendCenterY=-(bendLeftY+bendRightY) ; 
			return thrust::make_tuple(bendCenterX,bendCenterY,
					bendLeftX, bendLeftY, bendRightX, bendRightY);

	}
}; 




struct SumBendForce: public thrust::unary_function<IDD,DD> {

	double* _fBendLeftXAddr;
	double* _fBendLeftYAddr;
	double* _fBendRightXAddr;
	double* _fBendRightYAddr;
	int  	_numECMNodes ;

	__host__ __device__ SumBendForce (double* fBendLeftXAddr, double* fBendLeftYAddr, double *fBendRightXAddr, double *fBendRightYAddr, int numECMNodes) :
				_fBendLeftXAddr(fBendLeftXAddr), _fBendLeftYAddr(fBendLeftYAddr),_fBendRightXAddr(fBendRightXAddr),_fBendRightYAddr(fBendRightYAddr),_numECMNodes(numECMNodes){
	}
	
	__host__ __device__ DD operator() (const IDD & iDD) const {

		int   nodeRank = thrust::get<0>(iDD);
		double fBendCenterX = thrust::get<1>(iDD);
		double fBendCenterY = thrust::get<2>(iDD);

		int index_left = nodeRank - 1;
		if (index_left == -1) {
			index_left = _numECMNodes - 1;
		}
		
		int index_right = nodeRank + 1;
		if (index_right ==  _numECMNodes) {
			index_right = 0;
		}

	return thrust::make_tuple(fBendCenterX+_fBendLeftXAddr[index_right]+_fBendRightXAddr[index_left],
				  fBendCenterY+_fBendLeftYAddr[index_right]+_fBendRightYAddr[index_left]); 
	}
}; 



