#include "commonData.h"
#include "SceNodes.h" 
#include <string>
#include <sstream>
#include <fstream>  
//#include "SimulationDomainGPU.h"


typedef thrust ::tuple<int,double,double> IDD ; 
typedef thrust ::tuple<int,double,double,bool,MembraneType1> IDDBT ; 
typedef thrust ::tuple<double,double> DD ; 
typedef thrust ::tuple<double,double,double> DDD ; 
typedef thrust ::tuple<double,double,bool> DDB ; 
typedef thrust ::tuple<double,double,MembraneType1,int> DDTI ; 
typedef thrust ::tuple<double,double,double,double,double,double> DDDDDD ;
typedef thrust ::tuple<double,double,double,double,int> DDDDI ; 


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
double restLenECMSpring ;
double eCMLinSpringStiff ; 
double restLenECMAdhSpring ; 
double maxLenECMAdhSpring ; 
double kAdhECM ; 
int outputFrameECM ;  
int lastPrintECM ;
int numNodesECM ;
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

        void Initialize(uint maxAllNodePerCellECM, uint maxMembrNodePerCellECM, uint maxTotalNodesECM); 
};
 
__device__
double calMorse_ECM (const double & linkLength); 

__device__
double calWLC_ECM (const double & linkLength); 


__device__
bool IsValidAdhPair (const double & dist); 

__device__
bool IsValidAdhPairForNotInitPhase (const double & dist); 

__device__
double CalAdhECM (const double & dist);

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


struct MoveNodes2_Cell: public thrust::unary_function<IDDBT,DDTI> {
	 double  *_locXAddr_ECM; 
         double  *_locYAddr_ECM; 
        uint _maxMembrNodePerCell ; 
	 int _numNodes_ECM ;
	 double _dt ; 
	 double _Damp_Coef ;
	 bool _isInitPhase ;
	 EType*  _peripORexcmAddr ; 
	__host__ __device__ MoveNodes2_Cell (double * locXAddr_ECM, double * locYAddr_ECM, uint maxMembrNodePerCell, int numNodes_ECM, double dt, double Damp_Coef, bool isInitPhase, EType * peripORexcmAddr) :
				_locXAddr_ECM(locXAddr_ECM),_locYAddr_ECM(locYAddr_ECM),_maxMembrNodePerCell(maxMembrNodePerCell),_numNodes_ECM(numNodes_ECM),_dt(dt),
			    _Damp_Coef(Damp_Coef), _isInitPhase (isInitPhase), _peripORexcmAddr(peripORexcmAddr)	{
	}
	__device__ DDTI  operator()(const IDDBT & iDDBT) const {
	
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
	double distMin=10000 ; 
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
	 if ( nodeIsActive && nodeRankInOneCell<_maxMembrNodePerCell && (nodeType==basal1 || nodeType==apical1)) {
      //  if (_isInitPhase==true) {
			for (int i=0 ; i<_numNodes_ECM ; i++) {
				locX_ECM=_locXAddr_ECM[i]; 
				locY_ECM=_locYAddr_ECM[i];
				dist=sqrt((locX-locX_ECM)*(locX-locX_ECM)+(locY-locY_ECM)*(locY-locY_ECM)) ;
				if (dist < distMin) {
					distMin=dist ; 
					distMinX=(locX_ECM-locX) ;
					distMinY=(locY_ECM-locY) ; 
					iPair=i ; 
				}
				fMorse=calMorse_ECM(dist);
				//if (abs(fMorse)>smallNumber  && _peripORexcmAddr[i]==excm ){
				//	nodeType=basal1 ;
			//	}

			//	if (abs(fMorse)>smallNumber  && _peripORexcmAddr[i]==perip ){
			//		nodeType=apical1 ;
			//	}
				fTotalMorseX=fTotalMorseX+fMorse*(locX_ECM-locX)/dist ; 
				fTotalMorseY=fTotalMorseY+fMorse*(locY_ECM-locY)/dist ; 
				fTotalMorse=fTotalMorse+fMorse ; 
			}

			//if (fTotalMorse!=0.0) {
			//	nodeType=basal1 ; 
		//	}
			if (IsValidAdhPairForNotInitPhase(distMin)&& iPair!=-1) {

        		fAdhMemECM=CalAdhECM(distMin) ; 

				fAdhMemECMX=fAdhMemECM*distMinX/distMin ;  
				fAdhMemECMY=fAdhMemECM*distMinY/distMin ; 
				adhPairECM=iPair ; 
			}

	//	}

	/*	if (_isInitPhase==false) {
			for (int i=0 ; i<_numNodes_ECM ; i++) {
				locX_ECM=_locXAddr_ECM[i]; 
				locY_ECM=_locYAddr_ECM[i];
				dist=sqrt((locX-locX_ECM)*(locX-locX_ECM)+(locY-locY_ECM)*(locY-locY_ECM)) ;
				if (dist < distMin && nodeType==basal1) {
					distMin=dist ; 
					distMinX=(locX_ECM-locX) ;
					distMinY=(locY_ECM-locY) ; 
					iPair=i ; 
				}
				fMorse=calMorse_ECM(dist);  
				fTotalMorseX=fTotalMorseX+fMorse*(locX_ECM-locX)/dist ; 
				fTotalMorseY=fTotalMorseY+fMorse*(locY_ECM-locY)/dist ; 
				fTotalMorse=fTotalMorse+fMorse ; 
			}
			if (IsValidAdhPairForNotInitPhase(distMin)&& iPair!=-1) {

        		fAdhMemECM=CalAdhECM(distMin) ; 

				fAdhMemECMX=fAdhMemECM*distMinX/distMin ;  
				fAdhMemECMY=fAdhMemECM*distMinY/distMin ; 
				adhPairECM=iPair ; 
			}
		}
*/

	 }


    return thrust::make_tuple ((locX+(fTotalMorseX+fAdhMemECMX)*_dt/_Damp_Coef),(locY+(fTotalMorseY+fAdhMemECMY)*_dt/_Damp_Coef),nodeType,adhPairECM)  ; 
        	
		
}
	

	
} ;

        



struct LinSpringForceECM: public thrust::unary_function<IDD,DDD> {
         int   _numNodes ; 	
	 double  _restLen ; 
	 double  _linSpringStiff ;
         double  *_locXAddr; 
         double  *_locYAddr; 
	  

	__host__ __device__ LinSpringForceECM (double numNodes, double restLen, double linSpringStiff , double * locXAddr, double * locYAddr) :
	 _numNodes(numNodes),_restLen(restLen),_linSpringStiff(linSpringStiff),_locXAddr(locXAddr),_locYAddr(locYAddr) {
	}
	 __device__ DDD operator()(const IDD & iDD) const {
	
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
		}
	else {
		locXLeft=_locXAddr[_numNodes-1] ;
		locYLeft=_locYAddr[_numNodes-1] ;
	}

	distLeft=sqrt( ( locX-locXLeft )*(locX-locXLeft) +( locY-locYLeft )*(locY-locYLeft) ) ;
	//	forceLeft=calWLC_ECM(distLeft) ; 
		forceLeft=_linSpringStiff*(distLeft-_restLen) ; 
		forceLeftX=forceLeft*(locXLeft-locX)/distLeft ; 
		forceLeftY=forceLeft*(locYLeft-locY)/distLeft ; 

	if (index != _numNodes-1) {

		
		locXRight=_locXAddr[index+1] ; 
		locYRight=_locYAddr[index+1] ;
		}
	else {

		locXRight=_locXAddr[0] ; 
		locYRight=_locYAddr[0] ;
	}
	distRight=sqrt( ( locX-locXRight )*(locX-locXRight) +( locY-locYRight )*(locY-locYRight) ) ; 
   	    	forceRight=_linSpringStiff*(distRight-_restLen) ; 
        //  	forceRight=calWLC_ECM(distRight) ; 
		forceRightX=forceRight*(locXRight-locX) /distRight ; 
		forceRightY=forceRight*(locYRight-locY) /distRight ; 
		//for open ECM.
	if (index == 0 || index==int(_numNodes/2) ) {
		return thrust::make_tuple(forceRightX,forceRightY,forceRight) ;
   }
   else if (index ==_numNodes-1 || index==(int(_numNodes/2)-1) ) {
		return thrust::make_tuple(forceLeftX,forceLeftY,forceLeft) ;
	}
	else {
		return thrust::make_tuple(forceLeftX+forceRightX,forceLeftY+forceRightY,0.5*(forceLeft+forceRight)) ;
	}
        

}

	
} ;

 
 struct MorseAndAdhForceECM: public thrust::unary_function<IDD,DD> {
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
	DD operator()(const IDD & iDD) const {
	
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
			fTotalMorseX=fTotalMorseX+fMorse*(locX_C-locX)/dist ; 
			fTotalMorseY=fTotalMorseY+fMorse*(locY_C-locY)/dist ; 
			fTotalMorse=fTotalMorse+fMorse ;
			if (_adhPairECM_Cell[i]==index) {
				fAdh=CalAdhECM(dist) ; 
				fAdhX=fAdhX+fAdh*(locX_C-locX)/dist ; 
				fAdhY=fAdhY+fAdh*(locY_C-locY)/dist ; 

			}			 
		}
	}
	
	return thrust::make_tuple(fTotalMorseX+fAdhX,fTotalMorseY+fAdhY) ;

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



struct MoveNodeECM: public thrust::unary_function<DDDDI,DD> {

	double _dt ; 
	double _damp ; 
	int _numNodes ; 
	__host__ __device__ MoveNodeECM (double dt, double damp, int numNodes): _dt(dt), _damp(damp), _numNodes(numNodes) {
	}

	__host__ __device__ DD operator() (const DDDDI & dDDDI) const  {


	double locXOld= thrust:: get <0> (dDDDI) ;
	double locYOld= thrust:: get <1> (dDDDI) ;
	double fx= 	thrust:: get <2> (dDDDI) ;
	double fy= 	thrust:: get <3> (dDDDI) ;
	int    index=    thrust::get<4>(dDDDI) ; 

	if (index == 0 || index==_numNodes-1 || index==( int (_numNodes/2)-1) || index == int (_numNodes/2) ) {
		return thrust::make_tuple (locXOld, locYOld) ;
	}
	else {

		return thrust::make_tuple (locXOld+fx*_dt/_damp, locYOld+fy*_dt/_damp) ;
	}
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

		double bendCenterX = 0;
		double bendCenterY = 0;
		double bendLeftX = 0;
		double bendLeftY = 0;
		double bendRightX = 0;
		double bendRightY = 0;
		double PI=3.141592 ;
		double leftPosX,leftPosY;
		double leftDiffX,leftDiffY;
		double lenLeft;

		double rightPosX,rightPosY;
		double rightDiffX,rightDiffY;
		double lenRight;
			
			int index_left = nodeRank - 1;
			if (index_left == -1) {
				index_left = _numECMNodes - 1;
			}
		
				leftPosX = _locXAddr[index_left];
				leftPosY = _locYAddr[index_left];
				leftDiffX = leftPosX - locX;
				leftDiffY = leftPosY - locY;
				lenLeft = sqrt(leftDiffX * leftDiffX + leftDiffY * leftDiffY);


			int index_right = nodeRank + 1;
			if (index_right ==  _numECMNodes) {
				index_right = 0;
			}
				rightPosX = _locXAddr[index_right];
				rightPosY = _locYAddr[index_right];
				rightDiffX = rightPosX - locX;
				rightDiffY = rightPosY - locY;
				lenRight = sqrt(rightDiffX * rightDiffX + rightDiffY * rightDiffY);
				// if two nodes are extremely close no bending force is applied 
				if (lenLeft>1.0e-8 && lenRight>1.0e-8) {
					double dotP = -leftDiffX * rightDiffX -leftDiffY * rightDiffY;
					double vecP = dotP / (lenLeft * lenRight); //It is cose theta

					// because of numerical error, 1 - vecP*vecP could be less than 0, although it is not possible in mathematics.
					// sqrt(negative number) would cause term0 to be nan.
					// if an nan number is produced, it will not be accepted by bigEnough function.
					// this is OK, because we know at that time bending energy should be 0.
					double term0 = sqrt(1 - vecP * vecP);
					// this if statement is required for numerical purpose only.
					// Whole term would go to zero when term 0 close to zero, but the computation
					// would cause numerical errors, so need to make sure term0 is big enough.
					if (term0>1.0e-7) {
						double angle;
						// value of cross product in z direction: vecA_X * vecB_Y - vecA_Y * vecB_X
						double crossZ = leftDiffY * rightDiffX
								- leftDiffX * rightDiffY;
						if (crossZ > 0) {
							// means angle > PI (concave)
							angle = PI + acos(vecP);
						} else {
							// means angle < PI (convex)
							angle = PI - acos(vecP);
						}
						
						double term1x = -rightDiffX / (lenLeft * lenRight);
						double term2x = leftDiffX / (lenLeft * lenRight);
						double term3x = (dotP * leftDiffX)
								/ (lenLeft * lenLeft * lenLeft * lenRight);
						double term4x = (-dotP * rightDiffX)
								/ (lenLeft * lenRight * lenRight * lenRight);
						double term1y = -rightDiffY / (lenLeft * lenRight);
						double term2y = leftDiffY / (lenLeft * lenRight);
						double term3y = (dotP * leftDiffY)
								/ (lenLeft * lenLeft * lenLeft * lenRight);
						double term4y = (-dotP * rightDiffY)
								/ (lenLeft * lenRight * lenRight * lenRight);

						double bendMultiplier=_eCMBendStiff*(angle-(PI-PI/_numECMNodes)) ; // -calBendMulti_Mitotic(angle,
								//activeMembrCount, progress, _mitoticCri);//AAMIRI modified the arguments
						// because sign of angle formula would change if crossZ < 0
						if (crossZ > 0) {
							bendMultiplier = -bendMultiplier;
						}
						bendLeftX = bendMultiplier * (term1x - term3x) / term0;
                        bendCenterX= bendMultiplier* (term2x - term1x + term3x - term4x)/ term0 ;
						bendRightX = bendMultiplier * (term4x - term2x) / term0;
						bendLeftY = bendMultiplier * (term1y - term3y) / term0;
						bendCenterY=bendMultiplier* (term2y - term1y + term3y - term4y)/ term0;
						bendRightY = bendMultiplier * (term4y - term2y) / term0;

				}
			}
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



