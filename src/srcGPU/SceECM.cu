#include "SceECM.h"

/*
struct InitECMX 
{
  const double  a ; 
 
  InitECMX (double _a) : a(_a) {}

  __host__ __device__ 

   double operator() (const double & x , const double & y) const {
   return (x*a-50) ; 
  }
} ; 
*/
void SceECM::applyECMForce_M() {
double Ali=5 ; 
} 



void SceECM::Initialize() {

double eCMMinX= -50;
double eCMMaxX= 50;
// eCMY=24.0;
double eCMMinDist=0.04;
eCMY=23.7 ; 
int  numNodesECM= (eCMMaxX-eCMMinX)/eCMMinDist ; 
thrust:: device_vector<double> tmp ; 
cout<<" I am inside ECM initialization0  " << endl ; 
tmp.resize(numNodesECM,0.0) ;
nodeECMLocX.resize(numNodesECM,0.0) ;
nodeECMLocY.resize(numNodesECM,eCMY) ;

thrust::sequence (tmp.begin(),tmp.begin()+numNodesECM);

 
thrust::fill (nodeECMLocX.begin(),nodeECMLocX.begin()+numNodesECM,0.0); 
thrust::fill (nodeECMLocY.begin(),nodeECMLocY.begin()+numNodesECM,eCMY); 
thrust::transform(tmp.begin(),tmp.begin()+numNodesECM,nodeECMLocY.begin(),nodeECMLocX.begin(),InitECMLoc(eCMMinX,eCMMinDist));

for (int i=0;  i<nodeECMLocX.size() ;  i++) {
cout<< nodeECMLocX[i]<<endl; 
} 
cout<<" I am inside ECM initialization4  " << endl ; 

int t=0 ; 
std::string vtkFileName = "ECM_" + std::to_string(t) + ".vtk";
			ofstream ECMOut;
			ECMOut.open(vtkFileName.c_str());
			ECMOut<< "# vtk DataFile Version 3.0" << endl;
			ECMOut<< "Result for paraview 2d code" << endl;
			ECMOut << "ASCII" << endl;
			ECMOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
			ECMOut << "POINTS " << nodeECMLocX.size() << " float" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				fs << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
				<< 0.0 << std::endl;
			}
			ECMOut<< std::endl;
			 ECMOut<< "CELLS " << nodeECMLocX.size()-1<< " " << 3 *(nodeECMLocX.size()-1)<< std::endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << 2 << " " << tmp[i]+1 << " "
				<< tmp[i+1]+1 << std::endl;
			}

			ECMOut << "CELL_TYPES " << nodeECMLocX.size()-1<< endl;
			for (uint i = 0; i < (nodeECMLocX.size()-1); i++) {
				ECMOut << "3" << endl;
			}

			ECMOut.close();

}

void SceECM:: ApplyECMConstrain(int totalNodeCountForActiveCellsECM){   

thrust::counting_iterator<int> iBegin(0) ; 
nodeDeviceTmpLocX.resize(totalNodeCountForActiveCellsECM,0.0) ;
nodeDeviceTmpLocY.resize(totalNodeCountForActiveCellsECM,0.0) ;
 
thrust::copy(nodeDeviceLocX.begin(),nodeDeviceLocX.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocX.begin()) ; 
thrust::copy(nodeDeviceLocY.begin(),nodeDeviceLocY.begin()+totalNodeCountForActiveCellsECM,nodeDeviceTmpLocY.begin()) ; 


double yBound=10 ;
 eCMY=23.7 ;  
int maxAllNodePerCell=360 ;
int maxMembrNodePerCell=280 ;  

thrust:: transform (
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					nodeDeviceTmpLocX.begin(),
					nodeDeviceTmpLocY.begin())), 
		thrust::make_zip_iterator (
				thrust:: make_tuple (
					 make_transform_iterator (iBegin,
							ModuloFunctor2(maxAllNodePerCell)),
					 nodeDeviceTmpLocX.begin(),
                                         nodeDeviceTmpLocY.begin()))+totalNodeCountForActiveCellsECM,
		thrust::make_zip_iterator (
				thrust::make_tuple (
					nodeDeviceLocX.begin(),
					nodeDeviceLocY.begin())),
				MyFunctor2(eCMY,yBound,maxMembrNodePerCell)); 


/*
void print ECM
std::string vtkFileName = "DPP_" + std::to_string(t) + ".vtk";
			ofstream ECMOut;
			ECMOut.open(vtkFileName.c_str());
			ECMOut<< "# vtk DataFile Version 3.0" << endl;
			ECMOut<< "Result for paraview 2d code" << endl;
			ECMOut << "ASCII" << endl;
			ECMOut << "DATASET UNSTRUCTURED_GRID" << std::endl;
			ECMOut << "POINTS " << pointsAniData.size() << " float" << std::endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				fs << nodeECMLocX[i] << " " << nodeECMLocY[i] << " "
				<< 0.0 << std::endl;
			}
			ECMOut<< std::endl;
			 ECMOut<< "CELLS " << nodeECMLocX.size()-1<< " " << 3 *(nodeECMLocX.size()-1)<< std::endl;
			for (uint i = 0; i < nodeECMLocX.size()-1; i++) {
				ECMOut << 2 << " " << tmp[i]+1 << " "
				<< tmp[i+1]+1 << std::endl;
			}

			ECMOut << "CELL_TYPES " << linksAniData.size() << endl;
			for (uint i = 0; i < nodeECMLocX.size(); i++) {
				ECMOut << "3" << endl;
			}
std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << rank;
	std::string scriptNameRank = ss.str();
	std::string vtkFileName = scriptNameBase + scriptNameRank + ".vtk";
	std::cout << "start to create vtk file" << vtkFileName << std::endl;
	std::ofstream fs;
	fs.open(vtkFileName.c_str());
	fs << "# vtk DataFile Version 3.0" << std::endl;
	fs << "Lines and points representing subcelluar element cells "
			<< std::endl;
	fs << "ASCII" << std::endl;
	fs << std::endl;
	fs << "DATASET UNSTRUCTURED_GRID" << std::endl;
	fs << "POINTS " << pointsAniData.size() << " float" << std::endl;
	for (uint i = 0; i < pointsAniData.size(); i++) {
		fs << pointsAniData[i].pos.x << " " << pointsAniData[i].pos.y << " "
				<< pointsAniData[i].pos.z << std::endl;
	}
fs << std::endl;
	fs << "CELLS " << linksAniData.size() << " " << 3 * linksAniData.size()
			<< std::endl;
	for (uint i = 0; i < linksAniData.size(); i++) {
		fs << 2 << " " << linksAniData[i].node1Index << " "
				<< linksAniData[i].node2Index << std::endl;
	}
	fs << "CELL_TYPES " << linksAniData.size() << endl;
	for (uint i = 0; i < linksAniData.size(); i++) {
		fs << "3" << endl;
	}
	fs << "POINT_DATA " << pointsAniData.size() << endl;
	//fs << "SCALARS relative_tension float" << endl;  //Ali
	fs << "SCALARS Number_of_Neighbouring_Cells float" << endl;
	fs << "LOOKUP_TABLE default" << endl;

	for (uint i = 0; i < pointsAniData.size(); i++) {
		fs << pointsAniData[i].colorScale << endl;
	}

	fs << std::endl;

	//AAMIR wrote the curvature data here
	fs << "SCALARS Stress float" << endl;
	fs << "LOOKUP_TABLE default" << endl; //AliE

	for (uint i = 0; i < pointsAniData.size(); i++) {
		//std::cout << "****************      SIZE IS:     " << pointsAniData.size() << std::endl;
//		if (pointsAniData[i].colorScale2 < 0.000001 || pointsAniData[i].colorScale2>1.0){
//			pointsAniData[i].colorScale2 = 0.0;}
		fs << pointsAniData[i].F_MI_M_MagN_Int << endl;  //AliE
	}

	fs << std::endl;
	//AAMIRI finished writing the Node Curvature


	//AAMIR wrote the cell rank information data here
	fs << "SCALARS cellRank int" << endl;
	fs << "LOOKUP_TABLE default" << endl;
	for (uint i = 0; i < pointsAniData.size(); i++) {
fs << pointsAniData[i].rankScale << endl;
	}
	//AAMIRI finished writing the cell rank of each node


	//AAMIRI starts writing tension vector data
	fs << "VECTORS F_MI_M float" << endl;
		for (uint i = 0; i < pointsAniData.size(); i++) {

			fs << pointsAniData[i].F_MI_M.x << " " << pointsAniData[i].F_MI_M.y << " "
					<< pointsAniData[i].F_MI_M.z << endl;
		}
	//AAMIRI finished writing the node tension vector

	//AAMIRI starts writing external force vector
 	fs << "VECTORS ExternalForce float" << endl;
		for (uint i = 0; i < pointsAniData.size(); i++) {

			fs << pointsAniData[i].extForce.x << " " << pointsAniData[i].extForce.y << " "
					<< pointsAniData[i].extForce.z << endl;
		}
	//AAMIRI finished writing the node ext force vector


	if (isArrowIncluded) {
		fs << "VECTORS vectors float" << endl;
		for (uint i = 0; i < pointsAniData.size(); i++) {
			fs << pointsAniData[i].dir.x << " " << pointsAniData[i].dir.y << " "
					<< pointsAniData[i].dir.z << endl;
		}
	}
	fs.close();
}




}
*/
