#ifndef SimulationDomainGPU_H_
#define SimulationDomainGPU_H_

#include "SceNodes.h"
#include "SceCells.h"
#include "commonData.h"
#include "NetworkInfo.h"

#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

/**
 * This class will control the process of stabilizing cell centers.
 */
class StabPara {
public:
	bool isProcessStab;

	int outputFrameCount;
	int totalIterCount;
	double bdrySpacingRatio;
	double dt;
	std::string outputAniName;
};

/**
 * This class is responsible for domain-wise highest level logic, e.g. output animation.
 */
class SimulationDomainGPU {
	/**
	 * Variable that contains information for nodes.
	 * Handles node level interaction logic.
	 */
	SceNodes nodes;

	/**
	 * Variable that contains information for cells.
	 * Handles cell level logics like growth and division.
	 */
	SceCells cells;

	NetworkInfo netInfo;

	std::vector<std::vector<PreT1State> > preT1Vec;

	std::set<int> t1CellSet;
	std::vector<double> cellColorVec;

	/**
	 * memory related parameters.
	 */
	SceMemPara memPara;

	/**
	 * domain related parameters.
	 */
	SceDomainPara domainPara;

	/**
	 * chemical related parameters.
	 */
	SceChemPara chemPara;

	/**
	 * parameters used for stabilize the initial cell positions
	 */
	StabPara stabPara;

	/**
	 * reads memory related parameters.
	 */
	void readMemPara();

	/**
	 * reads domain related parameters.
	 */
	void readDomainPara();

	/**
	 * reads chemical related parameters.
	 */
	void readChemPara();

	/**
	 * reads all parameters by calling all other reading methods.
	 */
	void readAllParameters();

	/**
	 * Initializes data vectors by given vectors.
	 * improved from the previous version.
	 */
	void initializeNodes_M(std::vector<SceNodeType> &nodeTypes,
			std::vector<bool> &nodeIsActive, std::vector<CVector> &initNodesVec,
			std::vector<uint> &numOfInitActiveEpiNodeCounts,
			std::vector<uint> &numOfInitActiveInternalNodeCounts,
			std::vector<double> &initGrowProgVec);

	NetworkInfo buildNetInfo(CellsStatsData &polyData);
	std::set<int> findT1Transition();

	void outputVtkGivenCellColor(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri, std::vector<double>& cellColorVec);
	std::vector<double> processT1Color();

public:
	/**
	 * Default constructor.
	 * Reads all the configuration values from file.
	 */
	SimulationDomainGPU();

	/**
	 * Domain initialization.
	 * Assigns values to the data fields in simulation domain.
	 * @param initData initial data set for simulation domain
	 */
	void initialize_v2(SimulationInitData_V2 &initData);

	/**
	 * Domain initialization.
	 * Assigns values to the data fields in simulation domain.
	 * @param initData initial data set for simulation domain
	 */
	void initialize_v2_M(SimulationInitData_V2_M &initData);

	/**
	 * This method is used for producing more or less evenly distributed cell center positions.
	 * Because CGAL mesh cannot guarantee that mesh vertices are relatively evenly distributed,
	 * we would need this extra step in order to stabilize the cell center positions
	 * so that cells can be allocated evenly in our simulation domain.
	 */
	std::vector<CVector> stablizeCellCenters(SimulationInitData_V2 &initData);

	/**
	 * Checks if all data fields are valid.
	 * This methods only loosely checks the data validity.
	 */
	void printDomainInformation();

	/**
	 * Run one step of simulation in the domain.
	 * Contains cell level logics and node level logics.
	 * @param dt timestep
	 */
	void runAllLogic(double dt);

	/**
	 * Run one step of simulation in the domain.
	 * Contains cell level logics and node level logics.
	 * @param dt timestep
	 */
	void runAllLogic_M(double dt);

	bool isDividing_ForAni();

	/**
	 * Method that animates the domain to VTK format.
	 * @param scriptNameBase name of the vtk animation series.
	 * @param rank frame sequence in the vtk animation series.
	 * @param aniCri criteria for outputing animation.
	 */
	void outputVtkFilesWithCri(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);

	/**
	 * Method that animates the domain to VTK format.
	 * @param scriptNameBase name of the vtk animation series.
	 * @param rank frame sequence in the vtk animation series.
	 * @param aniCri criteria for outputing animation.
	 */
	void outputVtkFilesWithCri_M(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);

	void outputVtkColorByCell(std::string scriptNameBase, int rank,
			AnimationCriteria aniCri);
	/**
	 * Method that animates the domain to VTK format.
	 * @param resultNameBase name of the labelMatrix series.
	 * @param rank frame sequence in the labelMatrix series.
	 * @param pixelPara criteria for labelMatrix generation.
	 */
	vector<vector<int> > outputLabelMatrix(std::string resultNameBase, int rank,
			PixelizePara &pixelPara);

	/**
	 * method that prints out the growth progress vector of cells to a file.
	 */
	void outputGrowthProgressAuxFile(int step);

	/**
	 * Post processing for the label matrix.
	 */
	void analyzeLabelMatrix(vector<vector<int> > &labelMatrix, int step,
			std::string &imageFileNameBase, std::string &statFileName);

	void performAblation(AblationEvent &ablEvent);

	CellsStatsData outputPolyCountData();

	void processT1Info(int maxStepTraceBack, CellsStatsData &polyData);
};

#endif
