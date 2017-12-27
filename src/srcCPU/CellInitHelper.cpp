
 /* CellInitHelper.cpp
 *
 *  Created on: Sep 22, 2013
 *      Author: wsun2
 */
#include <fstream>
#include "CellInitHelper.h"
//Ali 
ForReadingData_M2 ReadFile_M2(std::string CellCentersFileName) {

          std::vector<GEOMETRY::Point2D> result;
          std::fstream inputc;
          ForReadingData_M2  ForReadingData1; 
          double TempPos_X,TempPos_Y,TempPos_Z ; 
		  string eCellTypeString ; 
          //inputc.open("./resources/CellCenters_General.txt");
          //inputc.open("./resources/CellCenters_General.txt");
//          inputc.open("./resources/CellCenters2.txt");
          inputc.open(CellCentersFileName.c_str());

          if (inputc.is_open())
          {
            cout << "File successfully open";
          }
         else
         {
          cout << "Error opening file";
          }
          inputc >> ForReadingData1.CellNumber ; 
          for (int i = 0; i <ForReadingData1.CellNumber; i = i + 1) {
	    cout << "i=" << i << endl;		
	    inputc >> TempPos_X >> TempPos_Y >> TempPos_Z>> eCellTypeString ;	
	    ForReadingData1.TempSX.push_back(TempPos_X);
	    ForReadingData1.TempSY.push_back(TempPos_Y);
	    ForReadingData1.TempSZ.push_back(TempPos_Z);
		ECellType eCellType=StringToECellTypeConvertor (eCellTypeString) ; 
		ForReadingData1.eCellTypeV.push_back (eCellType) ; 
            }     
         cout << "Cell center positions read successfully";
		
          for (int i = 0; i <ForReadingData1.CellNumber; i = i + 1) {
			  cout << ForReadingData1.eCellTypeV.at(i) << "cell type read" << endl ;				 }
return ForReadingData1;
}

ECellType StringToECellTypeConvertor ( const string & eCellTypeString) {
if (eCellTypeString=="bc") {return bc ;  }
else if (eCellTypeString=="peri") {return peri ;  }
else if (eCellTypeString=="pouch") {return pouch ;  }
//else { throw std:: invslid_argument ("recevied invalid cell type") ; } 

}
//Ali 


CellInitHelper::CellInitHelper() {
	int type = globalConfigVars.getConfigValue("SimulationType").toInt();
	simuType = parseTypeFromConfig(type);
	if (simuType == Beak) {
		initInternalBdry();
	}
}

CVector CellInitHelper::getPointGivenAngle(double currentAngle, double r,
		CVector centerPos) {
	double xPos = centerPos.x + r * cos(currentAngle);
	double yPos = centerPos.y + r * sin(currentAngle);
	return CVector(xPos, yPos, 0);
}

RawDataInput CellInitHelper::generateRawInputWithProfile(
		std::vector<CVector> &cellCenterPoss, bool isInnerBdryIncluded) {

	RawDataInput rawData;
	vector<CVector> outsideBdryNodePos;
	vector<CVector> outsideProfileNodePos;
	std::string bdryInputFileName = globalConfigVars.getConfigValue(
			"Bdry_InputFileName").toString();

	GEOMETRY::MeshGen meshGen;

	double genBdryRatio =
			globalConfigVars.getConfigValue("GenBdrySpacingRatio").toDouble();

	GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2DWithProfile(
			bdryInputFileName, genBdryRatio, isInnerBdryIncluded);

	std::vector<GEOMETRY::Point2D> bdryPoints = mesh.getFinalBdryPts();

	std::vector<GEOMETRY::Point2D> profilePoints = mesh.getFinalProfilePts();

	for (uint i = 0; i < bdryPoints.size(); i++) {
		outsideBdryNodePos.push_back(
				CVector(bdryPoints[i].getX(), bdryPoints[i].getY(), 0));
	}

	rawData.bdryNodes = outsideBdryNodePos;

	for (uint i = 0; i < profilePoints.size(); i++) {
		outsideProfileNodePos.push_back(
				CVector(profilePoints[i].getX(), profilePoints[i].getY(), 0));
	}

	rawData.profileNodes = outsideProfileNodePos;

	// calculate average length of profile links
	double sumLength = 0;
	for (uint i = 0; i < profilePoints.size() - 1; i++) {
		CVector tmpVec = outsideProfileNodePos[i]
				- outsideProfileNodePos[i + 1];
		sumLength += tmpVec.getModul();
	}

	for (unsigned int i = 0; i < cellCenterPoss.size(); i++) {
		CVector centerPos = cellCenterPoss[i];
		centerPos.Print();
		if (isMXType(centerPos)) {
			rawData.MXCellCenters.push_back(centerPos);
		} else {
			rawData.FNMCellCenters.push_back(centerPos);
		}

	}

	generateCellInitNodeInfo_v2(rawData.initCellNodePoss);

	return rawData;
}

RawDataInput CellInitHelper::generateRawInput_simu(
		std::vector<CVector>& cellCenterPoss) {
	if (simuType == Beak) {
		RawDataInput baseRawInput = generateRawInputWithProfile(cellCenterPoss,
				false);

		GEOMETRY::MeshGen meshGen;
		double genBdryRatio = globalConfigVars.getConfigValue(
				"GenBdrySpacingRatio").toDouble();
		std::string bdryInputFileName = globalConfigVars.getConfigValue(
				"Bdry_InputFileName").toString();
		GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2DWithProfile(
				bdryInputFileName, genBdryRatio);
		GEOMETRY::MeshInput input = meshGen.obtainMeshInput();
		baseRawInput.cartilageData = meshGen.obtainCartilageData(mesh, input);

		baseRawInput.isStab = false;
		baseRawInput.simuType = simuType;
		return baseRawInput;
	} else if (simuType == Disc) {
		RawDataInput rawInput;
		initializeRawInput(rawInput, cellCenterPoss);
		rawInput.isStab = false;
		rawInput.simuType = simuType;
		return rawInput;
	} else {
		throw SceException(
				"Simulation Type is not defined when trying to generate raw input",
				InvalidInput);
	}
}

void CellInitHelper::transformRawCartData(CartilageRawData& cartRawData,
		CartPara& cartPara, std::vector<CVector>& initNodePos) {
	// step 1, switch tip node1 to pos 0
	CVector tmpPos = cartRawData.tipVerticies[0];
	cartRawData.tipVerticies[0] =
			cartRawData.tipVerticies[cartRawData.growNode1Index_on_tip];
	cartRawData.tipVerticies[cartRawData.growNode1Index_on_tip] = tmpPos;
	cartPara.growNode1Index = 0;

	// step 2, switch tip node 2 to pos 1
	tmpPos = cartRawData.tipVerticies[1];
	cartRawData.tipVerticies[1] =
			cartRawData.tipVerticies[cartRawData.growNode2Index_on_tip];
	cartRawData.tipVerticies[cartRawData.growNode2Index_on_tip] = tmpPos;
	cartPara.growNode2Index = 1;
	cartPara.tipNodeStartPos = 2;
	cartPara.tipNodeIndexEnd = cartRawData.tipVerticies.size();

	// step 3, calculate size for tip nodes
	double tipMaxExpansionRatio = globalConfigVars.getConfigValue(
			"TipMaxExpansionRatio").toDouble();
	int maxTipSize = tipMaxExpansionRatio * cartRawData.tipVerticies.size();
	cartPara.nonTipNodeStartPos = maxTipSize;
	cartPara.nodeIndexEnd = cartPara.nonTipNodeStartPos
			+ cartRawData.nonTipVerticies.size();
	cartPara.pivotNode1Index = cartPara.nonTipNodeStartPos
			+ cartRawData.pivotNode1Index;
	cartPara.pivotNode2Index = cartPara.nonTipNodeStartPos
			+ cartRawData.pivotNode2Index;
	cartPara.growNodeBehind1Index = cartPara.nonTipNodeStartPos
			+ cartRawData.growNodeBehind1Index;
	cartPara.growNodeBehind2Index = cartPara.nonTipNodeStartPos
			+ cartRawData.growNodeBehind2Index;

	// step 4, calculate size for all nodes
	double cartmaxExpRatio = globalConfigVars.getConfigValue(
			"CartMaxExpansionRatio").toDouble();
	int maxCartNodeSize = (cartRawData.tipVerticies.size()
			+ cartRawData.nonTipVerticies.size()) * cartmaxExpRatio;
	cartPara.nodeIndexTotal = maxCartNodeSize;

	// step 5, initialize the first part of initNodePos
	initNodePos.resize(cartPara.nodeIndexTotal);
	for (uint i = 0; i < cartRawData.tipVerticies.size(); i++) {
		initNodePos[i] = cartRawData.tipVerticies[i];
	}

	// step 6, initialize the second part of initNodePos
	for (uint i = 0; i < cartRawData.nonTipVerticies.size(); i++) {
		initNodePos[i + cartPara.nonTipNodeStartPos] =
				cartRawData.nonTipVerticies[i];
	}

	for (uint i = 0; i < initNodePos.size(); i++) {
		initNodePos[i].Print();
	}

}

void CellInitHelper::initializeRawInput(RawDataInput& rawInput,
		std::vector<CVector>& cellCenterPoss) {

	for (unsigned int i = 0; i < cellCenterPoss.size(); i++) {
		CVector centerPos = cellCenterPoss[i];
		if (isMXType(centerPos)) {
			rawInput.MXCellCenters.push_back(centerPos);
		} else {
			rawInput.FNMCellCenters.push_back(centerPos);
		}
	}

	generateCellInitNodeInfo_v2(rawInput.initCellNodePoss);

}

/**
 * Initialize inputs for five different components.
 */

SimulationInitData CellInitHelper::initInputsV2(RawDataInput &rawData) {
	SimulationInitData initData;

	uint FnmCellCount = rawData.FNMCellCenters.size();
	uint MxCellCount = rawData.MXCellCenters.size();
	uint ECMCount = rawData.ECMCenters.size();

	uint maxNodePerCell =
			globalConfigVars.getConfigValue("MaxNodePerCell").toInt();
	uint maxNodePerECM = 0;
	if (rawData.simuType == Beak) {
		maxNodePerECM =
				globalConfigVars.getConfigValue("MaxNodePerECM").toInt();
	}

	uint initTotalCellCount = rawData.initCellNodePoss.size();
	initData.initBdryCellNodePosX.resize(rawData.bdryNodes.size(), 0.0);
	initData.initBdryCellNodePosY.resize(rawData.bdryNodes.size(), 0.0);
	initData.initProfileNodePosX.resize(rawData.profileNodes.size());
	initData.initProfileNodePosY.resize(rawData.profileNodes.size());
	initData.initECMNodePosX.resize(maxNodePerECM * ECMCount);
	initData.initECMNodePosY.resize(maxNodePerECM * ECMCount);
	initData.initFNMCellNodePosX.resize(maxNodePerCell * FnmCellCount, 0.0);
	initData.initFNMCellNodePosY.resize(maxNodePerCell * FnmCellCount, 0.0);
	initData.initMXCellNodePosX.resize(maxNodePerCell * MxCellCount, 0.0);
	initData.initMXCellNodePosY.resize(maxNodePerCell * MxCellCount, 0.0);

	for (uint i = 0; i < FnmCellCount; i++) {
		initData.cellTypes.push_back(FNM);
		initData.numOfInitActiveNodesOfCells.push_back(initTotalCellCount);
	}

	for (uint i = 0; i < MxCellCount; i++) {
		initData.cellTypes.push_back(MX);
		initData.numOfInitActiveNodesOfCells.push_back(initTotalCellCount);
	}

	for (uint i = 0; i < rawData.bdryNodes.size(); i++) {
		initData.initBdryCellNodePosX[i] = rawData.bdryNodes[i].x;
		initData.initBdryCellNodePosY[i] = rawData.bdryNodes[i].y;
	}

	for (uint i = 0; i < rawData.profileNodes.size(); i++) {
		initData.initProfileNodePosX[i] = rawData.profileNodes[i].x;
		initData.initProfileNodePosY[i] = rawData.profileNodes[i].y;
	}

	uint index;

	uint ECMInitNodeCount = rawData.initECMNodePoss.size();
	for (uint i = 0; i < ECMCount; i++) {
		vector<CVector> rotatedCoords = rotate2D(rawData.initECMNodePoss,
				rawData.ECMAngles[i]);
		for (uint j = 0; j < ECMInitNodeCount; j++) {
			index = i * maxNodePerECM + j;
			initData.initECMNodePosX[index] = rawData.ECMCenters[i].x
					+ rotatedCoords[j].x;
			initData.initECMNodePosY[index] = rawData.ECMCenters[i].y
					+ rotatedCoords[j].y;
		}
	}

	for (uint i = 0; i < FnmCellCount; i++) {
		for (uint j = 0; j < initTotalCellCount; j++) {
			index = i * maxNodePerCell + j;
			initData.initFNMCellNodePosX[index] = rawData.FNMCellCenters[i].x
					+ rawData.initCellNodePoss[j].x;
			initData.initFNMCellNodePosY[index] = rawData.FNMCellCenters[i].y
					+ rawData.initCellNodePoss[j].y;
		}
	}

	for (uint i = 0; i < MxCellCount; i++) {
		for (uint j = 0; j < initTotalCellCount; j++) {
			index = i * maxNodePerCell + j;
			initData.initMXCellNodePosX[index] = rawData.MXCellCenters[i].x
					+ rawData.initCellNodePoss[j].x;
			initData.initMXCellNodePosY[index] = rawData.MXCellCenters[i].y
					+ rawData.initCellNodePoss[j].y;
		}
	}

	return initData;
}

SimulationInitData_V2 CellInitHelper::initInputsV3(RawDataInput& rawData) {
	SimulationInitData_V2 initData;
	initData.isStab = rawData.isStab;
	initData.simuType = rawData.simuType;

	uint FnmCellCount = rawData.FNMCellCenters.size();
	uint MxCellCount = rawData.MXCellCenters.size();
	uint ECMCount = rawData.ECMCenters.size();

	uint maxNodePerCell =
			globalConfigVars.getConfigValue("MaxNodePerCell").toInt();
	uint maxNodePerECM = 0;
	if (initData.simuType == Beak) {
		maxNodePerECM =
				globalConfigVars.getConfigValue("MaxNodePerECM").toInt();
	}

	uint initTotalCellCount = rawData.initCellNodePoss.size();
	//uint initTotalECMCount = rawData.ECMCenters.size();
	initData.initBdryNodeVec.resize(rawData.bdryNodes.size());
	initData.initProfileNodeVec.resize(rawData.profileNodes.size());

	if (simuType == Beak && !initData.isStab) {
		transformRawCartData(rawData.cartilageData, initData.cartPara,
				initData.initCartNodeVec);
	}

	initData.initECMNodeVec.resize(maxNodePerECM * ECMCount);
	initData.initFNMNodeVec.resize(maxNodePerCell * FnmCellCount);
	initData.initMXNodeVec.resize(maxNodePerCell * MxCellCount);

	for (uint i = 0; i < FnmCellCount; i++) {
		initData.cellTypes.push_back(FNM);
		initData.numOfInitActiveNodesOfCells.push_back(initTotalCellCount);
	}

	for (uint i = 0; i < MxCellCount; i++) {
		initData.cellTypes.push_back(MX);
		initData.numOfInitActiveNodesOfCells.push_back(initTotalCellCount);
	}

	for (uint i = 0; i < rawData.bdryNodes.size(); i++) {
		initData.initBdryNodeVec[i] = rawData.bdryNodes[i];
	}

	for (uint i = 0; i < rawData.profileNodes.size(); i++) {
		initData.initProfileNodeVec[i] = rawData.profileNodes[i];
	}

	uint index;

	uint ECMInitNodeCount = rawData.initECMNodePoss.size();
	for (uint i = 0; i < ECMCount; i++) {
		vector<CVector> rotatedCoords = rotate2D(rawData.initECMNodePoss,
				rawData.ECMAngles[i]);
		for (uint j = 0; j < ECMInitNodeCount; j++) {
			index = i * maxNodePerECM + j;
			initData.initECMNodeVec[index] = rawData.ECMCenters[i]
					+ rotatedCoords[j];
		}
	}

	for (uint i = 0; i < FnmCellCount; i++) {
		for (uint j = 0; j < initTotalCellCount; j++) {
			index = i * maxNodePerCell + j;
			initData.initFNMNodeVec[index] = rawData.FNMCellCenters[i]
					+ rawData.initCellNodePoss[j];
		}
	}

	for (uint i = 0; i < MxCellCount; i++) {
		for (uint j = 0; j < initTotalCellCount; j++) {
			index = i * maxNodePerCell + j;
			initData.initMXNodeVec[index] = rawData.MXCellCenters[i]
					+ rawData.initCellNodePoss[j];
		}
	}

	return initData;
}

SimulationInitData_V2_M CellInitHelper::initInputsV3_M(
		RawDataInput_M& rawData_m) {

	if (rawData_m.simuType != Disc_M) {
		throw SceException("V2_M data can be used for Disc_M simulation only!",
				ConfigValueException);
	}

	uint maxMembrNodePerCell = globalConfigVars.getConfigValue(
			"MaxMembrNodeCountPerCell").toInt();
	uint maxAllNodeCountPerCell = globalConfigVars.getConfigValue(
			"MaxAllNodeCountPerCell").toInt();
	uint maxCellInDomain =
			globalConfigVars.getConfigValue("MaxCellInDomain").toInt();
	uint maxNodeInDomain = maxCellInDomain * maxAllNodeCountPerCell;
	uint initCellCount = rawData_m.initCellCenters.size();

	uint initMaxNodeCount = initCellCount * maxAllNodeCountPerCell;

	uint cellRank, nodeRank, activeMembrNodeCountThisCell,
			activeIntnlNodeCountThisCell;

	SimulationInitData_V2_M initData;
	initData.isStab = rawData_m.isStab;
	initData.simuType = rawData_m.simuType;

	initData.nodeTypes.resize(maxNodeInDomain);
	initData.initNodeVec.resize(initMaxNodeCount);
	initData.initIsActive.resize(initMaxNodeCount, false);
	//initData.initGrowProgVec.resize(initCellCount, 0);
    ECellType eCellTypeTmp2 ; // Ali
	for (uint i = 0; i < initCellCount; i++) {
		initData.initActiveMembrNodeCounts.push_back(
				rawData_m.initMembrNodePoss[i].size());
		initData.initActiveIntnlNodeCounts.push_back(
				rawData_m.initIntnlNodePoss[i].size());
		initData.initGrowProgVec.push_back(rawData_m.cellGrowProgVec[i]);
		eCellTypeTmp2=rawData_m.cellsTypeCPU.at(i) ;  // Ali
		initData.eCellTypeV1.push_back(eCellTypeTmp2); // Ali

	}

	for (uint i = 0; i < initCellCount; i++) {
		cout << "third check cell type" << initData.eCellTypeV1.at(i)<< endl ;  // Ali
	}
	for (uint i = 0; i < maxNodeInDomain; i++) {
		nodeRank = i % maxAllNodeCountPerCell;
		if (nodeRank < maxMembrNodePerCell) {
			initData.nodeTypes[i] = CellMembr;
		} else {
			initData.nodeTypes[i] = CellIntnl;
		}
	}

	for (uint i = 0; i < initMaxNodeCount; i++) {
		cellRank = i / maxAllNodeCountPerCell;
		nodeRank = i % maxAllNodeCountPerCell;
		activeMembrNodeCountThisCell =
				rawData_m.initMembrNodePoss[cellRank].size();
		activeIntnlNodeCountThisCell =
				rawData_m.initIntnlNodePoss[cellRank].size();
		if (nodeRank < maxMembrNodePerCell) {
			if (nodeRank < activeMembrNodeCountThisCell) {
				initData.initNodeVec[i] =
						rawData_m.initMembrNodePoss[cellRank][nodeRank];
				initData.initIsActive[i] = true;
			} else {
				initData.initIsActive[i] = false;
			}
		} else {
			uint intnlIndex = nodeRank - maxMembrNodePerCell;
			if (intnlIndex < activeIntnlNodeCountThisCell) {
				initData.initNodeVec[i] =
						rawData_m.initIntnlNodePoss[cellRank][intnlIndex];
				initData.initIsActive[i] = true;
			} else {
				initData.initIsActive[i] = false;
			}
		}
	}

	return initData;
}

vector<CVector> CellInitHelper::rotate2D(vector<CVector> &initECMNodePoss,
		double angle) {
	uint inputVectorSize = initECMNodePoss.size();
	CVector centerPosOfInitVector = CVector(0);
	for (uint i = 0; i < inputVectorSize; i++) {
		centerPosOfInitVector = centerPosOfInitVector + initECMNodePoss[i];
	}
	centerPosOfInitVector = centerPosOfInitVector / inputVectorSize;
	for (uint i = 0; i < inputVectorSize; i++) {
		initECMNodePoss[i] = initECMNodePoss[i] - centerPosOfInitVector;
	}
	vector<CVector> result;
	for (uint i = 0; i < inputVectorSize; i++) {
		CVector posNew;
		posNew.x = cos(angle) * initECMNodePoss[i].x
				- sin(angle) * initECMNodePoss[i].y;
		posNew.y = sin(angle) * initECMNodePoss[i].x
				+ cos(angle) * initECMNodePoss[i].y;
		result.push_back(posNew);
	}
	return result;
}

RawDataInput CellInitHelper::generateRawInput_stab() {
	RawDataInput rawData;
	rawData.simuType = simuType;
	vector<CVector> insideCellCenters;
	vector<CVector> outsideBdryNodePos;
	std::string bdryInputFileName = globalConfigVars.getConfigValue(
			"Bdry_InputFileName").toString();

	GEOMETRY::MeshGen meshGen;

	GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2DFromFile(
			bdryInputFileName);

	std::vector<GEOMETRY::Point2D> insideCenterPoints =
			mesh.getAllInsidePoints();

	double fine_Ratio =
			globalConfigVars.getConfigValue("StabBdrySpacingRatio").toDouble();

	for (uint i = 0; i < insideCenterPoints.size(); i++) {
		insideCellCenters.push_back(
				CVector(insideCenterPoints[i].getX(),
						insideCenterPoints[i].getY(), 0));
	}

	mesh = meshGen.generateMesh2DFromFile(bdryInputFileName, fine_Ratio);

	std::vector<GEOMETRY::Point2D> bdryPoints = mesh.getOrderedBdryPts();

	for (uint i = 0; i < bdryPoints.size(); i++) {
		outsideBdryNodePos.push_back(
				CVector(bdryPoints[i].getX(), bdryPoints[i].getY(), 0));
	}

	for (unsigned int i = 0; i < insideCellCenters.size(); i++) {
		CVector centerPos = insideCellCenters[i];
		rawData.MXCellCenters.push_back(centerPos);
		centerPos.Print();
	}

	for (uint i = 0; i < outsideBdryNodePos.size(); i++) {
		rawData.bdryNodes.push_back(outsideBdryNodePos[i]);
	}

	generateCellInitNodeInfo_v2(rawData.initCellNodePoss);

	rawData.isStab = true;
	return rawData;
}

RawDataInput_M CellInitHelper::generateRawInput_M() {
	RawDataInput_M rawData;

	rawData.simuType = simuType;
	vector<CVector> insideCellCenters;
	vector<CVector> outsideBdryNodePos;
	std::string bdryInputFileName = globalConfigVars.getConfigValue(
			"Bdry_InputFileName").toString();

	std::string CellCentersFileName = globalConfigVars.getConfigValue(
			"CellCenters_FileName").toString() ;
         //Ali 
        ForReadingData_M2 ForReadingData2 = ReadFile_M2(CellCentersFileName);
        GEOMETRY::Point2D Point2D1[ForReadingData2.CellNumber];
        //Ali 

	GEOMETRY::MeshGen meshGen;

	GEOMETRY::UnstructMesh2D mesh = meshGen.generateMesh2DFromFile(
			bdryInputFileName);


         //Ali
        std::vector<GEOMETRY::Point2D> insideCenterCenters ; 
        for (int ii = 0; ii <ForReadingData2.CellNumber; ii = ii + 1) {
		
		Point2D1[ii].Assign_M2(ForReadingData2.TempSX[ii], ForReadingData2.TempSY[ii]);
		cout << "x coordinate=" << Point2D1[ii].getX() << "y coordinate=" << Point2D1[ii].getY() << "Is on Boundary=" << Point2D1[ii].isIsOnBdry() << endl;
	insideCenterCenters.push_back(Point2D1[ii]); 
	}
         
        //Ali 


         //Ali comment
//	std::vector<GEOMETRY::Point2D> insideCenterCenters =
//			mesh.getAllInsidePoints();

       //Ali comment

	uint initCellCt = insideCenterCenters.size();

	for (uint i = 0; i < initCellCt; i++) {
		insideCellCenters.push_back(
				CVector(insideCenterCenters[i].getX(),
						insideCenterCenters[i].getY(), 0));
	}



	double randNum;
	double progDivStart =
			globalConfigVars.getConfigValue("GrowthPrgrCriVal").toDouble();
	for (uint i = 0; i < initCellCt; i++) {
		randNum = (double) rand() / ((double) RAND_MAX + 1) * progDivStart;
		//std::cout << "rand init growth progress = " << randNum << std::endl;
//Ali to make the initial progree of all nodes zero

 
		rawData.cellGrowProgVec.push_back(randNum);
		//rawData.cellGrowProgVec.push_back(0.75);
		ECellType eCellTypeTmp=ForReadingData2.eCellTypeV.at(i);  
		rawData.cellsTypeCPU.push_back(eCellTypeTmp);
	}
	for (uint i = 0; i < initCellCt; i++) {
    	cout << "second check for cell type" <<rawData.cellsTypeCPU.at(i) << endl ; 
	}

	std::cout << "Printing initial cell center positions ......" << std::endl;
	for (unsigned int i = 0; i < insideCellCenters.size(); i++) {
		CVector centerPos = insideCellCenters[i];
		rawData.initCellCenters.push_back(centerPos);
		std::cout << "    ";
		centerPos.Print();
	}

	generateCellInitNodeInfo_v3(rawData.initCellCenters,
			rawData.cellGrowProgVec, rawData.initMembrNodePoss,
			rawData.initIntnlNodePoss);

	//std::cout << "finished generate raw data" << std::endl;
	//std::cout.flush();

	rawData.isStab = true;
	return rawData;
}

void CellInitHelper::generateRandomAngles(vector<double> &randomAngles,
		int initProfileNodeSize) {
	static const double PI = acos(-1.0);
	randomAngles.clear();
	for (int i = 0; i < initProfileNodeSize; i++) {
		double randomNum = rand() / ((double) RAND_MAX + 1);
		randomAngles.push_back(randomNum * 2.0 * PI);
	}
}

/**
 * Initially, ECM nodes are alignd vertically.
 */
void CellInitHelper::generateECMInitNodeInfo(vector<CVector> &initECMNodePoss,
		int initNodeCountPerECM) {
	initECMNodePoss.clear();
	double ECMInitNodeInterval = globalConfigVars.getConfigValue(
			"ECM_Init_Node_Interval").toDouble();
	int numOfSegments = initNodeCountPerECM - 1;
//double totalLength = ECMInitNodeInterval * numOfSegments;
	if (numOfSegments % 2 == 0) {
		CVector initPt = CVector(0, 0, 0);
		initECMNodePoss.push_back(initPt);
		for (int i = 1; i <= numOfSegments / 2; i++) {
			CVector posSide = initPt + CVector(0, i * ECMInitNodeInterval, 0);
			CVector negSide = initPt - CVector(0, i * ECMInitNodeInterval, 0);
			initECMNodePoss.push_back(posSide);
			initECMNodePoss.push_back(negSide);
		}
	} else {
		CVector initPosPt = CVector(0, ECMInitNodeInterval / 2.0, 0);
		CVector initNegPt = CVector(0, -ECMInitNodeInterval / 2.0, 0);
		initECMNodePoss.push_back(initPosPt);
		initECMNodePoss.push_back(initNegPt);
		for (int i = 1; i <= numOfSegments / 2; i++) {
			CVector posSide = initPosPt
					+ CVector(0, i * ECMInitNodeInterval, 0);
			CVector negSide = initNegPt
					- CVector(0, i * ECMInitNodeInterval, 0);
			initECMNodePoss.push_back(posSide);
			initECMNodePoss.push_back(negSide);
		}
	}
}

void CellInitHelper::generateECMCenters(vector<CVector> &ECMCenters,
		vector<CVector> &CellCenters, vector<CVector> &bdryNodes) {
	ECMCenters.clear();
	vector<double> angles;
	vector<CVector> vecs;
	static const double PI = acos(-1.0);

	const int numberOfECMAroundCellCenter = globalConfigVars.getConfigValue(
			"ECM_Around_Cell_Center").toInt();
	const double distFromCellCenter = globalConfigVars.getConfigValue(
			"Dist_From_Cell_Center").toDouble();
	double unitAngle = 2 * PI / numberOfECMAroundCellCenter;
	for (int i = 0; i < numberOfECMAroundCellCenter; i++) {
		angles.push_back(i * unitAngle);
		CVector vec(sin(i * unitAngle), cos(i * unitAngle), 0);
		vec = vec * distFromCellCenter;
		vecs.push_back(vec);
	}
	for (uint i = 0; i < CellCenters.size(); i++) {
		for (int j = 0; j < numberOfECMAroundCellCenter; j++) {
			CVector pos = CellCenters[i] + vecs[j];
			if (anyCellCenterTooClose(CellCenters, pos)) {
				continue;
			}
			if (anyECMCenterTooClose(ECMCenters, pos)) {
				continue;
			}
			if (anyBoundaryNodeTooClose(bdryNodes, pos)) {
				continue;
			}
			ECMCenters.push_back(pos);
			if (std::isnan(pos.GetX())) {
				throw SceException("number is NAN!", InputInitException);
			}
		}
	}
}

bool CellInitHelper::anyCellCenterTooClose(vector<CVector> &cellCenters,
		CVector position) {
	double MinDistToOtherCellCenters = globalConfigVars.getConfigValue(
			"MinDistToCellCenter").toDouble();
	int size = cellCenters.size();
	for (int i = 0; i < size; i++) {
		if (Modul(cellCenters[i] - position) < MinDistToOtherCellCenters) {
			return true;
		}
	}
	return false;
}

bool CellInitHelper::anyECMCenterTooClose(vector<CVector> &ecmCenters,
		CVector position) {
	double MinDistToOtherECMCenters = globalConfigVars.getConfigValue(
			"MinDistToECMCenter").toDouble();
	int size = ecmCenters.size();
	for (int i = 0; i < size; i++) {
		if (Modul(ecmCenters[i] - position) < MinDistToOtherECMCenters) {
			return true;
		}
	}
	return false;
}

bool CellInitHelper::anyBoundaryNodeTooClose(vector<CVector> &bdryNodes,
		CVector position) {
	double MinDistToOtherBdryNodes = globalConfigVars.getConfigValue(
			"MinDistToBdryNodes").toDouble();
	int size = bdryNodes.size();
	for (int i = 0; i < size; i++) {
		if (Modul(bdryNodes[i] - position) < MinDistToOtherBdryNodes) {
			return true;
		}
	}
	return false;
}

CellInitHelper::~CellInitHelper() {
}

void CellInitHelper::generateCellInitNodeInfo_v2(vector<CVector>& initPos) {
	initPos = generateInitCellNodes();
}

void CellInitHelper::generateCellInitNodeInfo_v3(vector<CVector>& initCenters,
		vector<double>& initGrowProg, vector<vector<CVector> >& initMembrPos,
		vector<vector<CVector> >& initIntnlPos) {
	assert(initCenters.size() == initGrowProg.size());
	vector<CVector> initMembrPosTmp;
	vector<CVector> initIntnlPosTmp;
	for (uint i = 0; i < initCenters.size(); i++) {
		initMembrPosTmp = generateInitMembrNodes(initCenters[i],
				initGrowProg[i]);
		initIntnlPosTmp = generateInitIntnlNodes(initCenters[i],
				initGrowProg[i]);
		initMembrPos.push_back(initMembrPosTmp);
		initIntnlPos.push_back(initIntnlPosTmp);
	}
}

double CellInitHelper::getRandomNum(double min, double max) {
	double rand01 = rand() / ((double) RAND_MAX + 1);
	double randNum = min + (rand01 * (max - min));
	return randNum;
}

vector<CVector> CellInitHelper::generateInitCellNodes() {
	bool isSuccess = false;
	vector<CVector> attemptedPoss;
	while (!isSuccess) {
		attemptedPoss = tryGenInitCellNodes();
		if (isPosQualify(attemptedPoss)) {
			isSuccess = true;
		}
	}
	// also need to make sure center point is (0,0,0).
	CVector tmpSum(0, 0, 0);
	for (uint i = 0; i < attemptedPoss.size(); i++) {
		tmpSum = tmpSum + attemptedPoss[i];
	}
	tmpSum = tmpSum / (double) (attemptedPoss.size());
	for (uint i = 0; i < attemptedPoss.size(); i++) {
		attemptedPoss[i] = attemptedPoss[i] - tmpSum;
	}
	return attemptedPoss;
}

vector<CVector> CellInitHelper::generateInitIntnlNodes(CVector& center,
		double initProg) {
	bool isSuccess = false;

	uint minInitNodeCount =
			globalConfigVars.getConfigValue("InitCellNodeCount").toInt();
	uint maxInitNodeCount = globalConfigVars.getConfigValue(
			"MaxIntnlNodeCountPerCell").toInt();
//Ali

//	uint initIntnlNodeCt = minInitNodeCount ; 
//Ali
//Ali comment
	uint initIntnlNodeCt = minInitNodeCount
			+ (maxInitNodeCount - minInitNodeCount) * initProg;

	vector<CVector> attemptedPoss;
	while (!isSuccess) {
		attemptedPoss = tryGenInitCellNodes(initIntnlNodeCt);
		if (isPosQualify(attemptedPoss)) {
			isSuccess = true;
		}
	}
	/*
	 // also need to make sure center point is (0,0,0).
	 CVector tmpSum(0, 0, 0);
	 for (uint i = 0; i < attemptedPoss.size(); i++) {
	 tmpSum = tmpSum + attemptedPoss[i];
	 }
	 tmpSum = tmpSum / (double) (attemptedPoss.size());
	 for (uint i = 0; i < attemptedPoss.size(); i++) {
	 attemptedPoss[i] = attemptedPoss[i] - tmpSum;
	 }
	 */
	for (uint i = 0; i < attemptedPoss.size(); i++) {
		attemptedPoss[i] = attemptedPoss[i] + center;
	}
	return attemptedPoss;
}

vector<CVector> CellInitHelper::generateInitMembrNodes(CVector& center,
		double initProg) {
	double initRadius =
			globalConfigVars.getConfigValue("InitMembrRadius").toDouble();
	uint initMembrNodeCount = globalConfigVars.getConfigValue(
			"InitMembrNodeCount").toInt();
	vector<CVector> initMembrNodes;
	double unitAngle = 2 * acos(-1.0) / (double) (initMembrNodeCount);
	for (uint i = 0; i < initMembrNodeCount; i++) {
		CVector node;
		node.x = initRadius * cos(unitAngle * i) + center.x;
		node.y = initRadius * sin(unitAngle * i) + center.y;
		initMembrNodes.push_back(node);
	}
	return initMembrNodes;
}

vector<CVector> CellInitHelper::tryGenInitCellNodes() {
	double radius =
			globalConfigVars.getConfigValue("InitCellRadius").toDouble();
	//int initCellNodeCount =
	//    globalConfigVars.getConfigValue("InitCellNodeCount").toInt();
	// now we need non-uniform initial growth progress.
	int initCellNodeCount =
			globalConfigVars.getConfigValue("InitCellNodeCount").toInt();
	vector<CVector> poss;
	int foundCount = 0;
	double randX, randY;
	while (foundCount < initCellNodeCount) {
		bool isInCircle = false;
               //Ali
               cout << "I am in the wrong one" << endl ; 
               //Ali
		while (!isInCircle) {
			randX = getRandomNum(-radius, radius);
			randY = getRandomNum(-radius, radius);
			isInCircle = (sqrt(randX * randX + randY * randY) < radius);
		}
		poss.push_back(CVector(randX, randY, 0));
		foundCount++;
	}
	return poss;
}

vector<CVector> CellInitHelper::tryGenInitCellNodes(uint initNodeCt) {
	double radius =
			globalConfigVars.getConfigValue("InitCellRadius").toDouble();
	vector<CVector> poss;
	uint foundCount = 0;
	double randX, randY;

               //Ali
            //   cout << "I am in the right one" << endl ; 
             //  cout << "# of internal Nodes" << initNodeCt <<endl ; 
               //Ali
	while (foundCount < initNodeCt) {
		bool isInCircle = false;
		//while (!isInCircle) {
			randX = getRandomNum(-radius, radius);
			randY = getRandomNum(-radius, radius);
			isInCircle = (sqrt(randX * randX + randY * randY) < radius);
	//	}
                //Ali
                 if (isInCircle) {
                //Ali
		 poss.push_back(CVector(randX, randY, 0));
		 foundCount++;
         //      cout << "#internal nodes location" << foundCount<<"isInCircle"<<isInCircle <<endl ; 
               //Ali
                 }
               //Ali 
	}
	return poss;
}

bool CellInitHelper::isPosQualify(vector<CVector>& poss) {
	double minInitDist = globalConfigVars.getConfigValue(
			"MinInitDistToOtherNodes").toDouble();
	bool isQualify = true;
	for (uint i = 0; i < poss.size(); i++) {
		for (uint j = 0; j < poss.size(); j++) {
			if (i == j) {
				continue;
			} else {
				CVector distDir = poss[i] - poss[j];
				if (distDir.getModul() < minInitDist) {
					isQualify = false;
					return isQualify;
				}
			}
		}
	}
	return isQualify;
}

bool CellInitHelper::isMXType(CVector position) {
	if (simuType != Beak) {
		return true;
	}
	if (position.y >= internalBdryPts[0].y) {
		return false;
	}
	if (position.x >= internalBdryPts[2].x) {
		return false;
	}
	for (uint i = 0; i < internalBdryPts.size() - 1; i++) {
		CVector a = internalBdryPts[i + 1] - internalBdryPts[i];
		CVector b = position - internalBdryPts[i];
		CVector crossProduct = Cross(a, b);
		if (crossProduct.z > 0) {
			return false;
		}
	}
	return true;
}

void CellInitHelper::initInternalBdry() {
	GEOMETRY::MeshGen meshGen;
	GEOMETRY::MeshInput input = meshGen.obtainMeshInput();
	internalBdryPts = input.internalBdryPts;
}

SimulationInitData_V2 CellInitHelper::initStabInput() {
	RawDataInput rawInput = generateRawInput_stab();
	SimulationInitData_V2 initData = initInputsV3(rawInput);
	initData.isStab = true;
	return initData;
}

//RawDataInput rawInput = generateRawInput_stab();
SimulationInitData_V2_M CellInitHelper::initInput_M() {
	RawDataInput_M rawInput_m = generateRawInput_M();
	SimulationInitData_V2_M initData = initInputsV3_M(rawInput_m);
	initData.isStab = false;
	return initData;
}

SimulationInitData_V2 CellInitHelper::initSimuInput(
		std::vector<CVector> &cellCenterPoss) {
	RawDataInput rawInput = generateRawInput_simu(cellCenterPoss);
	SimulationInitData_V2 simuInitData = initInputsV3(rawInput);
	simuInitData.isStab = false;
	return simuInitData;
}

void SimulationGlobalParameter::initFromConfig() {
	int type = globalConfigVars.getConfigValue("SimulationType").toInt();
	SimulationType simuType = parseTypeFromConfig(type);

	animationNameBase =
			globalConfigVars.getConfigValue("AnimationFolder").toString()
					+ globalConfigVars.getConfigValue("AnimationName").toString()
					+ globalConfigVars.getConfigValue("UniqueSymbol").toString();
        //A & A 
	InitTimeStage=
			globalConfigVars.getConfigValue("InitTimeStage").toDouble();
        //A & A 
	totalSimuTime =
			globalConfigVars.getConfigValue("SimulationTotalTime").toDouble();

	dt = globalConfigVars.getConfigValue("SimulationTimeStep").toDouble();
	Damp_Coef= globalConfigVars.getConfigValue("DampingCoef").toDouble();

	totalTimeSteps = totalSimuTime / dt;

	totalFrameCount =
			globalConfigVars.getConfigValue("TotalNumOfOutputFrames").toInt();

	aniAuxVar = totalTimeSteps / totalFrameCount;

	aniCri.pairDisplayDist = globalConfigVars.getConfigValue(
			"IntraLinkDisplayRange").toDouble();

	aniCri.animationType = parseAniTpFromConfig(
			globalConfigVars.getConfigValue("AnimationType").toInt());

	aniCri.threshold = globalConfigVars.getConfigValue("DeltaValue").toDouble();
	if (simuType != Disc_M) {
		aniCri.arrowLength = globalConfigVars.getConfigValue(
				"DisplayArrowLength").toDouble();
	}

	if (simuType != Beak && simuType != Disc_M) {
		dataOutput =
				globalConfigVars.getConfigValue("PolygonStatFileName").toString()
						+ globalConfigVars.getConfigValue("UniqueSymbol").toString()
						+ ".txt";
		imgOutput =
				globalConfigVars.getConfigValue("DataOutputFolder").toString()
						+ globalConfigVars.getConfigValue("ImgSubFolder").toString()
						+ globalConfigVars.getConfigValue("ImgFileNameBase").toString();
		dataFolder = globalConfigVars.getConfigValue("DataFolder").toString();
		dataName = dataFolder
				+ globalConfigVars.getConfigValue("DataName").toString();
	}

}

RawDataInput CellInitHelper::generateRawInput_singleCell() {
	RawDataInput rawData;
	rawData.simuType = simuType;

	std::string initPosFileName = globalConfigVars.getConfigValue(
			"SingleCellCenterPos").toString();

	fstream fs(initPosFileName.c_str());
	vector<CVector> insideCellCenters = GEOMETRY::MeshInputReader::readPointVec(
			fs);
	fs.close();

	for (unsigned int i = 0; i < insideCellCenters.size(); i++) {
		CVector centerPos = insideCellCenters[i];
		rawData.MXCellCenters.push_back(centerPos);
	}

	generateCellInitNodeInfo_v2(rawData.initCellNodePoss);
	rawData.isStab = true;
	return rawData;
}

SimulationInitData_V2 CellInitHelper::initSingleCellTest() {
	RawDataInput rawInput = generateRawInput_singleCell();
	SimulationInitData_V2 initData = initInputsV3(rawInput);
	initData.isStab = true;
	return initData;
}
