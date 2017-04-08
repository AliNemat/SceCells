#include "commonData.h"

const std::string SceInputPoint::delimiter = " ";

std::string toString(SceExceptionType type) {
	std::string result("Undefined type");
	switch (type) {
	case BaseException:
		result = "Base exception (usage discouraged)";
		break;
	case InputInitException:
		result = "Initialization Exception";
		break;
	case ConfigFileNotFound:
		result = "Configuration file not found Exception";
		break;
	case ConfigValueException:
		result = "Configuration value Exception";
		break;
	case OutputAnalysisDataException:
		result = "Exception while outputting analysis result";
		break;
	case FileIOException:
		result = "Exception while file IO";
		break;
	case MemoryInvalidAccess:
		result = "Exception while executing memory operation";
		break;
	case InvalidInput:
		result = "Exception processing method input";
		break;
	case AlgorithmBug:
		result = "Possible bug in algorithm used";
		break;
	}
	return result;
}

std::string toString(SceNodeType type) {
	std::string result("Undefined type");
	switch (type) {
	case Boundary:
		result = "Boundary";
		break;
	case Profile:
		result = "Profile";
		break;
	case ECM:
		result = "ECM";
		break;
	case FNM:
		result = "FNM";
		break;
	case MX:
		result = "MX";
		break;
	case Cart:
		result = "Cartilage";
		break;
	case Base:
		result = "Base";
		break;
	case CellIntnl:
		result = "Base";
		break;
	case CellMembr:
		result = "Base";
		break;
	}
	return result;
}

double nodeTypeToScale(SceNodeType type) {
	double result = 0.0;
	switch (type) {
	case Boundary:
		result = 1;
		break;
	case Profile:
		result = 7;
		break;
	case ECM:
		result = 3;
		break;
	case FNM:
		result = 4;
		break;
	case MX:
		result = 5;
		break;
	case Cart:
		result = 2;
		break;
	case Base:
		result = 6;
		break;
	case CellIntnl:
		result = 7;
		break;
	case CellMembr:
		result = 8;
		break;
	}
	return result;
}

void inputInitialData::addNewPoints(std::vector<SceInputPoint>& newPoints) {
	for (unsigned int i = 0; i < newPoints.size(); i++) {
		inputPoints.push_back(newPoints[i]);
	}
}

bool AnimationCriteria::isPairQualify(uint seq1, uint seq2, double x1,
		double y1, double z1, SceNodeType t1, uint r1, double x2, double y2,
		double z2, SceNodeType t2, uint r2) {
	bool condi1 = false, condi2 = false;
	if (t1 == t2 && r1 == r2) {
		if (t1 == Boundary || t1 == ECM) {
			if (abs((int) seq1 - (int) seq2) == 1) {
				condi1 = true;
			}
		} else if (t1 == MX || t1 == FNM) {
			condi1 = true;
		}
		if (condi1) {
			double dist = compuDistHost(x1, y1, z1, x2, y2, z2);
			if (dist < pairDisplayDist) {
				condi2 = true;
			}
		}
	}
	return condi1 && condi2;
}

bool AnimationCriteria::isPairQualify_M(double x1, double y1, double x2,
		double y2) {
	double dummy = 0;
	double dist = compuDistHost(x1, y1, dummy, x2, y2, dummy);
	if (dist < pairDisplayDist) {
		return true;
	} else {
		return false;
	}
}

double compuDistHost(double &xPos, double &yPos, double &zPos, double &xPos2,
		double &yPos2, double &zPos2) {
	return sqrt(
			(xPos - xPos2) * (xPos - xPos2) + (yPos - yPos2) * (yPos - yPos2)
					+ (zPos - zPos2) * (zPos - zPos2));
}

void VtkAnimationData::outputVtkAni(std::string scriptNameBase, int rank) {
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
	fs << "SCALARS curvature float" << endl;
	fs << "LOOKUP_TABLE default" << endl;

	for (uint i = 0; i < pointsAniData.size(); i++) {
		std::cout << "****************      SIZE IS:     " << pointsAniData.size() << std::endl;
//		if (pointsAniData[i].colorScale2 < 0.000001 || pointsAniData[i].colorScale2>1.0){
//			pointsAniData[i].colorScale2 = 0.0;}
		fs << pointsAniData[i].colorScale2 << endl;
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

std::vector<double> getArrayXComp(std::vector<CVector>& nodePosVec) {
	std::vector<double> result;
	for (uint i = 0; i < nodePosVec.size(); i++) {
		result.push_back(nodePosVec[i].GetX());
	}
	return result;
}

std::vector<double> getArrayYComp(std::vector<CVector>& nodePosVec) {
	std::vector<double> result;
	for (uint i = 0; i < nodePosVec.size(); i++) {
		result.push_back(nodePosVec[i].GetY());
	}
	return result;
}

std::vector<double> getArrayZComp(std::vector<CVector>& nodePosVec) {
	std::vector<double> result;
	for (uint i = 0; i < nodePosVec.size(); i++) {
		result.push_back(nodePosVec[i].GetZ());
	}
	return result;
}

SimulationType parseTypeFromConfig(int configValue) {
	if (configValue == 0) {
		return Beak;
	} else if (configValue == 1) {
		return Disc;
	} else if (configValue == 2) {
		return SingleCellTest;
	} else if (configValue == 3) {
		return Disc_M;
	} else {
		throw SceException("Simulation type in config file is not defined",
				ConfigValueException);
	}
}

bool valueToType(int value) {
	if (value != 0) {
		return true;
	} else {
		return false;
	}
}

AniType parseAniTpFromConfig(int configValue) {
	if (configValue == 0) {
		return CellType;
	} else if (configValue == 1) {
		return ForceAbsVal;
	} else if (configValue == 2) {
		return Force;
	} else if (configValue == 3) {
		return Tension;
	} else if (configValue == 4) {
		return T1Tran;
	} else if (configValue == 5) {
		return PolySide;
	} else {
		throw SceException("Animation type in config file is not defined",
				ConfigValueException);
	}
}

uint findClosestArrIndexGivenPos(std::vector<CVector>& vecArr, CVector& pos) {
	if (vecArr.size() == 0) {
		std::string errorMsg =
				"while finding closest array index, the input array has zero elements";
		throw SceException(errorMsg, InvalidInput);
	} else {
		uint index = 0;
		CVector tmpVec = vecArr[0] - pos;
		double minDis = tmpVec.getModul();
		for (uint i = 1; i < vecArr.size(); i++) {
			tmpVec = vecArr[i] - pos;
			double dist = tmpVec.getModul();
			if (dist < minDis) {
				index = i;
				minDis = dist;
			}
		}
		return index;
	}
}

void AblationEvent::printInfo() {
	std::cout << "ablation experiment is scheduled on timestep " << timeStep
			<< ", " << ablationCells.size() << " cells will be cut" << endl;
	std::cout << "Detailed info: " << std::endl;
	for (uint i = 0; i < ablationCells.size(); i++) {
		std::cout << "cell rank: " << ablationCells[i].cellNum
				<< ", node ranks: [ ";
		for (uint j = 0; j < ablationCells[i].nodeNums.size(); j++) {
			std::cout << ablationCells[i].nodeNums[j] << " ";
		}
		std::cout << "]" << std::endl;
	}
}

AblationEvent readAblationEvent(std::string inputName) {
	AblationEvent ablaEvent;
	fstream fs(inputName.c_str());
	char specialChar;
	fs >> specialChar;
	assert(specialChar == '#');

	int timestep;
	fs >> timestep;
	ablaEvent.timeStep = timestep;

	fs >> specialChar;
	assert(specialChar == '{');

	fs >> specialChar;
	assert(specialChar == '$');

	int numOfCells;
	fs >> numOfCells;
	assert(numOfCells > 0);

	fs >> specialChar;
	assert(specialChar == '{');

	for (int i = 0; i < numOfCells; i++) {
		AblaInfo info;
		fs >> specialChar;
		assert(specialChar == '%');
		int cellRank;
		fs >> cellRank;
		info.cellNum = cellRank;
		fs >> specialChar;
		assert(specialChar == '{');
		int nodeCount;
		fs >> nodeCount;
		fs >> specialChar;
		assert(specialChar == '[');
		for (int j = 0; j < nodeCount; j++) {
			int nodeRank;
			fs >> nodeRank;
			info.nodeNums.push_back(nodeRank);
		}
		fs >> specialChar;
		assert(specialChar == ']');
		fs >> specialChar;
		assert(specialChar == '}');
		ablaEvent.ablationCells.push_back(info);
	}

	fs >> specialChar;
	assert(specialChar == '}');
	fs >> specialChar;
	assert(specialChar == '}');

	ablaEvent.printInfo();
	return ablaEvent;
}

std::vector<CVector> obtainPtsBetween(CVector& start, CVector& end,
		double& spacing, uint maxNodeCount) {
	std::vector<CVector> result;
	double spacingNew = spacing;
	CVector ptVec = end - start;
	CVector unitVec = ptVec.getUnitVector();
	double edgeLength = ptVec.getModul();
	uint numNewMemNodes = edgeLength / spacing;
	if (edgeLength - numNewMemNodes * spacing < 0.5 * spacing) {
		numNewMemNodes--;
	}
	numNewMemNodes--;
	if (numNewMemNodes > maxNodeCount) {
		spacingNew = edgeLength / (maxNodeCount + 1);
		numNewMemNodes = maxNodeCount;
	}
	CVector nodeOld = start, nodeNew;
	for (uint j = 0; j < numNewMemNodes; j++) {
		nodeNew = nodeOld + spacingNew * unitVec;
		result.push_back(nodeNew);
		nodeOld = nodeNew;
	}
	return result;
}

void CellsStatsData::printPolyCountToFile(std::string fileName,
		double divThreshold) {
	//std::remove(fileName.c_str());
	std::map<uint, uint> countBdry, countNormal, countDiv;
	for (uint i = 0; i < cellsStats.size(); i++) {
		if (cellsStats[i].isBdryCell == true) {
			insertCount(cellsStats[i].numNeighbors, countBdry);
		} else {
			if (cellsStats[i].cellGrowthProgress <= divThreshold) {
				insertCount(cellsStats[i].numNeighbors, countNormal);
			} else {
				insertCount(cellsStats[i].numNeighbors, countDiv);
			}
		}
	}
	printCountsToFile(fileName, countNormal, countDiv, countBdry);
}

void insertCount(uint numNeighbor, std::map<uint, uint>& count) {
	std::map<uint, uint>::iterator it = count.find(numNeighbor);
	if (it == count.end()) {
		count.insert(std::pair<uint, uint>(numNeighbor, 1));
	} else {
		it->second = it->second + 1;
	}
}

void printCountsToFile(std::string fileName, std::map<uint, uint>& countNormal,
		std::map<uint, uint>& countDiv, std::map<uint, uint>& countBdry) {
	ofstream ofs(fileName.c_str(), ios::app);
	std::vector<CountEntry> normalEntries = processCountMap(countNormal);
	printEntriesToFile(ofs, normalEntries);
	ofs << "# ";
	std::vector<CountEntry> divEntries = processCountMap(countDiv);
	printEntriesToFile(ofs, divEntries);
	ofs << "# ";
	std::vector<CountEntry> bdryEntries = processCountMap(countBdry);
	printEntriesToFile(ofs, bdryEntries);
	ofs << std::endl;
	ofs.close();
}

std::vector<CountEntry> processCountMap(std::map<uint, uint>& countMap) {
	std::vector<CountEntry> result;
	std::map<uint, uint>::iterator it;
	for (it = countMap.begin(); it != countMap.end(); ++it) {
		CountEntry tmpEntry;
		tmpEntry.numOfNeighbor = it->first;
		tmpEntry.count = it->second;
		result.push_back(tmpEntry);
	}
	sort(result.begin(), result.end());
	return result;
}

void printEntriesToFile(ofstream& fs, std::vector<CountEntry>& countEntries) {
	for (uint i = 0; i < countEntries.size(); i++) {
		fs << countEntries[i].numOfNeighbor << "," << countEntries[i].count
				<< " ";
	}
}

void CellsStatsData::printDetailStatsToFile(std::string fileNameBase,
		int timestep) {
	std::stringstream ss;
	ss << std::setw(5) << std::setfill('0') << timestep;
	std::string nameRank = ss.str();
	std::string fileName = fileNameBase + nameRank + ".txt";
	ofstream ofs(fileName.c_str(), ios::out);
	for (uint i = 0; i < cellsStats.size(); i++) {
		cellsStats[i].printToFile(ofs);
	}
	ofs.close();
}

void CellStats::printToFile(ofstream& ofs) {
	ofs << "CellRank:" << cellRank << std::endl;
	ofs << "    GrowthProgress:" << cellGrowthProgress << std::endl;
	ofs << "    MembrGrowthProgress:" << membrGrowthProgress << std::endl;
	ofs << "    IsBoundrayCell:" << isBdryCell << std::endl;
	ofs << "    NumOfNeighbors:" << numNeighbors << std::endl;
	ofs << "    CellArea:" << cellArea << std::endl;
	ofs << "    NeighborCells:{ ";
	for (std::set<int>::iterator it = neighborVec.begin();
			it != neighborVec.end(); ++it) {
		ofs << *it << " ";
	}
	ofs << "}" << std::endl;
	ofs << "    CurrentActiveIntnlNode:" << currentActiveIntnlNodes
			<< std::endl;
	ofs << "    CurrentActiveMembrNodes:" << currentActiveMembrNodes
			<< std::endl;
	ofs << "    CellCenter:" << cellCenter << std::endl;
	ofs << std::endl;
}

vector<double> CellsStatsData::outputPolySides() {
	vector<double> result;
	for (uint i = 0; i < cellsStats.size(); i++) {
		result.push_back(cellsStats[i].numNeighbors);
	}
	return result;
}

//Ali
void CellsStatsData::printStressStrain(std::string FileName1,double curTime,double Init_Displace) {
  ofstream ofs1(FileName1.c_str(),ios::app); 
//  double F_Ext=60*SceMechPara_M.F_Ext_Incline*curTime ; 
  ofs1 << curTime<<","<<50*F_Ext_Out/(Cells_Extrem_Loc[3]-Cells_Extrem_Loc[2])<<","
       <<((Cells_Extrem_Loc[1]-Cells_Extrem_Loc[0])-Init_Displace)/Init_Displace<<"," <<(MaxDistanceX-Init_Displace)/Init_Displace<<std::endl ; 

}
void CellsStatsData::printStressStrain_Ini(std::string FileName1) {
  ofstream ofs1(FileName1.c_str(),ios::out); 
  ofs1 << "Time"<<","<<"Stress"<<","<<"Strain_M"<<","<< "Strain_Center"<<std::endl ; 
}

//Ali     

