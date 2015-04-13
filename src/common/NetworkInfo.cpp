/*
 * NetworkInfo.cpp
 *
 *  Created on: Apr 11, 2015
 *      Author: wenzhao
 */

#include "NetworkInfo.h"

NetworkNode::NetworkNode() {
	_nodeRank = 0;
	_pos = CVector(0);
	_growP = 0;
}

NetworkNode::NetworkNode(int nodeRank, CVector pos, std::vector<int>& ngbrList,
		double growP) {
	_nodeRank = nodeRank;
	_pos = pos;
	_ngbrList = ngbrList;
	_growP = growP;
}

NetworkNode::~NetworkNode() {
}

NetworkEdge::NetworkEdge(int n1, int n2) {
	_n1 = n1;
	_n2 = n2;
}

bool NetworkEdge::isMatch(int n1, int n2) {
	if ((n1 == _n1 && n2 == _n2) || (n1 == _n2 && n2 == _n1)) {
		return true;
	} else {
		return false;
	}
}

NetworkInfo::NetworkInfo() {
// TODO Auto-generated constructor stub

}

NetworkInfo::NetworkInfo(std::vector<NetworkNode>& nodeList) {
	networkNodes = nodeList;
	for (uint i = 0; i < networkNodes.size(); i++) {
		for (uint j = 0; j < networkNodes[i].getNgbrList().size(); j++) {
			int n1 = networkNodes[i].getNodeRank();
			int n2 = networkNodes[i].getNgbrList()[j];
			if (n1 > n2) {
				int tmp = n2;
				n2 = n1;
				n1 = tmp;
			}
			ostringstream ss;
			ss << n1 << "_" << n2;
			std::string str = ss.str();
			edgeSet.insert(str);
		}
	}
}

bool NetworkInfo::isEdgePresent(int n1, int n2) {
	ostringstream ss;
	if (n1 < n2) {
		ss << n1 << "_" << n2;
	} else {
		ss << n2 << "_" << n1;
	}
	std::string str = ss.str();
	if (edgeSet.find(str) == edgeSet.end()) {
		return false;
	} else {
		return true;
	}
}

bool NetworkInfo::isT1Tran(PreT1State& preT1) {
	if (isEdgePresent(preT1.nodeRank, preT1.sideNgbrs[0])
			&& isEdgePresent(preT1.nodeRank, preT1.sideNgbrs[1])
			&& isEdgePresent(preT1.centerNgbr, preT1.sideNgbrs[0])
			&& isEdgePresent(preT1.centerNgbr, preT1.sideNgbrs[1])
			&& isEdgePresent(preT1.sideNgbrs[0], preT1.sideNgbrs[1])
			&& !isEdgePresent(preT1.nodeRank, preT1.centerNgbr)
			&& preT1.gp1 < networkNodes[preT1.nodeRank].getGrowP()
			&& preT1.gp2 < networkNodes[preT1.centerNgbr].getGrowP()
			&& preT1.gp3 < networkNodes[preT1.sideNgbrs[0]].getGrowP()
			&& preT1.gp4 < networkNodes[preT1.sideNgbrs[1]].getGrowP()) {
		return true;
	} else {
		return false;
	}
}

NetworkInfo::~NetworkInfo() {
// TODO Auto-generated destructor stub
}

PreT1State::PreT1State() {
	nodeRank = 0;
	centerNgbr = 0;
	sideNgbrs.resize(2);
	gp1 = 0;
	gp2 = 0;
	gp3 = 0;
	gp4 = 0;
}

vector<PreT1State> NetworkInfo::scanForPreT1States() {
	std::vector<PreT1State> result;
	for (uint i = 0; i < networkNodes.size(); i++) {
		int nodeRank = networkNodes[i].getNodeRank();
		double nodeProg = networkNodes[i].getGrowP();
		std::vector<int> ngbrList = networkNodes[i].getNgbrList();
		CVector nodeCenter = networkNodes[i].getPos();
		std::vector<RankWithCor> ngbrListWithCor;
		RankWithCor tmpData;
		for (uint j = 0; j < ngbrList.size(); j++) {
			CVector tmpPos = networkNodes[ngbrList[j]].getPos();
			tmpData.nodeRank = ngbrList[j];
			tmpData.pos = tmpPos - nodeCenter;
			tmpData.growP = networkNodes[ngbrList[j]].getGrowP();
			ngbrListWithCor.push_back(tmpData);
		}
		angleSortVec(ngbrListWithCor);
		vector<PreT1State> tmpPreT1s = extractPreT1States(ngbrListWithCor,
				nodeRank, nodeProg, this);
		result.insert(result.end(), tmpPreT1s.begin(), tmpPreT1s.end());
	}
	return result;
}

void angleSortVec(std::vector<RankWithCor>& ngbrListWithCor) {
	for (uint i = 0; i < ngbrListWithCor.size(); i++) {
		ngbrListWithCor[i].angle = ngbrListWithCor[i].pos.getAngle2DPlane();
	}
	std::sort(ngbrListWithCor.begin(), ngbrListWithCor.end(), AngleSortComp());
}

vector<PreT1State> extractPreT1States(std::vector<RankWithCor>& ngbrListWithCor,
		int nodeRank, double nodePg, NetworkInfo* netInfo) {
	vector<PreT1State> result;
	PreT1State preT1Tmp;
	for (uint i = 0; i < ngbrListWithCor.size() - 2; i++) {
		int side1Rank = ngbrListWithCor[i].nodeRank;
		int centerRank = ngbrListWithCor[i + 1].nodeRank;
		int side2Rank = ngbrListWithCor[i + 2].nodeRank;
		if (centerRank < nodeRank) {
			continue;
		}
		if (!netInfo->isEdgePresent(side1Rank, side2Rank)
				&& netInfo->isEdgePresent(side1Rank, centerRank)
				&& netInfo->isEdgePresent(side2Rank, centerRank)) {
			preT1Tmp.nodeRank = nodeRank;
			preT1Tmp.gp1 = nodePg;
			preT1Tmp.centerNgbr = centerRank;
			preT1Tmp.gp2 = ngbrListWithCor[i + 1].growP;
			preT1Tmp.sideNgbrs.clear();
			preT1Tmp.sideNgbrs.push_back(side1Rank);
			preT1Tmp.gp3 = ngbrListWithCor[i].growP;
			preT1Tmp.sideNgbrs.push_back(side2Rank);
			preT1Tmp.gp4 = ngbrListWithCor[i + 2].growP;
			result.push_back(preT1Tmp);
		}
	}
	return result;
}
