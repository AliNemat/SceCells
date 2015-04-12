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
}

NetworkNode::NetworkNode(int nodeRank, CVector pos,
		std::vector<int>& ngbrList) {
	_nodeRank = nodeRank;
	_pos = pos;
	_ngbrList = ngbrList;
}

NetworkNode::~NetworkNode() {
}

NetworkInfo::NetworkInfo() {
	// TODO Auto-generated constructor stub

}

NetworkInfo::NetworkInfo(std::vector<NetworkNode>& nodeList) {
}

NetworkInfo::~NetworkInfo() {
	// TODO Auto-generated destructor stub
}

PreT1State::PreT1State() {
	nodeRank = 0;
	centerNgbr = 0;
	sideNgbrs.resize(2);
}

vector<PreT1State> NetworkInfo::scanForPreT1States() {
	std::vector<PreT1State> result;
	for (uint i = 0; i < networkNodes.size(); i++) {
		int nodeRank = networkNodes[i].getNodeRank();
		std::vector<int> ngbrList = networkNodes[i].getNgbrList();
		CVector nodeCenter = networkNodes[i].getPos();
		std::vector<RankWithCor> ngbrListWithCor;
		RankWithCor tmpData;
		for (uint j = 0; j < ngbrList.size(); j++) {
			CVector tmpPos = networkNodes[ngbrList[j]].getPos();
			tmpData.nodeRank = ngbrList[j];
			tmpData.pos = tmpPos - nodeCenter;
			ngbrListWithCor.push_back(tmpData);
		}
		unifyVec(ngbrListWithCor);
		angleSortVec(ngbrListWithCor);
		vector<PreT1State> tmpPreT1s = extractPreT1States(ngbrListWithCor);
		result.insert(result.end(), tmpPreT1s.begin(), tmpPreT1s.end());
	}
	return result;
}

void unifyVec(std::vector<RankWithCor>& ngbrListWithCor) {
}

void angleSortVec(std::vector<RankWithCor>& ngbrListWithCor) {
}

vector<PreT1State> extractPreT1States(
		std::vector<RankWithCor>& ngbrListWithCor) {
	vector<PreT1State> result;
	return result;
}
