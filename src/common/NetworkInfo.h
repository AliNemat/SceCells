/*
 * NetworkInfo.h
 *
 *  Created on: Apr 11, 2015
 *      Author: wenzhao
 */

#ifndef NETWORKINFO_H_
#define NETWORKINFO_H_

#include <vector>
#include <GeoVector.h>
#include <commonData.h>

class RankWithCor {
public:
	int nodeRank;
	CVector pos;
};

class NetworkNode {
	int _nodeRank;
	CVector _pos;
	std::vector<int> _ngbrList;
public:
	NetworkNode();
	NetworkNode(int nodeRank, CVector pos, std::vector<int>& ngbrList);
	virtual ~NetworkNode();

	const std::vector<int>& getNgbrList() const {
		return _ngbrList;
	}

	void setNgbrList(const vector<int>& ngbrList) {
		_ngbrList = ngbrList;
	}

	int getNodeRank() const {
		return _nodeRank;
	}

	void setNodeRank(int nodeRank) {
		_nodeRank = nodeRank;
	}

	const CVector& getPos() const {
		return _pos;
	}

	void setPos(const CVector& pos) {
		_pos = pos;
	}
};

class PreT1State {
	int nodeRank;
	int centerNgbr;
	vector<int> sideNgbrs;
public:
	PreT1State();
};

class NetworkInfo {
	std::vector<NetworkNode> networkNodes;
public:
	NetworkInfo();
	NetworkInfo(std::vector<NetworkNode> &nodeList);
	vector<PreT1State> scanForPreT1States();
	virtual ~NetworkInfo();
};

void unifyVec(std::vector<RankWithCor> &ngbrListWithCor);
void angleSortVec(std::vector<RankWithCor> &ngbrListWithCor);
vector<PreT1State> extractPreT1States(
		std::vector<RankWithCor> &ngbrListWithCor);

#endif /* NETWORKINFO_H_ */
