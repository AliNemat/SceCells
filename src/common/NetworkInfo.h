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
#include <sstream>

class RankWithCor {
public:
	int nodeRank;
	CVector pos;
	double angle;
	double growP;
};

struct AngleSortComp {
	bool operator()(const RankWithCor& n1, const RankWithCor& n2) {
		return (n1.angle < n2.angle);
	}
};

class NetworkNode {
	int _nodeRank;
	CVector _pos;
	std::vector<int> _ngbrList;
	double _growP;
public:
	NetworkNode();
	NetworkNode(int nodeRank, CVector pos, std::vector<int>& ngbrList,
			double growP);
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

	double getGrowP() const {
		return _growP;
	}

	void setGrowP(double growP) {
		_growP = growP;
	}
};

/**
 * Undirected edge
 */
class NetworkEdge {
	int _n1;
	int _n2;
public:
	NetworkEdge(int n1, int n2);
	bool isMatch(int n1, int n2);
};

class PreT1State {
public:
	int nodeRank;
	int centerNgbr;
	vector<int> sideNgbrs;
	double gp1, gp2, gp3, gp4;
	PreT1State();
};

class NetworkInfo {
	std::vector<NetworkNode> networkNodes;
	//std::vector<NetworkEdge> networkEdges;
	std::set<std::string> edgeSet;
public:
	NetworkInfo();
	NetworkInfo(std::vector<NetworkNode> &nodeList);
	bool isEdgePresent(int n1, int n2);
	bool isT1Tran(PreT1State &preT1);
	vector<PreT1State> scanForPreT1States();
	virtual ~NetworkInfo();
};

void angleSortVec(std::vector<RankWithCor> &ngbrListWithCor);
vector<PreT1State> extractPreT1States(std::vector<RankWithCor> &ngbrListWithCor,
		int nodeRank, double nodePg, NetworkInfo* netInfo);

#endif /* NETWORKINFO_H_ */
