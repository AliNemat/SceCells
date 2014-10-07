#include <iostream>
#include "gtest/gtest.h"
#include <vector>
#include <algorithm>
#include "commonData.h"
#include "ConfigParser.h"
#include "ResAnalysisHelper.h"

using namespace std;

extern GlobalConfigVars globalConfigVars;

class ResAnalysisUtilTest: public ::testing::Test {
protected:
	virtual void SetUp() {
		ConfigParser parser;
		std::string configFileName = "./resources/unitTest.cfg";
		globalConfigVars = parser.parseConfigFile(configFileName);
	}
};

TEST_F(ResAnalysisUtilTest, emptyTest) {
	ResAnalysisHelper resHelper;
	vector<NodeWithLabel> nodeWithLabels;
	vector<vector<int> > result = resHelper.outputLabelMatrix(nodeWithLabels);
	EXPECT_EQ(result.size(), (uint )100);
	for (uint i = 0; i < 100; i++) {
		EXPECT_EQ(result[i].size(), (uint )200);
		for (uint j = 0; j < 200; j++) {
			EXPECT_EQ(result[i][j], -1);
		}
	}
}

TEST_F(ResAnalysisUtilTest, singleTest) {
	ResAnalysisHelper resHelper;
	vector<NodeWithLabel> nodeWithLabels;
	NodeWithLabel nodeLabel;
	nodeLabel.cellRank = 2;
	nodeLabel.position = CVector(26, 17.5, 0);
	nodeWithLabels.push_back(nodeLabel);

	vector<vector<int> > result = resHelper.outputLabelMatrix(nodeWithLabels);
	EXPECT_EQ(result.size(), (uint )100);
	for (uint i = 0; i < 100; i++) {
		EXPECT_EQ(result[i].size(), (uint )200);
		for (uint j = 0; j < 200; j++) {
			double yCoord = i * 0.25 + 5.0 + 0.125;
			double xCoord = j * 0.25 + 1.0 + 0.125;
			double xDiff = xCoord - nodeLabel.position.GetX();
			double yDiff = yCoord - nodeLabel.position.GetY();
			if (sqrt(xDiff * xDiff + yDiff * yDiff) < 0.1) {
				EXPECT_EQ(2, result[i][j]);
			} else {
				EXPECT_EQ(-1, result[i][j]);
			}
		}
	}
}

TEST_F(ResAnalysisUtilTest, realTest) {
	ResAnalysisHelper resHelper;
	vector<NodeWithLabel> nodeWithLabels;
	NodeWithLabel nodeLabel1;
	nodeLabel1.cellRank = 2;
	nodeLabel1.position = CVector(26.125, 17.5, 0);
	nodeWithLabels.push_back(nodeLabel1);
	NodeWithLabel nodeLabel2;
	nodeLabel2.cellRank = 7;
	nodeLabel2.position = CVector(26, 18.125, 0);
	nodeWithLabels.push_back(nodeLabel2);

	vector<vector<int> > result = resHelper.outputLabelMatrix(nodeWithLabels);
	EXPECT_EQ(result.size(), (uint )100);
	for (uint i = 0; i < 100; i++) {
		EXPECT_EQ(result[i].size(), (uint )200);
		for (uint j = 0; j < 200; j++) {
			double yCoord = i * 0.25 + 5.0 + 0.125;
			double xCoord = j * 0.25 + 1.0 + 0.125;

			double xDiff1 = xCoord - nodeLabel1.position.GetX();
			double yDiff1 = yCoord - nodeLabel1.position.GetY();
			double xDiff2 = xCoord - nodeLabel2.position.GetX();
			double yDiff2 = yCoord - nodeLabel2.position.GetY();

			double dist1 = sqrt(xDiff1 * xDiff1 + yDiff1 * yDiff1);
			double dist2 = sqrt(xDiff2 * xDiff2 + yDiff2 * yDiff2);

			if (dist1 > 0.5) {
				if (dist2 > 0.1) {
					EXPECT_EQ(-1, result[i][j]);
				} else {
					EXPECT_EQ(7, result[i][j]);
				}
			} else {
				if (dist2 > 0.5) {
					if (dist1 < 0.1) {
						EXPECT_EQ(2, result[i][j]);
					} else {
						EXPECT_EQ(-1, result[i][j]);
					}
				} else {
					if (dist1 < dist2) {
						EXPECT_EQ(2, result[i][j]);
					} else if (dist1 > dist2) {
						EXPECT_EQ(7, result[i][j]);
					}
				}
			}
		}
	}
}

