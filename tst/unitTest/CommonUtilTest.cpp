#include <iostream>
#include "gtest/gtest.h"
#include <vector>
#include <algorithm>
#include "commonData.h"
#include "ConfigParser.h"

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

}

TEST_F(ResAnalysisUtilTest, singleTest) {

}

TEST_F(ResAnalysisUtilTest, realTest) {

}

