#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include <tr1/unordered_map>
#include "commonData.h"

#ifndef CONFIGPARSER_H_
#define CONFIGPARSER_H_

/**
 * Config value represented by string but could be translated to other types.
 */
class ConfigVarValue {
private:
	std::string varibleValue;
public:
	void setValue(std::string &value);
	std::string toString() const;
	int toInt() const;
	double toDouble() const;
};

/**
 * A config variable is represented by name-value pair.
 */
class ConfigVar {
private:
	std::string varibleName;
	ConfigVarValue varibleValue;
public:
	ConfigVar(std::string varName);
	ConfigVarValue getValue() const {
		return varibleValue;
	}
	void setValue(std::string &value) {
		varibleValue.setValue(value);
	}
	std::string getVarName() const {
		return varibleName;
	}
};

/**
 * Global configuration variables.
 */
class GlobalConfigVars {
	std::vector<ConfigVar> configVars;
	std::tr1::unordered_map<std::string, ConfigVar> configMap;
public:
	std::vector<ConfigVar> getConfigVars() const {
		return configVars;
	}
	GlobalConfigVars();
	void insertData(std::string varName, std::string varValue);
	void updateData(std::string varName, std::string varValue);
	ConfigVarValue getConfigValue(std::string varName);
	void printAll();
};

/**
 * Parser for configuration file.
 */
class ConfigParser {
	std::string removeLeadingAndTrailingSpace(const std::string& str,
			const std::string& whitespace = " \t");
	std::vector<std::string> splitLineByEqualSign(std::string &inputString,
			const std::string& delimiters = "=");
	std::string removeTrailingSemicolon(const std::string& str);
public:
	GlobalConfigVars parseConfigFile(std::string configFileName);
	void updateConfigFile(GlobalConfigVars &configVar,
			std::string configFileName);
};

/**
 * All files that include this header file will aware of this global config variable.
 */
extern GlobalConfigVars globalConfigVars;

#endif /* CONFIGPARSER_H_ */
