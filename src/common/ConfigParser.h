#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include <tr1/unordered_map>

#ifndef CONFIGPARSER_H_
#define CONFIGPARSER_H_

class ConfigVarValue {
private:
	std::string varibleValue;
public:
	void setValue(std::string &value);
	std::string toString() const;
	int toInt() const;
	double toDouble() const;
};

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

class GlobalConfigVars {
	std::vector<ConfigVar> configVars;
	std::tr1::unordered_map<std::string, ConfigVar> configMap;
public:
	std::vector<ConfigVar> getConfigVars() const {
		return configVars;
	}
	GlobalConfigVars();
	void insertData(std::string varName, std::string varValue);
	ConfigVarValue getConfigValue(std::string varName);
	void printAll();
};

class ConfigParserException: public std::exception {
private:
	std::string _message;
public:
	ConfigParserException(const std::string& message) :
			_message(message) {
	}
	~ConfigParserException() throw () {
	}
	virtual const char* what() const throw () {
		return _message.c_str();
	}
};

class ConfigParserWarning: public std::exception {
private:
	std::string _message;
public:
	ConfigParserWarning(const std::string& message) :
			_message(message) {
	}
	~ConfigParserWarning() throw () {
	}
	virtual const char* what() const throw () {
		return _message.c_str();
	}
};

class ConfigParser {
	std::string removeLeadingAndTrailingSpace(const std::string& str,
			const std::string& whitespace = " \t");
	std::vector<std::string> splitLineByEqualSign(std::string &inputString,
			const std::string& delimiters = "=");
	std::string removeTrailingSemicolon(const std::string& str);
public:
	GlobalConfigVars parseConfigFile(std::string configFileName);
};

extern GlobalConfigVars globalConfigVars;

#endif /* CONFIGPARSER_H_ */
