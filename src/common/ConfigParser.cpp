#include "ConfigParser.h"
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>

std::string ConfigVarValue::toString() const {
	return this->varibleValue;
}
int ConfigVarValue::toInt() const {
	int result;
	std::stringstream strStream(this->varibleValue);
	strStream >> result;
	return result;
}

double ConfigVarValue::toDouble() const {
	double result;
	std::stringstream strStream(this->varibleValue);
	strStream >> result;
	return result;
}

void ConfigVarValue::setValue(std::string &value) {
	varibleValue = value;
}

ConfigVar::ConfigVar(std::string varName) {
	std::string emptyStr("");
	varibleName = varName;
	varibleValue.setValue(emptyStr);
}

GlobalConfigVars::GlobalConfigVars() {
}

/*
 * The insertion and lookup process could be done much more efficiently with std_map.
 * But it is not necessary because there are not so many entries in the config file
 */
void GlobalConfigVars::insertData(std::string varName, std::string varValue) {
	std::vector<ConfigVar>::iterator it = configVars.begin();
	//std::cout << "inserting name :" << varName << "inserting value:" << varValue
	//		<< std::endl;
	while (it != configVars.end()) {
		if (it->getVarName() == varName) {
			if (it->getValue().toString().length() != 0) {
				throw SceException(
						"Multiple definition of config value is not permitted",
						ConfigValueException);
			} else {
				it->setValue(varValue);
			}
			return;
		}
		++it;
	}
	ConfigVar dataEntry(varName);
	dataEntry.setValue(varValue);
	configMap.insert(std::pair<std::string, ConfigVar>(varName, dataEntry));
	configVars.push_back(dataEntry);
}

void GlobalConfigVars::updateData(std::string varName, std::string varValue) {
	std::vector<ConfigVar>::iterator it = configVars.begin();
	uint keyCount = 0;
	while (it != configVars.end()) {
		if (it->getVarName() == varName) {
			it->setValue(varValue);
			keyCount++;
		}
		++it;
	}
	if (keyCount > 1) {
		throw SceException(
				"Try to update, but Old config file is illegal because multiple definitions found",
				ConfigValueException);
	}
	ConfigVar dataEntry(varName);
	dataEntry.setValue(varValue);
	std::tr1::unordered_map<std::string, ConfigVar>::iterator it2 =
			configMap.find(varName);
	if (it2 != configMap.end()) {
		it2->second = dataEntry;
	} else {
		configMap.insert(std::pair<std::string, ConfigVar>(varName, dataEntry));
	}
	configVars.push_back(dataEntry);

}

ConfigVarValue GlobalConfigVars::getConfigValue(std::string varName) {
	std::tr1::unordered_map<std::string, ConfigVar>::iterator it =
			configMap.find(varName);
	if (it == configMap.end()) {
		throw SceException(
				"Error in lookup config varible, following key cannot be found:"
						+ varName, ConfigValueException);
	} else {
		//std::cout << "return value:" << (it->second).getValue().toString()
		//		<< std::endl;
		return (it->second.getValue());
	}
}

void GlobalConfigVars::printAll() {
	std::vector<ConfigVar>::iterator it = configVars.begin();
	while (it != configVars.end()) {
		std::cout << std::left << "Key: " << std::setw(25) << it->getVarName()
				<< " Value: " << std::setw(25) << it->getValue().toString()
				<< std::endl;
		++it;
	}
	//getchar();
}

std::string ConfigParser::removeLeadingAndTrailingSpace(const std::string& str,
		const std::string& whitespace) {
	const size_t strBegin = str.find_first_not_of(whitespace);
	if (strBegin == std::string::npos)
		return ""; // no content

	const size_t strEnd = str.find_last_not_of(whitespace);
	const size_t strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}

std::vector<std::string> ConfigParser::splitLineByEqualSign(
		std::string &inputString, const std::string& delimiters) {
	std::string s = inputString;
	std::vector<std::string> result;
	size_t pos = 0;
	std::string token;
	while ((pos = s.find(delimiters)) != std::string::npos) {
		token = s.substr(0, pos);
		result.push_back(token);
		s.erase(0, pos + delimiters.length());
	}
	result.push_back(s);
	return result;
}

GlobalConfigVars ConfigParser::parseConfigFile(std::string configFileName) {
	std::ifstream infile(configFileName.c_str());
	if (!infile.is_open()) {
		throw SceException(
				"Fatal error: Config file not found!" + configFileName,
				ConfigFileNotFound);
	}

	std::string line;
	std::vector<std::string> tmpReading;
	GlobalConfigVars result;
	while (std::getline(infile, line)) {
		std::string tmp = removeLeadingAndTrailingSpace(line);
		if (tmp.length() == 0) {
			continue;
		} else if (tmp[0] == '#') {
			continue;
		}
		tmpReading = splitLineByEqualSign(line);
		if (tmpReading.size() != 2) {
			throw SceException(
					"Error in Config file: More than one equal sign found in one line :"
							+ line, ConfigValueException);
		}
		std::string varName = removeLeadingAndTrailingSpace(tmpReading[0]);
		std::string varValue = removeLeadingAndTrailingSpace(tmpReading[1]);
		varValue = removeTrailingSemicolon(varValue);
		result.insertData(varName, varValue);
	}
	infile.close();
	std::vector<ConfigVar> configVaribles = result.getConfigVars();
	std::vector<ConfigVar>::iterator it = configVaribles.begin();
	while (it != configVaribles.end()) {
		if (it->getValue().toString().length() == 0) {
			throw SceException(
					"one or more config value is not defined in config file",
					ConfigValueException);
		}
		++it;
	}
	return result;
}

std::string ConfigParser::removeTrailingSemicolon(const std::string& str) {
	std::string trailingSignal = " ;";

	const size_t strBegin = 0;
	const size_t strEnd = str.find_last_not_of(trailingSignal);
	const size_t strRange = strEnd - strBegin + 1;

	return str.substr(strBegin, strRange);
}

void ConfigParser::updateConfigFile(GlobalConfigVars& configVar,
		std::string configFileName) {
	std::ifstream infile(configFileName.c_str());
	if (!infile.is_open()) {
		throw SceException(
				"Fatal error: Config file not found!" + configFileName,
				ConfigFileNotFound);
	}
	std::string line;
	std::vector<std::string> tmpReading;
	while (std::getline(infile, line)) {
		std::string tmp = removeLeadingAndTrailingSpace(line);
		if (tmp.length() == 0) {
			continue;
		} else if (tmp[0] == '#') {
			continue;
		}
		tmpReading = splitLineByEqualSign(line);
		if (tmpReading.size() != 2) {
			throw SceException(
					"Error in Config file: More than one equal sign found in one line :"
							+ line, ConfigValueException);
		}
		std::string varName = removeLeadingAndTrailingSpace(tmpReading[0]);
		std::string varValue = removeLeadingAndTrailingSpace(tmpReading[1]);
		varValue = removeTrailingSemicolon(varValue);
		configVar.updateData(varName, varValue);
	}
	infile.close();
	std::vector<ConfigVar> configVaribles = configVar.getConfigVars();
	std::vector<ConfigVar>::iterator it = configVaribles.begin();
	while (it != configVaribles.end()) {
		if (it->getValue().toString().length() == 0) {
			throw SceException(
					"one or more config value is not defined in config file",
					ConfigValueException);
		}
		++it;
	}
}

ConfigVarsCollection ConfigParser::parseConfigCollection(
		std::string configFileName) {
	ConfigVarsCollection result;
	uint totalConfigCount = 0;
	uint currentConfigSeq;
	std::tr1::unordered_map<int, uint> sequenceMap;
	std::ifstream infile(configFileName.c_str());
	if (!infile.is_open()) {
		throw SceException(
				"Fatal error: Config collection file not found!"
						+ configFileName, ConfigFileNotFound);
	}
	std::string line;
	std::vector<std::string> tmpReading;
	while (std::getline(infile, line)) {
		std::string tmp = removeLeadingAndTrailingSpace(line);
		if (tmp.length() == 0) {
			continue;
		} else if (tmp[0] == '#') {
			continue;
		} else if (tmp[0] == '$') {
			if (tmp.length() < 3) {
				continue;
			} else {
				uint i = 1;
				while (i < tmp.length() && tmp[i] != '$') {
					i++;
				}
				// need to make sure format of the line is correct.
				if (i >= tmp.length() - 1) {
					throw SceException(
							"Error in Config file: '$' symbol found but line is not complete:"
									+ line, ConfigValueException);
				} else {
					std::string sequenceStr = tmp.substr(1, i - 1);
					std::istringstream buffer1(sequenceStr);
					int seq;
					buffer1 >> seq;
					std::cout << "sequence is " << seq << std::endl;
					if (sequenceMap.find(seq) == sequenceMap.end()) {
						currentConfigSeq = totalConfigCount;
						sequenceMap.insert(
								std::pair<int, uint>(seq, currentConfigSeq));
						totalConfigCount++;
						result.configVarSets.push_back(GlobalConfigVars());
					} else {
						currentConfigSeq = sequenceMap.find(seq)->second;
					}
					std::string configPairStr = tmp.substr(i + 1);
					tmpReading = splitLineByEqualSign(configPairStr);
					if (tmpReading.size() != 2) {
						throw SceException(
								"Error in Config file: More than one equal sign found in one line :"
										+ line, ConfigValueException);
					}
					std::string varName = removeLeadingAndTrailingSpace(
							tmpReading[0]);
					std::string varValue = removeLeadingAndTrailingSpace(
							tmpReading[1]);
					varValue = removeTrailingSemicolon(varValue);
					if (varValue.length() == 0) {
						throw SceException(
								"empty config value is not allowed in config file",
								ConfigValueException);
					}
					result.configVarSets[currentConfigSeq].insertData(varName,
							varValue);
				}
			}
		}
	}
	infile.close();
	return result;
}

void GlobalConfigVars::updateFromConfig(GlobalConfigVars& otherConfigVar) {
	std::vector<ConfigVar> configVarVector = otherConfigVar.getConfigVars();
	std::vector<ConfigVar>::iterator it = configVarVector.begin();
	while (it != configVarVector.end()) {
		updateData(it->getVarName(), it->getValue().toString());
		++it;
	}
}
