#include "Profiler.h"

Profiler::Profiler(std::string identifier, bool rowEnd) {
	this->identifier = identifier;
	this->rowEnd = rowEnd;
}

Profiler::~Profiler() {}

