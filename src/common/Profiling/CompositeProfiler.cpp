#include "CompositeProfiler.h"

CompositeProfiler::CompositeProfiler(std::string identifier, bool rowEnd): Profiler(identifier, rowEnd) {}

CompositeProfiler::~CompositeProfiler() {
	for (unsigned i = 0; i < profilers.size(); i++)
		delete profilers.at(i);	
}

void CompositeProfiler::setStrategyOutputStream(OFStreamCSV* outputStream) {	
	for (unsigned i = 0; i < profilers.size(); i++)
		profilers.at(i)->setStrategyOutputStream(outputStream);

}

void CompositeProfiler::setSummingOutputStream(OFStreamCSV* outputStream) {
	for (unsigned i = 0; i < profilers.size(); i++)
		profilers.at(i)->setSummingOutputStream(outputStream);

}

void CompositeProfiler::setChildOutputStream() {
	for (unsigned i = 0; i < profilers.size(); i++)
		profilers.at(i)->setChildOutputStream();
}

void CompositeProfiler::addChild(Profiler* child) {
	profilers.push_back(child);
}

void CompositeProfiler::start() {
	for (unsigned i = 0; i < profilers.size(); i++)
		profilers.at(i)->start();
}

void CompositeProfiler::stop() {
	for (unsigned i = 0; i < profilers.size(); i++)
		profilers.at(i)->stop();
}

void CompositeProfiler::end() {
	for (unsigned i = 0; i < profilers.size(); i++)
		profilers.at(i)->end();
}	
