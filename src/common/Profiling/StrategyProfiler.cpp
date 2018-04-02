#include "StrategyProfiler.h"

StrategyProfiler::StrategyProfiler(std::string identifier,  
                                   ProfilingStrategy* strategy,
                                   bool rowEnd): Profiler(identifier, rowEnd) {
	this->strategy = strategy;
	strategy->setRowEnd(this->rowEnd);
}

StrategyProfiler::~StrategyProfiler() {
	//FIXME: strategy has virtual destructor
	delete strategy;
	delete eventPacket;
}

void StrategyProfiler::start() {
	eventPacket = timer.start();
}

void StrategyProfiler::stop() {
	if (eventPacket != NULL) {
		time = timer.getElapsedTime(eventPacket);
		eventPacket = NULL;
		
		strategy->stop(time);
	}	
}

void StrategyProfiler::end() {}

void StrategyProfiler::setStrategyOutputStream(OFStreamCSV* outputStream) {
	this->outputStream = outputStream;
	outputStream->write(identifier);

	if (rowEnd)
		outputStream->newRow();
	else
		outputStream->newColumn();

	setChildOutputStream();
}	

void StrategyProfiler::setSummingOutputStream(OFStreamCSV*) {}

void StrategyProfiler::setChildOutputStream() {
	strategy->setOutputStream(outputStream);
}
