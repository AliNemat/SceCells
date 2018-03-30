#include "SingleEventProfiler.h"

SingleEventProfiler::SingleEventProfiler(std::string identifier, bool rowEnd): Profiler(identifier, rowEnd) {}


SingleEventProfiler::~SingleEventProfiler() {
	delete eventPacket;
}

void SingleEventProfiler::start() {
	eventPacket = timer.start();	
}

void SingleEventProfiler::stop() {}

void SingleEventProfiler::end() {
	if (eventPacket != NULL) {
		time = timer.getElapsedTime(eventPacket);	
		outputStream->write(time);
 
		outputStream->newColumn();

		eventPacket = NULL;
	}
}

void SingleEventProfiler::setStrategyOutputStream(OFStreamCSV* outputStream) {}

void SingleEventProfiler::setSummingOutputStream(OFStreamCSV* outputStream) {
	this->outputStream = outputStream;
	outputStream->write(identifier);

	if (rowEnd)
		outputStream->newRow();
	else
		outputStream->newColumn();

	setChildOutputStream();
}	

void SingleEventProfiler::setChildOutputStream() {}
