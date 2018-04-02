#include "CompoundingEventProfiler.h"

CompoundingEventProfiler::CompoundingEventProfiler(std::string identifier, bool rowEnd): Profiler(identifier, rowEnd) {
	time = 0;
}

CompoundingEventProfiler::~CompoundingEventProfiler() { 
	delete eventPacket;
}

void CompoundingEventProfiler::start() {
	eventPacket = timer.start();
}

void CompoundingEventProfiler::stop() {	if (eventPacket != NULL) {	
		float lastTimeBlock = timer.getElapsedTime(eventPacket);	
		time += lastTimeBlock;	
		eventPacket = NULL;
	}
}

void CompoundingEventProfiler::end() {
	outputStream->write(time);
	outputStream->newColumn();
}

void CompoundingEventProfiler::setStrategyOutputStream(OFStreamCSV* outputStream) {}

void CompoundingEventProfiler::setSummingOutputStream(OFStreamCSV* outputStream) {
	this->outputStream = outputStream;
	outputStream->write(identifier);

	if (rowEnd)
		outputStream->newRow();
	else
		outputStream->newColumn();

	setChildOutputStream();
}	

void CompoundingEventProfiler::setChildOutputStream() {}
