#include "ProfilingStrategy.h"

void ProfilingStrategy::setOutputStream(OFStreamCSV* outputStream) {
	this->outputStream = outputStream;
}

void ProfilingStrategy::setRowEnd(bool rowEnd) {
	this->rowEnd = rowEnd;
}
