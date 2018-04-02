#include "OFStreamCSV.h"

OFStreamCSV::OFStreamCSV(PROFILE_TYPE type) {
	switch(type) {
		case STRATEGY: 
			streamCSV.open(OUTPUT_PATH_STRATEGY);
			break;
		case SUMMED:
			streamCSV.open(OUTPUT_PATH_SUMMING);
			break;
		default:
			break;
	}
}

void OFStreamCSV::write(std::string toWrite) {
	streamCSV << toWrite;
}

void OFStreamCSV::write(float toWrite) {
	streamCSV << toWrite;
}

void OFStreamCSV::newColumn() {
	streamCSV << COLUMN_DELIMITER_CSV << " ";
}

void OFStreamCSV::newRow() {
	streamCSV << NEWLINE_CSV;
}

void OFStreamCSV::close() {
	streamCSV.close();
}
