#ifndef OF_STREAM_CSV_H
#define OF_STREAM_CSV_H

#include <iostream>
#include <fstream>
#include <string>

#define OUTPUT_PATH_STRATEGY "dataOutput/strategy_profiling_output.csv"
#define OUTPUT_PATH_SUMMING  "dataOutput/summed_profiling_output.csv"
#define COLUMN_DELIMITER_CSV ", "
#define NEWLINE_CSV          ",,,,,,\n"

enum PROFILE_TYPE {
	STRATEGY, SUMMED
};

class OFStreamCSV {
	private:
		std::ofstream streamCSV;
	public:
		OFStreamCSV(PROFILE_TYPE);
		void write(std::string);
		void write(float);
		void newColumn();	
		void newRow();
		void close();	
};

#endif
