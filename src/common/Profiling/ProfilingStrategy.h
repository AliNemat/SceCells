#ifndef PROFILING_STRATEGY_H
#define PROFILING_STRATEGY_H

#include "OFStreamCSV.h"

class ProfilingStrategy {
	protected:
		bool rowEnd;
		OFStreamCSV* outputStream;
	public:
		virtual	void stop(float)=0;
		void setOutputStream(OFStreamCSV*);
		void setRowEnd(bool);
};

#endif
