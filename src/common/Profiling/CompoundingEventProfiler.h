#ifndef COMPOUNDING_EVENT_PROFILER
#define COMPOUNDING_EVENT_PROFILER

#include "Profiler.h"

class CompoundingEventProfiler: public Profiler {
	public:	
		CompoundingEventProfiler(std::string identifier, bool rowEnd = false);
		~CompoundingEventProfiler();
		void start();
		void stop();
		void end();
		void setStrategyOutputStream(OFStreamCSV*);
		void setSummingOutputStream(OFStreamCSV*);		
		void setChildOutputStream();
};

#endif
