#ifndef SINGLE_EVENT_PROFILER
#define SINGLE_EVENT_PROFILER

#include "Profiler.h"

class SingleEventProfiler: public Profiler {
	public:
		SingleEventProfiler(std::string identifier, bool rowEnd = false);
		~SingleEventProfiler();
		void start();
		void stop();
		void end();
		void setStrategyOutputStream(OFStreamCSV*);
		void setSummingOutputStream(OFStreamCSV*);
		void setChildOutputStream();
};

#endif
