#ifndef STRATEGY_PROFILER_H
#define STRATEGY_PROFILER_H

#include "Profiler.h"
#include "ProfilingStrategy.h"

class StrategyProfiler: public Profiler {
	protected:
		ProfilingStrategy* strategy;			
	public:
		StrategyProfiler(std::string identifier, ProfilingStrategy* strategy, bool rowEnd = false);
		~StrategyProfiler();
		void start();
		void stop();
		void end();
		void setStrategyOutputStream(OFStreamCSV*);
		void setSummingOutputStream(OFStreamCSV*);
		void setChildOutputStream();
};

#endif
