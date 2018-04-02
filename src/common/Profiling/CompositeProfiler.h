#ifndef COMPOSITE_PROFILER_H
#define COMPOSITE_PROFILER_H

#include <vector>
#include "Profiler.h"
#include "OFStreamCSV.h"

class CompositeProfiler:public Profiler {
	protected:
		std::vector<Profiler*> profilers;		
	public:
		CompositeProfiler(std::string identifier, bool rowEnd = false);
		~CompositeProfiler();
		void setStrategyOutputStream(OFStreamCSV*);
		void setSummingOutputStream(OFStreamCSV*);
		void setChildOutputStream();
		void addChild(Profiler*);
		void start();
		void stop();
		void end();
};

#endif
