#ifndef PROFILING_COORDINATOR_H
#define PROFILING_COORDINATOR_H

#include <vector>
#include "CompositeProfiler.h"
#include "OFStreamCSV.h"

class ProfilingCoordinator {
	private:
		static std::vector<Profiler*> profilers;
		static OFStreamCSV* strategyOutputStream;
		static OFStreamCSV* summingOutputStream;
	public:
		//DO NOT define destructor, otherwise, static member will be deallocated
		unsigned addProfiler(Profiler*) const;
		unsigned addCompProfiler(CompositeProfiler*) const;
		void startProfiler(unsigned) const;
		void stopProfiler(unsigned) const;
		void end();
};

#endif
