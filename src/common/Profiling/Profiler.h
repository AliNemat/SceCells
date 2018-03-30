#ifndef PROFILER_H
#define PROFILER_H

#include <iostream>
#include "CudaEventLifecycleHandler.h"
#include "OFStreamCSV.h"

class Profiler {
	protected:
		CudaEventLifecycleHandler timer;
		CudaEventPacket* eventPacket;
		OFStreamCSV* outputStream;
		std::string identifier;
		float time;
		bool rowEnd;

		Profiler(std::string identifier, bool rowEnd = false);

	public:
		virtual ~Profiler()=0;
		virtual void start()=0;
		virtual void stop()=0;
		virtual void end()=0;
		virtual void setStrategyOutputStream(OFStreamCSV*)=0;
		virtual void setSummingOutputStream(OFStreamCSV*)=0;
		virtual void setChildOutputStream()=0;
};

#endif

