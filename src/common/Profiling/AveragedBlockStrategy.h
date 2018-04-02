#ifndef AVERAGED_BLOCK_STRATEGY_H
#define AVERAGED_BLOCK_STRATEGY_H

#include "ProfilingStrategy.h"

#define BLOCK_SIZE 250

class AveragedBlockStrategy: public ProfilingStrategy {
	private:
		unsigned blockSize;
		unsigned currentCount;
		float time;
	public:
		AveragedBlockStrategy();
		void stop(float);
};

#endif
