#ifndef CUDA_EVENT_LIFECYCLE_HANDLER
#define CUDA_EVENT_LIFECYCLE_HANDLER

#include "cuda_runtime.h"

struct CudaEventPacket {
	cudaEvent_t* start;
	cudaEvent_t* stop;

	CudaEventPacket() {
		start = new cudaEvent_t();
		stop = new cudaEvent_t();
	}

	~CudaEventPacket() {
		delete start;
		delete stop;
	}
};

class CudaEventLifecycleHandler {
	public:
		CudaEventPacket* start();
		float getElapsedTime(CudaEventPacket*);
};

#endif
