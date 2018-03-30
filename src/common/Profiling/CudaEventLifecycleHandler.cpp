#include "CudaEventLifecycleHandler.h"

CudaEventPacket* CudaEventLifecycleHandler::start() {
	CudaEventPacket* eventPacket = new CudaEventPacket();	

	cudaEventCreate(eventPacket->start);
	cudaEventCreate(eventPacket->stop);		
	cudaEventRecord(*(eventPacket->start));

	return eventPacket;
}

float CudaEventLifecycleHandler::getElapsedTime(CudaEventPacket* eventPacket) {
	cudaEvent_t start = *(eventPacket->start);
	cudaEvent_t stop = *(eventPacket->stop);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	delete eventPacket;

	return elapsedTime; 
}


