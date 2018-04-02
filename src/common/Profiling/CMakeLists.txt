include_directories(${COMMON_INCLUDES})

add_library(SceProfiling
	Profiler.h
	Profiler.cpp
	SingleEventProfiler.h
	SingleEventProfiler.cpp
	CompoundingEventProfiler.h
	CompoundingEventProfiler.cpp
	StrategyProfiler.h
	StrategyProfiler.cpp
	CompositeProfiler.h
	CompositeProfiler.cpp
	ProfilingCoordinator.h
	ProfilingCoordinator.cpp
	CudaEventLifecycleHandler.h
	CudaEventLifecycleHandler.cpp
	OFStreamCSV.h
	OFStreamCSV.cpp
	ProfilingStrategy.h
	ProfilingStrategy.cpp
	AveragedBlockStrategy.h
	AveragedBlockStrategy.cpp
)

