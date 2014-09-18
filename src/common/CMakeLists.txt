include_directories(${COMMON_INCLUDES})

# For some reason I encountered strange linking problems with this lib
# therefore I switched to header file only instead of lib
# by using inline functions
#CUDA_ADD_LIBRARY(SceUtilsGPU
#      commonGPUData.h
#      #commonGPUData.cu
#)

add_library(SceUtils
      GeoVector.h
      GeoVector.cpp
      Point2D.h
      Point2D.cpp
      ConfigParser.h      
      ConfigParser.cpp
      commonData.h
      commonData.cpp
)
