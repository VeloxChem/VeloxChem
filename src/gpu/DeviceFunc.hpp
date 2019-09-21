//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DeviceFunc_hpp
#define DeviceFunc_hpp

#include <cstdint>
#include <tuple>

#include "MemBlock2D.hpp"
#include "ErrorHandler.hpp"

namespace gpu {

/**
 Sets up requested  CUDA compute capable device.

 @param iDevice the identifier of CUDA compute capable device.
 */
void setDevice(const int32_t iDevice);
    
/**
 Allocates CUDA compute device memory.

 @param pointer the pointer to device memory.
 @param dataPitch the pointer  to pitch of device memory.
 @param dataWidth the width of data.
 @param dataHeight the height  of data.
 */
void allocateDeviceMemory(void**  pointer,
                          size_t* dataPitch,
                          size_t  dataWidth,
                          size_t  dataHeight);

/**
 Frees CUDA compute device memory.

 @param pointer the pointer to device memory.
 */
void freeDeviceMemory(void* pointer);

/**
 Copies 2D data from host memory to device memory.

 @param destination the pointer to device memory.
 @param destinationPitch the pitch of device memory.
 @param source the pointer to host memory.
 @param sourcePitch the pitch of host memory.
 @param dataWidth the width of data.
 @param dataHeight the height of data.
 */
void copyToDeviceMemory(      void*  destination,
                              size_t destinationPitch,
                        const void*  source,
                              size_t sourcePitch,
                              size_t dataWidth,
                              size_t dataHeight);

/**
 Copies 2D data from device memory to host memory.

 @param destination the pointer to host memory.
 @param destinationPitch the pitch of host memory.
 @param source the pointer to device memory.
 @param sourcePitch the pitch of device memory.
 @param dataWidth the width of data.
 @param dataHeight the height of data.
 */
void copyFromDeviceMemory(      void*  destination,
                                size_t destinationPitch,
                          const void*  source,
                                size_t sourcePitch,
                                size_t dataWidth,
                                size_t dataHeight);
    
}

#endif
