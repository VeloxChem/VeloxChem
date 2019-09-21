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
 Allocates CUDA compute device memory for given 2D memory block.
 
 @param pointerToMemory the pointer to device memory block.
 @param pointerToPtich the pointer to pitch of memory block.
 @param memBlock2D the 2D memory block.
 */
template<class T>
void allocateDeviceMemory(      T**             pointerToMemory,
                                size_t*         pointerToPitch,
                          const CMemBlock2D<T>& memBlock2D)
{
#ifdef ENABLE_GPU
    auto cerr = cudaMallocPitch((void**)pointerToMemory, pointerToPitch, memBlock2D.size(0) * sizeof(T),
                                memBlock2D.blocks());
    
    errors::assertMsgCritical(cerr == cudaSuccess, {"allocateDeviceMemory"});
#endif
}

/**
 Frees CUDA compute device memory.

 @param pointerToMemory the pointer to device memory block.
 */
template<class T>
void freeDeviceMemory(T* pointerToMemory)
{
#ifdef ENABLE_GPU
    auto cerr = cudaFree(pointerToMemory);
    
    errors::assertMsgCritical(cerr == cudaSuccess, {"freeDeviceMemory"});
#endif
}

/**
 Copies data from 2D memory block to CUDA compute capable device memory.

 @param pointerToMemory the pointer to device memory.
 @param pitch the pitch of device memory.
 @param memBlock2D the 2D memory block.
 */
template<class T>
void copyToDeviceMemory(      T*              pointerToMemory,
                              size_t          pitch,
                        const CMemBlock2D<T>& memBlock2D)
{
#ifdef ENABLE_GPU
    
    auto cerr = cudaMemcpy2D(pointerToMemory, pitch, memBlock2D.data(), memBlock2D.pitched_size(0) * sizeof(T),
                             memBlock2D.size(0) * sizeof(T), memBlock2D.blocks(), cudaMemcpyHostToDevice);
        
    errors::assertMsgCritical(cerr == cudaSuccess, {"copyToDeviceMemory"});
#endif
}
    
/**
 Copies data from CUDA compute capable device memory to 2D memory block.

 @param pointerToMemory the pointer to device memory.
 @param pitch the pitch of device memory.
 @param memBlock2D the 2D memory block.
 */
template<class T>
void copyFromDeviceMemory(    T*          pointerToMemory,
                          size_t          pitch,
                          CMemBlock2D<T>& memBlock2D)
{
#ifdef ENABLE_GPU
    
    auto cerr = cudaMemcpy2D(memBlock2D.data(), memBlock2D.pitched_size(0) * sizeof(T), pointerToMemory, pitch,
                             memBlock2D.size(0) * sizeof(T), memBlock2D.blocks(), cudaMemcpyDeviceToHost);
    
    errors::assertMsgCritical(cerr == cudaSuccess, {"copyFromDeviceMemory"});
#endif
}

}

#endif
