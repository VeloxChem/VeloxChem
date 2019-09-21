//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef CudaDevices_hpp
#define CudaDevices_hpp

#include <cstdint>
#include <cstdlib>
#include <vector>
#include <string>

#include "DeviceFunc.hpp"

/**
 Class CudaDevices stores data about available CUDA devices and provides these devices managment functions.
 
 @author Z. Rinkevicius
 */
class CCudaDevices
{
    /**
     The vector of names of CUDA compute capable devices.
     */
    std::vector<std::string> _namesOfDevices;
    
    /**
     The vector of available global memory (in megabyte) on CUDA compute capable devices.
     */
    std::vector<int32_t> _globalMemoryOfDevices;
    
    /**
     The vector of compute major capabilities on CUDA compute capable devices.
     */
    std::vector<int32_t> _computeMajorCapabilityOfDevices;
    
    /**
     The vector of compute minor capabilities on CUDA compute capable devices.
     */
    std::vector<int32_t> _computeMinorCapabilityOfDevices;
    
public:
    /**
     Creates an CUDA devices object.
     */
    CCudaDevices();
    
    /**
     Destroys a CUDA devices object.
     */
    ~CCudaDevices();
    
    /**
     Sets up requested CUDA compute capable device.

     @param iDevice the identifier of CUDA compute capable device.
     */
    void setCudaDevice(const int32_t iDevice) const;
    
    /**
     Gets number of CUDA compute capable devices.

     @return the number of devices.
     */
    int32_t getNumberOfDevices() const;
    
    /**
     Allocates 2D memory block on CUDA compute capable device.

     @param pointer the pointer to 2D memory block.
     @param pitch the pointer to pitch of 2D  memory block.
     @param nElements the number of columns in 2D  memory block.
     @param nBlocks the number of rows in 2D memory block.
     */
    void allocate(double** pointer,
                  size_t*  pitch,
                  int32_t  nElements,
                  int32_t  nBlocks) const
    {
#ifdef ENABLE_GPU
        gpu::allocateDeviceMemory((void**)pointer, pitch, nElements * sizeof(double), static_cast<size_t>(nBlocks));
#endif
    }
    
    /**
     Deallocates device memory.
     
     @param pointer the pointer to device memory.
     */
    void free(double* pointer) const
    {
#ifdef ENABLE_GPU
        gpu::freeDeviceMemory((void*)pointer);
#endif
    }
    
    /**
     Gets string representation of CUDA devices object.

     @return the representation string.
     */
    std::string getString() const;
};

#endif /* CudaDevices_hpp */
