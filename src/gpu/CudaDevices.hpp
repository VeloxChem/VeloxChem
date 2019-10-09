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
#include "MemBlock2D.hpp"

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
     Synchorizes current CUDA compute capable device.
     */
    void synchronizeCudaDevice() const;
    
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
    template <class T>
    void allocate(T**     pointer,
                  size_t* pitch,
                  int32_t nElements,
                  int32_t nBlocks) const;
    
    /**
     Deallocates device memory.
     
     @param pointer the pointer to device memory.
     */
    template <class T>
    void free(T* pointer) const;
    
    /**
     Copies 2D data to CUDA device.

     @param pointer the pointer to device memory.
     @param pitch the pitch of device memory.
     @param memBlock2D the 2D memory block.
     */
    template <class T>
    void copyToDevice(      T*              pointer,
                            size_t          pitch,
                      const CMemBlock2D<T>& memBlock2D) const;
    
    /**
     Copies 2D data from CUDA device.
     
     @param pointer the pointer to device memory.
     @param pitch the pitch of device memory.
     @param memBlock2D the 2D memory block.
     */
    template <class T>
    void copyFromDevice(T*              pointer,
                        size_t          pitch,
                        CMemBlock2D<T>& memBlock2D) const;
    
    /**
     Gets string representation of CUDA devices object.

     @return the representation string.
     */
    std::string getString() const;
    
    /**
     Gets recommended size of grid block for specific CUDA compute capable device.

     @param iDevice the identifier of CUDA compute capable device.
     @return the recomended grid block size.
     */
    int32_t getGridBlockSize(const int32_t iDevice = 0) const;
};

template <class T>
void
CCudaDevices::allocate(T**     pointer,
                       size_t* pitch,
                       int32_t nElements,
                       int32_t nBlocks) const
{
#ifdef ENABLE_GPU
    gpu::allocate_device_memory((void**) pointer, pitch, nElements * sizeof(T), static_cast<size_t>(nBlocks));
#endif
}

template <class T>
void
CCudaDevices::free(T* pointer) const
{
#ifdef ENABLE_GPU
    gpu::free_device_memory((void*)pointer);
#endif
}

template <class T>
void
CCudaDevices::copyToDevice(T*                    pointer,
                           size_t                pitch,
                           const CMemBlock2D<T>& memBlock2D) const
{
#ifdef ENABLE_GPU
    gpu::copy_to_device_memory(pointer, pitch, memBlock2D.data(0), memBlock2D.pitched_size(0) * sizeof(T),
                               memBlock2D.size(0) * sizeof(T), static_cast<size_t>(memBlock2D.blocks()));
#endif
}

template <class T>
void
CCudaDevices::copyFromDevice(T*              pointer,
                             size_t          pitch,
                             CMemBlock2D<T>& memBlock2D) const
{
#ifdef ENABLE_GPU
    gpu::copy_from_device_memory(memBlock2D.data(0), memBlock2D.pitched_size(0) * sizeof(T), pointer, pitch,
                                 memBlock2D.size(0) * sizeof(T), static_cast<size_t>(memBlock2D.blocks()));
#endif
}

#endif /* CudaDevices_hpp */
