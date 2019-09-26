//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CudaDevices.hpp"

#include <sstream>
#include <cmath>

#include "StringFormat.hpp"

#ifdef ENABLE_GPU
#include "DeviceProp.hpp"
#include "DeviceFunc.hpp"
#endif

CCudaDevices::CCudaDevices()

    : _namesOfDevices(std::vector<std::string>())

    , _globalMemoryOfDevices(std::vector<int32_t>())

    , _computeMajorCapabilityOfDevices(std::vector<int32_t>())

    , _computeMinorCapabilityOfDevices(std::vector<int32_t>())
{
#ifdef ENABLE_GPU
    gpu::getDevicesProperty(_namesOfDevices, _globalMemoryOfDevices, _computeMajorCapabilityOfDevices, _computeMinorCapabilityOfDevices);
#endif
}

CCudaDevices::~CCudaDevices()
{
    
}

void
CCudaDevices::setCudaDevice(const int32_t iDevice) const
{
#ifdef ENABLE_GPU
    if (iDevice < getNumberOfDevices())
    {
        gpu::setDevice(iDevice); 
    }
#endif
}

void
CCudaDevices::synchronizeCudaDevice() const
{
#ifdef ENABLE_GPU
    gpu::synchronizeCudaDevice(); 
#endif
}

int32_t
CCudaDevices::getNumberOfDevices() const
{
    return static_cast<int32_t>(_namesOfDevices.size());
}

void
CCudaDevices::allocate(double** pointer,
                       size_t*  pitch,
                       int32_t  nElements,
                       int32_t  nBlocks) const
{
#ifdef ENABLE_GPU
    gpu::allocateDeviceMemory((void**) pointer, pitch, nElements * sizeof(double), static_cast<size_t>(nBlocks));
#endif
}

void
CCudaDevices::free(double* pointer) const
{
#ifdef ENABLE_GPU
    gpu::freeDeviceMemory((void*)pointer);
#endif
}

void
CCudaDevices::copyToDevice(      double*              pointer,
                                 size_t               pitch,
                           const CMemBlock2D<double>& memBlock2D) const
{
#ifdef ENABLE_GPU
    gpu::copyToDeviceMemory(pointer, pitch, memBlock2D.data(0), memBlock2D.pitched_size(0) * sizeof(double),
                            memBlock2D.size(0) * sizeof(double), static_cast<size_t>(memBlock2D.blocks()));
#endif
}

void
CCudaDevices::copyFromDevice(double*              pointer,
                             size_t               pitch,
                             CMemBlock2D<double>& memBlock2D) const
{
#ifdef ENABLE_GPU
    gpu::copyFromDeviceMemory(memBlock2D.data(0), memBlock2D.pitched_size(0) * sizeof(double), pointer, pitch,
                              memBlock2D.size(0) * sizeof(double), static_cast<size_t>(memBlock2D.blocks()));
#endif
}

std::string
CCudaDevices::getString() const
{
    std::stringstream ss;
    
    std::string str("CUDA Compute Capable Devices");
    
    ss << str << "\n";
    
    ss << std::string(str.size() + 2, '=') << "\n\n";
    
    const int32_t width = 50;
    
    for (size_t i = 0; i < _namesOfDevices.size(); i++)
    {
        str.assign("GPU device ID: ");
        
        str.append(std::to_string(i));
        
        ss << fstr::format(str, width, fmt::left) << "\n";
        
        str.assign("  Device name:             ");
        
        str.append(_namesOfDevices[i]);
        
        ss << fstr::format(str, width, fmt::left) << "\n";
        
        str.assign("  Compute capability:      ");
        
        str.append(std::to_string(_computeMajorCapabilityOfDevices[i]));
        
        str.append(".");
        
        str.append(std::to_string(_computeMinorCapabilityOfDevices[i]));
        
        ss << fstr::format(str, width, fmt::left) << "\n";
        
        str.assign("  Global memory:           ");
        
        auto glbmem = static_cast<double>(_globalMemoryOfDevices[i]) / 1024.0;
        
        str.append(fstr::to_string(glbmem, 0));
        
        str.append(" GB");
        
        ss << fstr::format(str, width, fmt::left) << "\n";
    }
    
    return ss.str();
}

int32_t
CCudaDevices::getGridBlockSize(const int32_t iDevice) const
{
    // FIX ME: set up optimal block size using device name
    
    return 256;
}
