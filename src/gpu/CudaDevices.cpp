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
#include "ErrorHandler.hpp"

#ifdef ENABLE_GPU
#include "DeviceProp.hpp"
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
        auto cerr = cudaSetDevice(iDevice);
        
        errors::assertMsgCritical(cerr == cudaSuccess, {"setCudaDevice"}); 
    }
#endif
}

int32_t
CCudaDevices::getNumberOfDevices() const
{
    return static_cast<int32_t>(_namesOfDevices.size());
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
