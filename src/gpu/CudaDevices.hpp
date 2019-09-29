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
     Gets string representation of CUDA devices object.

     @return the representation string.
     */
    std::string getString() const;
};

#endif /* CudaDevices_hpp */
