//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef DeviceProp_hpp
#define DeviceProp_hpp

#include <cstdint>
#include <string>
#include <vector>

namespace gpu {

/**
 Collects information about available CUDA compute capable devices.

 @param namesOfDevices the vector of names of CUDA compute capable devices.
 @param globalMemoryOfDevices the vector of global memory sizes of CUDA compute capable devices.
 @param computeMajorCapabilityOfDevices the vector of major compute capabilities of CUDA
        compute capable devices.
 @param computeMinorCapabilityOfDevices the vector of minor compute capabilities of CUDA
        compute capable devices.
 */
void getDevicesProperty(std::vector<std::string>& namesOfDevices,
                        std::vector<int32_t>&     globalMemoryOfDevices,
                        std::vector<int32_t>&     computeMajorCapabilityOfDevices,
                        std::vector<int32_t>&     computeMinorCapabilityOfDevices);
}

#endif
