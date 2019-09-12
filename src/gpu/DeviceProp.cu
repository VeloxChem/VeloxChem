//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include <cmath>
#include <sstream>
#include <string>

#include "DeviceProp.hpp"
#include "StringFormat.hpp"

namespace gpu {  // gpu namespace

void
getDeviceProperties(std::vector<std::string>& namesOfDevices,
                    std::vector<int32_t>&     globalMemoryOfDevices,
                    std::vector<int32_t>&     computeMajorCapabilityOfDevices,
                    std::vector<int32_t>&     computeMinorCapabilityOfDevices)
{
#ifdef ENABLE_GPU

    int devcnt = 0;

    cudaGetDeviceCount(&devcnt);

    for (int i = 0; i < devcnt; i++)
    {
        cudaDeviceProp prop;

        cudaGetDeviceProperties(&prop, i);

        namesOfDevices.push_back(std::string(prop.name));

        globalMemoryOfDevices.push_back(static_cast<int32_t>(prop.totalGlobalMem));

        computeMajorCapabilityOfDevices.push_back(static_cast<int32_t>(prop.major));

        computeMinorCapabilityOfDevices.push_back(static_cast<int32_t>(prop.minor));
    }
#endif
}

}  // namespace gpu
