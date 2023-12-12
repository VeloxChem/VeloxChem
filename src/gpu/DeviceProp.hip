//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2023 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include <hip/hip_runtime.h>

#include <iostream>

#include "DeviceProp.hpp"

#define hipSafe(e)                                                                                                       \
    {                                                                                                                    \
        hipError_t err = (e);                                                                                            \
        if (err != hipSuccess)                                                                                           \
        {                                                                                                                \
            std::cerr << "HIP error in " << __FILE__ << ":" << __LINE__ << ": " << hipGetErrorString(err) << std::endl;  \
            std::exit(EXIT_FAILURE);                                                                                     \
        }                                                                                                                \
    }

namespace gpu {  // gpu namespace

void
getDevicesProperty(std::vector<std::string>& namesOfDevices,
                   std::vector<int64_t>&     globalMemoryOfDevices,
                   std::vector<int64_t>&     computeMajorCapabilityOfDevices,
                   std::vector<int64_t>&     computeMinorCapabilityOfDevices)
{
    int device_count = 0;

    hipSafe(hipGetDeviceCount(&device_count));

    for (int i = 0; i < device_count; i++)
    {
        hipDeviceProp_t prop;

        hipSafe(hipGetDeviceProperties(&prop, i));

        namesOfDevices.push_back(std::string(prop.name));

        globalMemoryOfDevices.push_back(static_cast<int64_t>(prop.totalGlobalMem));

        computeMajorCapabilityOfDevices.push_back(static_cast<int64_t>(prop.major));

        computeMinorCapabilityOfDevices.push_back(static_cast<int64_t>(prop.minor));
    }
}

}  // namespace gpu
