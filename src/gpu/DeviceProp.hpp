//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
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

#ifndef DeviceProp_hpp
#define DeviceProp_hpp

#include <cstdint>
#include <string>
#include <vector>

namespace gpu {

/**
 Collects information about available GPU devices.

 @param namesOfDevices the vector of names of GPU devices.
 @param globalMemoryOfDevices the vector of global memory sizes of GPU devices.
 @param computeMajorCapabilityOfDevices the vector of major compute
        capabilities of GPU devices.
 @param computeMinorCapabilityOfDevices the vector of minor compute
        capabilities of GPU devices.
 */
void getDevicesProperty(std::vector<std::string>& namesOfDevices,
                        std::vector<int64_t>&     globalMemoryOfDevices,
                        std::vector<int64_t>&     computeMajorCapabilityOfDevices,
                        std::vector<int64_t>&     computeMinorCapabilityOfDevices);
}  // namespace gpu

#endif
