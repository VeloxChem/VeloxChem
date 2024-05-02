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

#ifndef GpuDevices_hpp
#define GpuDevices_hpp

#include <cstdint>
#include <string>
#include <vector>

/**
 Class GpuDevices stores data about available GPU devices and provides these
 devices managment functions.

 @author Z. Rinkevicius
 */
class CGpuDevices
{
    /**
     The vector of names of GPU devices.
     */
    std::vector<std::string> _namesOfDevices;

    /**
     The vector of available global memory on GPU devices.
     */
    std::vector<int64_t> _globalMemoryOfDevices;

    /**
     The vector of compute major capabilities on GPU devices.
     */
    std::vector<int64_t> _computeMajorCapabilityOfDevices;

    /**
     The vector of compute minor capabilities on GPU devices.
     */
    std::vector<int64_t> _computeMinorCapabilityOfDevices;

   public:
    /**
     Creates an GPU devices object.
     */
    CGpuDevices();

    /**
     Gets number of GPU devices.

     @return the number of devices.
     */
    int64_t getNumberOfDevices() const;

    /**
     Gets string representation of GPU devices object.

     @return the representation string.
     */
    std::string getString() const;
};

#endif /* GpuDevices_hpp */
