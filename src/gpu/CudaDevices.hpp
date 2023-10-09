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

#ifndef CudaDevices_hpp
#define CudaDevices_hpp

#include <cstdint>
#include <string>
#include <vector>

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
     The vector of available global memory on CUDA compute capable devices.
     */
    std::vector<int64_t> _globalMemoryOfDevices;

    /**
     The vector of compute major capabilities on CUDA compute capable devices.
     */
    std::vector<int64_t> _computeMajorCapabilityOfDevices;

    /**
     The vector of compute minor capabilities on CUDA compute capable devices.
     */
    std::vector<int64_t> _computeMinorCapabilityOfDevices;

   public:
    /**
     Creates an CUDA devices object.
     */
    CCudaDevices();

    /**
     Gets number of CUDA compute capable devices.

     @return the number of devices.
     */
    int64_t getNumberOfDevices() const;

    /**
     Gets string representation of CUDA devices object.

     @return the representation string.
     */
    std::string getString() const;
};

#endif /* CudaDevices_hpp */
