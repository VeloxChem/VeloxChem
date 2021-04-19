//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
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

#ifndef NumaPolicy_hpp
#define NumaPolicy_hpp

#include <cstdint>
#include <string>

/**
 Enumerate class numa:
 
 Defines all allowed key values for data policy
 numa::serial - the serial NUMA policy
 numa::parallel - the parallel NUMA policy
 */

enum class numa
{
    serial,
    parallel
};

/**
 Converts key value of numa to integer number.
 
 @param numaKey the key value of numa.
 @return the integer number.
 */
inline int32_t
to_int(const numa numaKey)
{
    return static_cast<int32_t>(numaKey);
}

/**
 Converts integer key value to numa type.
 
 @param keyValue the integer key value.
 @return the numa type.
 */
inline numa
to_numa(const int32_t keyValue)
{
    if (keyValue == to_int(numa::serial)) return numa::serial;
    
    if (keyValue == to_int(numa::parallel)) return numa::parallel;
    
    return numa::serial;
}

/**
 Converts enumerate class value to it's string label.
 
 @param numaKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const numa numaKey)
{
    if (numaKey == numa::serial)
    {
        return std::string("NUMA Policy: Serial");
    }
    
    if (numaKey == numa::parallel)
    {
        return std::string("NUMA Policy: Parallel");
    }
   
    return std::string("UNKNOWN");
}

#endif /* NumaPolicy_hpp */
