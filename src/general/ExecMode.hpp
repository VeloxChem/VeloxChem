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

#ifndef ExecMode_hpp
#define ExecMode_hpp

#include <cstdint>
#include <string>

/**
 Enumerate class execmode:

 Defines all allowed key values for jobs execution mode:
 execmode::cpu     - the CPU jobs execution mode
 execmode::cpu_gpu - the hybrid CPU/GPU jobs execution mode
 */

enum class execmode : int32_t
{
    cpu,
    cpu_gpu
};

/**
 Converts key value of job execution mode to integer number.

 @param execModeKey the key value of job execution mode.
 @return the integer number.
 */
inline int32_t
to_int(const execmode execModeKey)
{
    return static_cast<int32_t>(execModeKey);
}

/**
 Converts enumerate class value to it's string label.

 @param execModeKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const execmode execModeKey)
{
    if (execModeKey == execmode::cpu)
    {
        return std::string("Execution Mode: CPU");
    }

    if (execModeKey == execmode::cpu_gpu)
    {
        return std::string("Execution Mode: CPU-GPU");
    }

    return std::string("UNKNOWN");
}

#endif /* ExecMode_hpp */
