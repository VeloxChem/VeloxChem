//
//                           VELOXCHEM 1.0-RC2
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

#ifndef RecursionBlock_hpp
#define RecursionBlock_hpp

#include <cstdint>
#include <string>

/**
 Enumerate class recblock:

 Defines all allowed key values for angular form of recursion term:
 szblock::cc - the <cart|cart> form
 szblock::sc - the <spher|cart> form
 szblock::cs - the <cart|spher> from
 szblock::ss - the <spher|spher> form
 */

enum class recblock
{
    cc,
    sc,
    cs,
    ss
};

/**
 Converts key value of  recursion block to integer number.

 @param recBlockKey the key value of recursion block.
 @return the integer number.
 */
inline int32_t
to_int(const recblock recBlockKey)
{
    return static_cast<int32_t>(recBlockKey);
}

/**
 Converts enumerate class value to it's string label.

 @param recBlockKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const recblock recBlockKey)
{
    if (recBlockKey == recblock::cc)
    {
        return std::string("Recursion Block: Cartesian-Cartesian");
    }

    if (recBlockKey == recblock::sc)
    {
        return std::string("Recursion Block: Spherical-Cartesian");
    }

    if (recBlockKey == recblock::cs)
    {
        return std::string("Recursion Block: Cartesian-Spherical");
    }

    if (recBlockKey == recblock::ss)
    {
        return std::string("Recursion Block: Spherical-Spherical");
    }

    return std::string("UNKNOWN");
}

#endif /* RecursionBlock_hpp */
