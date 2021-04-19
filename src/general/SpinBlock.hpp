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

#ifndef SpinBlock_hpp
#define SpinBlock_hpp

#include <cstdint>

/**
 Enumerate class szblock:

 Defines all allowed key values for spin block:
 szblock::aa - the alpha-alpha spin block
 szblock::bb - the beta-beta spin block
 szblock::ab - the alpha-beta spin block
 szblock::ba - the beta-alpha spin block
 */

enum class szblock
{
    aa,
    bb,
    ab,
    ba
};

/**
 Converts key value of spin block to integer number.

 @param szBlockKey the key value of spin block.
 @return the integer number.
 */
inline int32_t
to_int(const szblock szBlockKey)
{
    return static_cast<int32_t>(szBlockKey);
}

/**
 Converts enumerate class value to it's string label.

 @param szBlockKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const szblock szBlockKey)
{
    if (szBlockKey == szblock::aa)
    {
        return std::string("Spin Block: Alpha-Alpha");
    }

    if (szBlockKey == szblock::bb)
    {
        return std::string("Spin Block: Beta-Beta");
    }

    if (szBlockKey == szblock::ab)
    {
        return std::string("Spin Block: Alpha-Beta");
    }

    if (szBlockKey == szblock::ba)
    {
        return std::string("Spin Block: Beta-Alpha");
    }

    return std::string("UNKNOWN");
}

#endif /* SpinBlock_hpp */
