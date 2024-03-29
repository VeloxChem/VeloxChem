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

#ifndef FmtType_hpp
#define FmtType_hpp

/**
 Enumerate class fmt:

 Defines supported formatting keys:
 fmt::left   - the left alignment of output line
 fmt::center - the center alignment of output line
 fmt::right  - the right alignment of output line
 */
enum class fmt
{
    left,
    center,
    right
};

/**
 Converts enumerate class value to it's string label.

 @param formatKey the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const fmt formatKey)
{
    if (formatKey == fmt::left)
    {
        return std::string("Format Key: Left");
    }

    if (formatKey == fmt::center)
    {
        return std::string("Format Key: Center");
    }

    if (formatKey == fmt::right)
    {
        return std::string("Format Key: Right");
    }

    return std::string("UNKNOWN");
}

#endif /* FmtType_hpp */
