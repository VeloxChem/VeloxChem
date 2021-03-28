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

#ifndef OneIntsDistType_hpp
#define OneIntsDistType_hpp

#include <string>

/**
 Enumerate class dist1e:

 Defines supported one electron integrals distribution keys:
 dist1e::symsq  - the symmetric square matrix
 dist1e::antisq - the antisymmetric square matrix
 dist1e::rect   - the general rectangular matrix
 dist1e::batch  - the batch with natural order of data
 */
enum class dist1e
{
    symsq,
    antisq,
    rect,
    batch
};

/**
 Converts enumerate class value to it's string label.

 @param distPattern the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const dist1e distPattern)
{
    if (distPattern == dist1e::symsq)
    {
        return std::string("Symmetric Square Matrix");
    }

    if (distPattern == dist1e::antisq)
    {
        return std::string("Anti-symmetric Square Matrix");
    }

    if (distPattern == dist1e::rect)
    {
        return std::string("Rectangular Matrix");
    }

    if (distPattern == dist1e::batch)
    {
        return std::string("Raw Integrals Batch");
    }

    return std::string("UNKNOWN");
}

#endif /* OneIntsDistType_hpp */
