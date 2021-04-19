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

#ifndef TwoIntsDistType_hpp
#define TwoIntsDistType_hpp

#include <string>

/**
 Enumerate class dist2e:
 
 Defines supported two electron integrals distribution keys:
 dist2e::batch   - the batch with natural order of data
 dist2e::fock    - the Fock matrix
 dist2e::qvalues - the Q values vectors for bra and ket sides
 */
enum class dist2e
{
    batch,
    fock,
    qvalues 
};

/**
 Converts enumerate class value to it's string label.
 
 @param distPattern the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string to_string(const dist2e distPattern)
{
    if (distPattern == dist2e::batch)
    {
        return std::string("Raw Integrals Batch");
    }
    
    if (distPattern == dist2e::fock)
    {
        return std::string("Fock Matrix");
    }
    
    if (distPattern == dist2e::qvalues)
    {
        return std::string("Q Values Vector");
    }
    
    return std::string("UNKNOWN");
}

#endif /* TwoIntsDistType_hpp */
