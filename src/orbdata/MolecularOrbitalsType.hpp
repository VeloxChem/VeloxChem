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

#ifndef MolecularOrbitalsType_hpp
#define MolecularOrbitalsType_hpp

#include <string>

/**
 Enumerate class molorb:

 Defines supported molecular orbital types:
 morb::rest   - the spin restricted molecular orbitals.
 morb::unrest - the spin unrestricted molecular orbitals.
 */
enum class molorb
{
    rest,
    unrest,
    restopen
};

/**
 Converts key value of molecular orbitals type to integer number.

 @param motyp the enumerate class value.
 @return the integer number.
 */
inline int32_t
to_int(const molorb motyp)
{
    return static_cast<int32_t>(motyp);
}

/**
 Converts integer key value to molecular orbitals type.

 @param keyValue the integer key value.
 @return the molecular orbital matrix type.
 */
inline molorb
to_molorb(const int32_t keyValue)
{
    if (keyValue == to_int(molorb::rest)) return molorb::rest;

    if (keyValue == to_int(molorb::unrest)) return molorb::unrest;

    if (keyValue == to_int(molorb::restopen)) return molorb::restopen;

    return molorb::rest;
}

/**
 Converts enumerate class value to it's string label.

 @param molecularOrbitals the enumerate class value.
 @return the label of enumerate class value.
 */
inline std::string
to_string(const molorb molecularOrbitals)
{
    if (molecularOrbitals == molorb::rest)
    {
        return std::string("Spin Restricted Molecular Orbitals");
    }

    if (molecularOrbitals == molorb::unrest)
    {
        return std::string("Spin Unrestricted Molecular Orbitals");
    }

    if (molecularOrbitals == molorb::restopen)
    {
        return std::string("Spin Restricted Open-Shell Molecular Orbitals");
    }

    return std::string("UNKNOWN");
}

#endif /* MolecularOrbitalsType_hpp */
