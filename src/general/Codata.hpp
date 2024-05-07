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

#ifndef Codata_hpp
#define Codata_hpp

namespace units {  // units namespace

// CODATA 2018
// https://physics.nist.gov/cuu/Constants/Table/allascii.txt

/**
 Gets Bohr value in Angstroms.

 @return the conversion factor.
 */
inline auto
getBohrValueInAngstroms() -> double
{
    // Bohr radius: 0.5291 772 109 03 e-10 [m]

    return 0.529177210903;
};

/**
 Gets Hartree value in electronvolts.

 @return the conversion factor.
 */
inline auto
getHartreeValueInElectronVolts() -> double
{
    // hartree-electron volt relationship: 27.211 386 245 988

    return 27.211386245988;
};

}  // namespace units

#endif /* Codata_hpp */
