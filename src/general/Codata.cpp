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

#include "Codata.hpp"

namespace units {  // units namespace

// CODATA 2010
// https://physics.nist.gov/cuu/Constants/allascii_2010.txt

double
getBohrValueInAngstroms()
{
    // Bohr radius: 0.529 177 210 92 e-10 [m]

    return 0.52917721092;
}

double
getHartreeValueInElectronVolts()
{
    // hartree-electron volt relationship: 27.211 385 05

    return 27.21138505;
}

double
getHartreeValueInKiloCaloriePerMole()
{
    // hartree-joule relationship: 4.359 744 34 e-18
    // Avogadro constant: 6.022 141 29 e23 [mol^-1]

    // hartree-kcal/mol relationship:
    // 4.35974434e-18 * 1e-3 * 6.02214129e23 / 4.184

    return 627.50947428;
}

double
getHartreeValueInInverseNanometer()
{
    // hartree-inverse meter relationship: 2.194 746 313 708 e7 m^-1
    //                                     2.194 746 313 708 e-2 nm^-1

    return 2.194746313708e-2;
}

double
getDipoleInDebye()
{
    // atomic unit of electric dipole mom.: 8.478 353 26 e-30 [C m]
    // ea0 = 8.47835326e-30 * 299792458*10 * 100 [statC cm]
    // Debye = 1e-18 [statC cm]

    // 1 [a.u.] = 2.54174636 Debye

    return 2.54174636;
}

double
getRotatoryStrengthInCGS()
{
    // atomic unit of electric dipole mom.: 8.478 353 26 e-30 [C m]
    // speed of light in vacuum: 299 792 458 [m s^-1]
    // ea0 = 8.47835326e-30 * 299792458*10 * 100 [statC cm]

    // Bohr magneton: 927.400 968 e-26 [J T^-1]
    // mu_B = 927.400968e-23 [erg G^-1]

    // 1 [a.u.] = 2 ea0 mu_B = 471.44360760 [10**(-40) cgs unit]

    return 471.44360760;
}

}  // namespace units
