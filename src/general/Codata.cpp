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

// CODATA 2018
// https://physics.nist.gov/cuu/Constants/Table/allascii.txt

double
getBohrValueInAngstroms()
{
    // Bohr radius: 0.5291 772 109 03 e-10 [m]

    return 0.529177210903;
}

double
getHartreeValueInElectronVolts()
{
    // hartree-electron volt relationship: 27.211 386 245 988

    return 27.211386245988;
}

double
getHartreeValueInKiloCaloriePerMole()
{
    // hartree-joule relationship: 4.359 744 722 2071 e-18
    // Avogadro constant: 6.022 140 76 e23 [mol^-1]

    // hartree-kcal/mol relationship:
    // 4.3597447222071e-18 * 1e-3 * 6.02214076e23 / 4.184

    return 627.509474063;
}

double
getHartreeValueInInverseNanometer()
{
    // hartree-inverse meter relationship: 2.194 746 313 6320 e7 m^-1
    //                                     2.194 746 313 6320 e-2 nm^-1

    return 2.1947463136320e-2;
}

double
getHartreeValueInWavenumbers()
{
    // hartree-inverse meter relationship: 2.194 746 313 6320 e7 m^-1
    //                                     2.194 746 313 6320 e5 cm^-1

    return 2.1947463136320e+5;
}

double
getAtomicMassUnitInElectronMasses()
{
    // electron mass in u: 5.485 799 090 65 e-4
    // u in electron masses: 1.0 / 5.485 799 090 65 e-4

    return 1822.888486209;
}

double
getDipoleInDebye()
{
    // atomic unit of electric dipole mom.: 8.478 353 6255 e-30 [C m]
    // speed of light in vacuum: 299 792 458 [m s^-1]
    // ea0 = 8.4783536255e-30 * 299792458*10 * 100 [statC cm]
    // Debye = 1e-18 [statC cm]

    // 1 [a.u.] = 2.541746473 Debye

    return 2.541746473;
}

double
getRotatoryStrengthInCGS()
{
    // atomic unit of electric dipole mom.: 8.478 353 6255 e-30 [C m]
    // speed of light in vacuum: 299 792 458 [m s^-1]
    // ea0 = 8.4783536255e-30 * 299792458*10 * 100 [statC cm]

    // Bohr magneton: 927.4 010 0783 e-26 [J T^-1]
    // mu_B = 927.40100783e-23 [erg G^-1]

    // 1 [a.u.] = 2 ea0 mu_B = 471.44360760 [10**(-40) cgs unit]

    return 471.443648175;
}

double
getBoltzmannConstantInElectronVoltsPerKelvin()
{
    // Boltzmann constant: 8.617 333 262... e-5 eV K^-1

    return 8.617333262e-5;
}

double
getBoltzmannConstantInHartreePerKelvin()
{
    // Boltzmann constant: 8.617 333 262... e-5 eV K^-1
    // hartree-electron volt relationship: 27.211 386 245 988

    // Boltzmann constant in hartree K^-1: 3.166 811 563 e-6

    return 3.166811563e-6;
}

double
getExtinctionCoefficientFromBeta()
{
    // Avogadro constant: N_A = 6.022 140 76 e23 [mol^-1]
    // inverse of fine-structure constant: c = 137.035 999 084
    // Bohr radius: a0 = 0.5291 772 109 03 e-10 [m]
    // Extinction coefficient in [L mol^-1 cm^-1]

    // factor = 16 * pi * N_A * 10 * a0^2 / (ln(10) * c^2)

    return 19.603697575813566;
}

double
getFineStructureConstant()
{
    // fine-structure constant: 7.297 352 5693 e-3

    return 7.2973525693e-3;
}

}  // namespace units
