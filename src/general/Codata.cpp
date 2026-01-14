//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Codata.hpp"

namespace units {  // units namespace

// CODATA 2018
// https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2018.txt

double
bohr_in_angstrom()
{
    // Bohr radius: 0.5291 772 109 03 e-10 [m]

    return 0.529177210903;
}

double
hartree_in_ev()
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
getHartreeValueInKiloJoulePerMole()
{
    // hartree-joule relationship: 4.359 744 722 2071 e-18
    // Avogadro constant: 6.022 140 76 e23 [mol^-1]

    // hartree-kcal/mol relationship:
    // 4.3597447222071e-18 * 1e-3 * 6.02214076e23

    return 2625.4996394798;
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
getElectronMassInAtomicMassUnit()
{
    // electron mass in u: 5.485 799 090 65 e-4

    return 5.48579909065e-4;
}

double
getAtomicMassUnitInElectronMasses()
{
    // electron mass in u: 5.485 799 090 65 e-4
    // u in electron masses: 1.0 / 5.485 799 090 65 e-4

    return 1.0 / getElectronMassInAtomicMassUnit();
}

double
getAtomicMassUnitInKg()
{
    // atomic mass constant 1.660 539 066 60 e-27 kg

    return 1.66053906660e-27;
}

double
getSpeedOfLightInVacuumInSI()
{
    // speed of light in vacuum: 299 792 458 [m s^-1]

    return 299792458.0;
}

double
getAvogadroConstant()
{
    // Avogadro constant: N_A = 6.022 140 76 e23 [mol^-1]

    return 6.02214076e+23;
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
