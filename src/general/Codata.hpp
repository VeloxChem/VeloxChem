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

#ifndef Codata_hpp
#define Codata_hpp

namespace units {  // units namespace

/**
 Gets Bohr value in Angstroms.

 @return the conversion factor.
 */
double bohr_in_angstrom();

/**
 Gets Hartree value in electronvolts.

 @return the conversion factor.
 */
double hartree_in_ev();

/**
 Gets Hartree value in kcal/mol.

 @return the conversion factor.
 */
double getHartreeValueInKiloCaloriePerMole();

/**
 Gets Hartree value in kJ/mol.

 @return the conversion factor.
 */
double getHartreeValueInKiloJoulePerMole();

/**
 Gets a Hartree value in inverse nanometer.

 @return the conversion factor.
 */
double getHartreeValueInInverseNanometer();

/**
 Gets Hartree value in reciprocal cm.

 @return the conversion factor.
 */
double getHartreeValueInWavenumbers();

/**
 Gets electron mass in amu.

 @return electron mass in amu.
 */
double getElectronMassInAtomicMassUnit();

/**
 Gets amu value in electron masses.

 @return the conversion factor.
 */
double getAtomicMassUnitInElectronMasses();

/**
 Gets amu value in kg.

 @return the conversion factor.
 */
double getAtomicMassUnitInKg();

/**
 Gets speed of light in vacuum in SI.

 @return the speed of light in vacuum.
 */
double getSpeedOfLightInVacuumInSI();

/**
 Gets Avogadro constant.

 @return Avogadro constant.
 */
double getAvogadroConstant();

/**
 Gets convertion factor for dipole moment (a.u. -> Debye).

 @return the conversion factor.
 */
double getDipoleInDebye();

/**
 Gets convertion factor for rotatory strength (a.u. -> 10^-40 cgs).

 @return the conversion factor.
 */
double getRotatoryStrengthInCGS();

/**
 Gets Boltzmann constant in eV/K.

 @return the conversion factor.
 */
double getBoltzmannConstantInElectronVoltsPerKelvin();

/**
 Gets Boltzmann constant in hartree/K.

 @return the conversion factor.
 */
double getBoltzmannConstantInHartreePerKelvin();

/**
 Gets factor needed for the calculation of extinction coefficient from the
 electric-dipole magnetic-dipole polarizability tensor beta.

 @return the factor.
 */
double getExtinctionCoefficientFromBeta();

/**
 Gets fine-structure constant.

 @return the fine-structure constant.
 */
double getFineStructureConstant();

}  // namespace units

#endif /* Codata_hpp */
