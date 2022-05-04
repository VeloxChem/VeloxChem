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

#ifndef Codata_hpp
#define Codata_hpp

namespace units {  // units namespace

/**
 Gets Bohr value in Angstroms.

 @return the conversion factor.
 */
double getBohrValueInAngstroms();

/**
 Gets Hartree value in electronvolts.

 @return the conversion factor.
 */
double getHartreeValueInElectronVolts();

/**
 Gets Hartree value in kcal/mol.

 @return the conversion factor.
 */
double getHartreeValueInKiloCaloriePerMole();

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
