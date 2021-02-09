//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef Codata_hpp
#define Codata_hpp

namespace units {  // units namespace

/**
 Gets a Bohr value in Angstroms.

 @return the conversion factor.
 */
double getBohrValueInAngstroms();

/**
 Gets a Hartree value in electronvolts.

 @return the conversion factor.
 */
double getHartreeValueInElectronVolts();

/**
 Gets a Hartree value in kcal/mol.

 @return the conversion factor.
 */
double getHartreeValueInKiloCaloriePerMole();

/**
 Gets a Hartree value in inverse nanometer.

 @return the conversion factor.
 */
double getHartreeValueInInverseNanometer();

/**
 Gets convertion factor for dipole moment (a.u.->Debye)

 @return the conversion factor.
 */
double getDipoleInDebye();

/**
 Gets convertion factor for rotatory strength (a.u.->cgs)

 @return the conversion factor.
 */
double getRotatoryStrengthInCGS();

}  // namespace units

#endif /* Codata_hpp */
