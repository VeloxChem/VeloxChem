//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef Codata_hpp
#define Codata_hpp

namespace units { // units namespace

/**
 Gets a Bohr value in Angstroms.

 @return the conversion factor.
 */
double getBohrValueInAngstroms();

/**
 Gets a Hatree value in electronvolts.

 @return the conversion factor.
 */
double getHatreeValueInElectronVolts();

} // units namespace

#endif /* Codata_hpp */
