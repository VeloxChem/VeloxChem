//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
