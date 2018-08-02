//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef MolecularBasisSetter_hpp
#define MolecularBasisSetter_hpp

#include "MolecularBasis.hpp"

namespace vlxbas { // vlxbas namespace
    
CMolecularBasis getMolecularBasisEmpty();

CMolecularBasis getMolecularBasisForLiH();
    
CMolecularBasis getTestBasisForLiH(const int32_t angularMomentumA,
                                   const int32_t angularMomentumB);
    
} // vlxbas namespace

#endif /* MolecularBasisSetter_hpp */
