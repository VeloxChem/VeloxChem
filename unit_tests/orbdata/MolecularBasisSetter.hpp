//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularBasisSetter_hpp
#define MolecularBasisSetter_hpp

#include "MolecularBasis.hpp"

namespace vlxbas { // vlxbas namespace
    
CMolecularBasis getMolecularBasisEmpty();

CMolecularBasis getMolecularBasisForLiH();
    
CMolecularBasis getMolecularBasisForHeAtom();

CMolecularBasis getMolecularBasisForH2O();

CMolecularBasis getMolecularBasisForH2Se();

CMolecularBasis getMinimalBasisForH2O();

CMolecularBasis getMinimalBasisForNH3CH4();
    
CMolecularBasis getTestBasisForLiH();
    
CMolecularBasis getReducedTestBasisForLiH();

} // vlxbas namespace

#endif /* MolecularBasisSetter_hpp */
