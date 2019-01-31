//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AtomBasisSetter_hpp
#define AtomBasisSetter_hpp

#include <cstdint>

#include "AtomBasis.hpp"

namespace vlxbas { // vlxbas namespace
    
CAtomBasis getAtomBasisEmpty();

CAtomBasis getAtomBasisForH();

CAtomBasis getAtomBasisForLi();
    
CAtomBasis getAtomBasisForLiX();

CAtomBasis getNormalizedAtomBasisForH();

CAtomBasis getNormalizedAtomBasisForHe();
    
CAtomBasis getAtomBasisSPDForHe();

CAtomBasis getNormalizedAtomBasisForO();

CAtomBasis getNormalizedAtomBasisForSe();

CAtomBasis getMinimalBasisForH();

CAtomBasis getMinimalBasisForHe();

CAtomBasis getMinimalBasisForC();

CAtomBasis getMinimalBasisForN();

CAtomBasis getMinimalBasisForO();
    
CAtomBasis getTestBasisForH();
    
CAtomBasis getTestBasisForLi();

} // vlxbas namespace

#endif /* AtomBasisSetter_hpp */
