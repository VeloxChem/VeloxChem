//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef AtomBasisSetter_hpp
#define AtomBasisSetter_hpp

#include <cstdint>

#include "AtomBasis.hpp"

namespace vlxbas { // vlxbas namespace
    
CAtomBasis getAtomBasisEmpty();

CAtomBasis getAtomBasisForH();

CAtomBasis getAtomBasisForLi();

CAtomBasis getAtomBasisForO();

CAtomBasis getMinimalBasisForH();

CAtomBasis getMinimalBasisForC();

CAtomBasis getMinimalBasisForN();

CAtomBasis getMinimalBasisForO();
    
CAtomBasis getTestBasisForH();
    
CAtomBasis getTestBasisForLi();

} // vlxbas namespace

#endif /* AtomBasisSetter_hpp */
