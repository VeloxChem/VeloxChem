//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef AtomBasisSetter_hpp
#define AtomBasisSetter_hpp

#include <cstdint>

#include "AtomBasis.hpp"

namespace vlxbas {  // vlxbas namespace

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

}  // namespace vlxbas

#endif /* AtomBasisSetter_hpp */
