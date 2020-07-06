//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
