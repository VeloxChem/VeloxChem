//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularBasisSetter_hpp
#define MolecularBasisSetter_hpp

#include "MolecularBasis.hpp"

namespace vlxbas {  // vlxbas namespace

CMolecularBasis getMolecularBasisEmpty();

CMolecularBasis getMolecularBasisForLiH();

CMolecularBasis getMolecularBasisForLiHX();

CMolecularBasis getMolecularBasisForHeAtom();

CMolecularBasis getMolecularBasisSPDForHeAtom();

CMolecularBasis getMolecularBasisForH2O();

CMolecularBasis getMolecularBasisForH2Se();

CMolecularBasis getMinimalBasisForHeAtom();

CMolecularBasis getMinimalBasisForH2O();

CMolecularBasis getMinimalBasisForNH3CH4();

CMolecularBasis getTestBasisForLiH();

CMolecularBasis getReducedTestBasisForLiH();

}  // namespace vlxbas

#endif /* MolecularBasisSetter_hpp */
