//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "MolecularBasisSetter.hpp"

#include "AtomBasisSetter.hpp"

namespace vlxbas { // vlxbas namespace
    
CMolecularBasis
getMolecularBasisEmpty()
{
    CMolecularBasis mbas;

    mbas.setMaxAngularMomentum(-1);

    return mbas;
}

CMolecularBasis
getMolecularBasisForLiH()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getAtomBasisForLi());

    mbas.addAtomBasis(getAtomBasisForH());

    return mbas;
}

CMolecularBasis
getTestBasisForLiH(const int32_t angularMomentumA,
                   const int32_t angularMomentumB)
{
    CMolecularBasis mbas;
    
    mbas.setLabel({"Test-Basis"});
    
    mbas.addAtomBasis(getTestBasisForLi(angularMomentumA));
    
    mbas.addAtomBasis(getTestBasisForH(angularMomentumB));
    
    return mbas;
}
    
    
} // vlxbas namespace

