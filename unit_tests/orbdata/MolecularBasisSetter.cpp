//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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
getMolecularBasisForHeAtom()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getNormalizedAtomBasisForHe());

    return mbas;
}

CMolecularBasis
getMolecularBasisForH2O()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getNormalizedAtomBasisForO());

    mbas.addAtomBasis(getNormalizedAtomBasisForH());

    return mbas;
}

CMolecularBasis
getMolecularBasisForH2Se()
{
    CMolecularBasis mbas;

    mbas.setLabel({"def2-SVP"});

    mbas.addAtomBasis(getNormalizedAtomBasisForSe());

    mbas.addAtomBasis(getNormalizedAtomBasisForH());

    return mbas;
}

CMolecularBasis
getMinimalBasisForH2O()
{
    CMolecularBasis mbas;

    mbas.setLabel({"MIN-CC-PVDZ"});

    mbas.addAtomBasis(getMinimalBasisForO());

    mbas.addAtomBasis(getMinimalBasisForH());

    return mbas;
}

CMolecularBasis
getMinimalBasisForNH3CH4()
{
    CMolecularBasis mbas;

    mbas.setLabel({"MIN-CC-PVDZ"});

    mbas.addAtomBasis(getMinimalBasisForC());

    mbas.addAtomBasis(getMinimalBasisForN());

    mbas.addAtomBasis(getMinimalBasisForH());

    return mbas;
}
    
CMolecularBasis
getTestBasisForLiH()
{
    CMolecularBasis mbas;
    
    mbas.setLabel({"Test-Basis"});
    
    mbas.addAtomBasis(getTestBasisForLi());
    
    mbas.addAtomBasis(getTestBasisForH());
    
    return mbas;
}
    
CMolecularBasis
getReducedTestBasisForLiH()
{
    CMolecularBasis mbas;
    
    mbas.setLabel({"Reduced-Basis"});
    
    mbas.addAtomBasis(getAtomBasisForLi());
    
    mbas.addAtomBasis(getTestBasisForH());
    
    return mbas;
}

} // vlxbas namespace

