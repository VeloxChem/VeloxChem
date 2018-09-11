//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem X developers. All rights reserved.

#ifndef MoleculeSetter_hpp
#define MoleculeSetter_hpp

#include "Molecule.hpp"

namespace vlxmol { // vlxmol namespace
    
    CMolecule getMoleculeEmpty();
    
    CMolecule getMoleculeLiH();
    
    CMolecule getTestLiH();
    
    CMolecule getMoleculeLiHCation();
    
    CMolecule getMoleculeH2O();

    CMolecule getMoleculeH2ODimer();

    CMolecule getMoleculeNH3CH4();
    
} // vlxmol namespace

#endif /* MoleculeSetter_hpp */
