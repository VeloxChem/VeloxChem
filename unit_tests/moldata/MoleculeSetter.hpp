//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef MoleculeSetter_hpp
#define MoleculeSetter_hpp

#include "Molecule.hpp"

namespace vlxmol {  // vlxmol namespace

CMolecule getMoleculeEmpty();

CMolecule getMoleculeLiH();

CMolecule getTestLiH();

CMolecule getMoleculeLiHCation();

CMolecule getMoleculeHeAtom();

CMolecule getMoleculeH2O();

CMolecule getMoleculeH2Se();

CMolecule getMoleculeH2ODimer();

CMolecule getMoleculeNH3CH4();

CMolecule getTestLiH2();

}  // namespace vlxmol

#endif /* MoleculeSetter_hpp */
