//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem X developers. All rights reserved.

#include "MoleculeSetter.hpp"

namespace vlxmol { // vlxmol namespace

CMolecule
getMoleculeEmpty()
{
    return CMolecule();
}

CMolecule
getMoleculeLiH()
{
    std::vector<double> coords({0.0, 0.0, 0.0, 0.0, 0.0, 1.2});
    
    std::vector<double> charges({3.0, 1.0});
    
    std::vector<double> masses({7.016005, 1.007825});
    
    std::vector<std::string> labels({{"Li"}, {"H"}});
    
    std::vector<int32_t> idselem({3, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);

    mol.setCharge(0.0);

    mol.setMultiplicity(1);

    return mol;
}

CMolecule
getTestLiH()
{
    std::vector<double> coords({0.0, 0.4, 0.0, 0.6, 0.0, 1.1});
    
    std::vector<double> charges({3.0, 1.0});
    
    std::vector<double> masses({7.016005, 1.007825});
    
    std::vector<std::string> labels({{"Li"}, {"H"}});
    
    std::vector<int32_t> idselem({3, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(0.0);
    
    mol.setMultiplicity(1);
    
    return mol;
}
    
CMolecule
getMoleculeLiHCation()
{
    std::vector<double> coords({0.0, 0.0, 0.0, 0.0, 0.0, 1.2});
    
    std::vector<double> charges({3.0, 1.0});
    
    std::vector<double> masses({7.016005, 1.007825});
    
    std::vector<std::string> labels({{"Li"}, {"H"}});
    
    std::vector<int32_t> idselem({3, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(1.0);

    mol.setMultiplicity(2);

    return mol;
}

CMolecule
getMoleculeH2O()
{
    std::vector<double> coords({0.0, 0.0,  0.0,
                                0.0, 1.4, -1.4,
                                0.0, 1.1,  1.1});
    
    std::vector<double> charges({8.0, 1.0, 1.0});
    
    std::vector<double> masses({15.994915, 1.007825, 1.007825});
    
    std::vector<std::string> labels({{"O"}, {"H"}, {"H"}});
    
    std::vector<int32_t> idselem({8, 1, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(0.0);

    mol.setMultiplicity(1);

    return mol;
}

CMolecule
getMoleculeH2ODimer()
{
    std::vector<double> coords({0.0, 0.0,  0.0, 3.0, 3.0,  3.0,
                                0.0, 1.4, -1.4, 0.0, 1.4, -1.4,
                                0.0, 1.1,  1.1, 0.0, 1.1,  1.1});
    
    std::vector<double> charges({8.0, 1.0, 1.0, 8.0, 1.0, 1.0});
    
    std::vector<double> masses({15.994915, 1.007825, 1.007825,
                                15.994915, 1.007825, 1.007825});
    
    std::vector<std::string> labels({"O", "H", "H", "O", "H", "H"});
    
    std::vector<int32_t> idselem({8, 1, 1, 8, 1, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(0.0);

    mol.setMultiplicity(1);

    return mol;
}

} // vlxmol namespace
