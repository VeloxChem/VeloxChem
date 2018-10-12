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
getMoleculeHeAtom()
{
    std::vector<double> coords({0.0, 0.0, 0.0});
    
    std::vector<double> charges({2.0});
    
    std::vector<double> masses({4.0026});
    
    std::vector<std::string> labels({{"He"}});
    
    std::vector<int32_t> idselem({2});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(0.0);

    mol.setMultiplicity(1);

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

CMolecule
getMoleculeNH3CH4()
{
    std::vector<double> coords({ -3.710, -3.702, -4.704, -4.780,
                                 -1.621, -0.819, -3.412, -0.381, -1.872,
                                  3.019,  4.942,  2.415,  2.569,
                                 -5.080, -6.698, -4.654, -3.498, -5.468,
                                 -0.037,  0.059,  1.497, -1.573,
                                  0.444, -0.465, -0.393,  0.222,  2.413});
    
    std::vector<double> charges({7.0, 1.0, 1.0, 1.0,
                                 6.0, 1.0, 1.0, 1.0, 1.0});
    
    std::vector<double> masses({14.007, 1.008, 1.008, 1.008,
                                12.011, 1.008, 1.008, 1.008, 1.008});
    
    std::vector<std::string> labels({"N", "H", "H", "H",
                                     "C", "H", "H", "H", "H"});
    
    std::vector<int32_t> idselem({7, 1, 1, 1,
                                  6, 1, 1, 1, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(0.0);

    mol.setMultiplicity(1);

    return mol;
}

CMolecule
getTestLiH2()
{
    std::vector<double> coords({0.00, 0.40, -0.20, 0.00, 0.60, 0.38,
                                0.00, 1.10, -1.50});
    
    std::vector<double> charges({3.0, 1.0, 1.0});
    
    std::vector<double> masses({7.016005, 1.007825, 1.007825});
    
    std::vector<std::string> labels({{"Li"}, {"H"}, {"H"}});
    
    std::vector<int32_t> idselem({3, 1, 1});
    
    CMolecule mol(coords, charges, masses, labels, idselem);
    
    mol.setCharge(1.0);
    
    mol.setMultiplicity(1);
    
    return mol;
}
    
} // vlxmol namespace
