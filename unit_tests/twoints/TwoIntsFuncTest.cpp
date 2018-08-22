//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "TwoIntsFuncTest.hpp"

#include "TwoIntsFunc.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CTwoIntsFuncTest, CompDistancesAQ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> raq(11, 3);

    twointsfunc::compDistancesAQ(raq, pbra, ppairs, 2);
    
    CMemBlock2D<double> refaq({0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               1.2000, 1.2000, 1.2000, 1.2000, 1.2000, 1.2000,
                               1.7400 / 2.250, 0.360 / 1.100,
                               1.2000, 0.0984 / 0.882, 0.0000},  11, 3);
    
    ASSERT_EQ(raq, refaq);
}

TEST_F(CTwoIntsFuncTest, CompFactorsForThreeCenterElectronRepulsion)
{
    CMolecularBasis bas = vlxbas::getReducedTestBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 2);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(7, 10);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 0);
    
    CMemBlock2D<double> reffacts({1.0000 / 9.0000, 1.0000 / 8.0000, 1.0000 / 8.0000,
                                  1.0000 / 7.0000, 1.0000 / 6.8000, 1.0000 / 5.8000,
                                  1.0000 / 4.6000,
                                  18.000 / 9.0000, 15.000 / 8.0000, 15.000 / 8.0000,
                                  12.000 / 7.0000, 11.400 / 6.8000, 8.4000 / 5.8000,
                                  4.8000 / 4.6000,
                                  1.0000 / 3.0000, 1.0000 / 3.0000, 1.0000 / 3.0000,
                                  1.0000 / 3.0000, 1.0000 / 3.0000, 1.0000 / 3.0000,
                                  1.0000 / 3.0000,
                                  6.0000 / 9.0000, 5.0000 / 8.0000, 5.0000 / 8.0000,
                                  4.0000 / 7.0000, 3.8000 / 6.8000, 2.8000 / 5.8000,
                                  1.6000 / 4.6000,
                                  3.0000 / 9.0000, 3.0000 / 8.0000, 3.0000 / 8.0000,
                                  3.0000 / 7.0000, 3.0000 / 6.8000, 3.0000 / 5.8000,
                                  3.0000 / 4.6000,
                                  1.0000 / 8.0000, 1.0000 / 7.0000, 1.0000 / 7.0000,
                                  1.0000 / 6.0000, 1.0000 / 5.8000, 1.0000 / 4.8000,
                                  1.0000 / 3.6000,
                                  12.000 / 8.0000, 10.000 / 7.0000, 10.000 / 7.0000,
                                  8.0000 / 6.0000, 7.6000 / 5.8000, 5.6000 / 4.8000,
                                  3.2000 / 3.6000,
                                  1.0000 / 2.0000, 1.0000 / 2.0000, 1.0000 / 2.0000,
                                  1.0000 / 2.0000, 1.0000 / 2.0000, 1.0000 / 2.0000,
                                  1.0000 / 2.0000,
                                  6.0000 / 8.0000, 5.0000 / 7.0000, 5.0000 / 7.0000,
                                  4.0000 / 6.0000, 3.8000 / 5.8000, 2.8000 / 4.8000,
                                  1.6000 / 3.6000,
                                  2.0000 / 8.0000, 2.0000 / 7.0000, 2.0000 / 7.0000,
                                  2.0000 / 6.0000, 2.0000 / 5.8000, 2.0000 / 4.8000,
                                  2.0000 / 3.6000}, 7, 10);
    
    ASSERT_EQ(facts, reffacts);
}

TEST_F(CTwoIntsFuncTest, CompCoordinatesForW)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(11, 10);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 0);
    
    CMemBlock2D<double> rw(11, 6);
    
    twointsfunc::compCoordinatesForW(rw, facts, 5, pbra, ppairs, 0);
    
    CMemBlock2D<double> refrw({0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.000 / 4.350, 0.000 / 3.200, 0.000 / 3.200, 0.000 / 2.050,
                               0.000 / 2.982, 0.000 / 1.832, 0.960 / 3.700, 0.960 / 2.550,
                               0.000 / 1.614, 0.960 / 2.332, 1.920 / 3.050,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.000 / 3.200, 0.000 / 2.050, 0.000 / 2.050, 0.000 / 0.900,
                               0.000 / 1.832, 0.000 / 0.682, 0.960 / 2.550, 0.960 / 1.400,
                               0.000 / 0.464, 0.960 / 1.182, 1.920 / 1.900}, 11, 6);
    
    ASSERT_EQ(rw, refrw);
}

TEST_F(CTwoIntsFuncTest, CompDistancesWA)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoPairsBlock ppairs(pbra, 1.0e-13);
    
    CMemBlock2D<double> facts(11, 5);
    
    twointsfunc::compFactorsForThreeCenterElectronRepulsion(facts, pbra, ppairs, 2);
    
    CMemBlock2D<double> rw(11, 3);
    
    twointsfunc::compCoordinatesForW(rw, facts, 5, pbra, ppairs, 2);
    
    CMemBlock2D<double> rwa(11, 3);
    
    twointsfunc::compDistancesWA(rwa, rw, pbra, ppairs, 2);
    
    CMemBlock2D<double> refrwa({0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000,
                                0.0000,  0.0000,  0.0000,  0.0000,
                                0.0000,  0.0000,  0.0000,  0.0000,  0.0000, 0.0000, 0.0000,
                                0.0000,  0.0000,  0.0000,  0.0000,
                               -3.4800 / 3.700, -2.1000 / 2.550, -2.1000 / 2.550, -0.7200 / 1.400,
                               -1.8384 / 2.332, -0.4584 / 1.182, -1.7400 / 3.050, -0.3600 / 1.900,
                               -0.1968 / 0.964, -0.0984 / 1.682,  0.0000 / 2.400}, 11, 3);
    
    ASSERT_EQ(rwa, refrwa);
}
