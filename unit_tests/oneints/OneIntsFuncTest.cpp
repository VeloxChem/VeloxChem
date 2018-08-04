//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OneIntsFuncTest.hpp"

#include "OneIntsFunc.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(COneIntsFuncTest, CompDistancesAB)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoBlock pket(lih, bas, 1);
    
    CMemBlock2D<double> rab(4, 3);
    
    intsfunc::compDistancesAB(rab, pbra, pket, 0);
    
    CMemBlock2D<double> refab({0.0, 0.0, 0.0,  0.0,
                               0.0, 0.0, 0.0,  0.0,
                               0.0, 0.0, 0.0, -1.2},
                               4, 3);
    
    ASSERT_EQ(rab, refab);
    
    intsfunc::compDistancesAB(rab, pbra, pket, 2);
    
    CMemBlock2D<double> refba({0.0, 0.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 0.0,
                               1.2, 1.2, 1.2, 0.0},
                              4, 3);
    
    ASSERT_EQ(rab, refba);
}

TEST_F(COneIntsFuncTest, compFactorsForOverlap)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoBlock pket(lih, bas, 1);
    
    CMemBlock2D<double> facts(4, 4);
    
    intsfunc::compFactorsForOverlap(facts, pbra, pket, 0);
    
    CMemBlock2D<double> refxz({1.0000/2.9000, 1.0000/1.7500, 1.0000/1.5320,
                               1.0000/2.2500,
                               2.1025/2.9000, 0.4350/1.7500, 0.1189/1.5320,
                               1.1600/2.2500,
                               1.0000/1.7500, 1.0000/0.6000, 1.0000/0.3820,
                               1.0000/1.1000,
                               0.4350/1.7500, 0.0900/0.6000, 0.0246/0.3820,
                               0.2400/1.1000}, 4, 4);
    
    ASSERT_EQ(facts, refxz);
}

TEST_F(COneIntsFuncTest, compFactorsForKineticEnergy)
{
    CMolecularBasis bas = vlxbas::getTestBasisForLiH(2, 2);
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 2);
    
    CGtoBlock pket(lih, bas, 2);
    
    CMemBlock2D<double> facts(6, 8);
    
    intsfunc::compFactorsForKineticEnergy(facts, pbra, pket, 0);
    
    CMemBlock2D<double> refxz({1.0000/5.6000, 1.0000/4.3000, 1.0000/4.0000,
                               1.0000/5.8000, 1.0000/4.8000, 1.0000/3.6000,
                               7.8400/5.6000, 4.2000/4.3000, 3.3600/4.0000,
                               8.4000/5.8000, 5.6000/4.8000, 2.2400/3.6000,
                               1.0000/2.8000, 1.0000/2.8000, 1.0000/2.8000,
                               1.0000/2.8000, 1.0000/2.8000, 1.0000/2.8000,
                               1.0000/2.8000, 1.0000/1.5000, 1.0000/1.2000,
                               1.0000/3.0000, 1.0000/2.0000, 1.0000/0.8000,
                               1.0000/4.3000, 1.0000/3.0000, 1.0000/2.7000,
                               1.0000/4.5000, 1.0000/3.5000, 1.0000/2.3000,
                               4.2000/4.3000, 2.2500/3.0000, 1.8000/2.7000,
                               4.5000/4.5000, 3.0000/3.5000, 1.2000/2.3000,
                               1.0000/1.5000, 1.0000/1.5000, 1.0000/1.5000,
                               1.0000/1.5000, 1.0000/1.5000, 1.0000/1.5000,
                               1.0000/2.8000, 1.0000/1.5000, 1.0000/1.2000,
                               1.0000/3.0000, 1.0000/2.0000, 1.0000/0.8000},
                              6, 8);
    
    ASSERT_EQ(facts, refxz);
}

TEST_F(COneIntsFuncTest, CompDistancesPA)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoBlock pket(lih, bas, 1);
    
    CMemBlock2D<double> rab(4, 3);
    
    intsfunc::compDistancesAB(rab, pbra, pket, 0);
    
    CMemBlock2D<double> facts(4, 4);
    
    intsfunc::compFactorsForOverlap(facts, pbra, pket, 0);
    
    CMemBlock2D<double> rpa(4, 6);
    
    intsfunc::compDistancesPA(rpa, rab, facts, 2, pbra, pket, 0);
    
    CMemBlock2D<double> refpa({0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.9600/2.2500,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.0000,
                               0.0000, 0.0000, 0.0000, 0.9600/1.1000},
                               4, 6);
    
    ASSERT_EQ(rpa, refpa);
}

TEST_F(COneIntsFuncTest, CompDistancesPB)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock pbra(lih, bas, 1);
    
    CGtoBlock pket(lih, bas, 1);
    
    CMemBlock2D<double> rab(4, 3);
    
    intsfunc::compDistancesAB(rab, pbra, pket, 0);
    
    CMemBlock2D<double> facts(4, 4);
    
    intsfunc::compFactorsForOverlap(facts, pbra, pket, 0);
    
    CMemBlock2D<double> rpb(4, 6);
    
    intsfunc::compDistancesPB(rpb, rab, facts, 2, pbra, pket, 0);
    
    CMemBlock2D<double> refpb({0.0000, 0.0000, 0.0000,  0.0000,
                               0.0000, 0.0000, 0.0000,  0.0000,
                               0.0000, 0.0000, 0.0000, -1.7400/2.2500,
                               0.0000, 0.0000, 0.0000,  0.0000,
                               0.0000, 0.0000, 0.0000,  0.0000,
                               0.0000, 0.0000, 0.0000, -0.3600/1.1000},
                              4, 6);
    
    ASSERT_EQ(rpb, refpb);
}


