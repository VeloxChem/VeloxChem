//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoPairsBlockTest.hpp"

#include "GtoPairsBlock.hpp"
#include "GtoBlock.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CGtoPairsBlockTest, ConstructorForSingleGtoBlock)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoPairsBlock apairs(agto, agto, 1.0e-13);
    
    CGtoPairsBlock bpairs(agto, 1.0e-13);
    
    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CGtoPairsBlockTest, CopyConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    CGtoPairsBlock bpairs(apairs);
    
    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CGtoPairsBlockTest, MoveConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    CGtoPairsBlock bpairs(CGtoPairsBlock(agto, bgto, 1.0e-13));
    
    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CGtoPairsBlockTest, CopyAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    CGtoPairsBlock bpairs = apairs;
    
    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CGtoPairsBlockTest, MoveAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    CGtoPairsBlock bpairs = CGtoPairsBlock(agto, bgto, 1.0e-13);
    
    ASSERT_EQ(apairs, bpairs);
}

TEST_F(CGtoPairsBlockTest, GetBraAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    ASSERT_EQ(0, apairs.getBraAngularMomentum());
    
    CGtoPairsBlock bpairs(bgto, 1.0e-13);
    
    ASSERT_EQ(1, bpairs.getBraAngularMomentum());
}

TEST_F(CGtoPairsBlockTest, GetKetAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    ASSERT_EQ(1, apairs.getKetAngularMomentum());
    
    CGtoPairsBlock bpairs(bgto, 1.0e-13);
    
    ASSERT_EQ(1, bpairs.getKetAngularMomentum());
}

TEST_F(CGtoPairsBlockTest, Empty)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    ASSERT_FALSE(apairs.empty());
    
    CGtoPairsBlock bpairs(bgto, 1.0e-13);
    
    ASSERT_FALSE(bpairs.empty());
    
    CGtoPairsBlock cpairs(bgto, 1.0e5);
    
    ASSERT_TRUE(cpairs.empty());
}

TEST_F(CGtoPairsBlockTest, getFactorsXi)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({2.900, 1.750, 1.750, 0.600, 1.532, 0.382, 2.250,
                      1.100, 0.164, 0.882, 1.600}, apairs.getFactorsXi());
}

TEST_F(CGtoPairsBlockTest, getFactorsOneOverXi)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({1.000 / 2.900, 1.000 / 1.750, 1.000 / 1.750, 1.000 / 0.600,
                      1.000 / 1.532, 1.000 / 0.382, 1.000 / 2.250, 1.000 / 1.100,
                      1.000 / 0.164, 1.000 / 0.882, 1.000 / 1.600},
                     apairs.getFactorsOneOverXi());
}

TEST_F(CGtoPairsBlockTest, getFactorsZeta)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({2.102500 / 2.900, 0.4350 / 1.750, 0.435 / 1.750, 0.090 / 0.600,
                      0.118900 / 1.532, 0.0246 / 0.382, 1.160 / 2.250, 0.240 / 1.100,
                      0.006724 / 0.164, 0.0656 / 0.882, 0.640 / 1.600},
                     apairs.getFactorsZeta());
}

TEST_F(CGtoPairsBlockTest, getOverlaps)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.0754023451246241, 0.622008409788908, 0.622008409788908,
                      11.981134221083200, 0.759390538190526, 23.58466786896330,
                      0.2030763406751850, 3.525237256775160, 83.84149948663410,
                      6.0395990205739400, 2.751343629511100},
                     apairs.getOverlaps());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesPX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getCoordinatesPX());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesPY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getCoordinatesPY());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesPZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.000 / 2.900, 0.000 / 1.750, 0.000 / 1.750, 0.000 / 0.600,
                      0.000 / 1.532, 0.000 / 0.382, 0.960 / 2.250, 0.960 / 1.100,
                      0.000 / 0.164, 0.960 / 0.882, 1.920 / 1.600},
                     apairs.getCoordinatesPZ());
}

TEST_F(CGtoPairsBlockTest, GetDistancesPAX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getDistancesPAX());
}

TEST_F(CGtoPairsBlockTest, GetDistancesPAY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getDistancesPAY());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesPAZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.000 / 2.900, 0.000 / 1.750, 0.000 / 1.750, 0.000 / 0.600,
                      0.000 / 1.532, 0.000 / 0.382, 0.960 / 2.250, 0.960 / 1.100,
                      0.000 / 0.164, 0.960 / 0.882, 0.000 / 1.600},
                     apairs.getDistancesPAZ());
}

TEST_F(CGtoPairsBlockTest, GetDistancesPBX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getDistancesPBX());
}

TEST_F(CGtoPairsBlockTest, GetDistancesPBY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getDistancesPBY());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesPBZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.000 / 2.900,  0.0000 / 1.750,  0.000 / 1.750,  0.000 / 0.600,
                      0.000 / 1.532,  0.0000 / 0.382, -1.740 / 2.250, -0.360 / 1.100,
                      0.000 / 0.164, -0.0984 / 0.882,  0.000 / 1.600},
                     apairs.getDistancesPBZ());
}

TEST_F(CGtoPairsBlockTest, GetStartPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0, 4, 6, 8, 9, 10}, apairs.getStartPositions());
}

TEST_F(CGtoPairsBlockTest, GetEndPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({4, 6, 8, 9, 10, 11}, apairs.getEndPositions());
}

TEST_F(CGtoPairsBlockTest, GetBraIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({5, 5, 5, 6, 6, 7}, apairs.getBraIdentifiers(0));
    
    vlxtest::compare({8, 8, 8, 9, 9, 10}, apairs.getBraIdentifiers(1));
    
    vlxtest::compare({11, 11, 11, 12, 12, 13}, apairs.getBraIdentifiers(2));
}

TEST_F(CGtoPairsBlockTest, GetKetIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({5, 6, 7, 6, 7, 7}, apairs.getKetIdentifiers(0));
    
    vlxtest::compare({8, 9, 10, 9, 10, 10}, apairs.getKetIdentifiers(1));
    
    vlxtest::compare({11, 12, 13, 12, 13, 13}, apairs.getKetIdentifiers(2));
}

TEST_F(CGtoPairsBlockTest, GetNumberOfOriginalPrimPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    CGtoPairsBlock bpairs(agto, 10.0);
    
    ASSERT_EQ(11, apairs.getNumberOfOriginalPrimPairs());
    
    ASSERT_EQ(11, bpairs.getNumberOfOriginalPrimPairs());
}

TEST_F(CGtoPairsBlockTest, GetNumberOfScreenedPrimPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    CGtoPairsBlock bpairs(agto, 10.0);
    
    ASSERT_EQ(11, apairs.getNumberOfScreenedPrimPairs());
    
    ASSERT_EQ(3, bpairs.getNumberOfScreenedPrimPairs());
}

TEST_F(CGtoPairsBlockTest, GetNumberOfOriginalContrPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    CGtoPairsBlock bpairs(agto, 10.0);
    
    ASSERT_EQ(6, apairs.getNumberOfOriginalContrPairs());
    
    ASSERT_EQ(6, bpairs.getNumberOfOriginalContrPairs());
}

TEST_F(CGtoPairsBlockTest, GetNumberOfScreenedContrPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    CGtoPairsBlock bpairs(agto, 10.0);
    
    ASSERT_EQ(6, apairs.getNumberOfScreenedContrPairs());
    
    ASSERT_EQ(3, bpairs.getNumberOfScreenedContrPairs());
}


