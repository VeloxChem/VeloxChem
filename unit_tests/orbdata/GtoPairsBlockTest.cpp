//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

TEST_F(CGtoPairsBlockTest, Split)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    auto ppvec = apairs.split(4);
    
    ASSERT_EQ(2u, ppvec.size());
    
    CMemBlock2D<int32_t> cpat00({ 0,  4,  6,
                                  4,  6,  8,
                                  5,  5,  5,
                                  8,  8,  8,
                                 11, 11, 11,
                                  5,  6,  7,
                                  8,  9, 10,
                                 11, 12, 13},
                                3, 8);
    
    CMemBlock2D<int32_t> cpat01({ 0,  1,  2,
                                  1,  2,  3,
                                  6,  6,  7,
                                  9,  9, 10,
                                 12, 12, 13,
                                  6,  7,  7,
                                  9, 10, 10,
                                 12, 13, 13},
                                3, 8);
    
    CMemBlock2D<double> ppat00({2.900, 1.750, 1.750, 0.600, 1.532, 0.382, 2.250, 1.100,
                                1.000 / 2.900, 1.000 / 1.750, 1.000 / 1.750, 1.000 / 0.600,
                                1.000 / 1.532, 1.000 / 0.382, 1.000 / 2.250, 1.000 / 1.100,
                                2.102500 / 2.900, 0.4350 / 1.750, 0.435 / 1.750, 0.090 / 0.600,
                                0.118900 / 1.532, 0.0246 / 0.382, 1.160 / 2.250, 0.240 / 1.100,
                                0.0754023451246241, 0.622008409788908, 0.622008409788908,
                                11.981134221083200, 0.759390538190526, 23.58466786896330,
                                0.2030763406751850, 3.525237256775160,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000 / 2.900, 0.000 / 1.750, 0.000 / 1.750, 0.000 / 0.600,
                                0.000 / 1.532, 0.000 / 0.382, 0.960 / 2.250, 0.960 / 1.100,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000 / 2.900, 0.000 / 1.750, 0.000 / 1.750, 0.000 / 0.600,
                                0.000 / 1.532, 0.000 / 0.382, 0.960 / 2.250, 0.960 / 1.100,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000 / 2.900,  0.0000 / 1.750,  0.000 / 1.750,  0.000 / 0.600,
                                0.000 / 1.532,  0.0000 / 0.382, -1.740 / 2.250, -0.360 / 1.100,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 1.200, 1.200,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000, 0.000, 0.000, -1.200, -1.200},
                               8, 22);
    
    CMemBlock2D<double> ppat01({0.164, 0.882, 1.600,
                                1.000 / 0.164, 1.000 / 0.882, 1.000 / 1.600,
                                0.006724 / 0.164, 0.0656 / 0.882, 0.640 / 1.600,
                                83.84149948663410, 6.0395990205739400, 2.751343629511100,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000 / 0.164, 0.960 / 0.882, 1.920 / 1.600,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000 / 0.164, 0.960 / 0.882, 0.000 / 1.600,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000 / 0.164, -0.0984 / 0.882,  0.000 / 1.600,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 1.200,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, 1.200, 1.200,
                                0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000,
                                0.000, -1.200, 0.000},
                               3, 22);
    
    ASSERT_EQ(ppvec[0], CGtoPairsBlock(cpat00, ppat00, 1, 1, 1.0e-13));
    
    ASSERT_EQ(ppvec[1], CGtoPairsBlock(cpat01, ppat01, 1, 1, 1.0e-13));
}


TEST_F(CGtoPairsBlockTest, Pick)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    auto p0pair = apairs.pick(0);
    
    auto p1pair = apairs.pick(1);
    
    auto p2pair = apairs.pick(2);
    
    auto p5pair = apairs.pick(5);
    
    CMemBlock2D<int32_t> cpat00({ 0, 4, 5, 8, 11, 5, 8, 11}, 1, 8);
    
    CMemBlock2D<int32_t> cpat01({ 0, 2, 5, 8, 11, 6, 9, 12}, 1, 8);
    
    CMemBlock2D<int32_t> cpat02({ 0, 2, 5, 8, 11, 7, 10, 13}, 1, 8);
    
    CMemBlock2D<int32_t> cpat05({ 0, 1, 7, 10, 13, 7, 10, 13}, 1, 8);
    
    CMemBlock2D<double> ppat00({2.900, 1.750, 1.750, 0.600,
                                1.000 / 2.900, 1.000 / 1.750, 1.000 / 1.750, 1.000 / 0.600,
                                2.102500 / 2.900, 0.4350 / 1.750, 0.435 / 1.750, 0.090 / 0.600,
                                0.0754023451246241, 0.622008409788908, 0.622008409788908,
                                11.981134221083200,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000 / 2.900, 0.000 / 1.750, 0.000 / 1.750, 0.000 / 0.600,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000 / 2.900, 0.000 / 1.750, 0.000 / 1.750, 0.000 / 0.600,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000 / 2.900,  0.0000 / 1.750,  0.000 / 1.750,  0.000 / 0.600,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000,
                                0.000, 0.000, 0.000, 0.000},
                               4, 22);
    
    CMemBlock2D<double> ppat01({ 1.532, 0.382,
                                 1.000 / 1.532, 1.000 / 0.382,
                                 0.118900 / 1.532, 0.0246 / 0.382,
                                 0.759390538190526, 23.58466786896330,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000 / 1.532, 0.000 / 0.382,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000 / 1.532, 0.000 / 0.382,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000 / 1.532,  0.0000 / 0.382,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000,
                                 0.000, 0.000},
                                2, 22);
    
    CMemBlock2D<double> ppat02({2.250, 1.100,
                                1.000 / 2.250, 1.000 / 1.100,
                                1.160 / 2.250, 0.240 / 1.100,
                                0.2030763406751850, 3.525237256775160,
                                0.000, 0.000,
                                0.000, 0.000,
                                0.960 / 2.250, 0.960 / 1.100,
                                0.000, 0.000,
                                0.000, 0.000,
                                0.960 / 2.250, 0.960 / 1.100,
                                0.000, 0.000,
                                0.000, 0.000,
                               -1.740 / 2.250, -0.360 / 1.100,
                                0.000, 0.000,
                                0.000, 0.000,
                                0.000, 0.000,
                                0.000, 0.000,
                                0.000, 0.000,
                                1.200, 1.200,
                                0.000, 0.000,
                                0.000, 0.000,
                               -1.200, -1.200},
                               2, 22);
    
    CMemBlock2D<double> ppat05({ 1.600,
                                 1.000 / 1.600,
                                 0.640 / 1.600,
                                 2.751343629511100,
                                 0.000,
                                 0.000,
                                 1.920 / 1.600,
                                 0.000,
                                 0.000,
                                 0.000 / 1.600,
                                 0.000,
                                 0.000,
                                 0.000 / 1.600,
                                 0.000,
                                 0.000,
                                 1.200,
                                 0.000,
                                 0.000,
                                 1.200,
                                 0.000,
                                 0.000,
                                 0.000},
                                1, 22);
    
    ASSERT_EQ(p0pair, CGtoPairsBlock(cpat00, ppat00, 1, 1, 1.0e-13));
    
    ASSERT_EQ(p1pair, CGtoPairsBlock(cpat01, ppat01, 1, 1, 1.0e-13));
    
    ASSERT_EQ(p2pair, CGtoPairsBlock(cpat02, ppat02, 1, 1, 1.0e-13));
   
    ASSERT_EQ(p5pair, CGtoPairsBlock(cpat05, ppat05, 1, 1, 1.0e-13));
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

TEST_F(CGtoPairsBlockTest, GetDistancesPAZ)
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

TEST_F(CGtoPairsBlockTest, GetDistancesPBZ)
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

TEST_F(CGtoPairsBlockTest, GetCoordinatesAX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getCoordinatesAX());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesAY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getCoordinatesAY());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesAZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2},
                     apairs.getCoordinatesAZ());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesBX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getCoordinatesBX());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesBY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getCoordinatesBY());
}

TEST_F(CGtoPairsBlockTest, GetCoordinatesBZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.2, 1.2, 0.0, 1.2, 1.2},
                     apairs.getCoordinatesBZ());
}

TEST_F(CGtoPairsBlockTest, GetDistancesABX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getDistancesABX());
}

TEST_F(CGtoPairsBlockTest, GetDistancesABY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare(std::vector<double>(11, 0.0), apairs.getDistancesABY());
}

TEST_F(CGtoPairsBlockTest, GetDistancesABZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    vlxtest::compare({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -1.2, -1.2, 0.0, -1.2, 0.0},
                     apairs.getDistancesABZ());
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

TEST_F(CGtoPairsBlockTest, GetDistancesAB)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    CMemBlock2D<double> rab({0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
                             0.0, 0.0,  0.0, 0.0,  0.0, 0.0,
                             0.0, 0.0, -1.2, 0.0, -1.2, 0.0},
                             6, 3);
    
    ASSERT_EQ(apairs.getDistancesAB(), rab); 
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

TEST_F(CGtoPairsBlockTest, GetMaxContractionDepth)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    ASSERT_EQ(4, apairs.getMaxContractionDepth());
}

TEST_F(CGtoPairsBlockTest, GetNumberOfPrimPairs)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, 1.0e-13);
    
    ASSERT_EQ(4, apairs.getNumberOfPrimPairs(0));
    
    ASSERT_EQ(6, apairs.getNumberOfPrimPairs(1));
    
    ASSERT_EQ(8, apairs.getNumberOfPrimPairs(2));
    
    ASSERT_EQ(9, apairs.getNumberOfPrimPairs(3));
    
    ASSERT_EQ(10, apairs.getNumberOfPrimPairs(4));
    
    ASSERT_EQ(11, apairs.getNumberOfPrimPairs(5));
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

TEST_F(CGtoPairsBlockTest, GetPairType)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    CGtoPairsBlock apairs(agto, bgto, 1.0e-13);
    
    ASSERT_EQ(std::string("(S,P)"), apairs.getPairType());
    
    CGtoPairsBlock bpairs(bgto, agto, 1.0e-13);
    
    ASSERT_EQ(std::string("(P,S)"), bpairs.getPairType());
    
    CGtoPairsBlock cpairs(agto, 1.0e-13);
    
    ASSERT_EQ(std::string("(S,S)"), cpairs.getPairType());
    
    CGtoPairsBlock dpairs(bgto, 1.0e-13);
    
    ASSERT_EQ(std::string("(P,P)"), dpairs.getPairType());
}

