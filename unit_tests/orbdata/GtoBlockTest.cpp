//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoBlockTest.hpp"

#include "GtoBlock.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CGtoBlockTest, ConstructorWithMoleculeForS)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock sorb(lih, bas, 0);

    CMemBlock2D<int32_t> sidx({0, 5, 6,  7, 10,
                               5, 6, 7, 10, 11,
                               0, 1, 2,  3, 4},
                               5, 3);

    CMemBlock2D<double> sprim({2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00,
                               2.450300905100e+00, 7.220957185500e-01, 5.281088472100e-02,
                               2.096094879800e-02, 1.301070100000e+01, 1.962257200000e+00,
                               4.445379600000e-01, 1.219496200000e-01,
                               6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01,
                               4.860657481700e-01, 4.362697795500e-01, 1.000000000000e+00,
                               1.000000000000e+00, 1.968215800000e-02, 1.379652400000e-01,
                               4.783193500000e-01, 1.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 1.200000000000e+00, 1.200000000000e+00,
                               1.200000000000e+00, 1.200000000000e+00},
                               11, 5);

    CGtoBlock sdat(sprim, sidx, 0);

    ASSERT_EQ(sorb, sdat);
}

TEST_F(CGtoBlockTest, ConstructorWithMoleculeForP)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock porb(lih, bas, 1);

    CMemBlock2D<int32_t> pidx({ 0,  2,  3,
                                2,  3,  4,
                                5,  6,  7,
                                8,  9, 10,
                               11, 12, 13},
                               3, 5);

    CMemBlock2D<double> pprim({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                               8.000000000000e-01,
                               2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                               1.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               1.200000000000e+00},
                               4, 5);

    CGtoBlock pdat(pprim, pidx, 1);

    ASSERT_EQ(porb, pdat);
}

TEST_F(CGtoBlockTest, ConstructorWithMoleculeForD)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock dorb(lih, bas, 2);

    CMemBlock2D<int32_t> didx;

    CMemBlock2D<double> dprim;

    CGtoBlock ddat(dprim, didx, 2);

    ASSERT_EQ(dorb, ddat);
}

TEST_F(CGtoBlockTest, CopyConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock agto(lih, bas, 1);

    CGtoBlock bgto(agto);

    ASSERT_EQ(agto, bgto);
}

TEST_F(CGtoBlockTest, MoveConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock agto(lih, bas, 1);

    CGtoBlock bgto(CGtoBlock(lih, bas, 1));

    ASSERT_EQ(agto, bgto);
}

TEST_F(CGtoBlockTest, CopyAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock agto(lih, bas, 1);

    CGtoBlock bgto = agto;

    ASSERT_EQ(agto, bgto);
}

TEST_F(CGtoBlockTest, MoveAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock agto(lih, bas, 1); 

    CGtoBlock bgto = CGtoBlock(lih, bas, 1);

    ASSERT_EQ(agto, bgto);
}

TEST_F(CGtoBlockTest, GetAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(0, agto.getAngularMomentum());
    
    ASSERT_EQ(1, bgto.getAngularMomentum());
}

TEST_F(CGtoBlockTest, Empty)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_FALSE(agto.empty());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_FALSE(bgto.empty());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_TRUE(cgto.empty());
}

TEST_F(CGtoBlockTest, GetNumberOfPrimGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(11, agto.getNumberOfPrimGtos());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(4, bgto.getNumberOfPrimGtos());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfPrimGtos());
}

TEST_F(CGtoBlockTest, GetNumberOfContrGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(5, agto.getNumberOfContrGtos());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(3, bgto.getNumberOfContrGtos());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfContrGtos());
}

TEST_F(CGtoBlockTest, GetStartPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0, 2, 3}, agto.getStartPositions());
}

TEST_F(CGtoBlockTest, GetEndPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({2, 3, 4}, agto.getEndPositions());
}

TEST_F(CGtoBlockTest, GetIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({5, 6, 7}, agto.getIdentifiers(0));
    
    vlxtest::compare({8, 9, 10}, agto.getIdentifiers(1));
    
    vlxtest::compare({11, 12, 13}, agto.getIdentifiers(2));
    
    EXPECT_TRUE(agto.getIdentifiers(3) == nullptr);
}

TEST_F(CGtoBlockTest, GetExponents)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock agto(lih, bas, 1);

    vlxtest::compare({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                      8.000000000000e-01}, agto.getExponents());
}

TEST_F(CGtoBlockTest, GetNormFactors)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                      1.000000000000e+00},  agto.getNormFactors());
}

TEST_F(CGtoBlockTest, getCoordinatesX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesX());
}

TEST_F(CGtoBlockTest, getCoordinatesY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesY());
}

TEST_F(CGtoBlockTest, getCoordinatesZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      1.200000000000e+00}, agto.getCoordinatesZ());
}

