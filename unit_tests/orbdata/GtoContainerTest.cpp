//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoContainerTest.hpp"

#include <vector>

#include "GtoBlock.hpp"
#include "GtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"
#include "CheckFunctions.hpp"

TEST_F(CGtoContainerTest, DefaultConstructor)
{
    std::vector<CGtoBlock> gbloks;
    
    CGtoContainer acont;
    
    CGtoContainer bcont(gbloks);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoContainerTest, ConstructorWithMolecule)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    CGtoBlock porb(lih, bas, 1);
    
    CGtoContainer acont({sorb, porb});
    
    CGtoContainer bcont(lih, bas);

    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoContainerTest, CopyConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    CGtoContainer bcont(acont);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoContainerTest, MoveConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    CGtoContainer bcont(CGtoContainer(lih, bas));
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoContainerTest, CopyAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    CGtoContainer bcont = acont;
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoContainerTest, MoveAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    CGtoContainer bcont = CGtoContainer(lih, bas);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoContainerTest, GetMaxAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(1, acont.getMaxAngularMomentum());
}

TEST_F(CGtoContainerTest, GetAngularMomentum)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(0, acont.getAngularMomentum(0));
    
    ASSERT_EQ(1, acont.getAngularMomentum(1));
}

TEST_F(CGtoContainerTest, GetNumberOfGtoBlocks)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(2, acont.getNumberOfGtoBlocks());
}

TEST_F(CGtoContainerTest, GetMaxNumberOfPrimGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(11, acont.getMaxNumberOfPrimGtos());
}

TEST_F(CGtoContainerTest, GetMaxNumberOfContrGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(5, acont.getMaxNumberOfContrGtos());
}

TEST_F(CGtoContainerTest, getNumberOfPrimGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(11, acont.getNumberOfPrimGtos(0));
    
    ASSERT_EQ(4, acont.getNumberOfPrimGtos(1));
}

TEST_F(CGtoContainerTest, getNumberOfContrGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(5, acont.getNumberOfContrGtos(0));
    
    ASSERT_EQ(3, acont.getNumberOfContrGtos(1));
}

TEST_F(CGtoContainerTest, GetStartPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({0, 5, 6, 7, 10}, acont.getStartPositions(0));
    
    vlxtest::compare({0, 2, 3}, acont.getStartPositions(1));
}

TEST_F(CGtoContainerTest, GetEndPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({5, 6, 7, 10, 11}, acont.getEndPositions(0));
    
    vlxtest::compare({2, 3, 4}, acont.getEndPositions(1));
}

TEST_F(CGtoContainerTest, GetIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({0, 1, 2,  3, 4}, acont.getIdentifiers(0, 0));
    
    vlxtest::compare({5, 6, 7}, acont.getIdentifiers(1, 0));
    
    vlxtest::compare({8, 9, 10}, acont.getIdentifiers(1, 1));
    
    vlxtest::compare({11, 12, 13}, acont.getIdentifiers(1, 2));
}

TEST_F(CGtoContainerTest, GetExponents)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
   
    vlxtest::compare({2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00,
                      2.450300905100e+00, 7.220957185500e-01, 5.281088472100e-02,
                      2.096094879800e-02, 1.301070100000e+01, 1.962257200000e+00,
                      4.445379600000e-01, 1.219496200000e-01},
                     acont.getExponents(0));
    
    vlxtest::compare({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                      8.000000000000e-01}, acont.getExponents(1));
}

TEST_F(CGtoContainerTest, GetNormFactors)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01,
                      4.860657481700e-01, 4.362697795500e-01, 1.000000000000e+00,
                      1.000000000000e+00, 1.968215800000e-02, 1.379652400000e-01,
                      4.783193500000e-01, 1.000000000000e+00},
                      acont.getNormFactors(0));
    
    vlxtest::compare({2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                      1.000000000000e+00}, acont.getNormFactors(1));
}

TEST_F(CGtoContainerTest, GetCoordinatesX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                      acont.getCoordinatesX(0));
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, acont.getCoordinatesX(1));
}

TEST_F(CGtoContainerTest, GetCoordinatesY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                      acont.getCoordinatesY(0));
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, acont.getCoordinatesY(1));
}

TEST_F(CGtoContainerTest, GetCoordinatesZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 1.200000000000e+00, 1.200000000000e+00,
                      1.200000000000e+00, 1.200000000000e+00},
                     acont.getCoordinatesZ(0));
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      1.200000000000e+00}, acont.getCoordinatesZ(1));
}
