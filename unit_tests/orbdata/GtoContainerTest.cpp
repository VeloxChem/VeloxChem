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
#include "VecMemBlocks.hpp"
#include "SphericalMomentum.hpp"

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

TEST_F(CGtoContainerTest, ConstructorWithMoleculeAndAtomsList)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 1, 1, 0);
    
    CGtoBlock porb(lih, bas, 1, 1, 1);
    
    CGtoContainer acont({sorb, porb});
    
    CGtoContainer bcont(lih, bas, 1, 1);
    
    ASSERT_EQ(acont, bcont);
    
    CGtoContainer ccont(lih, bas);
    
    CGtoContainer dcont(lih, bas, 0, 2);
    
    ASSERT_EQ(ccont, dcont);
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

TEST_F(CGtoContainerTest, Compress)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    CGtoContainer bcont(acont);
    
    // screening factors
    
    CVecMemBlock<double> tscreen;
    
    CMemBlock<double> sscreen({0.001, 0.005, 1.000, 3.000, 0.040, 2.000, 0.002,
                               2.000, 4.000, 0.001, 4.000});
    
    tscreen.push_back(sscreen);
    
    CMemBlock<double> pscreen({0.300, 0.000, 0.400, 0.000});
    
    tscreen.push_back(pscreen);
    
    // reduced dimensions
    
    CMemBlock2D<int32_t> redim(2, 2);
    
    bcont.compress(acont, redim, tscreen, 0.1);
    
    // check reduced dimensions
    
    vlxtest::compare({6, 2}, redim.data(0));
    
    vlxtest::compare({4, 2}, redim.data(1));
    
    // check number of GTOs blocks
    
    ASSERT_EQ(2, bcont.getNumberOfGtoBlocks());
    
    // check angular momentum
    
    ASSERT_EQ(0, bcont.getAngularMomentum(0));
    
    ASSERT_EQ(1, bcont.getAngularMomentum(1));
    
    // check original primitive GTOs dimensions
    
    ASSERT_EQ(11, bcont.getNumberOfPrimGtos(0));
    
    ASSERT_EQ(4, bcont.getNumberOfPrimGtos(1));
    
    // check original contracted GTOs dimensions
    
    ASSERT_EQ(5, bcont.getNumberOfContrGtos(0));
    
    ASSERT_EQ(3, bcont.getNumberOfContrGtos(1));
    
    // check start positions
    
    vlxtest::compare({0, 2, 3, 5, 0}, bcont.getStartPositions(0));
    
    vlxtest::compare({0, 1, 0}, bcont.getStartPositions(1));
    
    // check end positions
    
    vlxtest::compare({2, 3, 5, 6, 0}, bcont.getEndPositions(0));
    
    vlxtest::compare({1, 2, 0}, bcont.getEndPositions(1));
    
    // check identifiers
    
    vlxtest::compare({0, 1, 3, 4, 0}, bcont.getIdentifiers(0, 0));
    
    vlxtest::compare({5, 6, 0}, bcont.getIdentifiers(1, 0));
    
    vlxtest::compare({8, 9, 0}, bcont.getIdentifiers(1, 1));
    
    vlxtest::compare({11, 12, 0}, bcont.getIdentifiers(1, 2));
    
    // check exponents
    
    vlxtest::compare({9.055994438900e+00, 2.450300905100e+00, 5.281088472100e-02,
                      1.301070100000e+01, 1.962257200000e+00, 1.219496200000e-01,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                      bcont.getExponents(0));
    
    vlxtest::compare({1.450000000000e+00, 8.200000000000e-02, 0.000000000000e+00,
                      0.000000000000e+00}, bcont.getExponents(1));
    
    // check normalization factors
    
    vlxtest::compare({2.026879611100e-01, 4.860657481700e-01, 1.000000000000e+00,
                      1.968215800000e-02, 1.379652400000e-01, 1.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                     bcont.getNormFactors(0));
    
    vlxtest::compare({2.586000000000e-01, 1.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, bcont.getNormFactors(1));
    
    // check coordinates X
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                     bcont.getCoordinatesX(0));
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, bcont.getCoordinatesX(1));
    
    // check coordinates Y
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                     bcont.getCoordinatesY(0));
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, bcont.getCoordinatesY(1));
    
    // check coordinates Z
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      1.200000000000e+00, 1.200000000000e+00, 1.200000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00, 0.000000000000e+00},
                      bcont.getCoordinatesZ(0));
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, bcont.getCoordinatesZ(1));
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

TEST_F(CGtoContainerTest, GetNumberOfAtomicOrbitals)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    ASSERT_EQ(14, acont.getNumberOfAtomicOrbitals());
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

TEST_F(CGtoContainerTest, GetPrimBuffer)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    auto pbuff = acont.getPrimBuffer();
    
    ASSERT_EQ(2, pbuff.size());
    
    ASSERT_EQ(pbuff[0].size(), 11);
    
    ASSERT_EQ(pbuff[1].size(), 4);
}

TEST_F(CGtoContainerTest, GetPrimAngBuffer)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    auto abuff = acont.getPrimAngBuffer(2);
    
    ASSERT_EQ(2, abuff.size());
    
    ASSERT_EQ(abuff[0].size(0), 11);
    
    ASSERT_EQ(abuff[0].size(1), 11);
    
    ASSERT_EQ(abuff[0].blocks(), 2);
    
    ASSERT_EQ(abuff[1].size(0), 4);
    
    ASSERT_EQ(abuff[1].size(1), 4);
    
    ASSERT_EQ(abuff[1].size(2), 4);
    
    ASSERT_EQ(abuff[1].size(3), 4);
    
    ASSERT_EQ(abuff[1].size(4), 4);
    
    ASSERT_EQ(abuff[1].size(5), 4);
    
    ASSERT_EQ(abuff[1].size(6), 4);
    
    ASSERT_EQ(abuff[1].size(7), 4);
    
    ASSERT_EQ(abuff[1].blocks(), 8);
}

TEST_F(CGtoContainerTest, GetCartesianBuffer)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    auto abuff = acont.getCartesianBuffer(2);
    
    ASSERT_EQ(2, abuff.size());
    
    ASSERT_EQ(abuff[0].size(0), 5);
    
    ASSERT_EQ(abuff[0].size(1), 5);
    
    ASSERT_EQ(abuff[0].blocks(), 2);
    
    ASSERT_EQ(abuff[1].size(0), 3);
    
    ASSERT_EQ(abuff[1].size(1), 3);
    
    ASSERT_EQ(abuff[1].size(2), 3);
    
    ASSERT_EQ(abuff[1].size(3), 3);
    
    ASSERT_EQ(abuff[1].size(4), 3);
    
    ASSERT_EQ(abuff[1].size(5), 3);
    
    ASSERT_EQ(abuff[1].blocks(), 6);
}

TEST_F(CGtoContainerTest, GetSphericalBuffer)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    auto abuff = acont.getSphericalBuffer(2);
    
    ASSERT_EQ(2, abuff.size());
    
    ASSERT_EQ(abuff[0].size(0), 5);
    
    ASSERT_EQ(abuff[0].size(1), 5);
    
    ASSERT_EQ(abuff[0].blocks(), 2);
    
    ASSERT_EQ(abuff[1].size(0), 3);
    
    ASSERT_EQ(abuff[1].size(1), 3);
    
    ASSERT_EQ(abuff[1].size(2), 3);
    
    ASSERT_EQ(abuff[1].size(3), 3);
    
    ASSERT_EQ(abuff[1].size(4), 3);
    
    ASSERT_EQ(abuff[1].size(5), 3);
    
    ASSERT_EQ(abuff[1].blocks(), 6);
}

TEST_F(CGtoContainerTest, GetSphericalMomentumVector)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);
    
    auto momvec = acont.getSphericalMomentumVector();
    
    ASSERT_EQ(2, momvec.size());
    
    ASSERT_EQ(momvec[0], CSphericalMomentum(0));
    
    ASSERT_EQ(momvec[1], CSphericalMomentum(1));
}

TEST_F(CGtoContainerTest, GetGtoBlock)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoContainer acont(lih, bas);

    ASSERT_EQ(acont.getGtoBlock(0), CGtoBlock(lih, bas, 0));
    
    ASSERT_EQ(acont.getGtoBlock(1), CGtoBlock(lih, bas, 1));
}
