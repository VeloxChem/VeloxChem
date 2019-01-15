//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

    CMemBlock2D<int32_t> scntr({0, 5, 6,  7, 10,
                               5, 6, 7, 10, 11,
                               0, 0, 0,  1,  1,
                               0, 1, 2,  3,  4,
                               1, 2, 3,  4,  5},
                               5, 5);
    
    CMemBlock2D<int32_t> sidx({0, 5, 6,  7, 10,
                               5, 6, 7, 10, 11,
                               0, 1, 2,  3, 4},
                              5, 3);

    CMemBlock2D<double> sprim({2.662778551600e+02, 4.006978344700e+01, 9.055994438900e+00,
                               2.450300905100e+00, 7.220957185500e-01, 5.281088472100e-02,
                               2.096094879800e-02, 1.301070100000e+01, 1.962257200000e+00,
                               4.445379600000e-01, 1.219496200000e-01,
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
                               11, 4);
    
    CMemBlock<double> sfacts({6.492015032500e-03, 4.774786321500e-02, 2.026879611100e-01,
                              4.860657481700e-01, 4.362697795500e-01, 1.000000000000e+00,
                              1.000000000000e+00, 1.968215800000e-02, 1.379652400000e-01,
                              4.783193500000e-01, 1.000000000000e+00});

    CGtoBlock sdat(sprim, sfacts, scntr, sidx, 0);

    ASSERT_EQ(sorb, sdat);
}

TEST_F(CGtoBlockTest, ConstructorWithMoleculeAndGeneralBasisForS)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    CMemBlock2D<int32_t> scntr({0, 3, 4, 7,
                                3, 4, 7, 8,
                                0, 0, 1, 1,
                                0, 2, 3, 4,
                                2, 3, 4, 5},
                                4, 5);
    
    CMemBlock2D<int32_t> sidx({0, 3, 6,  7, 10,
                               3, 6, 7, 10, 11,
                               0, 1, 2,  3, 4},
                              5, 3);
    
    CMemBlock2D<double> sprim({1.400000000000e+00, 1.200000000000e+00, 0.300000000000e+00,
                               1.500000000000e+00, 1.301070100000e+01, 1.962257200000e+00,
                               4.445379600000e-01, 1.219496200000e-01,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 1.200000000000e+00, 1.200000000000e+00,
                               1.200000000000e+00, 1.200000000000e+00},
                              8, 4);
    
    CMemBlock<double> sfacts({0.800000000000e+00,  0.200000000000e+00, 0.300000000000e+00,
                              0.400000000000e+00, -0.200000000000e+00, 0.900000000000e+00,
                              1.000000000000e+00,  1.968215800000e-02, 1.379652400000e-01,
                              4.783193500000e-01, 1.000000000000e+00});
    
    CGtoBlock sdat(sprim, sfacts, scntr, sidx, 0);
    
    ASSERT_EQ(sorb, sdat);
}

TEST_F(CGtoBlockTest, ConstructorWithMoleculeForP)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock porb(lih, bas, 1);

    CMemBlock2D<int32_t> pcntr({0,  2,  3,
                                2,  3,  4,
                                0,  0,  1,
                                0,  1,  2,
                                1,  2,  3},
                               3, 5);
    
    CMemBlock2D<int32_t> pidx({0,  2,  3,
                               2,  3,  4,
                               5,  6,  7,
                               8,  9, 10,
                               11, 12, 13},
                              3, 5);

    CMemBlock2D<double> pprim({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                               8.000000000000e-01,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               1.200000000000e+00},
                               4, 4);
    
    CMemBlock<double> pfacts({2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                              1.000000000000e+00});

    CGtoBlock pdat(pprim, pfacts, pcntr, pidx, 1);

    ASSERT_EQ(porb, pdat);
}

TEST_F(CGtoBlockTest, ConstructorWithMoleculeAndAtomlistForP)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock aporb(lih, bas, 0, 1, 1);
    
    CMemBlock2D<int32_t> apcntr({ 0,  2,
                                 2,  3,
                                 0,  0,
                                 0,  1,
                                 1,  2},
                                2, 5);
    
    CMemBlock2D<int32_t> apidx({0,  2,
                                2,  3,
                                5,  6,
                                8,  9,
                               11, 12},
                                2, 5);
    
    CMemBlock2D<double> apprim({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                                0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                                0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                                0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00},
                               3, 4);
    
    CMemBlock<double> apfacts({2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00});
    
    CGtoBlock apdat(apprim, apfacts, apcntr, apidx, 1);
    
    ASSERT_EQ(aporb, apdat);
    
    CGtoBlock bporb(lih, bas, 1, 1, 1);
    
    CMemBlock2D<int32_t> bpcntr({ 0,
                                 1,
                                 1,
                                 0,
                                 1},
                                1, 5);
    
    CMemBlock2D<int32_t> bpidx({ 0,
                                 1,
                                 7,
                                10,
                                13},
                                1, 5);
    
    CMemBlock2D<double> bpprim({8.000000000000e-01,
                                0.000000000000e+00,
                                0.000000000000e+00,
                                1.200000000000e+00},
                               1, 4);
    
    CMemBlock<double> bpfacts(std::vector<double>({1.000000000000e+00}));
    
    CGtoBlock bpdat(bpprim, bpfacts, bpcntr, bpidx, 1);
    
    ASSERT_EQ(bporb, bpdat);
    
    CGtoBlock cporb(lih, bas, 0, 2, 1);
    
    CMemBlock2D<int32_t> cpcntr({ 0,  2,  3,
                                 2,  3,  4,
                                 0,  0,  1,
                                 0,  1,  2,
                                 1,  2,  3},
                                3, 5);
    
    CMemBlock2D<int32_t> cpidx({0,  2,  3,
                                2,  3,  4,
                                5,  6,  7,
                                8,  9, 10,
                               11, 12, 13},
                               3, 5);
    
    CMemBlock2D<double> cpprim({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                                8.000000000000e-01,
                                0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                                0.000000000000e+00,
                                0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                                0.000000000000e+00,
                                0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                                1.200000000000e+00},
                              4, 4);
    
    CMemBlock<double> cpfacts({ 2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                                1.000000000000e+00});
    
    CGtoBlock cpdat(cpprim, cpfacts, cpcntr, cpidx, 1);
    
    ASSERT_EQ(cporb, cpdat);
}

TEST_F(CGtoBlockTest, ConstructorWithMoleculeForD)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    CGtoBlock dorb(lih, bas, 2);

    CMemBlock2D<int32_t> dcntr;
    
    CMemBlock2D<int32_t> didx;

    CMemBlock2D<double> dprim;
    
    CMemBlock<double> dfacts;

    CGtoBlock ddat(dprim, dfacts, dcntr, didx, 2);

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

TEST_F(CGtoBlockTest, GetNumberOfPrimGtosWithGeneralContraction)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(8, agto.getNumberOfPrimGtos());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(2, bgto.getNumberOfPrimGtos());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfPrimGtos());
}

TEST_F(CGtoBlockTest, GetNumberOfNormFactors)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(11, agto.getNumberOfNormFactors());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(4, bgto.getNumberOfNormFactors());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfNormFactors());
}

TEST_F(CGtoBlockTest, GetNumberOfNormFactorsWithGeneralContraction)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(11, agto.getNumberOfNormFactors());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(2, bgto.getNumberOfNormFactors());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfNormFactors());
}


TEST_F(CGtoBlockTest, GetNumberOfRedContrGtos)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(5, agto.getNumberOfRedContrGtos());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(3, bgto.getNumberOfRedContrGtos());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfRedContrGtos());
}

TEST_F(CGtoBlockTest, GetNumberOfRedContrGtosWithGeneralContraction)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(4, agto.getNumberOfRedContrGtos());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(2, bgto.getNumberOfRedContrGtos());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfRedContrGtos());
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

TEST_F(CGtoBlockTest, GetNumberOfContrGtosWithGeneralContraction)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(5, agto.getNumberOfContrGtos());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(2, bgto.getNumberOfContrGtos());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getNumberOfContrGtos());
}

TEST_F(CGtoBlockTest, GetStartPositionsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0, 2, 3}, agto.getStartPositions());
}

TEST_F(CGtoBlockTest, GetStartPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto spos = agto.getStartPositions();
    
    vlxtest::compare({0, 2, 3}, spos);
    
    spos[0] = 4;
    
    vlxtest::compare({4, 2, 3}, agto.getStartPositions());
}

TEST_F(CGtoBlockTest, GetEndPositionsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({2, 3, 4}, agto.getEndPositions());
}

TEST_F(CGtoBlockTest, GetEndPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto epos = agto.getEndPositions();
    
    vlxtest::compare({2, 3, 4}, epos);
    
    epos[1] = 0;
    
    vlxtest::compare({2, 0, 4}, agto.getEndPositions());
}

TEST_F(CGtoBlockTest, GetAtomicIdentifiersConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0, 0, 1}, agto.getAtomicIdentifiers());
}

TEST_F(CGtoBlockTest, GetAtomicIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto idxatm = agto.getAtomicIdentifiers();
    
    vlxtest::compare({0, 0, 1}, idxatm);
    
    idxatm[0] = 1;
    
    idxatm[2] = 0;
    
    vlxtest::compare({1, 0, 0}, agto.getAtomicIdentifiers());
}

TEST_F(CGtoBlockTest, GetContrStartPositionsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0, 1, 2}, agto.getContrStartPositions());
}

TEST_F(CGtoBlockTest, GetContrStartPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto spos = agto.getContrStartPositions();
    
    vlxtest::compare({0, 1, 2}, spos);
    
    spos[0] = 4;
    
    vlxtest::compare({4, 1, 2}, agto.getContrStartPositions());
}

TEST_F(CGtoBlockTest, GetContrEndPositionsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({1, 2, 3}, agto.getContrEndPositions());
}

TEST_F(CGtoBlockTest, GetContrEndPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto epos = agto.getContrEndPositions();
    
    vlxtest::compare({1, 2, 3}, epos);
    
    epos[1] = 0;
    
    vlxtest::compare({1, 0, 3}, agto.getContrEndPositions());
}

TEST_F(CGtoBlockTest, GetNormFactorsStartPositionsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0, 2, 3}, agto.getNormFactorsStartPositions());
}

TEST_F(CGtoBlockTest, GetNormFactorsStartPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto spos = agto.getNormFactorsStartPositions();
    
    vlxtest::compare({0, 2, 3}, spos);
    
    spos[0] = 4;
    
    vlxtest::compare({4, 2, 3}, agto.getNormFactorsStartPositions());
}

TEST_F(CGtoBlockTest, GetNormFactorsEndPositionsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({2, 3, 4}, agto.getNormFactorsEndPositions());
}

TEST_F(CGtoBlockTest, GetNormFactorsEndPositions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto epos = agto.getNormFactorsEndPositions();
    
    vlxtest::compare({2, 3, 4}, epos);
    
    epos[1] = 0;
    
    vlxtest::compare({2, 0, 4}, agto.getNormFactorsEndPositions());
}

TEST_F(CGtoBlockTest, GetIdentifiersConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({5, 6, 7}, agto.getIdentifiers(0));
    
    vlxtest::compare({8, 9, 10}, agto.getIdentifiers(1));
    
    vlxtest::compare({11, 12, 13}, agto.getIdentifiers(2));
    
    EXPECT_TRUE(agto.getIdentifiers(3) == nullptr);
}

TEST_F(CGtoBlockTest, GetIdentifiers)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto idpx = agto.getIdentifiers(0);
    
    vlxtest::compare({5, 6, 7}, idpx);
    
    auto idpy = agto.getIdentifiers(1);
    
    vlxtest::compare({8, 9, 10}, idpy);
    
    auto idpz = agto.getIdentifiers(2);
    
    vlxtest::compare({11, 12, 13}, idpz);
    
    idpx[0] = 2;
    
    idpy[1] = 3;
    
    idpz[2] = 4;
    
    vlxtest::compare({2, 6, 7}, agto.getIdentifiers(0));
    
    vlxtest::compare({8, 3, 10}, agto.getIdentifiers(1));
    
    vlxtest::compare({11, 12, 4}, agto.getIdentifiers(2));
    
    EXPECT_TRUE(agto.getIdentifiers(3) == nullptr);
}

TEST_F(CGtoBlockTest, GetExponentsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();

    auto lih = vlxmol::getMoleculeLiH();

    const CGtoBlock agto(lih, bas, 1);

    vlxtest::compare({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                      8.000000000000e-01}, agto.getExponents());
}

TEST_F(CGtoBlockTest, GetExponents)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto pexp = agto.getExponents();
    
    vlxtest::compare({1.450000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                      8.000000000000e-01}, pexp);
    
    pexp[0] = 1.0;
    
    pexp[3] = 3.0;
    
    vlxtest::compare({1.000000000000e+00, 3.000000000000e-01, 8.200000000000e-02,
                      3.000000000000e+00}, agto.getExponents());
}

TEST_F(CGtoBlockTest, GetNormFactorsConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                      1.000000000000e+00},  agto.getNormFactors());
}

TEST_F(CGtoBlockTest, GetNormFactors)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto pnorm = agto.getNormFactors();
    
    vlxtest::compare({2.586000000000e-01, 1.000000000000e+00, 1.000000000000e+00,
                      1.000000000000e+00}, pnorm );
    
    pnorm[1] = 0.55;
    
    pnorm[3] = 0.00;
    
    vlxtest::compare({2.586000000000e-01, 5.500000000000e-01, 1.000000000000e+00,
                      0.000000000000e+00},  agto.getNormFactors());
}

TEST_F(CGtoBlockTest, getCoordinatesXConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesX());
}

TEST_F(CGtoBlockTest, getCoordinatesX)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto coordx = agto.getCoordinatesX();
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, coordx);
    
    coordx[0] = 1.00;
    
    coordx[2] = 2.00;
    
    vlxtest::compare({1.000000000000e+00, 0.000000000000e+00, 2.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesX());
}

TEST_F(CGtoBlockTest, getCoordinatesYConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);

    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesY());
}

TEST_F(CGtoBlockTest, getCoordinatesY)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto coordy = agto.getCoordinatesY();
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, coordy);
    
    coordy[1] = 3.00;
    
    coordy[2] = 2.00;
    
    vlxtest::compare({0.000000000000e+00, 3.000000000000e+00, 2.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesY());
}

TEST_F(CGtoBlockTest, getCoordinatesZConstant)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    const CGtoBlock agto(lih, bas, 1);
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      1.200000000000e+00}, agto.getCoordinatesZ());
}

TEST_F(CGtoBlockTest, getCoordinatesZ)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 1);
    
    auto coordz = agto.getCoordinatesZ();
    
    vlxtest::compare({0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                      1.200000000000e+00}, coordz);
    
    coordz[1] = 2.10;
    
    coordz[3] = 0.00;
    
    vlxtest::compare({0.000000000000e+00, 2.100000000000e+00, 0.000000000000e+00,
                      0.000000000000e+00}, agto.getCoordinatesZ());
}

TEST_F(CGtoBlockTest, Compress)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    CGtoBlock redorb(sorb);
    
    CMemBlock<double> tscreen({0.001, 0.005, 1.000, 3.000, 0.040, 2.000, 0.002,
                               2.000, 4.000, 0.001, 4.000});
    
    // compress GTOs block with 0.1 threshold
    
    auto pdim = redorb.compress(sorb, tscreen, 0.1);
    
    ASSERT_EQ(std::get<0>(pdim), 6);
    
    ASSERT_EQ(std::get<1>(pdim), 4);
    
    CMemBlock2D<int32_t> scntr({0, 2, 3, 5, 0,
                                2, 3, 5, 6, 0,
                                0, 0, 1, 1, 0,
                                0, 1, 2, 3, 0,
                                1, 2, 3, 4, 0},
                                5, 5);
    
    CMemBlock2D<int32_t> sidx({0, 2, 3, 5, 0,
                               2, 3, 5, 6, 0,
                               0, 1, 3, 4, 0},
                               5, 3);
    
    CMemBlock2D<double> sprim({9.055994438900e+00, 2.450300905100e+00, 5.281088472100e-02,
                               1.301070100000e+01, 1.962257200000e+00, 1.219496200000e-01,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               1.200000000000e+00, 1.200000000000e+00, 1.200000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                               0.000000000000e+00, 0.000000000000e+00},
                               11, 4);
    
    CMemBlock<double> cfacts({2.026879611100e-01, 4.860657481700e-01, 1.000000000000e+00,
                              1.968215800000e-02, 1.379652400000e-01, 1.000000000000e+00,
                              0.000000000000e+00, 0.000000000000e+00, 0.000000000000e+00,
                              0.000000000000e+00, 0.000000000000e+00});
    
    CGtoBlock sdat(sprim, cfacts, scntr, sidx, 0);
    
    ASSERT_EQ(redorb, sdat);
    
    // compress GTOs block with 1.0e-6 threshold
    
    pdim = redorb.compress(sorb, tscreen, 1.0e-6);
    
    ASSERT_EQ(std::get<0>(pdim), 11);
    
    ASSERT_EQ(std::get<1>(pdim), 5);
    
    ASSERT_EQ(sorb, redorb);
    
    // compress GTOs block with 50.0 threshold
    
    pdim = redorb.compress(sorb, tscreen, 50.0);
    
    ASSERT_EQ(std::get<0>(pdim), 0);
    
    ASSERT_EQ(std::get<1>(pdim), 0);
    
    sdat = CGtoBlock(CMemBlock2D<double>(11, 4), CMemBlock<double>(11),
                     CMemBlock2D<int32_t>(5, 5), CMemBlock2D<int32_t>(5, 3),
                     0);
    
    ASSERT_EQ(sdat, redorb);
}

TEST_F(CGtoBlockTest, GetMaxContractionDepth)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(5, agto.getMaxContractionDepth());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(2, bgto.getMaxContractionDepth());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getMaxContractionDepth());
}

TEST_F(CGtoBlockTest, GetMaxNumberContrFunctions)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(1, agto.getMaxNumberContrFunctions());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(1, bgto.getMaxNumberContrFunctions());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getMaxNumberContrFunctions());
}

TEST_F(CGtoBlockTest, GetMaxNumberContrFunctionsWithGeneralContraction)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock agto(lih, bas, 0);
    
    ASSERT_EQ(2, agto.getMaxNumberContrFunctions());
    
    CGtoBlock bgto(lih, bas, 1);
    
    ASSERT_EQ(1, bgto.getMaxNumberContrFunctions());
    
    CGtoBlock cgto(lih, bas, 2);
    
    ASSERT_EQ(0, cgto.getMaxNumberContrFunctions());
}

TEST_F(CGtoBlockTest, GetMaxNormFactor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 0), 6.492015032500e-03, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 1), 4.774786321500e-02, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 2), 2.026879611100e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 3), 4.860657481700e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 4), 4.362697795500e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(1, 0), 1.000000000000e+00, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(2, 0), 1.000000000000e+00, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(3, 0), 1.968215800000e-02, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(3, 1), 1.379652400000e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(3, 2), 4.783193500000e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(4, 0), 1.000000000000e+00, 1.0e-13);
}

TEST_F(CGtoBlockTest, GetMaxNormFactorWithGeneralContraction)
{
    CMolecularBasis bas = vlxbas::getGenContrBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 0), 0.800000000000e+00, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 1), 0.200000000000e+00, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(0, 2), 0.900000000000e+00, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(1, 0), 1.000000000000e+00, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(2, 0), 1.968215800000e-02, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(2, 1), 1.379652400000e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(2, 2), 4.783193500000e-01, 1.0e-13);
    
    ASSERT_NEAR(sorb.getMaxNormFactor(3, 0), 1.000000000000e+00, 1.0e-13);
}
