//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "GtoPairsContainerTest.hpp"

#include "GtoPairsBlock.hpp"
#include "GtoPairsContainer.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"

TEST_F(CGtoPairsContainerTest, DefaultConstructor)
{
    std::vector<CGtoPairsBlock> ppbloks;
    
    CGtoPairsContainer acont;
    
    CGtoPairsContainer bcont(ppbloks);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, ConstructorWithMolecule)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock sorb(lih, bas, 0);
    
    CGtoBlock porb(lih, bas, 1);
    
    CGtoPairsBlock bss(sorb, 1.0e-13);
    
    CGtoPairsBlock bsp(sorb, porb, 1.0e-13);
    
    CGtoPairsBlock bpp(porb, 1.0e-13);
    
    CGtoPairsContainer acont({bss, bsp, bpp});
    
    CGtoPairsContainer bcont(lih, bas, 1.0e-13);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, CopyConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont(acont);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, MoveConstructor)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont(CGtoPairsContainer(lih, bas, 1.0e-13));
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, CopyAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont = acont;
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, MoveAssignment)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont = CGtoPairsContainer(lih, bas, 1.0e-13);
    
    ASSERT_EQ(acont, bcont);
}

TEST_F(CGtoPairsContainerTest, Split)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoBlock porb(lih, bas, 1);
    
    CGtoPairsBlock bpp(porb, 1.0e-13);
    
    CGtoPairsContainer acont({bpp});
    
    auto bcont = acont.split(4);
    
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
    
    auto pp0 = CGtoPairsBlock(cpat00, ppat00, 1, 1, 1.0e-13);
    
    auto pp1 = CGtoPairsBlock(cpat01, ppat01, 1, 1, 1.0e-13);
    
    ASSERT_EQ(bcont, CGtoPairsContainer({pp0, pp1}));
    
    ASSERT_EQ(acont, acont.split(1000));
}

TEST_F(CGtoPairsContainerTest, GetNumberOfGtoPairsBlocks)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoPairsContainer bcont(lih, bas, 1.0e5);
    
    ASSERT_EQ(3, acont.getNumberOfGtoPairsBlocks());
    
    ASSERT_EQ(0, bcont.getNumberOfGtoPairsBlocks());
}

TEST_F(CGtoPairsContainerTest, GetGtoPairsBlock)
{
    CMolecularBasis bas = vlxbas::getMolecularBasisForLiH();
    
    auto lih = vlxmol::getMoleculeLiH();
    
    CGtoPairsContainer acont(lih, bas, 1.0e-13);
    
    CGtoBlock sorb(lih, bas, 0);
    
    CGtoBlock porb(lih, bas, 1);
    
    CGtoPairsBlock bss(sorb, 1.0e-13);
    
    ASSERT_EQ(bss, acont.getGtoPairsBlock(0));
    
    CGtoPairsBlock bsp(sorb, porb, 1.0e-13);
    
    ASSERT_EQ(bsp, acont.getGtoPairsBlock(1));
    
    CGtoPairsBlock bpp(porb, 1.0e-13);
    
    ASSERT_EQ(bpp, acont.getGtoPairsBlock(2));
}

