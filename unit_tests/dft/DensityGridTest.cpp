//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DensityGridTest.hpp"

#include "DensityGrid.hpp"
#include "CheckFunctions.hpp"

TEST_F(CDensityGridTest, DefaultConstructor)
{
    CDensityGrid dgrida;
    
    CMemBlock2D<double> mblock;
    
    CDensityGrid dgridb(0, 0, xcfun::lda, dengrid::undefined);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CDensityGridTest, ConstructorWithFunctionalType)
{
    CDensityGrid dgrida(5, 1, xcfun::lda, dengrid::ab);
    
    ASSERT_EQ(dgrida, CDensityGrid(CMemBlock2D<double>(5, 2), dengrid::ab));
    
    CDensityGrid dgridb(5, 1, xcfun::lda, dengrid::lima);
    
    ASSERT_EQ(dgridb, CDensityGrid(CMemBlock2D<double>(5, 1), dengrid::lima));
    
    CDensityGrid dgridc(5, 1, xcfun::lda, dengrid::limb);
    
    ASSERT_EQ(dgridc, CDensityGrid(CMemBlock2D<double>(5, 1), dengrid::limb));
    
    CDensityGrid dgride(5, 1, xcfun::gga, dengrid::ab);
    
    ASSERT_EQ(dgride, CDensityGrid(CMemBlock2D<double>(5, 5), dengrid::ab));
    
    CDensityGrid dgridf(5, 1, xcfun::gga, dengrid::lima);
    
    ASSERT_EQ(dgridf, CDensityGrid(CMemBlock2D<double>(5, 2), dengrid::lima));
    
    CDensityGrid dgridg(5, 1, xcfun::gga, dengrid::limb);
    
    ASSERT_EQ(dgridg, CDensityGrid(CMemBlock2D<double>(5, 2), dengrid::limb));
    
    CDensityGrid dgridk(5, 1, xcfun::mgga, dengrid::ab);
    
    ASSERT_EQ(dgridk, CDensityGrid(CMemBlock2D<double>(5, 7), dengrid::ab));
    
    CDensityGrid dgridl(5, 1, xcfun::mgga, dengrid::lima);
    
    ASSERT_EQ(dgridl, CDensityGrid(CMemBlock2D<double>(5, 3), dengrid::lima));
    
    CDensityGrid dgridm(5, 1, xcfun::mgga, dengrid::limb);
    
    ASSERT_EQ(dgridm, CDensityGrid(CMemBlock2D<double>(5, 3), dengrid::limb));
}

TEST_F(CDensityGridTest, CopyConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::lima);
    
    CDensityGrid dgridb(dgrida);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CDensityGridTest, MoveConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    CDensityGrid dgridb((CDensityGrid(mblock, dengrid::ab)));
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CDensityGridTest, CopyAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::limb);
    
    CDensityGrid dgridb = dgrida;
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CDensityGridTest, MoveAssignment)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    CDensityGrid dgridb = CDensityGrid(mblock, dengrid::ab);
    
    ASSERT_EQ(dgrida, dgridb);
}

TEST_F(CDensityGridTest, GetNumberOfGridPoints)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    ASSERT_EQ(2, dgrida.getNumberOfGridPoints());
}

TEST_F(CDensityGridTest, AlphaDensityConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    CMemBlock2D<double> mblocka({3.0, 6.0}, 2, 1);
    
    CMemBlock2D<double> mblockb({1.0, 2.0}, 2, 1);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    vlxtest::compare({1.0, 2.0}, dgridab.alphaDensity(0));
    
    vlxtest::compare({3.0, 6.0}, dgridlb.alphaDensity(0));
    
    ASSERT_TRUE(dgridla.alphaDensity(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensity)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    CMemBlock2D<double> mblocka({3.0, 6.0}, 2, 1);
    
    CMemBlock2D<double> mblockb({1.0, 2.0}, 2, 1);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto rhoa_ab = dgridab.alphaDensity(0);
    
    auto rhoa_lb = dgridlb.alphaDensity(0);
    
    vlxtest::compare({1.0, 2.0}, rhoa_ab);
    
    vlxtest::compare({3.0, 6.0}, rhoa_lb);
    
    rhoa_ab[1] = 7.0;
    
    rhoa_lb[0] = 4.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 7.0, 3.0, 6.0},  2, 2), dengrid::ab));
    
    ASSERT_EQ(dgridlb, CDensityGrid(CMemBlock2D<double>({4.0, 6.0},  2, 1), dengrid::limb));
    
    ASSERT_TRUE(dgridla.alphaDensity(0) == nullptr);
}


TEST_F(CDensityGridTest, BetaDensityConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    CMemBlock2D<double> mblocka({3.0, 6.0}, 2, 1);
    
    CMemBlock2D<double> mblockb({1.0, 2.0}, 2, 1);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    vlxtest::compare({3.0, 6.0}, dgridab.betaDensity(0));
    
    vlxtest::compare({1.0, 2.0}, dgridla.betaDensity(0));
    
    ASSERT_TRUE(dgridlb.betaDensity(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensity)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    CMemBlock2D<double> mblocka({3.0, 6.0}, 2, 1);
    
    CMemBlock2D<double> mblockb({1.0, 2.0}, 2, 1);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
   
    auto rhob_ab = dgridab.betaDensity(0);
    
    auto rhob_la = dgridla.betaDensity(0);
    
    vlxtest::compare({3.0, 6.0}, rhob_ab);
    
    vlxtest::compare({1.0, 2.0}, rhob_la);
    
    rhob_ab[0] = 8.0;
    
    rhob_la[0] = 9.0;
   
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 8.0, 6.0},  2, 2), dengrid::ab));
    
    ASSERT_EQ(dgridla, CDensityGrid(CMemBlock2D<double>({9.0, 2.0},  2, 1), dengrid::lima));
    
    ASSERT_TRUE(dgridlb.betaDensity(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGraddientConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0}, 2, 2);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    vlxtest::compare({4.0, 5.0}, dgridab.alphaDensityGradient(0));
    
    vlxtest::compare({7.0, 2.0}, dgridlb.alphaDensityGradient(0));
    
    ASSERT_TRUE(dgridla.alphaDensityGradient(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGraddient)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0}, 2, 2);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_ab = dgridab.alphaDensityGradient(0);
    
    auto grada_lb = dgridlb.alphaDensityGradient(0);
    
    vlxtest::compare({4.0, 5.0}, grada_ab);
    
    vlxtest::compare({7.0, 2.0}, grada_lb);
    
    grada_ab[1] = 2.0;
    
    grada_lb[0] = 3.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 4.0, 2.0, 7.0, 9.0, 2.0, 4.0}, 2, 5), dengrid::ab));
    
    ASSERT_EQ(dgridlb, CDensityGrid(CMemBlock2D<double>({3.0, 6.0, 3.0, 2.0}, 2, 2), dengrid::limb));
    
    ASSERT_TRUE(dgridla.alphaDensityGradient(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGraddientConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0}, 2, 2);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    vlxtest::compare({7.0, 9.0}, dgridab.betaDensityGradient(0));
    
    vlxtest::compare({3.0, 5.0}, dgridla.betaDensityGradient(0));
    
    ASSERT_TRUE(dgridlb.betaDensityGradient(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGraddient)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0}, 2, 2);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_ab = dgridab.betaDensityGradient(0);
    
    auto gradb_la = dgridla.betaDensityGradient(0);
    
    vlxtest::compare({7.0, 9.0}, gradb_ab);
    
    vlxtest::compare({3.0, 5.0}, gradb_la);
    
    gradb_ab[0] = 6.0;
    
    gradb_la[0] = 7.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 6.0, 9.0, 2.0, 4.0}, 2, 5), dengrid::ab));
    
    ASSERT_EQ(dgridla, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 7.0, 5.0}, 2, 2), dengrid::lima));
    
    ASSERT_TRUE(dgridlb.betaDensityGradient(0) == nullptr);
}

TEST_F(CDensityGridTest, MixedDensityGraddientConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0}, 2, 2);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    vlxtest::compare({2.0, 4.0}, dgridab.mixedDensityGradient(0));
    
    ASSERT_TRUE(dgridla.mixedDensityGradient(0) == nullptr);
    
    ASSERT_TRUE(dgridlb.mixedDensityGradient(0) == nullptr);
}

TEST_F(CDensityGridTest, MixedDensityGraddient)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0}, 2, 2);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0}, 2, 2);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradab_ab = dgridab.mixedDensityGradient(0);
    
    vlxtest::compare({2.0, 4.0}, gradab_ab);
    
    gradab_ab[0] = 3.0;
    
    gradab_ab[1] = 9.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 3.0, 9.0}, 2, 5), dengrid::ab));
    
    ASSERT_TRUE(dgridla.mixedDensityGradient(0) == nullptr);
    
    ASSERT_TRUE(dgridlb.mixedDensityGradient(0) == nullptr);
}
