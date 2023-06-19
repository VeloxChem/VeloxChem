//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
    
    ASSERT_EQ(dgride, CDensityGrid(CMemBlock2D<double>(5, 11), dengrid::ab));
    
    CDensityGrid dgridf(5, 1, xcfun::gga, dengrid::lima);
    
    ASSERT_EQ(dgridf, CDensityGrid(CMemBlock2D<double>(5, 5), dengrid::lima));
    
    CDensityGrid dgridg(5, 1, xcfun::gga, dengrid::limb);
    
    ASSERT_EQ(dgridg, CDensityGrid(CMemBlock2D<double>(5, 5), dengrid::limb));
    
    CDensityGrid dgridk(5, 1, xcfun::mgga, dengrid::ab);
    
    ASSERT_EQ(dgridk, CDensityGrid(CMemBlock2D<double>(5, 15), dengrid::ab));
    
    CDensityGrid dgridl(5, 1, xcfun::mgga, dengrid::lima);
    
    ASSERT_EQ(dgridl, CDensityGrid(CMemBlock2D<double>(5, 7), dengrid::lima));
    
    CDensityGrid dgridm(5, 1, xcfun::mgga, dengrid::limb);
    
    ASSERT_EQ(dgridm, CDensityGrid(CMemBlock2D<double>(5, 7), dengrid::limb));
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

TEST_F(CDensityGridTest, Zero)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    dgrida.zero();
    
    ASSERT_EQ(dgrida, CDensityGrid(CMemBlock2D<double>({0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, 2, 4), dengrid::ab));
}

TEST_F(CDensityGridTest, Slice)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    dgrida.slice(1);
    
    ASSERT_EQ(dgrida, CDensityGrid(CMemBlock2D<double>({1.0, 3.0, 2.0, 5.0}, 1, 4), dengrid::ab));
}

TEST_F(CDensityGridTest, GetNumberOfGridPoints)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    ASSERT_EQ(2, dgrida.getNumberOfGridPoints());
}

TEST_F(CDensityGridTest, GetNumberOfDensityMatrices)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    ASSERT_EQ(1, dgrida.getNumberOfDensityMatrices());
    
    CDensityGrid dgridb(5, 3, xcfun::lda, dengrid::lima);
    
    ASSERT_EQ(3, dgridb.getNumberOfDensityMatrices()); 
}

TEST_F(CDensityGridTest, GetDensityGridType)
{
    CDensityGrid dgrida(5, 1, xcfun::lda, dengrid::ab);
    
    ASSERT_EQ(dgrida.getDensityGridType(), dengrid::ab);
    
    CDensityGrid dgridb(5, 1, xcfun::lda, dengrid::lima);
    
    ASSERT_EQ(dgridb.getDensityGridType(), dengrid::lima);
    
    CDensityGrid dgridc(5, 1, xcfun::lda, dengrid::limb);
    
    ASSERT_EQ(dgridc.getDensityGridType(), dengrid::limb);
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

TEST_F(CDensityGridTest, AlphaDensityGradientConstant)
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

TEST_F(CDensityGridTest, AlphaDensityGradient)
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

TEST_F(CDensityGridTest, BetaDensityGradientConstant)
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

TEST_F(CDensityGridTest, BetaDensityGradient)
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

TEST_F(CDensityGridTest, MixedDensityGradientConstant)
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

TEST_F(CDensityGridTest, MixedDensityGradient)
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

TEST_F(CDensityGridTest, AlphaDensityGradientXConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_x_ab = dgridab.alphaDensityGradientX(0);
    
    auto grada_x_lb = dgridlb.alphaDensityGradientX(0);
    
    vlxtest::compare({1.0, 3.0}, grada_x_ab);
    
    vlxtest::compare({1.0, 3.0}, grada_x_lb);
    
    ASSERT_TRUE(dgridla.alphaDensityGradientX(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGradientX)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_x_ab = dgridab.alphaDensityGradientX(0);
    
    auto grada_x_lb = dgridlb.alphaDensityGradientX(0);
    
    vlxtest::compare({1.0, 3.0}, grada_x_ab);
    
    vlxtest::compare({1.0, 3.0}, grada_x_lb);
    
    grada_x_ab[1] = 4.0;
    
    grada_x_lb[0] = 2.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 4.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> nblocka({3.0, 6.0, 7.0, 2.0,
                                 2.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
    
    ASSERT_EQ(dgridlb, CDensityGrid(nblocka, dengrid::limb));
    
    ASSERT_TRUE(dgridla.alphaDensityGradientX(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGradientYConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_y_ab = dgridab.alphaDensityGradientY(0);
    
    auto grada_y_lb = dgridlb.alphaDensityGradientY(0);
    
    vlxtest::compare({3.0, 7.0}, grada_y_ab);
    
    vlxtest::compare({3.0, 7.0}, grada_y_lb);
    
    ASSERT_TRUE(dgridla.alphaDensityGradientY(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGradientY)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_y_ab = dgridab.alphaDensityGradientY(0);
    
    auto grada_y_lb = dgridlb.alphaDensityGradientY(0);
    
    vlxtest::compare({3.0, 7.0}, grada_y_ab);
    
    vlxtest::compare({3.0, 7.0}, grada_y_lb);
    
    grada_y_ab[1] = 4.0;
    
    grada_y_lb[0] = 2.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 4.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> nblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 2.0, 7.0, 8.0, 2.0}, 2, 5);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
    
    ASSERT_EQ(dgridlb, CDensityGrid(nblocka, dengrid::limb));
    
    ASSERT_TRUE(dgridla.alphaDensityGradientY(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGradientZConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_z_ab = dgridab.alphaDensityGradientZ(0);
    
    auto grada_z_lb = dgridlb.alphaDensityGradientZ(0);
    
    vlxtest::compare({8.0, 2.0}, grada_z_ab);
    
    vlxtest::compare({8.0, 2.0}, grada_z_lb);
    
    ASSERT_TRUE(dgridla.alphaDensityGradientZ(0) == nullptr);
}

TEST_F(CDensityGridTest, AlphaDensityGradientZ)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto grada_z_ab = dgridab.alphaDensityGradientZ(0);
    
    auto grada_z_lb = dgridlb.alphaDensityGradientZ(0);
    
    vlxtest::compare({8.0, 2.0}, grada_z_ab);
    
    vlxtest::compare({8.0, 2.0}, grada_z_lb);
    
    grada_z_ab[1] = 4.0;
    
    grada_z_lb[0] = 2.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 4.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> nblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 2.0, 2.0}, 2, 5);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
    
    ASSERT_EQ(dgridlb, CDensityGrid(nblocka, dengrid::limb));
    
    ASSERT_TRUE(dgridla.alphaDensityGradientX(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGradientXConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_x_ab = dgridab.betaDensityGradientX(0);
    
    auto gradb_x_la = dgridla.betaDensityGradientX(0);
    
    vlxtest::compare({5.0, 6.0}, gradb_x_ab);
    
    vlxtest::compare({5.0, 6.0}, gradb_x_la);
    
    ASSERT_TRUE(dgridlb.betaDensityGradientX(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGradientX)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_x_ab = dgridab.betaDensityGradientX(0);
    
    auto gradb_x_la = dgridla.betaDensityGradientX(0);
    
    vlxtest::compare({5.0, 6.0}, gradb_x_ab);
    
    vlxtest::compare({5.0, 6.0}, gradb_x_la);
    
    gradb_x_ab[1] = 4.0;
    
    gradb_x_la[0] = 2.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 4.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> nblockb({1.0, 2.0, 3.0, 5.0,
                                 2.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
    
    ASSERT_EQ(dgridla, CDensityGrid(nblockb, dengrid::lima));
    
    ASSERT_TRUE(dgridlb.betaDensityGradientX(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGradientYConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_y_ab = dgridab.betaDensityGradientY(0);
    
    auto gradb_y_la = dgridla.betaDensityGradientY(0);
    
    vlxtest::compare({2.0, 3.0}, gradb_y_ab);
    
    vlxtest::compare({2.0, 3.0}, gradb_y_la);
    
    ASSERT_TRUE(dgridlb.betaDensityGradientY(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGradientY)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_y_ab = dgridab.betaDensityGradientY(0);
    
    auto gradb_y_la = dgridla.betaDensityGradientY(0);
    
    vlxtest::compare({2.0, 3.0}, gradb_y_ab);
    
    vlxtest::compare({2.0, 3.0}, gradb_y_la);
    
    gradb_y_ab[1] = 4.0;
    
    gradb_y_la[0] = 2.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 4.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> nblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
    
    ASSERT_EQ(dgridla, CDensityGrid(nblockb, dengrid::lima));
    
    ASSERT_TRUE(dgridlb.betaDensityGradientY(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGradientZConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    const CDensityGrid dgridla(mblockb, dengrid::lima);
    
    const CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_z_ab = dgridab.betaDensityGradientZ(0);
    
    auto gradb_z_la = dgridla.betaDensityGradientZ(0);
    
    vlxtest::compare({1.0, 6.0}, gradb_z_ab);
    
    vlxtest::compare({1.0, 6.0}, gradb_z_la);
    
    ASSERT_TRUE(dgridlb.betaDensityGradientZ(0) == nullptr);
}

TEST_F(CDensityGridTest, BetaDensityGradientZ)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CMemBlock2D<double> mblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    CDensityGrid dgridla(mblockb, dengrid::lima);
    
    CDensityGrid dgridlb(mblocka, dengrid::limb);
    
    auto gradb_z_ab = dgridab.betaDensityGradientZ(0);
    
    auto gradb_z_la = dgridla.betaDensityGradientZ(0);
    
    vlxtest::compare({1.0, 6.0}, gradb_z_ab);
    
    vlxtest::compare({1.0, 6.0}, gradb_z_la);
    
    gradb_z_ab[1] = 4.0;
    
    gradb_z_la[0] = 2.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 4.0}, 2, 11);
    
    CMemBlock2D<double> nblockb({1.0, 2.0, 3.0, 5.0,
                                 5.0, 6.0, 2.0, 3.0, 2.0, 6.0}, 2, 5);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
    
    ASSERT_EQ(dgridla, CDensityGrid(nblockb, dengrid::lima));
    
    ASSERT_TRUE(dgridlb.betaDensityGradientZ(0) == nullptr);
}

TEST_F(CDensityGridTest, GetComponent)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto cgrid0 = dgridab.getComponent(0); cgrid0[0] = 3.0;
    
    auto cgrid1 = dgridab.getComponent(1); cgrid1[1] = 9.0;
    
    auto cgrid2 = dgridab.getComponent(2); cgrid2[1] = 0.0;
    
    auto cgrid3 = dgridab.getComponent(3); cgrid3[0] = 1.0;
    
    auto cgrid4 = dgridab.getComponent(4); cgrid4[0] = 3.0;
    
    CMemBlock2D<double> rblockab({3.0, 2.0, 3.0, 9.0, 4.0, 0.0, 1.0, 9.0, 3.0, 4.0}, 2, 5);
    
    CDensityGrid rgridab(rblockab, dengrid::ab);
    
    ASSERT_EQ(dgridab, rgridab);
}

TEST_F(CDensityGridTest, GetComponentConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    vlxtest::compare({1.0, 2.0}, dgridab.getComponent(0));
    
    vlxtest::compare({3.0, 6.0}, dgridab.getComponent(1));
    
    vlxtest::compare({4.0, 5.0}, dgridab.getComponent(2));
    
    vlxtest::compare({7.0, 9.0}, dgridab.getComponent(3));
    
    vlxtest::compare({2.0, 4.0}, dgridab.getComponent(4));
}
