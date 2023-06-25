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
    
    CDensityGrid dgride(5, 1, xcfun::gga, dengrid::ab);
    
    ASSERT_EQ(dgride, CDensityGrid(CMemBlock2D<double>(5, 11), dengrid::ab));
    
    CDensityGrid dgridk(5, 1, xcfun::mgga, dengrid::ab);
    
    ASSERT_EQ(dgridk, CDensityGrid(CMemBlock2D<double>(5, 15), dengrid::ab));
}

TEST_F(CDensityGridTest, MoveConstructor)
{
    CMemBlock2D<double> mblock({1.0, 2.0, 3.0, 6.0, 2.0, 4.0, 5.0, 9.0}, 2, 4);
    
    CDensityGrid dgrida(mblock, dengrid::ab);
    
    CDensityGrid dgridb((CDensityGrid(mblock, dengrid::ab)));
    
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
}

TEST_F(CDensityGridTest, GetDensityGridType)
{
    CDensityGrid dgrida(5, 1, xcfun::lda, dengrid::ab);
    
    ASSERT_EQ(dgrida.getDensityGridType(), dengrid::ab);
}

TEST_F(CDensityGridTest, AlphaDensityConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    vlxtest::compare({1.0, 2.0}, dgridab.alphaDensity(0));
}

TEST_F(CDensityGridTest, AlphaDensity)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto rhoa_ab = dgridab.alphaDensity(0);
    
    vlxtest::compare({1.0, 2.0}, rhoa_ab);
    
    rhoa_ab[1] = 7.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 7.0, 3.0, 6.0},  2, 2), dengrid::ab));
}

TEST_F(CDensityGridTest, BetaDensityConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    vlxtest::compare({3.0, 6.0}, dgridab.betaDensity(0));
}

TEST_F(CDensityGridTest, BetaDensity)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0}, 2, 2);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto rhob_ab = dgridab.betaDensity(0);
    
    vlxtest::compare({3.0, 6.0}, rhob_ab);
    
    rhob_ab[0] = 8.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 8.0, 6.0},  2, 2), dengrid::ab));
}

TEST_F(CDensityGridTest, AlphaDensityGradientConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    vlxtest::compare({4.0, 5.0}, dgridab.alphaDensityGradient(0));
}

TEST_F(CDensityGridTest, AlphaDensityGradient)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_ab = dgridab.alphaDensityGradient(0);
    
    vlxtest::compare({4.0, 5.0}, grada_ab);
    
    grada_ab[1] = 2.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 4.0, 2.0, 7.0, 9.0, 2.0, 4.0}, 2, 5), dengrid::ab));
}

TEST_F(CDensityGridTest, BetaDensityGradientConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    vlxtest::compare({7.0, 9.0}, dgridab.betaDensityGradient(0));
}

TEST_F(CDensityGridTest, BetaDensityGradient)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_ab = dgridab.betaDensityGradient(0);
    
    vlxtest::compare({7.0, 9.0}, gradb_ab);
    
    gradb_ab[0] = 6.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 6.0, 9.0, 2.0, 4.0}, 2, 5), dengrid::ab));
}

TEST_F(CDensityGridTest, MixedDensityGradientConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    vlxtest::compare({2.0, 4.0}, dgridab.mixedDensityGradient(0));
}

TEST_F(CDensityGridTest, MixedDensityGradient)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradab_ab = dgridab.mixedDensityGradient(0);
    
    vlxtest::compare({2.0, 4.0}, gradab_ab);
    
    gradab_ab[0] = 3.0;
    
    gradab_ab[1] = 9.0;
    
    ASSERT_EQ(dgridab, CDensityGrid(CMemBlock2D<double>({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 3.0, 9.0}, 2, 5), dengrid::ab));
}

TEST_F(CDensityGridTest, AlphaDensityGradientXConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_x_ab = dgridab.alphaDensityGradientX(0);
    
    vlxtest::compare({1.0, 3.0}, grada_x_ab);
}

TEST_F(CDensityGridTest, AlphaDensityGradientX)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_x_ab = dgridab.alphaDensityGradientX(0);
    
    vlxtest::compare({1.0, 3.0}, grada_x_ab);
    
    grada_x_ab[1] = 4.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 4.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
}

TEST_F(CDensityGridTest, AlphaDensityGradientYConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_y_ab = dgridab.alphaDensityGradientY(0);
    
    vlxtest::compare({3.0, 7.0}, grada_y_ab);
}

TEST_F(CDensityGridTest, AlphaDensityGradientY)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_y_ab = dgridab.alphaDensityGradientY(0);
    
    vlxtest::compare({3.0, 7.0}, grada_y_ab);
    
    grada_y_ab[1] = 4.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 4.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
}

TEST_F(CDensityGridTest, AlphaDensityGradientZConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_z_ab = dgridab.alphaDensityGradientZ(0);
    
    vlxtest::compare({8.0, 2.0}, grada_z_ab);
}

TEST_F(CDensityGridTest, AlphaDensityGradientZ)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto grada_z_ab = dgridab.alphaDensityGradientZ(0);
    
    vlxtest::compare({8.0, 2.0}, grada_z_ab);
    
    grada_z_ab[1] = 4.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 4.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
}

TEST_F(CDensityGridTest, BetaDensityGradientXConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_x_ab = dgridab.betaDensityGradientX(0);
    
    vlxtest::compare({5.0, 6.0}, gradb_x_ab);
}

TEST_F(CDensityGridTest, BetaDensityGradientX)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_x_ab = dgridab.betaDensityGradientX(0);
    
    vlxtest::compare({5.0, 6.0}, gradb_x_ab);
    
    gradb_x_ab[1] = 4.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 4.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
}

TEST_F(CDensityGridTest, BetaDensityGradientYConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_y_ab = dgridab.betaDensityGradientY(0);
    
    vlxtest::compare({2.0, 3.0}, gradb_y_ab);
}

TEST_F(CDensityGridTest, BetaDensityGradientY)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CMemBlock2D<double> mblocka({3.0, 6.0, 7.0, 2.0,
                                 1.0, 3.0, 3.0, 7.0, 8.0, 2.0}, 2, 5);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_y_ab = dgridab.betaDensityGradientY(0);
    
    vlxtest::compare({2.0, 3.0}, gradb_y_ab);
    
    gradb_y_ab[1] = 4.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 4.0, 1.0, 6.0}, 2, 11);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
}

TEST_F(CDensityGridTest, BetaDensityGradientZConstant)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    const CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_z_ab = dgridab.betaDensityGradientZ(0);
    
    vlxtest::compare({1.0, 6.0}, gradb_z_ab);
}

TEST_F(CDensityGridTest, BetaDensityGradientZ)
{
    CMemBlock2D<double> mblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 6.0}, 2, 11);
    
    CDensityGrid dgridab(mblockab, dengrid::ab);
    
    auto gradb_z_ab = dgridab.betaDensityGradientZ(0);
    
    vlxtest::compare({1.0, 6.0}, gradb_z_ab);
    
    gradb_z_ab[1] = 4.0;
    
    CMemBlock2D<double> nblockab({1.0, 2.0, 3.0, 6.0, 4.0, 5.0, 7.0, 9.0, 2.0, 4.0,
                                  1.0, 3.0, 3.0, 7.0, 8.0, 2.0,
                                  5.0, 6.0, 2.0, 3.0, 1.0, 4.0}, 2, 11);
    
    ASSERT_EQ(dgridab, CDensityGrid(nblockab, dengrid::ab));
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
