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

#include "AngularMomentumMatrixTest.hpp"

#include <array>

#include "AngularMomentumMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAngularMomentumMatrixTest, DefaultConstructor)
{
    CAngularMomentumMatrix smata;

    CAngularMomentumMatrix smatb(CDenseMatrix(0, 0), CDenseMatrix(0, 0), CDenseMatrix(0, 0), 0.0, 0.0, 0.0);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, CopyConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.0, 2.0, 4.0);

    CAngularMomentumMatrix smatb(smata);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, MoveConstructor)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 2.0, 3.0, 5.0);

    CAngularMomentumMatrix smatb(CAngularMomentumMatrix(mx, my, mz, 2.0, 3.0, 5.0));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, CopyAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 2.9);

    CAngularMomentumMatrix smatb = smata;

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, MoveAssignment)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 0.3, -2.1, 1.8);

    CAngularMomentumMatrix smatb = CAngularMomentumMatrix(mx, my, mz, 0.3, -2.1, 1.8);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, SetOriginCoordinates)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    smata.setOriginCoordinates(std::array<double, 3>({2.3, 1.6, 8.1}));

    CAngularMomentumMatrix smatb(mx, my, mz, 2.3, 1.6, 8.1);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CAngularMomentumMatrixTest, GetStringForComponentX)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 2.0, 8.0);

    ASSERT_EQ(mx.getString(), smata.getStringForComponentX());
}

TEST_F(CAngularMomentumMatrixTest, GetStringForComponentY)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 2.0, 8.0);

    ASSERT_EQ(my.getString(), smata.getStringForComponentY());
}

TEST_F(CAngularMomentumMatrixTest, GetStringForComponentZ)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    CAngularMomentumMatrix smata(mx, my, mz, 1.2, 2.0, 8.0);

    ASSERT_EQ(mz.getString(), smata.getStringForComponentZ());
}

TEST_F(CAngularMomentumMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAngularMomentumMatrix smata(ma, ma, ma, 1.2, 3.2, 8.9);

    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CAngularMomentumMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAngularMomentumMatrix smata(ma, ma, ma, 0.2, 3.1, 4.5);

    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CAngularMomentumMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAngularMomentumMatrix smata(ma, ma, ma, 1.2, 6.7, 8.1);

    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CAngularMomentumMatrixTest, XValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, smata.xvalues());
}

TEST_F(CAngularMomentumMatrixTest, YValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    vlxtest::compare({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, smata.yvalues());
}

TEST_F(CAngularMomentumMatrixTest, ZValues)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    vlxtest::compare({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, smata.zvalues());
}

TEST_F(CAngularMomentumMatrixTest, GetOriginCoordinates)
{
    CDenseMatrix mx({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix my({1.0, -2.0, -1.0, -6.0, 6.0, 1.0, 3.0, 4.0, -8.0}, 3, 3);

    CDenseMatrix mz({7.0, -1.0, -2.0, -9.0, 5.0, 8.0, 1.0, 3.0, -2.0}, 3, 3);

    const CAngularMomentumMatrix smata(mx, my, mz, 1.2, 3.4, 6.1);

    ASSERT_NEAR(1.2, smata.getOriginCoordinates()[0], 1.0e-13);

    ASSERT_NEAR(3.4, smata.getOriginCoordinates()[1], 1.0e-13);

    ASSERT_NEAR(6.1, smata.getOriginCoordinates()[2], 1.0e-13);
}
