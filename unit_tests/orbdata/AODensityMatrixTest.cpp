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

#include "AODensityMatrixTest.hpp"

#include "AODensityMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAODensityMatrixTest, DefaultConstructor)
{
    CAODensityMatrix dmata;

    CAODensityMatrix dmatb({}, denmat::rest);

    ASSERT_EQ(dmata, dmatb);
}

TEST_F(CAODensityMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CAODensityMatrix dmatb(dmata);

    ASSERT_EQ(dmata, dmatb);
}

TEST_F(CAODensityMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CAODensityMatrix dmatb(CAODensityMatrix({ma, ma, mb, mb}, denmat::unrest));

    ASSERT_EQ(dmata, dmatb);
}

TEST_F(CAODensityMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma}, denmat::rest);

    CAODensityMatrix dmatb = dmata;

    ASSERT_EQ(dmata, dmatb);
}

TEST_F(CAODensityMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma}, denmat::rest);

    CAODensityMatrix dmatb = CAODensityMatrix({ma}, denmat::rest);

    ASSERT_EQ(dmata, dmatb);
}

TEST_F(CAODensityMatrixTest, SetDensityType)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma}, denmat::rest);

    dmata.setDensityType(denmat::rmoij);

    CAODensityMatrix dmatb = CAODensityMatrix({ma}, denmat::rmoij);

    ASSERT_EQ(dmata, dmatb);
}

TEST_F(CAODensityMatrixTest, Append)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma}, denmat::rest);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAODensityMatrix dmatb({ma, ma, mb, mb}, denmat::rgen);

    dmata.append(dmatb);

    ASSERT_EQ(dmata, CAODensityMatrix({ma}, denmat::rest));

    dmata.setDensityType(denmat::rgen);

    dmata.append(dmatb);

    ASSERT_EQ(dmata, CAODensityMatrix({ma, ma, ma, mb, mb}, denmat::rgen));
}

TEST_F(CAODensityMatrixTest, Sub)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CDenseMatrix mc({0.0, -2.0, 1.0, 4.0, 1.0, 2.0, -6.0, 1.0, -7.0}, 3, 3);

    CDenseMatrix md({1.0, -2.0, -1.0, -8.0, 5.0, 2.0}, 3, 2);

    CAODensityMatrix dmatb({mc, mc, md, md}, denmat::unrest);

    CDenseMatrix me({-1.0, -1.0, 4.0, 6.0, -4.0, -2.0, -12.0, -3.0, -3.0}, 3, 3);

    CDenseMatrix mf({0.0, -1.0, 2.0, -6.0, 0.0, -2.0}, 3, 2);

    CAODensityMatrix dmatc({me, me, mf, mf}, denmat::unrest);

    ASSERT_EQ(dmatc, dmatb.sub(dmata));
}

TEST_F(CAODensityMatrixTest, GetNumberOfDensityMatrices)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma, ma}, denmat::unrest);

    CAODensityMatrix dmatb({ma, ma}, denmat::rest);

    ASSERT_EQ(1, dmata.getNumberOfDensityMatrices());

    ASSERT_EQ(2, dmatb.getNumberOfDensityMatrices());
}

TEST_F(CAODensityMatrixTest, GetNumberOfMatrices)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma, ma}, denmat::unrest);

    CAODensityMatrix dmatb({ma, ma, ma}, denmat::rest);

    ASSERT_EQ(2, dmata.getNumberOfMatrices());

    ASSERT_EQ(3, dmatb.getNumberOfMatrices());
}

TEST_F(CAODensityMatrixTest, GetDensityType)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CAODensityMatrix dmata({ma, ma}, denmat::unrest);

    CAODensityMatrix dmatb({ma, ma}, denmat::rest);

    ASSERT_TRUE(dmata.getDensityType() == denmat::unrest);

    ASSERT_TRUE(dmatb.getDensityType() == denmat::rest);
}

TEST_F(CAODensityMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0}, 2, 2);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CAODensityMatrix dmatb({ma}, denmat::rest);

    ASSERT_EQ(3, dmata.getNumberOfRows(0));

    ASSERT_EQ(2, dmata.getNumberOfRows(1));

    ASSERT_EQ(3, dmatb.getNumberOfRows(0));
}

TEST_F(CAODensityMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CAODensityMatrix dmatb({mb}, denmat::rest);

    ASSERT_EQ(3, dmata.getNumberOfColumns(0));

    ASSERT_EQ(4, dmata.getNumberOfColumns(1));

    ASSERT_EQ(4, dmatb.getNumberOfColumns(0));
}

TEST_F(CAODensityMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CAODensityMatrix dmatb({mb}, denmat::rest);

    ASSERT_EQ(9, dmata.getNumberOfElements(0));

    ASSERT_EQ(8, dmata.getNumberOfElements(1));

    ASSERT_EQ(8, dmatb.getNumberOfElements(0));
}

TEST_F(CAODensityMatrixTest, RestrictedDensity)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    const CAODensityMatrix dmat({ma, mb}, denmat::rest);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, dmat.alphaDensity(0));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, dmat.betaDensity(0));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, dmat.alphaDensity(1));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, dmat.betaDensity(1));

    ASSERT_TRUE(dmat.alphaDensity(2) == nullptr);

    ASSERT_TRUE(dmat.betaDensity(2) == nullptr);
}

TEST_F(CAODensityMatrixTest, UnrestrictedAlphaDensity)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmat({ma, ma, mb, mb}, denmat::unrest);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, dmat.alphaDensity(0));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, dmat.alphaDensity(1));

    ASSERT_TRUE(dmat.alphaDensity(2) == nullptr);
}

TEST_F(CAODensityMatrixTest, UnrestrictedBetaDensity)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmat({mb, ma, ma, mb}, denmat::unrest);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, dmat.betaDensity(0));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, dmat.betaDensity(1));

    ASSERT_TRUE(dmat.betaDensity(2) == nullptr);
}

TEST_F(CAODensityMatrixTest, GetDensity)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmat({ma, ma, mb, mb}, denmat::unrest);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, dmat.getDensity(0));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, dmat.getDensity(1));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, dmat.getDensity(2));

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, dmat.getDensity(3));
}

TEST_F(CAODensityMatrixTest, GetReferenceToDensity)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmat({ma, ma, mb, mb}, denmat::unrest);

    ASSERT_EQ(ma, dmat.getReferenceToDensity(0));

    ASSERT_EQ(ma, dmat.getReferenceToDensity(1));

    ASSERT_EQ(mb, dmat.getReferenceToDensity(2));

    ASSERT_EQ(mb, dmat.getReferenceToDensity(3));
}

TEST_F(CAODensityMatrixTest, IsRestricted)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0}, 2, 4);

    CAODensityMatrix dmata({ma, ma, mb, mb}, denmat::unrest);

    CAODensityMatrix dmatb({mb}, denmat::rest);

    ASSERT_FALSE(dmata.isClosedShell());

    ASSERT_TRUE(dmatb.isClosedShell());
}
