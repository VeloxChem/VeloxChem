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

#include "NuclearPotentialMatrixTest.hpp"

#include "CheckFunctions.hpp"
#include "NuclearPotentialMatrix.hpp"

TEST_F(CNuclearPotentialMatrixTest, DefaultConstructor)
{
    CNuclearPotentialMatrix smata;

    CNuclearPotentialMatrix smatb(CDenseMatrix(0, 0));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CNuclearPotentialMatrix smata(ma);

    CNuclearPotentialMatrix smatb(smata);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CNuclearPotentialMatrix smata(ma);

    CNuclearPotentialMatrix smatb(CNuclearPotentialMatrix({ma}));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CNuclearPotentialMatrix smata(ma);

    CNuclearPotentialMatrix smatb = smata;

    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CNuclearPotentialMatrix smata(ma);

    CNuclearPotentialMatrix smatb = CNuclearPotentialMatrix({ma});

    ASSERT_EQ(smata, smatb);
}

TEST_F(CNuclearPotentialMatrixTest, GetString)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CNuclearPotentialMatrix smata(ma);

    ASSERT_EQ(ma.getString(), smata.getString());
}

TEST_F(CNuclearPotentialMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CNuclearPotentialMatrix smata(ma);

    ASSERT_EQ(3, smata.getNumberOfRows());
}

TEST_F(CNuclearPotentialMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CNuclearPotentialMatrix smata(ma);

    ASSERT_EQ(2, smata.getNumberOfColumns());
}

TEST_F(CNuclearPotentialMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    CNuclearPotentialMatrix smata(ma);

    ASSERT_EQ(6, smata.getNumberOfElements());
}

TEST_F(CNuclearPotentialMatrixTest, Values)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    const CNuclearPotentialMatrix smata(ma);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, smata.values());
}

TEST_F(CNuclearPotentialMatrixTest, GetKineticEnergy)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CNuclearPotentialMatrix smat(ma);

    CDenseMatrix da({1.0, -1.0, -3.0, -2.0, 2.0, 3.0, 6.0, 3.0, -4.0}, 3, 3);

    CDenseMatrix db({1.0, -1.0, -3.0, -2.0, 1.0, 4.0, 3.0, 7.0, 3.0}, 3, 3);

    CAODensityMatrix dmat({da, db}, denmat::rest);

    ASSERT_NEAR(19.0, smat.getNuclearPotentialEnergy(dmat, 0), 1.0e-13);

    ASSERT_NEAR(15.0, smat.getNuclearPotentialEnergy(dmat, 1), 1.0e-13);

    ASSERT_NEAR(0.0, smat.getNuclearPotentialEnergy(dmat, 2), 1.0e-13);
}
