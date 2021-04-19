//
//                           VELOXCHEM 1.0-RC
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

#include "SparseVectorTest.hpp"

#include "CheckFunctions.hpp"
#include "SparseVector.hpp"

TEST_F(CSparseVectorTest, DefaultConstructor)
{
    CSparseVector veca;

    CSparseVector vecb(std::vector<double>{}, std::vector<int32_t>{}, 0, 1.0e-13);

    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, CopyConstructor)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);

    CSparseVector vecb(veca);

    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, MoveConstructor)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);

    CSparseVector vecb(CSparseVector({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13));

    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, CopyAssignment)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);

    CSparseVector vecb = veca;

    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, MoveAssignment)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);

    CSparseVector vecb = CSparseVector({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);

    ASSERT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, ValuesConstant)
{
    const CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0}, veca.values());
}

TEST_F(CSparseVectorTest, Values)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0}, veca.values());

    auto vdat = veca.values();

    vdat[2] = 9.0;

    CSparseVector vecb({1.0, -1.0, 9.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);

    EXPECT_EQ(veca, vecb);
}

TEST_F(CSparseVectorTest, IndexesConstant)
{
    const CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);

    vlxtest::compare({0, 1, 2, 4}, veca.indexes());
}

TEST_F(CSparseVectorTest, GetThreshold)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 0.28);

    ASSERT_NEAR(veca.getThreshold(), 0.28, 1.0e-13);
}

TEST_F(CSparseVectorTest, GetSparsity)
{
    CSparseVector veca({1.0, -1.0, -3.0, -2.0}, {0, 1, 2, 4}, 10, 1.0e-13);

    ASSERT_NEAR(veca.getSparsity(), 0.40, 1.0e-13);
}
