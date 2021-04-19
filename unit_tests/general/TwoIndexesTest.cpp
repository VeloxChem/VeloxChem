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

#include "TwoIndexesTest.hpp"

#include "TwoIndexes.hpp"

TEST_F(CTwoIndexesTest, DefaultConstructor)
{
    CTwoIndexes ida(-1, -1);

    CTwoIndexes idb;

    ASSERT_EQ(ida, idb);
}

TEST_F(CTwoIndexesTest, CopyConstructor)
{
    CTwoIndexes ida(2, 7);

    CTwoIndexes idb(ida);

    ASSERT_EQ(ida, idb);
}

TEST_F(CTwoIndexesTest, MoveConstructor)
{
    CTwoIndexes ida(2, 7);

    CTwoIndexes idb(CTwoIndexes(2, 7));

    ASSERT_EQ(ida, idb);
}

TEST_F(CTwoIndexesTest, CopyAssignment)
{
    CTwoIndexes ida(2, 7);

    CTwoIndexes idb = ida;

    ASSERT_EQ(ida, idb);
}

TEST_F(CTwoIndexesTest, MoveAssignment)
{
    CTwoIndexes ida(2, 7);

    CTwoIndexes idb = CTwoIndexes(2, 7);

    ASSERT_EQ(ida, idb);
}

TEST_F(CTwoIndexesTest, First)
{
    CTwoIndexes ida(2, 7);

    ASSERT_EQ(2, ida.first());

    CTwoIndexes idb(3, 3);

    ASSERT_EQ(3, idb.first());
}

TEST_F(CTwoIndexesTest, Second)
{
    CTwoIndexes ida(2, 7);

    ASSERT_EQ(7, ida.second());

    CTwoIndexes idb(3, 3);

    ASSERT_EQ(3, idb.second());
}

TEST_F(CTwoIndexesTest, IsValidPair)
{
    CTwoIndexes ida(2, 7);

    ASSERT_TRUE(ida.isValidPair());

    CTwoIndexes idb(3, -3);

    ASSERT_FALSE(idb.isValidPair());

    CTwoIndexes idc(-3, 0);

    ASSERT_FALSE(idc.isValidPair());

    CTwoIndexes idd(0, 0);

    ASSERT_TRUE(idd.isValidPair());
}
