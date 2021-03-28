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

#include "ThreeIndexesTest.hpp"

#include "ThreeIndexes.hpp"

TEST_F(CThreeIndexesTest, DefaultConstructor)
{
    CThreeIndexes ida(-1, -1, -1);

    CThreeIndexes idb;

    ASSERT_EQ(ida, idb);
}

TEST_F(CThreeIndexesTest, CopyConstructor)
{
    CThreeIndexes ida(2, 7, 3);

    CThreeIndexes idb(ida);

    ASSERT_EQ(ida, idb);
}

TEST_F(CThreeIndexesTest, MoveConstructor)
{
    CThreeIndexes ida(2, 7, 3);

    CThreeIndexes idb(CThreeIndexes(2, 7, 3));

    ASSERT_EQ(ida, idb);
}

TEST_F(CThreeIndexesTest, CopyAssignment)
{
    CThreeIndexes ida(2, 7, 3);

    CThreeIndexes idb = ida;

    ASSERT_EQ(ida, idb);
}

TEST_F(CThreeIndexesTest, MoveAssignment)
{
    CThreeIndexes ida(2, 7, 3);

    CThreeIndexes idb = CThreeIndexes(2, 7, 3);

    ASSERT_EQ(ida, idb);
}

TEST_F(CThreeIndexesTest, First)
{
    CThreeIndexes ida(2, 7, 1);

    ASSERT_EQ(2, ida.first());

    CThreeIndexes idb(3, 2, 4);

    ASSERT_EQ(3, idb.first());
}

TEST_F(CThreeIndexesTest, Second)
{
    CThreeIndexes ida(2, 7, 1);

    ASSERT_EQ(7, ida.second());

    CThreeIndexes idb(3, 2, 4);

    ASSERT_EQ(2, idb.second());
}

TEST_F(CThreeIndexesTest, Third)
{
    CThreeIndexes ida(2, 7, 1);

    ASSERT_EQ(1, ida.third());

    CThreeIndexes idb(3, 2, 4);

    ASSERT_EQ(4, idb.third());
}

TEST_F(CThreeIndexesTest, IsValidTriple)
{
    CThreeIndexes ida(2, 7, 0);

    ASSERT_TRUE(ida.isValidTriple());

    CThreeIndexes idb(3, -3, 1);

    ASSERT_FALSE(idb.isValidTriple());

    CThreeIndexes idc(-3, 0, 1);

    ASSERT_FALSE(idc.isValidTriple());

    CThreeIndexes idd(3, 0, -1);

    ASSERT_FALSE(idd.isValidTriple());

    CThreeIndexes ide(0, 0, 0);

    ASSERT_TRUE(ide.isValidTriple());
}
