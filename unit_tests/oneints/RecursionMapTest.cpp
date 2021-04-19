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

#include "RecursionMapTest.hpp"

#include "RecursionMap.hpp"

TEST_F(CRecursionMapTest, DefaultConstructor)
{
    CRecursionMap rma;

    CRecursionMap rmb({}, {}, recblock::cc, 0);

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, AlternativeConstructor)
{
    CRecursionMap rma(recblock::ss, 2);

    CRecursionMap rmb({}, {}, recblock::ss, 2);

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, CopyConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3}, 2, 3, 0);

    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc, 3);

    CRecursionMap rmb(rma);

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, MoveConstructor)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3}, 2, 3, 0);

    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc, 3);

    CRecursionMap rmb(CRecursionMap({rta, rtb}, {0, 4}, recblock::cc, 3));

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, CopyAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3}, 2, 3, 0);

    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc, 4);

    CRecursionMap rmb = rma;

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, MoveAssignment)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {6, 1, 7, 3}, 2, 3, 0);

    CRecursionMap rma({rta, rtb}, {0, 4}, recblock::cc, 4);

    CRecursionMap rmb = CRecursionMap({rta, rtb}, {0, 4}, recblock::cc, 4);

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, Add)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, -1, 7, 3}, 1, 2, 0);

    CRecursionMap rma(recblock::cc, 1);

    rma.add(rta);

    rma.add(rtb);

    rma.add(rta);

    rma.add(rta);

    rma.add(rtb);

    rma.add(rtc);

    CRecursionMap rmb({rta, rtb}, {0, 18}, recblock::cc, 1);

    ASSERT_EQ(rma, rmb);
}

TEST_F(CRecursionMapTest, Append)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, -1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtd({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb}, {0, 36}, recblock::cc, 2);

    CRecursionMap rmb({rtc, rtd, rta}, {0, 0, 0}, recblock::cc, 2);

    rma.append(rmb);

    CRecursionMap rmc({rta, rtb, rtd}, {0, 36, 1656}, recblock::cc, 2);

    ASSERT_EQ(rma, rmc);
}

TEST_F(CRecursionMapTest, AppendWithVector)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, -1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtd({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb}, {0, 18}, recblock::cc, 1);

    std::vector<CRecursionTerm> rvec({rtc, rtd, rta});

    rma.append(rvec);

    CRecursionMap rmc({rta, rtb, rtd}, {0, 18, 828}, recblock::cc, 1);

    ASSERT_EQ(rma, rmc);
}

TEST_F(CRecursionMapTest, Find)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {3, 3, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb}, {0, 18}, recblock::cc, 1);

    ASSERT_TRUE(rma.find(rta));

    ASSERT_TRUE(rma.find(rtb));

    ASSERT_FALSE(rma.find(rtc));
}

TEST_F(CRecursionMapTest, GetNumberOfComponents)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb, rtc}, {0, 18, 828}, recblock::cc, 1);

    ASSERT_EQ(1152, rma.getNumberOfComponents());
}

TEST_F(CRecursionMapTest, GetNumberOfTerms)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb, rtc}, {0, 18, 828}, recblock::cc, 1);

    ASSERT_EQ(3, rma.getNumberOfTerms());
}

TEST_F(CRecursionMapTest, GetTerm)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb, rtc}, {0, 18, 828}, recblock::cc, 1);

    ASSERT_EQ(rta, rma.getTerm(0));

    ASSERT_EQ(rtb, rma.getTerm(1));

    ASSERT_EQ(rtc, rma.getTerm(2));
}

TEST_F(CRecursionMapTest, GetIndexOfTerm)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionTerm rtd({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 0);

    CRecursionMap rma({rta, rtb, rtc}, {0, 18, 828}, recblock::cc, 1);

    ASSERT_EQ(0, rma.getIndexOfTerm(rta));

    ASSERT_EQ(18, rma.getIndexOfTerm(rtb));

    ASSERT_EQ(828, rma.getIndexOfTerm(rtc));

    ASSERT_EQ(-1, rma.getIndexOfTerm(rtd));
}

TEST_F(CRecursionMapTest, GetMaxOrder)
{
    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2, 0);

    CRecursionTerm rtc({"Kinetic Energy"}, 1, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 3);

    CRecursionTerm rtd({"Kinetic Energy"}, 2, false, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2, 8);

    CRecursionMap rma({rta, rtb, rtc, rtd}, {0, 18, 828, 2300}, recblock::cc, 1);

    ASSERT_EQ(5, rma.getMaxOrder({"Overlap"}, {2, 3, 4, 5}, {1, 0, 2, 3}, 1, 2));

    ASSERT_EQ(0, rma.getMaxOrder({"Kinetic Energy"}, {1, 2, 6, 3}, {3, 1, 7, 3}, 1, 2));

    ASSERT_EQ(8, rma.getMaxOrder({"Kinetic Energy"}, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 2));

    ASSERT_EQ(-1, rma.getMaxOrder({"Kinetic Energy"}, {1, 2, 6, 3}, {2, 2, 7, 3}, 1, 1));
}
