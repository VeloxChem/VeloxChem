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

#include "RecursionFunctionTest.hpp"

#include "DummyFunctions.hpp"
#include "RecursionFunction.hpp"

TEST_F(CRecursionFunctionTest, DefaultConstructor)
{
    CRecursionFunction rfa({}, nullptr);

    CRecursionFunction rfb;

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionTest, CopyConstructor)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction rfb(rfa);

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionTest, MoveConstructor)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction rfb(CRecursionFunction({"Overlap"}, &vlxtest::dummy_func_11));

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionTest, CopyAssignment)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction rfb = rfa;

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionTest, MoveAssignment)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction rfb = CRecursionFunction({"Overlap"}, &vlxtest::dummy_func_11);

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionTest, GetLabel)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    ASSERT_EQ(std::string("Overlap"), rfa.getLabel());
}

TEST_F(CRecursionFunctionTest, Compute)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Overlap"}, 0, true, {1, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    auto rtvecs = rfa.compute(rta);

    ASSERT_EQ(2, rtvecs.size());

    ASSERT_EQ(rta, rtvecs[0]);

    ASSERT_EQ(rtb, rtvecs[1]);
}

TEST_F(CRecursionFunctionTest, IsMatch)
{
    CRecursionFunction rfa({"Overlap"}, &vlxtest::dummy_func_11);

    ASSERT_FALSE(rfa.isMatch({"Kinetic Energy"}));

    ASSERT_TRUE(rfa.isMatch({"Overlap"}));
}
