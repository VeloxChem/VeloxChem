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

#include "RecursionFunctionsListTest.hpp"

#include "RecursionFunctionsList.hpp"

#include "DummyFunctions.hpp"

TEST_F(CRecursionFunctionsListTest, DefaultConstructor)
{
    CRecursionFunctionsList rfa(std::vector<CRecursionFunction>({}));

    CRecursionFunctionsList rfb;

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionsListTest, CopyConstructor)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa({fa, fb, fc});

    CRecursionFunctionsList rfb(rfa);

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionsListTest, MoveConstructor)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa({fa, fb, fc});

    CRecursionFunctionsList rfb(CRecursionFunctionsList({fa, fb, fc}));

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionsListTest, CopyAssignment)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa({fa, fb, fc});

    CRecursionFunctionsList rfb = rfa;

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionsListTest, MoveAssignment)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa({fa, fb, fc});

    CRecursionFunctionsList rfb = CRecursionFunctionsList({fa, fb, fc});

    ASSERT_EQ(rfa, rfb);
}

TEST_F(CRecursionFunctionsListTest, Add)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa;

    rfa.add(fa);

    ASSERT_EQ(rfa, CRecursionFunctionsList({fa}));

    rfa.add(fa);

    ASSERT_EQ(rfa, CRecursionFunctionsList({fa}));

    rfa.add(fb);

    ASSERT_EQ(rfa, CRecursionFunctionsList({fa, fb}));

    rfa.add(fa);

    ASSERT_EQ(rfa, CRecursionFunctionsList({fa, fb}));

    rfa.add(fc);

    ASSERT_EQ(rfa, CRecursionFunctionsList({fa, fb, fc}));
}

TEST_F(CRecursionFunctionsListTest, Compute)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa({fa, fb, fc});

    CRecursionTerm rta({"Overlap"}, 0, true, {2, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    CRecursionTerm rtb({"Overlap"}, 0, true, {1, 3, 4, 5}, {1, 7, 2, 3}, 1, 2, 5);

    auto rtvecs = rfa.compute(rta);

    ASSERT_EQ(2, rtvecs.size());

    ASSERT_EQ(rta, rtvecs[0]);

    ASSERT_EQ(rtb, rtvecs[1]);
}

TEST_F(CRecursionFunctionsListTest, Find)
{
    CRecursionFunction fa({"Overlap"}, &vlxtest::dummy_func_11);

    CRecursionFunction fb({"Kinetic Energy"}, &vlxtest::dummy_func_01);

    CRecursionFunction fc({"Nuclear Potential"}, &vlxtest::dummy_func_10);

    CRecursionFunctionsList rfa({fa, fb, fc});

    ASSERT_EQ(-1, rfa.find({"Kinetic Energy Gradient"}));

    ASSERT_EQ(0, rfa.find({"Overlap"}));

    ASSERT_EQ(1, rfa.find({"Kinetic Energy"}));

    ASSERT_EQ(2, rfa.find({"Nuclear Potential"}));
}
