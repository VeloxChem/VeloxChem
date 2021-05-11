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

#include "ChemicalElementTest.hpp"

#include "ChemicalElement.hpp"

TEST_F(CChemicalElementTest, DefaultConstructor)
{
    CChemicalElement achem;

    CChemicalElement bchem(std::string(), 0.0, 0.0, -1);

    ASSERT_EQ(achem, bchem);
}

TEST_F(CChemicalElementTest, SetAtomType)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType({"CU"}));

    CChemicalElement bchem({"Cu"}, 29.0, 62.929598, 29);

    ASSERT_EQ(achem, bchem);

    ASSERT_FALSE(achem.setAtomType({"XZ"}));

    ASSERT_EQ(achem, bchem);
}

TEST_F(CChemicalElementTest, SetAtomTypeWithElementNumber)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType(29));

    CChemicalElement bchem({"Cu"}, 29.0, 62.929598, 29);

    ASSERT_EQ(achem, bchem);

    ASSERT_FALSE(achem.setAtomType(99));

    ASSERT_EQ(achem, bchem);
}

TEST_F(CChemicalElementTest, SetIsotope)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType(29));

    CChemicalElement bchem({"Cu"}, 29.0, 62.929598, 29);

    ASSERT_EQ(achem, bchem);

    ASSERT_TRUE(achem.setIsotope(65));

    CChemicalElement cchem({"Cu"}, 29.0, 64.927790, 29);

    ASSERT_EQ(achem, cchem);

    ASSERT_TRUE(achem.setIsotope(0));

    ASSERT_EQ(achem, bchem);

    ASSERT_FALSE(achem.setIsotope(64));

    ASSERT_EQ(achem, bchem);
}

TEST_F(CChemicalElementTest, GetName)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType(6));

    ASSERT_EQ(std::string("C"), achem.getName());
}

TEST_F(CChemicalElementTest, GetIdentifier)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType({"SC"}));

    ASSERT_EQ(21, achem.getIdentifier());
}

TEST_F(CChemicalElementTest, GetAtomicMass)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType({"H"}));

    ASSERT_NEAR(1.007825, achem.getAtomicMass(), 1.0e-13);
}

TEST_F(CChemicalElementTest, GetAtomicCharge)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType({"N"}));

    ASSERT_NEAR(7.0, achem.getAtomicCharge(), 1.0e-13);
}

TEST_F(CChemicalElementTest, GetMaxAngularMomentum)
{
    CChemicalElement achem;

    ASSERT_TRUE(achem.setAtomType({"N"}));

    ASSERT_EQ(1, achem.getMaxAngularMomentum());

    ASSERT_TRUE(achem.setAtomType({"H"}));

    ASSERT_EQ(0, achem.getMaxAngularMomentum());

    ASSERT_TRUE(achem.setAtomType({"CU"}));

    ASSERT_EQ(2, achem.getMaxAngularMomentum());

    ASSERT_TRUE(achem.setAtomType({"AU"}));

    ASSERT_EQ(3, achem.getMaxAngularMomentum());
}
