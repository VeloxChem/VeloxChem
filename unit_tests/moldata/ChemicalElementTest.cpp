//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

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

