//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem X developers. All rights reserved.

#include "MoleculeTest.hpp"

#include "Molecule.hpp"
#include "MoleculeSetter.hpp"
#include "OutputStream.hpp"
#include "CheckFunctions.hpp"

TEST_F(CMoleculeTest, DefaultConstructor)
{
    CMolecule mol;

    ASSERT_EQ(mol, vlxmol::getMoleculeEmpty());
}

TEST_F(CMoleculeTest, CopyConstructor)
{
    CMolecule mola = vlxmol::getMoleculeLiHCation();

    CMolecule molb(mola);

    ASSERT_EQ(mola, molb);
}

TEST_F(CMoleculeTest, MoveConstructor)
{
    CMolecule mola = vlxmol::getMoleculeLiHCation();

    CMolecule molb(vlxmol::getMoleculeLiHCation());

    ASSERT_EQ(mola, molb);
}

TEST_F(CMoleculeTest, CopyAssignment)
{
    CMolecule mola = vlxmol::getMoleculeLiHCation();

    CMolecule molb = mola;

    ASSERT_EQ(mola, molb);
}

TEST_F(CMoleculeTest, MoveAssignment)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    ASSERT_EQ(mol, vlxmol::getMoleculeLiHCation());
}

TEST_F(CMoleculeTest, SetChargeAndSetMultiplicity)
{
    CMolecule mol = vlxmol::getMoleculeLiH();

    mol.setCharge(1.0);

    mol.setMultiplicity(2);

    ASSERT_EQ(mol, vlxmol::getMoleculeLiHCation());
}

TEST_F(CMoleculeTest, GetCharge)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();
    
    ASSERT_NEAR(1.0, mol.getCharge(), 1.0e-13);
}

TEST_F(CMoleculeTest, GetMultiplicity)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();
    
    ASSERT_EQ(2, mol.getMultiplicity());
}

TEST_F(CMoleculeTest, GetNumberOfAtoms)
{
    CMolecule mol = vlxmol::getMoleculeH2O();
    
    ASSERT_EQ(3, mol.getNumberOfAtoms());
}

TEST_F(CMoleculeTest, GetNumberOfElectrons)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();
    
    ASSERT_EQ(3, mol.getNumberOfElectrons());
}
