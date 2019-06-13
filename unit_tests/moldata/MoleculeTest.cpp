//
//                       V.E.L.O.X. C.H.E.M. X
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2019 by Velox Chem X developers. All rights reserved.

#include "MoleculeTest.hpp"

#include "CheckFunctions.hpp"
#include "Molecule.hpp"
#include "MoleculeSetter.hpp"

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

TEST_F(CMoleculeTest, GetNumberOfAtomsWithIdElemental)
{
    CMolecule mol = vlxmol::getMoleculeLiH();

    ASSERT_EQ(1, mol.getNumberOfAtoms(1));

    ASSERT_EQ(1, mol.getNumberOfAtoms(3));

    ASSERT_EQ(0, mol.getNumberOfAtoms(2));
}

TEST_F(CMoleculeTest, GetNumberOfAtomsWithIdElementalAndAtomsList)
{
    CMolecule mol = vlxmol::getMoleculeLiH();

    ASSERT_EQ(0, mol.getNumberOfAtoms(0, 1, 1));

    ASSERT_EQ(1, mol.getNumberOfAtoms(0, 1, 3));

    ASSERT_EQ(0, mol.getNumberOfAtoms(0, 1, 2));

    ASSERT_EQ(1, mol.getNumberOfAtoms(1, 1, 1));

    ASSERT_EQ(0, mol.getNumberOfAtoms(1, 1, 3));

    ASSERT_EQ(0, mol.getNumberOfAtoms(1, 1, 2));

    ASSERT_EQ(1, mol.getNumberOfAtoms(0, 2, 1));

    ASSERT_EQ(1, mol.getNumberOfAtoms(0, 2, 3));

    ASSERT_EQ(0, mol.getNumberOfAtoms(0, 2, 2));
}

TEST_F(CMoleculeTest, GetElementalComposition)
{
    CMolecule mol = vlxmol::getMoleculeLiH();

    ASSERT_EQ(std::set<int32_t>({1, 3}), mol.getElementalComposition());

    mol = vlxmol::getMoleculeEmpty();

    ASSERT_EQ(std::set<int32_t>(), mol.getElementalComposition());
}

TEST_F(CMoleculeTest, GetNumberOfElectrons)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    ASSERT_EQ(3, mol.getNumberOfElectrons());
}

TEST_F(CMoleculeTest, GetIdsElemental)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    vlxtest::compare({3, 1}, mol.getIdsElemental());
}

TEST_F(CMoleculeTest, GetCoordinatesX)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    vlxtest::compare({0.0, 0.0}, mol.getCoordinatesX());
}

TEST_F(CMoleculeTest, GetCoordinatesY)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    vlxtest::compare({0.0, 0.0}, mol.getCoordinatesY());
}

TEST_F(CMoleculeTest, GetCoordinatesZ)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    vlxtest::compare({0.0, 1.2}, mol.getCoordinatesZ());
}

TEST_F(CMoleculeTest, GetCoordinates)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    CMemBlock2D<double> coords({0.0, 0.0, 0.0, 0.0, 0.0, 1.2}, 2, 3);

    ASSERT_EQ(mol.getCoordinates(), coords);
}

TEST_F(CMoleculeTest, GetCharges)
{
    CMolecule mol = vlxmol::getMoleculeLiHCation();

    CMemBlock<double> chrg({3.0, 1.0});

    ASSERT_EQ(mol.getCharges(), chrg);
}

TEST_F(CMoleculeTest, GetMinDistances)
{
    CMolecule mol = vlxmol::getMoleculeH2O();

    CMemBlock<double> mdist({1.78044938147649, 1.78044938147649, 1.78044938147649});

    ASSERT_EQ(mdist, mol.getMinDistances());
}

TEST_F(CMoleculeTest, GetNuclearRepulsionEnergy)
{
    CMolecule mol = vlxmol::getMoleculeH2O();

    ASSERT_NEAR(9.34363815797054450919, mol.getNuclearRepulsionEnergy(), 1.0e-13);
}

TEST_F(CMoleculeTest, GetLabel)
{
    CMolecule mol = vlxmol::getMoleculeH2O();

    ASSERT_EQ(std::string("O"), mol.getLabel(0));

    ASSERT_EQ(std::string("H"), mol.getLabel(1));

    ASSERT_EQ(std::string("H"), mol.getLabel(2));

    ASSERT_EQ(std::string(), mol.getLabel(3));
}

TEST_F(CMoleculeTest, GetSubMolecule)
{
    CMolecule empty_molecule;

    CMolecule water = vlxmol::getMoleculeH2O();

    CMolecule dimer = vlxmol::getMoleculeH2ODimer();

    ASSERT_EQ(dimer.getSubMolecule(-1, 3), empty_molecule);

    ASSERT_EQ(dimer.getSubMolecule(0, -1), empty_molecule);

    ASSERT_EQ(dimer.getSubMolecule(0, 0), empty_molecule);

    ASSERT_EQ(dimer.getSubMolecule(0, 7), empty_molecule);

    ASSERT_EQ(dimer.getSubMolecule(0, 3), water);

    ASSERT_EQ(dimer.getSubMolecule(0, 6), dimer);
}

TEST_F(CMoleculeTest, CombineMolecule)
{
    CMolecule empty_molecule;

    CMolecule dimer = vlxmol::getMoleculeH2ODimer();

    CMolecule h2o_1 = dimer.getSubMolecule(0, 3);

    CMolecule h2o_2 = dimer.getSubMolecule(3, 3);

    CMolecule h2o_11(h2o_1, empty_molecule);

    ASSERT_EQ(h2o_11, h2o_1);

    CMolecule h2o_dimer(h2o_1, h2o_2);

    ASSERT_EQ(h2o_dimer, dimer);
}
