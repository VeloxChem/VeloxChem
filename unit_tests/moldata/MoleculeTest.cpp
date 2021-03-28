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

#include "MoleculeTest.hpp"

#include "CheckFunctions.hpp"
#include "CoordinationNumber.hpp"
#include "DispersionModel.hpp"
#include "Molecule.hpp"
#include "MoleculeSetter.hpp"
#include "PartialCharges.hpp"

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

TEST_F(CMoleculeTest, GetCoordinationNumber)
{
    CMolecule mol = vlxmol::getMoleculeNH3CH4();

    std::vector<double> refcn({2.984, 0.996, 0.996, 0.996, 3.966, 0.995, 0.995, 0.995, 0.995});

    auto cn = coordnum::getCoordinationNumber(mol);

    vlxtest::compare(refcn, cn.data(), 1.0e-3);
}

TEST_F(CMoleculeTest, GetPartialCharges)
{
    CMolecule mol = vlxmol::getMoleculeNH3CH4();

    std::vector<double> refchg({-0.835, 0.257, 0.255, 0.255, -0.320, 0.093, 0.097, 0.103, 0.094});

    auto chg = parchg::getPartialCharges(mol, 0.0);

    vlxtest::compare(refchg, chg.data(), 1.0e-3);
}

TEST_F(CMoleculeTest, DispersionModel)
{
    CMolecule mol = vlxmol::getMoleculeNH3CH4();

    CDispersionModel disp;

    disp.compute(mol, "B3LYP");

    auto e = disp.getEnergy();

    auto g = disp.getGradient();

    double refEnergy = -0.00242153;

    std::vector<double> refGradient({-0.1299360785248e-03, 0.2173451050590e-03, -0.3709704540840e-05, 0.3994969804870e-05,  -0.4278600323727e-04,
                                     -0.3004853785695e-05, 0.2248231831000e-04, 0.4826552264307e-04,  -0.4026908304668e-04, 0.2585427033048e-04,
                                     0.3687322138623e-04,  0.3605346888461e-04, 0.3668558637179e-04,  -0.1301671081015e-03, 0.5463254511935e-05,
                                     -0.3229412701673e-05, 0.4922085484071e-05, 0.5884251321327e-05,  0.1936253825266e-04,  -0.5468305617267e-04,
                                     0.4693862097277e-05,  0.1839250629302e-04, -0.6545014186048e-04, 0.3934710919238e-05,  0.6393301863664e-05,
                                     -0.1431962520046e-04, -0.9045906361170e-05});

    ASSERT_NEAR(refEnergy, e, 1.0e-8);

    vlxtest::compare(refGradient, g.values());
}
