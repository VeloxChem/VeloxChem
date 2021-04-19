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

#include "MolecularOrbitalsTest.hpp"

#include "CheckFunctions.hpp"
#include "MolecularBasisSetter.hpp"
#include "MolecularOrbitals.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CMolecularOrbitalsTest, DefaultConstructor)
{
    CMolecularOrbitals moa;

    CMolecularOrbitals mob({}, {}, molorb::rest);

    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    CMolecularOrbitals mob(moa);

    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    CMolecularOrbitals mob(CMolecularOrbitals({ma, mb}, {ea, eb}, molorb::unrest));

    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    std::vector<double> ea({1.0, 2.0, 4.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    CMolecularOrbitals mob = moa;

    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    std::vector<double> ea({1.0, 2.0, 4.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    CMolecularOrbitals mob = CMolecularOrbitals({ma}, {ea}, molorb::rest);

    ASSERT_EQ(moa, mob);
}

TEST_F(CMolecularOrbitalsTest, Insert)
{
    auto mbas = vlxbas::getMolecularBasisForLiH();

    auto vbas = mbas.reduceToValenceBasis();

    auto mlih = vlxmol::getTestLiH();

    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 2.0, 3.0, 1.0, 3.0, -6.0, -8.0}, 5, 3);

    std::vector<double> ea({1.0, 2.0, 4.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    auto xbas = vlxbas::getMolecularBasisForLiHX();

    auto mob = moa.insert(mlih, xbas, vbas);

    CDenseMatrix mc({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0, 0.0, 0.0, 0.0, 2.0, 3.0, 1.0, 3.0, -6.0, -8.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                     0.0, 0.0,  0.0,  0.0,  0.0, 0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,  0.0,  0.0, 0.0, 0.0, 0.0},
                    15,
                    3);

    CMolecularOrbitals moc({mc}, {ea}, molorb::rest);

    ASSERT_EQ(moc, mob);
}

TEST_F(CMolecularOrbitalsTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    ASSERT_EQ(3, moa.getNumberOfRows());
}

TEST_F(CMolecularOrbitalsTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    ASSERT_EQ(2, moa.getNumberOfColumns());
}

TEST_F(CMolecularOrbitalsTest, AlphaOrbitals)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    const CMolecularOrbitals mob({ma, mb}, {ea, eb}, molorb::unrest);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, moa.alphaOrbitals());

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, mob.alphaOrbitals());
}

TEST_F(CMolecularOrbitalsTest, BetaOrbitals)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    const CMolecularOrbitals mob({ma, mb}, {ea, eb}, molorb::unrest);

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, moa.betaOrbitals());

    vlxtest::compare({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, mob.betaOrbitals());
}

TEST_F(CMolecularOrbitalsTest, AlphaEnergies)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    const CMolecularOrbitals mob({mb}, {eb}, molorb::rest);

    vlxtest::compare(ea, moa.alphaEnergies());

    vlxtest::compare(eb, mob.alphaEnergies());
}

TEST_F(CMolecularOrbitalsTest, BetaEnergies)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    const CMolecularOrbitals mob({mb}, {eb}, molorb::rest);

    vlxtest::compare(eb, moa.betaEnergies());

    vlxtest::compare(eb, mob.betaEnergies());
}

TEST_F(CMolecularOrbitalsTest, AlphaOrbitalsWithRange)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    const CMolecularOrbitals mob({ma}, {ea}, molorb::rest);

    ASSERT_EQ(CDenseMatrix({1.0, -2.0, 6.0}, 3, 1), moa.alphaOrbitals(0, 1));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3), moa.alphaOrbitals(0, 3));

    ASSERT_EQ(CDenseMatrix({-1.0, -3.0, 5.0, 4.0, 4.0, -4.0}, 3, 2), moa.alphaOrbitals(1, 2));

    ASSERT_EQ(CDenseMatrix({1.0, -2.0, 6.0}, 3, 1), mob.alphaOrbitals(0, 1));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3), mob.alphaOrbitals(0, 3));

    ASSERT_EQ(CDenseMatrix({-1.0, -3.0, 5.0, 4.0, 4.0, -4.0}, 3, 2), mob.alphaOrbitals(1, 2));
}

TEST_F(CMolecularOrbitalsTest, BetaOrbitalsWithRange)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    const CMolecularOrbitals mob({mb}, {eb}, molorb::rest);

    ASSERT_EQ(CDenseMatrix({-1.0, -2.0, 4.0}, 3, 1), moa.betaOrbitals(1, 1));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2), moa.betaOrbitals(0, 2));

    ASSERT_EQ(CDenseMatrix({1.0, -3.0, 5.0}, 3, 1), moa.betaOrbitals(0, 1));

    ASSERT_EQ(CDenseMatrix({-1.0, -2.0, 4.0}, 3, 1), mob.betaOrbitals(1, 1));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2), mob.betaOrbitals(0, 2));

    ASSERT_EQ(CDenseMatrix({1.0, -3.0, 5.0}, 3, 1), mob.betaOrbitals(0, 1));
}

TEST_F(CMolecularOrbitalsTest, AlphaOrbitalsWithList)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    const CMolecularOrbitals mob({ma}, {ea}, molorb::rest);

    ASSERT_EQ(CDenseMatrix({1.0, -2.0, 6.0}, 3, 1), moa.alphaOrbitals(0, 1));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3), moa.alphaOrbitals({0, 1, 2}));

    ASSERT_EQ(CDenseMatrix({-1.0, -3.0, 5.0, 4.0, 4.0, -4.0}, 3, 2), moa.alphaOrbitals({1, 2}));

    ASSERT_EQ(CDenseMatrix({1.0, -2.0, 6.0}, 3, 1), mob.alphaOrbitals(0, 1));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3), mob.alphaOrbitals({0, 1, 2}));

    ASSERT_EQ(CDenseMatrix({-1.0, -3.0, 5.0, 4.0, 4.0, -4.0}, 3, 2), mob.alphaOrbitals({1, 2}));
}

TEST_F(CMolecularOrbitalsTest, BetaOrbitalsWithList)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    const CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    const CMolecularOrbitals mob({mb}, {eb}, molorb::rest);

    ASSERT_EQ(CDenseMatrix({-1.0, -2.0, 4.0}, 3, 1), moa.betaOrbitals({1}));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2), moa.betaOrbitals({0, 1}));

    ASSERT_EQ(CDenseMatrix({1.0, -3.0, 5.0}, 3, 1), moa.betaOrbitals({0}));

    ASSERT_EQ(CDenseMatrix({-1.0, -2.0, 4.0}, 3, 1), mob.betaOrbitals({1}));

    ASSERT_EQ(CDenseMatrix({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2), mob.betaOrbitals({0, 1}));

    ASSERT_EQ(CDenseMatrix({1.0, -3.0, 5.0}, 3, 1), mob.betaOrbitals({0}));
}

TEST_F(CMolecularOrbitalsTest, GetAODensityForRestrictedCase)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    std::vector<double> ea({1.0, 2.0, 4.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    CDenseMatrix ref01mat({1.0, -2.0, 6.0, -2.0, 4.0, -12.0, 6.0, -12.0, 36.0}, 3, 3);

    ASSERT_EQ(moa.getAODensity(2), CAODensityMatrix({ref01mat}, denmat::rest));

    CDenseMatrix ref02mat({2.0, -7.0, 2.0, -7.0, 29.0, 8.0, 2.0, 8.0, 52.0}, 3, 3);

    ASSERT_EQ(moa.getAODensity(4), CAODensityMatrix({ref02mat}, denmat::rest));

    CDenseMatrix ref03mat({11.0, -19.0, 14.0, -19.0, 45.0, -8.0, 14.0, -8.0, 68}, 3, 3);

    ASSERT_EQ(moa.getAODensity(6), CAODensityMatrix({ref03mat}, denmat::rest));
}

TEST_F(CMolecularOrbitalsTest, GetAODensityForUnrestrictedCase)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CDenseMatrix mb({1.0, -1.0, -3.0, -2.0, 5.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0, 4.0});

    std::vector<double> eb({3.0, 5.0});

    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    CDenseMatrix ref01mata({1.0, -2.0, 6.0, -2.0, 4.0, -12.0, 6.0, -12.0, 36.0}, 3, 3);

    CDenseMatrix ref02mata({2.0, -7.0, 2.0, -7.0, 29.0, 8.0, 2.0, 8.0, 52.0}, 3, 3);

    CDenseMatrix ref03mata({11.0, -19.0, 14.0, -19.0, 45.0, -8.0, 14.0, -8.0, 68}, 3, 3);

    CDenseMatrix ref01matb({1.0, -3.0, 5.0, -3.0, 9.0, -15.0, 5.0, -15.0, 25.0}, 3, 3);

    CDenseMatrix ref02matb({2.0, -1.0, 1.0, -1.0, 13.0, -23.0, 1.0, -23.0, 41.0}, 3, 3);

    ASSERT_EQ(moa.getAODensity(2, 1), CAODensityMatrix({ref02mata, ref01matb}, denmat::unrest));

    ASSERT_EQ(moa.getAODensity(1, 1), CAODensityMatrix({ref01mata, ref01matb}, denmat::unrest));

    ASSERT_EQ(moa.getAODensity(2, 2), CAODensityMatrix({ref02mata, ref02matb}, denmat::unrest));

    ASSERT_EQ(moa.getAODensity(3, 2), CAODensityMatrix({ref03mata, ref02matb}, denmat::unrest));
}

TEST_F(CMolecularOrbitalsTest, getRestrictedPairDensitySinglePair)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    std::vector<double> ea({1.0, 2.0, 4.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    CDenseMatrix ref00mat({1.0, -2.0, 6.0, -2.0, 4.0, -12.0, 6.0, -12.0, 36.0}, 3, 3);

    ASSERT_EQ(moa.getRestrictedPairDensity(0, 0), CAODensityMatrix({ref00mat}, denmat::rmoij));

    CDenseMatrix ref11mat({1.0, -5.0, -4.0, -5.0, 25.0, 20.0, -4.0, 20.0, 16.0}, 3, 3);

    ASSERT_EQ(moa.getRestrictedPairDensity(1, 1), CAODensityMatrix({ref11mat}, denmat::rmoij));

    CDenseMatrix ref01mat({-1.0, 5.0, 4.0, 2.0, -10.0, -8.0, -6.0, 30.0, 24.0}, 3, 3);

    ASSERT_EQ(moa.getRestrictedPairDensity(0, 1), CAODensityMatrix({ref01mat}, denmat::rmoij));

    CDenseMatrix ref10mat({-1.0, 2.0, -6.0, 5.0, -10.0, 30.0, 4.0, -8.0, 24.0}, 3, 3);

    ASSERT_EQ(moa.getRestrictedPairDensity(1, 0), CAODensityMatrix({ref10mat}, denmat::rmoij));
}

TEST_F(CMolecularOrbitalsTest, getRestrictedPairDensityVectorOfPairs)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    std::vector<double> ea({1.0, 2.0, 4.0});

    CMolecularOrbitals moa({ma}, {ea}, molorb::rest);

    CDenseMatrix ref00mat({1.0, -2.0, 6.0, -2.0, 4.0, -12.0, 6.0, -12.0, 36.0}, 3, 3);

    CDenseMatrix ref11mat({1.0, -5.0, -4.0, -5.0, 25.0, 20.0, -4.0, 20.0, 16.0}, 3, 3);

    CDenseMatrix ref01mat({-1.0, 5.0, 4.0, 2.0, -10.0, -8.0, -6.0, 30.0, 24.0}, 3, 3);

    CDenseMatrix ref10mat({-1.0, 2.0, -6.0, 5.0, -10.0, 30.0, 4.0, -8.0, 24.0}, 3, 3);

    ASSERT_EQ(moa.getRestrictedPairDensity({0, 1, 0, 1}, {0, 1, 1, 0}), CAODensityMatrix({ref00mat, ref11mat, ref01mat, ref10mat}, denmat::rmoij));
}

TEST_F(CMolecularOrbitalsTest, Transform)
{
    CDenseMatrix ma({1.0, -1.0, -2.0, 5.0, 6.0, 4.0}, 3, 2);

    std::vector<double> ea({1.0, 2.0});

    CDenseMatrix mb({2.0, 4.0, -1.0, 2.0, 3.0, 5.0}, 3, 2);

    std::vector<double> eb({1.0, 1.0});

    CMolecularOrbitals moa({ma, mb}, {ea, eb}, molorb::unrest);

    CDenseMatrix aomat({2.0, 3.0, -1.0, 4.0, 1.0, 2.0, 5.0, -2.0, 4.0}, 3, 3);

    CDenseMatrix refmaa({160.0, -3.0, 212.0, 40.0}, 2, 2);

    ASSERT_EQ(refmaa, moa.transform(aomat, szblock::aa));

    CDenseMatrix refmab({116.0, 169.0, 163.0, 275.0}, 2, 2);

    ASSERT_EQ(refmab, moa.transform(aomat, szblock::ab));

    CDenseMatrix refmba({65.0, 12.0, 153.0, 59.0}, 2, 2);

    ASSERT_EQ(refmba, moa.transform(aomat, szblock::ba));

    CDenseMatrix refmbb({55.0, 98.0, 138.0, 272.0}, 2, 2);

    ASSERT_EQ(refmbb, moa.transform(aomat, szblock::bb));
}
