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

#include "AOFockMatrixTest.hpp"

#include "AOFockMatrix.hpp"
#include "CheckFunctions.hpp"

TEST_F(CAOFockMatrixTest, DefaultConstructor)
{
    CAOFockMatrix fmata;

    CAOFockMatrix fmatb({}, {}, {}, {});

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, ConstructorWithDensity)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAODensityMatrix dmata({ma, mb}, denmat::rest);

    CAOFockMatrix fmata(dmata);

    ASSERT_EQ(2, fmata.getNumberOfFockMatrices());

    ASSERT_EQ(2, fmata.getNumberOfRows(0));

    ASSERT_EQ(2, fmata.getNumberOfRows(1));

    ASSERT_EQ(2, fmata.getNumberOfColumns(0));

    ASSERT_EQ(2, fmata.getNumberOfColumns(1));

    ASSERT_EQ(0, fmata.getDensityIdentifier(0));

    ASSERT_EQ(1, fmata.getDensityIdentifier(1));

    ASSERT_TRUE(fockmat::restjk == fmata.getFockType(0));

    ASSERT_TRUE(fockmat::restjk == fmata.getFockType(1));

    CAODensityMatrix dmatb({ma, mb}, denmat::unrest);

    CAOFockMatrix fmatb(dmatb);

    ASSERT_EQ(1, fmatb.getNumberOfFockMatrices());

    ASSERT_EQ(2, fmatb.getNumberOfRows(0));

    ASSERT_EQ(2, fmatb.getNumberOfColumns(0));

    ASSERT_EQ(0, fmatb.getDensityIdentifier(0));

    ASSERT_TRUE(fockmat::unrestjk == fmatb.getFockType(0));
}

TEST_F(CAOFockMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1});

    CAOFockMatrix fmatb(fmata);

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1});

    CAOFockMatrix fmatb(CAOFockMatrix({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1}));

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1});

    CAOFockMatrix fmatb = fmata;

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1});

    CAOFockMatrix fmatb = CAOFockMatrix({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, SetFockType)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restj, fockmat::restj}, {1.0, 1.0}, {0, 1});

    fmata.setFockType(fockmat::restjk, 0);

    fmata.setFockType(fockmat::restjk, 1);

    CAOFockMatrix fmatb({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);

    CAOFockMatrix fmatc({ma, mb}, {fockmat::unrestj, fockmat::unrestj}, {1.0, 1.0}, {0, 0});

    fmatc.setFockType(fockmat::unrestjk, 0, "alpha");

    fmatc.setFockType(fockmat::unrestjk, 0, "beta");

    CAOFockMatrix fmatd({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(fmatc, fmatd);
}

TEST_F(CAOFockMatrixTest, SetFockScaleFactor)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    fmata.setFockScaleFactor(2.0, 0);

    fmata.setFockScaleFactor(2.0, 1);

    CAOFockMatrix fmatb({ma, mb}, {fockmat::restjk, fockmat::restjk}, {2.0, 2.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);

    CAOFockMatrix fmatc({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    fmatc.setFockScaleFactor(0.5, 0, "alpha");

    fmatc.setFockScaleFactor(0.5, 0, "beta");

    CAOFockMatrix fmatd({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {0.5, 0.5}, {0, 0});

    ASSERT_EQ(fmatc, fmatd);
}

TEST_F(CAOFockMatrixTest, Zero)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    fmata.zero();

    ma.zero();

    mb.zero();

    CAOFockMatrix fmatb({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, Symmetrize)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    fmata.symmetrize();

    CDenseMatrix mc({2.0, 0.3, 0.3, 4.0}, 2, 2);

    CDenseMatrix md({4.0, 0.7, 0.7, 6.0}, 2, 2);

    CAOFockMatrix fmatb({mc, md}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);
}

TEST_F(CAOFockMatrixTest, Add)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    CAOFockMatrix fmatb({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    fmatb.add(fmata);

    CDenseMatrix mc({2.0, 0.2, 0.4, 4.0}, 2, 2);

    CDenseMatrix md({4.0, 0.6, 0.8, 6.0}, 2, 2);

    CAOFockMatrix fmatc({mc, md}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmatb, fmatc);
}

TEST_F(CAOFockMatrixTest, Scale)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    fmata.scale(2.0, 0);

    fmata.scale(3.0, 1);

    CDenseMatrix mc({2.0, 0.2, 0.4, 4.0}, 2, 2);

    CDenseMatrix md({6.0, 0.9, 1.2, 9.0}, 2, 2);

    CAOFockMatrix fmatb({mc, md}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);

    CAOFockMatrix fmatc({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    fmatc.scale(2.0, 0, "alpha");

    fmatc.scale(3.0, 0, "beta");

    CAOFockMatrix fmatd({mc, md}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(fmatc, fmatd);
}

TEST_F(CAOFockMatrixTest, GetNumberOfFockMatrices)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(2, fmata.getNumberOfFockMatrices());

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(1, fmatb.getNumberOfFockMatrices());
}

TEST_F(CAOFockMatrixTest, GetNumberOfRows)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(2, fmata.getNumberOfRows(0));

    ASSERT_EQ(2, fmata.getNumberOfRows(1));

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(2, fmatb.getNumberOfRows(0));
}

TEST_F(CAOFockMatrixTest, GetNumberOfColumns)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(2, fmata.getNumberOfColumns(0));

    ASSERT_EQ(2, fmata.getNumberOfColumns(1));

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(2, fmatb.getNumberOfColumns(0));
}

TEST_F(CAOFockMatrixTest, GetNumberOfElements)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(4, fmata.getNumberOfElements(0));

    ASSERT_EQ(4, fmata.getNumberOfElements(1));

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(4, fmatb.getNumberOfElements(0));
}

TEST_F(CAOFockMatrixTest, GetFockConstant)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    const CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    vlxtest::compare({1.0, 0.1, 0.2, 2.0}, fmata.getFock(0));

    vlxtest::compare({2.0, 0.3, 0.4, 3.0}, fmata.getFock(1));

    const CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    vlxtest::compare({1.0, 0.1, 0.2, 2.0}, fmatb.getFock(0, "alpha"));

    vlxtest::compare({2.0, 0.3, 0.4, 3.0}, fmatb.getFock(0, "beta"));
}

TEST_F(CAOFockMatrixTest, GetFock)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    vlxtest::compare({1.0, 0.1, 0.2, 2.0}, fmata.getFock(0));

    vlxtest::compare({2.0, 0.3, 0.4, 3.0}, fmata.getFock(1));

    auto pdat = fmata.getFock(1);

    pdat[1] = 0.5;

    pdat[3] = 2.0;

    CDenseMatrix mc({2.0, 0.5, 0.4, 2.0}, 2, 2);

    ASSERT_EQ(fmata, CAOFockMatrix({ma, mc}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1}));

    CAOFockMatrix fmatc({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    vlxtest::compare({1.0, 0.1, 0.2, 2.0}, fmatc.getFock(0, "alpha"));

    vlxtest::compare({2.0, 0.3, 0.4, 3.0}, fmatc.getFock(0, "beta"));

    auto pdatb = fmatc.getFock(0, "beta");

    pdatb[1] = 0.5;

    pdatb[3] = 2.0;

    ASSERT_EQ(fmatc, CAOFockMatrix({ma, mc}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0}));
}

TEST_F(CAOFockMatrixTest, GetReferenceToFock)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    const CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(ma, fmata.getReferenceToFock(0));

    ASSERT_EQ(mb, fmata.getReferenceToFock(1));

    const CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(ma, fmatb.getReferenceToFock(0, "alpha"));

    ASSERT_EQ(mb, fmatb.getReferenceToFock(0, "beta"));
}

TEST_F(CAOFockMatrixTest, GetFockType)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_TRUE(fmata.getFockType(0) == fockmat::restjk);

    ASSERT_TRUE(fmata.getFockType(1) == fockmat::restjk);

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestj, fockmat::unrestj}, {1.0, 1.0}, {0, 0});

    ASSERT_TRUE(fmatb.getFockType(0, "alpha") == fockmat::unrestj);

    ASSERT_TRUE(fmatb.getFockType(0, "beta") == fockmat::unrestj);
}

TEST_F(CAOFockMatrixTest, GetScaleFactor)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_NEAR(1.0, fmata.getScaleFactor(0), 1.0e-13);

    ASSERT_NEAR(1.0, fmata.getScaleFactor(1), 1.0e-13);

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestj, fockmat::unrestj}, {2.0, 2.0}, {0, 0});

    ASSERT_NEAR(2.0, fmatb.getScaleFactor(0, "alpha"), 1.0e-13);

    ASSERT_NEAR(2.0, fmatb.getScaleFactor(0, "beta"), 1.0e-13);
}

TEST_F(CAOFockMatrixTest, GetDensityIdentifier)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(0, fmata.getDensityIdentifier(0));

    ASSERT_EQ(1, fmata.getDensityIdentifier(1));

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(0, fmatb.getDensityIdentifier(0, "alpha"));

    ASSERT_EQ(0, fmatb.getDensityIdentifier(0, "beta"));
}

TEST_F(CAOFockMatrixTest, AddOneElectronProperty)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    fmata.addOneElectronMatrix(ma, 0);

    CDenseMatrix mc({2.0, 0.2, 0.4, 4.0}, 2, 2);

    CAOFockMatrix fmatb({mc, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_EQ(fmata, fmatb);

    CAOFockMatrix fmatc({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    fmatc.addOneElectronMatrix(ma, 0, "alpha");

    fmatc.addOneElectronMatrix(mb, 0, "beta");

    CDenseMatrix md({4.0, 0.6, 0.8, 6.0}, 2, 2);

    CAOFockMatrix fmatd({mc, md}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_EQ(fmatc, fmatd);
}

TEST_F(CAOFockMatrixTest, IsSymmetric)
{
    CDenseMatrix ma({1.0, 0.1, 0.2, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.3, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmata({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    ASSERT_TRUE(fmata.isSymmetric(0));

    ASSERT_TRUE(fmata.isSymmetric(1));

    CAOFockMatrix fmatb({ma, mb}, {fockmat::rgenjk, fockmat::rgenjk}, {1.0, 1.0}, {0, 1});

    ASSERT_FALSE(fmatb.isSymmetric(0));

    ASSERT_FALSE(fmatb.isSymmetric(1));

    CAOFockMatrix fmatc({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    ASSERT_TRUE(fmatc.isSymmetric(0));
}

TEST_F(CAOFockMatrixTest, GetElectronicEnergy)
{
    CDenseMatrix ma({1.0, 0.3, 0.3, 2.0}, 2, 2);

    CDenseMatrix mb({2.0, 0.4, 0.4, 3.0}, 2, 2);

    CAOFockMatrix fmat({ma, mb}, {fockmat::restjk, fockmat::restjk}, {1.0, 1.0}, {0, 1});

    CDenseMatrix da({1.0, 0.5, 0.5, 3.0}, 2, 2);

    CDenseMatrix db({2.0, 0.6, 0.6, 5.0}, 2, 2);

    CAODensityMatrix dmat({da, db}, denmat::rest);

    ASSERT_NEAR(7.3, fmat.getElectronicEnergy(0, dmat, 0), 1.0e-13);

    ASSERT_NEAR(12.36, fmat.getElectronicEnergy(0, dmat, 1), 1.0e-13);

    ASSERT_NEAR(11.4, fmat.getElectronicEnergy(1, dmat, 0), 1.0e-13);

    ASSERT_NEAR(19.48, fmat.getElectronicEnergy(1, dmat, 1), 1.0e-13);

    CAOFockMatrix fmatb({ma, mb}, {fockmat::unrestjk, fockmat::unrestjk}, {1.0, 1.0}, {0, 0});

    CAODensityMatrix dmatb({da, db}, denmat::unrest);

    ASSERT_NEAR(13.39, fmatb.getElectronicEnergy(0, dmatb, 0), 1.0e-13);
}
