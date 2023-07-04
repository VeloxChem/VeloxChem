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

#include "GridDriverTest.hpp"

#include <cmath>

#include <mpi.h>

#include "GridDriver.hpp"
#include "MathConst.hpp"
#include "MoleculeSetter.hpp"
#include "MpiFunc.hpp"
#include "MolecularGrid.hpp"

TEST_F(CGridDriverTest, DefaultConstructor)
{
    CGridDriver gdrv(MPI_COMM_WORLD);

    gdrv.setLevel(8);

    auto mlih = vlxmol::getMoleculeLiH();

    auto mgrid = gdrv.generate(mlih);

    auto npnt = mgrid.getNumberOfGridPoints();

    EXPECT_EQ(741140, npnt);

    auto rx = mgrid.getCoordinatesX();

    auto ry = mgrid.getCoordinatesY();

    auto rz = mgrid.getCoordinatesZ();

    auto w = mgrid.getWeights();

    double fa = 0.0;

    double fb = 0.0;

    double fab = 0.0;

    for (int32_t i = 0; i < npnt; i++)
    {
        auto r2a = rx[i] * rx[i] + ry[i] * ry[i] + rz[i] * rz[i];

        fa += w[i] * std::exp(-2.3 * r2a);

        auto r2b = rx[i] * rx[i] + ry[i] * ry[i]

                   + (rz[i] - 1.20) * (rz[i] - 1.20);

        fb += w[i] * std::exp(-0.5 * r2b);

        fab += w[i] * std::exp(-2.3 * r2a) * std::exp(-0.5 * r2b);
    }

    EXPECT_NEAR(fa, 1.5963681525241624, 1.0e-10);

    EXPECT_NEAR(fb, 15.749609945385632, 1.0e-10);

    EXPECT_NEAR(fab, 0.65786017622805693, 1.0e-10);
}
