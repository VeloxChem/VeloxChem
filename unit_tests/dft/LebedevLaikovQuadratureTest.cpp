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

#include "LebedevLaikovQuadratureTest.hpp"

#include <vector>

#include "CheckFunctions.hpp"
#include "LebedevLaikovQuadrature.hpp"

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith6Points)
{
    CLebedevLaikovQuadrature lquad(6);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(6, qpoints.size(0));

    ASSERT_EQ(6, qpoints.size(1));

    ASSERT_EQ(6, qpoints.size(2));

    ASSERT_EQ(6, qpoints.size(3));

    std::vector<double> coordx({1.0, -1.0, 0.0, 0.0, 0.0, 0.0});

    vlxtest::compare(coordx, qpoints.data(0));

    std::vector<double> coordy({0.0, 0.0, 1.0, -1.0, 0.0, 0.0});

    vlxtest::compare(coordy, qpoints.data(1));

    std::vector<double> coordz({0.0, 0.0, 0.0, 0.0, 1.0, -1.0});

    vlxtest::compare(coordz, qpoints.data(2));

    std::vector<double> weights(
        {0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667, 0.1666666666666667});

    vlxtest::compare(weights, qpoints.data(3));

    vlxtest::checkNorm(qpoints.data(3), 6);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith50Points)
{
    CLebedevLaikovQuadrature lquad(50);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(50, qpoints.size(0));

    ASSERT_EQ(50, qpoints.size(1));

    ASSERT_EQ(50, qpoints.size(2));

    ASSERT_EQ(50, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 50);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith110Points)
{
    CLebedevLaikovQuadrature lquad(110);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(110, qpoints.size(0));

    ASSERT_EQ(110, qpoints.size(1));

    ASSERT_EQ(110, qpoints.size(2));

    ASSERT_EQ(110, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 110);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith194Points)
{
    CLebedevLaikovQuadrature lquad(194);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(194, qpoints.size(0));

    ASSERT_EQ(194, qpoints.size(1));

    ASSERT_EQ(194, qpoints.size(2));

    ASSERT_EQ(194, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 194);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith302Points)
{
    CLebedevLaikovQuadrature lquad(302);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(302, qpoints.size(0));

    ASSERT_EQ(302, qpoints.size(1));

    ASSERT_EQ(302, qpoints.size(2));

    ASSERT_EQ(302, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 302);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith434Points)
{
    CLebedevLaikovQuadrature lquad(434);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(434, qpoints.size(0));

    ASSERT_EQ(434, qpoints.size(1));

    ASSERT_EQ(434, qpoints.size(2));

    ASSERT_EQ(434, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 434);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith590Points)
{
    CLebedevLaikovQuadrature lquad(590);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(590, qpoints.size(0));

    ASSERT_EQ(590, qpoints.size(1));

    ASSERT_EQ(590, qpoints.size(2));

    ASSERT_EQ(590, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 590);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith770Points)
{
    CLebedevLaikovQuadrature lquad(770);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(770, qpoints.size(0));

    ASSERT_EQ(770, qpoints.size(1));

    ASSERT_EQ(770, qpoints.size(2));

    ASSERT_EQ(770, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 770);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith974Points)
{
    CLebedevLaikovQuadrature lquad(974);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(974, qpoints.size(0));

    ASSERT_EQ(974, qpoints.size(1));

    ASSERT_EQ(974, qpoints.size(2));

    ASSERT_EQ(974, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 974);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith2030Points)
{
    CLebedevLaikovQuadrature lquad(2030);

    auto qpoints = lquad.generate();

    ASSERT_EQ(4, qpoints.blocks());

    ASSERT_EQ(2030, qpoints.size(0));

    ASSERT_EQ(2030, qpoints.size(1));

    ASSERT_EQ(2030, qpoints.size(2));

    ASSERT_EQ(2030, qpoints.size(3));

    vlxtest::checkNorm(qpoints.data(3), 2030);
}

TEST_F(CLebedevLaikovQuadratureTest, ConstructorWith100Points)
{
    CLebedevLaikovQuadrature lquad(100);

    auto qpoints = lquad.generate();

    ASSERT_EQ(0, qpoints.blocks());

    ASSERT_EQ(0, qpoints.size(0));

    ASSERT_EQ(0, qpoints.size(1));

    ASSERT_EQ(0, qpoints.size(2));

    ASSERT_EQ(0, qpoints.size(3));
}
