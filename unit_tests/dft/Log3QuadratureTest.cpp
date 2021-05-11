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

#include "Log3QuadratureTest.hpp"

#include "CheckFunctions.hpp"
#include "Log3Quadrature.hpp"

TEST_F(CLog3QuadratureTest, ConstructorWith5Points)
{
    CLog3Quadrature rquad(5, 1);

    auto qpoints = rquad.generate();

    ASSERT_EQ(2, qpoints.blocks());

    ASSERT_EQ(5, qpoints.size(0));

    ASSERT_EQ(5, qpoints.size(1));

    std::vector<double> rad({4.536298573026920000, 2.040679201001270000, 0.800000000000003000, 0.219058105426335000, 0.023957400390054500});

    vlxtest::compare(rad, qpoints.data(0));

    std::vector<double> wgt({3.661009388224700000, 1.705129857584120000, 0.855642097864144000, 0.349387213052902000, 0.076565310395870600});

    vlxtest::compare(wgt, qpoints.data(1));
}
