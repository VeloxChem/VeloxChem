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

#include "ElectronicPotentialMatrixTest.hpp"

#include "ElectronicPotentialMatrix.hpp"

TEST_F(CElectronicPotentialMatrixTest, DefaultConstructor)
{
    CElectronicPotentialMatrix smata;

    CElectronicPotentialMatrix smatb(CDenseMatrix(0, 0));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb(smata);

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveConstructor)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb(CElectronicPotentialMatrix({ma}));

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, CopyAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb = smata;

    ASSERT_EQ(smata, smatb);
}

TEST_F(CElectronicPotentialMatrixTest, MoveAssignment)
{
    CDenseMatrix ma({1.0, -1.0, -3.0, -2.0, 5.0, 4.0, 6.0, 4.0, -4.0}, 3, 3);

    CElectronicPotentialMatrix smata(ma);

    CElectronicPotentialMatrix smatb = CElectronicPotentialMatrix({ma});

    ASSERT_EQ(smata, smatb);
}
