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

#include "GtoTransformTest.hpp"

#include "GtoTransform.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CGtoTransformTest, To_VeloxChem)
{
    auto ambas = vlxbas::getMolecularBasisForHeAtom();
    
    auto he = vlxmol::getMoleculeHeAtom();
    
    CDenseMatrix dalmat({1.0, 2.0, 3.0, 4.0, 5.0,
                         1.1, 2.1, 3.1, 4.1, 5.1,
                         1.2, 2.2, 3.2, 4.2, 5.2,
                         1.3, 2.3, 3.3, 4.3, 5.3,
                         1.4, 2.4, 3.4, 4.4, 5.4},
                        5, 5);
    
    auto vlxmat = gtotra::to_veloxchem(dalmat, ambas, he);
    
    CDenseMatrix refmat({1.0, 2.0, 4.0, 5.0, 3.0,
                         1.1, 2.1, 4.1, 5.1, 3.1,
                         1.3, 2.3, 4.3, 5.3, 3.3,
                         1.4, 2.4, 4.4, 5.4, 3.4,
                         1.2, 2.2, 4.2, 5.2, 3.2},
                        5, 5);
    
    ASSERT_EQ(vlxmat, refmat);
}

TEST_F(CGtoTransformTest, To_Dalton)
{
    auto ambas = vlxbas::getMolecularBasisForHeAtom();
    
    auto he = vlxmol::getMoleculeHeAtom();
    
    CDenseMatrix vlxmat({1.0, 2.0, 4.0, 5.0, 3.0,
                         1.1, 2.1, 4.1, 5.1, 3.1,
                         1.3, 2.3, 4.3, 5.3, 3.3,
                         1.4, 2.4, 4.4, 5.4, 3.4,
                         1.2, 2.2, 4.2, 5.2, 3.2},
                        5, 5);
    
    auto dalmat = gtotra::to_dalton(vlxmat, ambas, he);
    
    CDenseMatrix refmat({1.0, 2.0, 3.0, 4.0, 5.0,
                         1.1, 2.1, 3.1, 4.1, 5.1,
                         1.2, 2.2, 3.2, 4.2, 5.2,
                         1.3, 2.3, 3.3, 4.3, 5.3,
                         1.4, 2.4, 3.4, 4.4, 5.4},
                        5, 5);
    
    ASSERT_EQ(dalmat, refmat); 
}
