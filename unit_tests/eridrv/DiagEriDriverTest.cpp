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

#include "DiagEriDriverTest.hpp"

#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "DiagEriDriver.hpp"

TEST_F(CDiagEriDriverTest, Compute)
{
    const auto mlih = vlxmol::getTestLiH();
    
    const auto mbas = vlxbas::getMolecularBasisForLiHX();
    
    CDiagEriDriver<double> eridrv;
    
    eridrv.compute(mlih, mbas);
}
