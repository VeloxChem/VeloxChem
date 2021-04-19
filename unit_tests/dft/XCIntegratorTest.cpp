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

#include "XCIntegratorTest.hpp"

#include <mpi.h>

#include "GridDriver.hpp"
#include "XCIntegrator.hpp"
#include "MoleculeSetter.hpp"
#include "AOFockMatrix.hpp"
#include "MolecularBasisSetter.hpp"
#include "AODensityMatrixSetter.hpp"
#include "DenseLinearAlgebra.hpp"

TEST_F(CXCIntegratorTest, IntegrateKohnSham)
{
    CGridDriver gdrv(MPI_COMM_WORLD);
    
    gdrv.setLevel(6);
    
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mgrid = gdrv.generate(mh2o);
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    auto dmat = vlxden::getRestDensityMatrixForH2O();
    
    CXCIntegrator xcdrv(MPI_COMM_WORLD);
    
    auto ksmat = xcdrv.integrate(dmat, mh2o, mbas, mgrid, std::string("LYP"));
    
    ASSERT_NEAR(ksmat.getNumberOfElectrons(), 10.0, 1.0e-6);
}
