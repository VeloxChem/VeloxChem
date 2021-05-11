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

#include "DensityGridDriverTest.hpp"

#include <mpi.h>

#include "AODensityMatrixSetter.hpp"
#include "DensityGrid.hpp"
#include "DensityGridDriver.hpp"
#include "GridDriver.hpp"
#include "MolecularBasisSetter.hpp"
#include "MolecularGrid.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CDensityGridDriverTest, Generate)
{
    CGridDriver gdrv(MPI_COMM_WORLD);
    
    gdrv.setLevel(6);
    
    auto mh2o = vlxmol::getMoleculeH2O();
    
    auto mgrid = gdrv.generate(mh2o);
    
    auto mbas = vlxbas::getMolecularBasisForH2O();
    
    auto dmat = vlxden::getRestDensityMatrixForH2O();
    
    CDensityGridDriver ddrv(MPI_COMM_WORLD);
    
    auto dgrid = ddrv.generate(dmat, mh2o, mbas, mgrid, xcfun::lda);
    
    // compute number of electrons
    
    auto gw = mgrid.getWeights();
    
    auto rhoa = dgrid.alphaDensity(0);
    
    auto rhob = dgrid.betaDensity(0);
    
    double nelec = 0.0;
    
    for (int32_t i = 0; i < dgrid.getNumberOfGridPoints(); i++)
    {
        nelec += (rhoa[i] + rhob[i]) * gw[i]; 
    }
    
    ASSERT_NEAR(nelec, 10.0, 1.0e-6); 
}
