//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

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
