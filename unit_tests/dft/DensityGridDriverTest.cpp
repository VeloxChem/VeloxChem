//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "DensityGridDriverTest.hpp"

#include "GridDriver.hpp"
#include "DensityGridDriver.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"
#include "AODensityMatrixSetter.hpp"

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
