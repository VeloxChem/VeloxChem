//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "XCIntegratorTest.hpp"

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
