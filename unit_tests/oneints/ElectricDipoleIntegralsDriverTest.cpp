//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.


#include "ElectricDipoleIntegralsDriverTest.hpp"

#include "ElectricDipoleIntegralsDriver.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"

TEST_F(CElectricDipoleIntegralsDriverTest, ComputeSSForLiH)
{
    CElectricDipoleIntegralsDriver dipdrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                          MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoBlock kgtos(mlih, mbas, 0);
    
    CDenseMatrix dipxmat(4, 4);
    
    CDenseMatrix dipymat(4, 4);
    
    CDenseMatrix dipzmat(4, 4);
    
    dipdrv.compute(dipxmat.values(), dipymat.values(), dipzmat.values(),
                   bgtos, kgtos);
    
    std::vector<double> intxvals{ 0.000000000000000,  0.000000000000000,  0.034080416572876,
                                  0.038381743586952,  0.000000000000000,  0.000000000000000,
                                  0.060847045273318,  0.067637145335229,  0.034080416572876,
                                  0.060847045273318,  0.400000000000000,  0.321445559231693,
                                  0.038381743586952,  0.067637145335229,  0.321445559231693,
                                  0.400000000000000};
    
    ASSERT_EQ(dipxmat, CDenseMatrix(intxvals, 4, 4));
    
    std::vector<double> intyvals{ 0.000000000000000,  0.000000000000000,  0.051120624859315,
                                  0.057572615380427,  0.000000000000000,  0.000000000000000,
                                  0.091270567909978,  0.101455718002844,  0.051120624859315,
                                  0.091270567909978,  0.600000000000000,  0.482168338847539,
                                  0.057572615380427,  0.101455718002844,  0.482168338847539,
                                  0.600000000000000};
    
    ASSERT_EQ(dipymat, CDenseMatrix(intyvals, 4, 4));
    
    std::vector<double> intzvals{ 0.000000000000000,  0.000000000000000,  0.093721145575410,
                                  0.105549794864117,  0.000000000000000,  0.000000000000000,
                                  0.167329374501625,  0.186002149671880,  0.093721145575410,
                                  0.167329374501625,  1.100000000000001,  0.883975287887155,
                                  0.105549794864117,  0.186002149671880,  0.883975287887155,
                                  1.100000000000001};
    
    ASSERT_EQ(dipzmat, CDenseMatrix(intzvals, 4, 4));
}
