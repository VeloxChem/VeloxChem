//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronRepulsionIntegralsDriverTest.hpp"

#include "mpi.h"

#include "ElectronRepulsionIntegralsDriver.hpp"
#include "MolecularBasisSetter.hpp"
#include "MoleculeSetter.hpp"

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSSForLiH2)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih2 = vlxmol::getTestLiH2();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih2, mbas, 0);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> fints(231);
    
    eridrv.compute(fints.data(), bpairs, bpairs);
    
//    auto idxa = bpairs.getBraIdentifiers(0);
//    
//    auto idxb = bpairs.getKetIdentifiers(0);
//    
//    int32_t idx = 0;
//    
//    for (int32_t i = 0; i < bpairs.getNumberOfScreenedContrPairs(); i++)
//    {
//        for (int32_t j = 0; j < i + 1; j++)
//        {
//            printf("(%i,%i|%i,%i) = %lf\n", idxa[i], idxb[i], idxa[j], idxb[j], fints.at(idx));
//            
//            idx++;
//        }
//    }
//    
//    printf("Number of pairs: %i Number of integrals: %i\n",
//           bpairs.getNumberOfScreenedContrPairs(), idx); 
}
