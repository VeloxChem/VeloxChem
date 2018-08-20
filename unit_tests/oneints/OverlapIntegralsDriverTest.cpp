//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "OverlapIntegralsDriverTest.hpp"

#include "OverlapIntegralsDriver.hpp"
#include "MoleculeSetter.hpp"
#include "MolecularBasisSetter.hpp"

TEST_F(COverlapIntegralsDriverTest, ComputeSSForLiH)
{
    COverlapIntegralsDriver ovldrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                   MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH(0, 0);
    
//    auto ovlmat = ovldrv.compute(mlih, mbas, MPI_COMM_WORLD);
//
//    std::vector<double> intvals{ 1.000000000000000,  0.950707129188576,  0.153183585372556,
//                                 0.323354079718979,  0.950707129188576,  1.000000000000001,
//                                 0.229357720620491,  0.422732158345182,  0.153183585372556,
//                                 0.229357720620491,  1.000000000000000,  0.803613898079232,
//                                 0.323354079718979,  0.422732158345182,  0.803613898079232,
//                                 1.000000000000000};
//
//    std::vector<int32_t> rows{0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3};
//
//    std::vector<int32_t> cols{0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3};
//
//    ASSERT_EQ(ovlmat, CSparseMatrix(intvals, rows, cols, 4, 4, 1.0e-13));
}


