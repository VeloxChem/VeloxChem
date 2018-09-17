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
#include "CheckFunctions.hpp"

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSSForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock bgtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(bgtos, 1.0e-13);
    
    CMemBlock<double> fints(1);
    
    eridrv.compute(fints.data(), bpairs.pick(0), bpairs.pick(0));
    
    std::vector<double> r0000vals{ 1.629232440467330};
    
    vlxtest::compare(r0000vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(1));
    
    std::vector<double> r0101vals{ 1.296490921648330};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), bpairs.pick(8));
    
    std::vector<double> r0123vals{ 0.569573705357336};
    
    vlxtest::compare(r0123vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(3), bpairs.pick(5));
    
    std::vector<double> r0312vals{ 0.091957062417462};
    
    vlxtest::compare(r0312vals, fints.data());
}

TEST_F(CElectronRepulsionIntegralsDriverTest, ComputeSSSPForLiH)
{
    CElectronRepulsionIntegralsDriver eridrv(mpi::master(), mpi::nodes(MPI_COMM_WORLD),
                                             MPI_COMM_WORLD);
    
    auto mlih = vlxmol::getTestLiH();
    
    auto mbas = vlxbas::getTestBasisForLiH();
    
    CGtoBlock agtos(mlih, mbas, 0);
    
    CGtoPairsBlock bpairs(agtos, 1.0e-13);
    
    CGtoBlock cgtos(mlih, mbas, 0);
    
    CGtoBlock dgtos(mlih, mbas, 1);
    
    CGtoPairsBlock kpairs(cgtos, dgtos, 1.0e-13);
    
    CMemBlock<double> fints(3);
    
    eridrv.compute(fints.data(), bpairs.pick(0), kpairs.pick(0));
    
    std::vector<double> r0000vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0000vals,fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(1));
    
    std::vector<double> r0101vals{ 0.000000000000000,  0.000000000000000,  0.000000000000000};
    
    vlxtest::compare(r0101vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(1), kpairs.pick(11));
    
    std::vector<double> r0123vals{-0.048499014664650, -0.088914860218526, -0.032332676443100};
    
    vlxtest::compare(r0123vals, fints.data());
    
    eridrv.compute(fints.data(), bpairs.pick(3), kpairs.pick(6));
    
    std::vector<double> r0312vals{-0.065342004894742, -0.119793675640361, -0.043561336596495};
    
    vlxtest::compare(r0312vals, fints.data());
}
