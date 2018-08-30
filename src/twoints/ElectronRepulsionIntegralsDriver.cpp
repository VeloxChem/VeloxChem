//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ElectronRepulsionIntegralsDriver.hpp"

#include "SystemClock.hpp"
#include "GtoPairsContainer.hpp"

CElectronRepulsionIntegralsDriver::CElectronRepulsionIntegralsDriver(const int32_t  globRank,
                                                                     const int32_t  globNodes,
                                                                           MPI_Comm comm)

    : _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CElectronRepulsionIntegralsDriver::~CElectronRepulsionIntegralsDriver()
{
    
}

void
CElectronRepulsionIntegralsDriver::compute(const CMolecule&       molecule,
                                           const CMolecularBasis& aoBasis,
                                           const double           threshold,
                                                 COutputStream&   oStream,
                                                 MPI_Comm         comm) const
{
    CSystemClock eritim;
    
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-13);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split(10000);
    
    // FIX ME: ....
}
