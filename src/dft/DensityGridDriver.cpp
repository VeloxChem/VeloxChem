//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "DensityGridDriver.hpp"

#include "SystemClock.hpp"
#include "GtoContainer.hpp"
#include "XCFuncType.hpp"
#include "OMPTasks.hpp"

CDensityGridDriver::CDensityGridDriver(const int32_t  globRank,
                                       const int32_t  globNodes,
                                             MPI_Comm comm)
    :_globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)

    , _thresholdOfDensity(1.0e-12)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CDensityGridDriver::~CDensityGridDriver()
{
    
}

void
CDensityGridDriver::generate(const CMolecule&       molecule,
                             const CMolecularBasis& basis,
                             const CMolecularGrid&  molGrid,
                             const xcfun            xcFunctional, 
                                   COutputStream&   oStream,
                                   MPI_Comm         comm)
{
    CSystemClock tim;
    
    COMPTasks omptaks(5); 
    
    printf("Density Grid: rank %i points: %i time %lf s.\n",
           _locRank, molGrid.getNumberOfGridPoints(), tim.getElapsedTimeInSeconds());
    
    return;
}
