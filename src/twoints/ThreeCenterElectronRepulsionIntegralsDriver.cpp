//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ThreeCenterElectronRepulsionIntegralsDriver.hpp"

#ifdef MAC_OS_OMP
#include "/opt/intel/compilers_and_libraries/mac/include/omp.h"
#else
#include "omp.h"
#endif

#include "MpiFunc.hpp"
#include "SystemClock.hpp"
#include "GtoPairsContainer.hpp"
#include "GtoContainer.hpp"
#include "MathFunc.hpp"

CThreeCenterElectronRepulsionIntegralsDriver::CThreeCenterElectronRepulsionIntegralsDriver(const int32_t  globRank,
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

CThreeCenterElectronRepulsionIntegralsDriver::~CThreeCenterElectronRepulsionIntegralsDriver()
{
    
}

void
CThreeCenterElectronRepulsionIntegralsDriver::compute(const CMolecule&       molecule,
                                                      const CMolecularBasis& aoBasis,
                                                      const CMolecularBasis& riBasis,
                                                      const double           threshold, 
                                                            COutputStream&   oStream,
                                                            MPI_Comm         comm) const
{
    CSystemClock eritim;
    
    // generate GTOs pairs blocks for AO basis
    
    CGtoPairsContainer kgtopairs(molecule, aoBasis, 1.0e-13);
    
    // set up GTOs splitting pattern for RI basis
    
    auto gtopat = _getBatchesOfGtoBlocks(molecule, riBasis, kgtopairs);
    
    // generate RI gtos for on each MPI node
    
    // CGtoContainer bgtos(molecule, riBasis, nodoff, nodatm);
    
    // print start header
    
    if (_globRank == mpi::master()) _startHeader(oStream);
}

void
CThreeCenterElectronRepulsionIntegralsDriver::_startHeader(COutputStream& oStream) const
{
    oStream << fmt::header << "Three-Center Electron Repulsion Integrals" << fmt::end;
    
    oStream << std::string(43, '=') << fmt::end << fmt::blank;
}


CMemBlock2D<int32_t>
CThreeCenterElectronRepulsionIntegralsDriver::_getBatchesOfGtoBlocks(const CMolecule&          molecule,
                                                                     const CMolecularBasis&    riBasis,
                                                                     const CGtoPairsContainer& gtoPairs) const
{
    // determine dimensions of atoms batch for each MPI node
    
    auto natoms = molecule.getNumberOfAtoms();
    
    auto nodatm = mpi::batch_size(natoms, _locRank, _locNodes);
    
    auto nodoff = mpi::batch_offset(natoms, _locRank, _locNodes);
    
    // determine max. task vector dimensions for single atoms batch
    
    auto mtasks = (riBasis.getMaxAngularMomentum() + 1)
    
                * gtoPairs.getNumberOfGtoPairsBlocks();
    
    // determine max. number of threads
    
    int32_t mthreads = 10 * omp_get_max_threads();
    
    // determine number of atomic batches
    
    auto nblocks = mthreads / mtasks + 1;
    
    CMemBlock2D<int32_t> batches(nblocks, 2);
    
    auto bpos = batches.data(0);
   
    auto bdim = batches.data(1);
    
    for (int32_t i = 0; i < nblocks; i++)
    {
        bdim[i] = mpi::batch_size(nodatm, i, nblocks);
    }
    
    mathfunc::indexes(bpos, bdim, nodoff, nblocks);
    
    return batches;
}
