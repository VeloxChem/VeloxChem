//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#include "ThreeCenterElectronRepulsionIntegralsDriver.hpp"

#include "MpiFunc.hpp"
#include "SystemClock.hpp"
#include "GtoPairsContainer.hpp"
#include "GtoContainer.hpp"

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
    
    // determine dimensions of atoms batch for each MPI node
    
    auto natoms = molecule.getNumberOfAtoms();
    
    auto nodatm = mpi::batch_size(natoms, _locRank, _locNodes);
    
    auto nodoff = mpi::batch_offset(natoms, _locRank, _locNodes);
    
    printf("Local rank: %i Atoms: %i Position: %i\n", _locRank, nodatm, nodoff);
    
    // generate AO pairs
    
    CGtoPairsContainer kgtopairs(molecule, aoBasis, 1.0e-13);
    
    // generate RI gtos for on each MPI node
    
    CGtoContainer bgtos(molecule, riBasis, nodoff, nodatm);
    
    // print start header
    
    if (_globRank == mpi::master()) _startHeader(oStream);
    
    

    

}

void
CThreeCenterElectronRepulsionIntegralsDriver::_startHeader(COutputStream& oStream) const
{
    oStream << fmt::header << "Three-Center Electron Repulsion Integrals" << fmt::end;
    
    oStream << std::string(43, '=') << fmt::end << fmt::blank;
}
