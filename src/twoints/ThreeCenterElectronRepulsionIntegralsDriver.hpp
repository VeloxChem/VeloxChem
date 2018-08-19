//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ThreeCenterElectronRepulsionIntegralsDriver_hpp
#define ThreeCenterElectronRepulsionIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "OutputStream.hpp"
#include "MemBlock2D.hpp"
#include "GtoPairsContainer.hpp"
#include "GtoContainer.hpp"

/**
 Class CThreeCenterElectronicRepulsionIntegralsDriver computes electronic potential
 <f(r)| 1/|r-r'||g(r') h(r')> integrals.
 
 @author Z. Rinkevicius
 */
class CThreeCenterElectronRepulsionIntegralsDriver
{
    /**
     The rank of associated global MPI process.
     */
    int32_t _globRank;
    
    /**
     The total number of global MPI processes.
     */
    int32_t _globNodes;
    
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;
    
    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;
    
    /**
     The flag for local execution mode.
     */
    bool _isLocalMode;
    
    /**
     Prints start header for computation of three-center electron repulsion
     integrals.
     
     @param gtoPairs the GTOs pairs container on ket side.
     @param oStream the output stream.
     */
    void _startHeader(const CGtoPairsContainer& gtoPairs,
                            COutputStream&      oStream) const;
    
    /**
     Creates atoms list splitting pattern for generation of GTO blocks on each
     MPI process. The splitting pattern enables effective task based computation
     of integrals on systems with variable number of CPU cores.

     @param molecule the molecule.
     @param riBasis the molecular RI basis.
     @param gtoPairs the GTOs pairs container.
     @return the splitting pattern for atoms list.
     */
    CMemBlock2D<int32_t> _getBatchesOfGtoBlocks(const CMolecule&          molecule,
                                                const CMolecularBasis&    riBasis,
                                                const CGtoPairsContainer& gtoPairs) const;
    
public:
    
    /**
     Creates a three center electron repulsion integrals driver object using
     MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CThreeCenterElectronRepulsionIntegralsDriver(const int32_t  globRank,
                                                 const int32_t  globNodes,
                                                       MPI_Comm comm);
    
    /**
     Destroys a three center electron repulsion integrals driver object.
     */
    ~CThreeCenterElectronRepulsionIntegralsDriver();
    
    /**
     Computes three center electron repulsion integrals for molecule with
     specific combination of AO and RI basis sets and stores results in
     compressed format on each MPI process within domain of MPI communicator.
     
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @param riBasis the molecular RI basis.
     @param threshold the integrals cut-off threshold. 
     @param oStream the output stream.
     @param comm the MPI communicator.
     */
    void  compute(const CMolecule&       molecule,
                  const CMolecularBasis& aoBasis,
                  const CMolecularBasis& riBasis,
                  const double           threshold,
                        COutputStream&   oStream,
                        MPI_Comm         comm) const;
};

#endif /* ThreeCenterElectronRepulsionIntegralsDriver_hpp */
