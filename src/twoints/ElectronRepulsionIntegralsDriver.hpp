//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ElectronRepulsionIntegralsDriver_hpp
#define ElectronRepulsionIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "OutputStream.hpp"

/**
 Class CElectronicRepulsionIntegralsDriver computes electron repulsion
 <f(r)g(r')| 1/|r-r'||h(r') i(r')> integrals.
 
 @author Z. Rinkevicius
 */
class CElectronRepulsionIntegralsDriver
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
    
    
public:
    
    /**
     Creates an electron repulsion integrals driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CElectronRepulsionIntegralsDriver(const int32_t  globRank,
                                      const int32_t  globNodes,
                                            MPI_Comm comm);
    
    /**
     Destroys an electron repulsion integrals driver object.
     */
    ~CElectronRepulsionIntegralsDriver();
    
    /**
     Computes electron repulsion integrals for molecule with specific AO basis
     set and process results according to provided distribution function.
     
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @param threshold the integrals cut-off threshold.
     @param oStream the output stream.
     @param comm the MPI communicator.
     */
    void  compute(const CMolecule&       molecule,
                  const CMolecularBasis& aoBasis,
                  const double           threshold,
                        COutputStream&   oStream,
                        MPI_Comm         comm) const;
};


#endif /* ElectronRepulsionIntegralsDriver_hpp */
