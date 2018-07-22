//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef OverlapIntegralsDriver_hpp
#define OverlapIntegralsDriver_hpp

#include "mpi.h"

#include "OutputStream.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "OutputStream.hpp"
#include "ExecMode.hpp"
#include "OverlapMatrix.hpp"

/**
 Class COverlapIntegralsDriver computes one-electron integrals.
 
 @author Z. Rinkevicius
 */
class COverlapIntegralsDriver
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
     The execution mode of grid driver object.
     */
    execmode _runMode;
    
public:
    
    /**
     Creates a overlap integrals driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param runMode the execution mode.
     @param comm the MPI communicator.
     */
    COverlapIntegralsDriver(const int32_t  globRank,
                            const int32_t  globNodes,
                            const execmode runMode,
                                  MPI_Comm comm);
    
    /**
     Destroys a overlap integrals driver object.
     */
    ~COverlapIntegralsDriver();
    
    /**
     Computes overlap integrals for molecule in specific basis set and stores
     results in overlap matrix object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param comm the MPI communicator.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       molecule,
                           const CMolecularBasis& basis,
                                 MPI_Comm         comm) const;
    
    /**
     Computes overlap integrals for molecule in two basis sets and stores
     results in overlap matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @param comm the MPI communicator.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       molecule,
                           const CMolecularBasis& braBasis,
                           const CMolecularBasis& ketBasis,
                                 MPI_Comm         comm) const;
    
    /**
     Computes overlap integrals for two molecules in basis set and stores
     results in overlap matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param basis the molecular basis.
     @param comm the MPI communicator.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       braMolecule,
                           const CMolecule&       ketMolecule,
                           const CMolecularBasis& basis,
                                 MPI_Comm         comm) const;
    
    /**
     Computes overlap integrals for two molecules in different basis sets and
     stores results in overlap matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @param comm the MPI communicator.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       braMolecule,
                           const CMolecule&       ketMolecule,
                           const CMolecularBasis& braBasis,
                           const CMolecularBasis& ketBasis,
                                 MPI_Comm         comm) const;
};

#endif /* OverlapIntegralsDriver_hpp */
