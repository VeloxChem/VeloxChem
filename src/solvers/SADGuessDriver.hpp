//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef SADGuessDriver_hpp
#define SADGuessDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "OutputStream.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"
#include "SystemClock.hpp"

/**
 Class CSADGuessDriver computes SAD initial guess.
 
 @author 
 */
class CSADGuessDriver
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
     Sets occupation numbers.

     @return a vector of vector with the first dimension for nuclear charge and
             the second dimension for occupation numbers.
     */
    std::vector< std::vector<double> >
    _buildQocc() const;
    
    /**
     Computes SAD initial guess.

     @param molecule the molecule.
     @param basis_1 the minimal (smaller) basis set.
     @param basis_2 the molecular (larger) basis set.
     @param S12 the crossing overlap matrix between basis_1 and basis_2.
     @param S22 the overlap matrix computed from basis_2.
     @return the density matrix of SAD guess.
     */
    CDenseMatrix
    _compSADGuess(const CMolecule&       molecule,
                  const CMolecularBasis& basis_1,
                  const CMolecularBasis& basis_2,
                  const COverlapMatrix&  S12,
                  const COverlapMatrix&  S22) const;
    
    /**
     Prints overlap integrals computation time to output stream.

     @param timer the system clock timer.
     @param oStream the output stream.
     */
    void _printComputationTime(const CSystemClock&  timer,
                                     COutputStream& oStream) const;

public:
    
    /**
     Creates a SAD guess driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CSADGuessDriver(const int32_t  globRank,
                    const int32_t  globNodes,
                          MPI_Comm comm);
    
    /**
     Destroys a SAD guess driver object.
     */
    ~CSADGuessDriver();
    
    /**
     Computes SAD initial guess.

     @param molecule the molecule.
     @param basis_1 the minimal (smaller) basis set.
     @param basis_2 the molecular (larger) basis set.
     @param S12 the crossing overlap matrix between basis_1 and basis_2.
     @param S22 the overlap matrix computed from basis_2.
     @param oStream the output stream.
     @param comm the MPI communicator.
     @return the density matrix of SAD guess.
     */
    CDenseMatrix
    compute(const CMolecule&       molecule,
            const CMolecularBasis& basis_1,
            const CMolecularBasis& basis_2,
            const COverlapMatrix&  S12,
            const COverlapMatrix&  S22,
                  COutputStream&   oStream,
                  MPI_Comm         comm) const;
    
    /**
     Computes indicies of atomic orbitals that are located on each atom.
     
     @param molecule the molecule.
     @param basis the molecular basis set.
     @return a vector of vector with the first dimension for atoms and the
             second dimension for atomic orbitals.
     */
    std::vector< std::vector<int32_t> >
    getAOIndicesOfAtoms(const CMolecule&       molecule,
                        const CMolecularBasis& basis) const;
};

#endif /* SADGuessDriver_hpp */
