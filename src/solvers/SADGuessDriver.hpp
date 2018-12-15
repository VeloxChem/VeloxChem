//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef SADGuessDriver_hpp
#define SADGuessDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "DenseMatrix.hpp"
#include "OverlapMatrix.hpp"
#include "AODensityMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "SystemClock.hpp"

/**
 Class CSADGuessDriver computes SAD initial guess.
 
 @author X. Li
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
     Gets occupation numbers for 1s elements.

     @param occupation number of 1s orbital.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc1s(double occ) const;

    /**
     Gets occupation numbers for 2s elements.

     @param occupation number of 2s orbital.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc2s(double occ) const;

    /**
     Gets occupation numbers for 2p elements.

     @param occupation number of 2s2p orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc2s2p(double occ) const;

    /**
     Gets occupation numbers for 3s elements.

     @param occupation number of 3s orbital.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc3s(double occ) const;

    /**
     Gets occupation numbers for 3p elements.

     @param occupation number of 3s3p orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc3s3p(double occ) const;

    /**
     Gets occupation numbers for 4s elements.

     @param occupation number of 4s orbital.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc4s(double occ) const;

    /**
     Gets occupation numbers for 3d elements.

     @param occupation number of 3d orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc3d(double occ) const;

    /**
     Gets occupation numbers for 4p elements.

     @param occupation number of 4s4p orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc4s4p(double occ) const;

    /**
     Sets occupation numbers.

     @return a vector of vector with the first dimension for nuclear charge and
             the second dimension for occupation numbers.
     */
    std::vector<std::vector<double>> _buildQocc() const;
    
    /**
     Computes SAD initial guess.

     @param molecule the molecule.
     @param basis_1 the minimal (smaller) basis set.
     @param basis_2 the molecular (larger) basis set.
     @param S12 the crossing overlap matrix between basis_1 and basis_2.
     @param S22 the overlap matrix computed from basis_2.
     @return the density matrix of SAD guess.
     */
    CAODensityMatrix
    _compSADGuess(const CMolecule&       molecule,
                  const CMolecularBasis& basis_1,
                  const CMolecularBasis& basis_2,
                  const COverlapMatrix&  S12,
                  const COverlapMatrix&  S22) const;
    
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
     @param comm the MPI communicator.
     @return the density matrix of SAD guess.
     */
    CAODensityMatrix
    compute(const CMolecule&       molecule,
            const CMolecularBasis& basis_1,
            const CMolecularBasis& basis_2,
            const COverlapMatrix&  S12,
            const COverlapMatrix&  S22,
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
