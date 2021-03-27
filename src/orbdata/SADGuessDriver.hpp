//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#ifndef SADGuessDriver_hpp
#define SADGuessDriver_hpp

#include <cstdint>

#include <mpi.h>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"
#include "OverlapMatrix.hpp"

/**
 Class CSADGuessDriver computes SAD initial guess.

 @author X. Li
 */
class CSADGuessDriver
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The local MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     Gets occupation numbers for 1s elements.

     @param nocc number of 1s orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc1s(const double nocc) const;

    /**
     Gets occupation numbers for 2s elements.

     @param nocc number of 2s orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc2s(const double nocc) const;

    /**
     Gets occupation numbers for 2p elements.

     @param nocc number of 2s2p orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc2s2p(const double nocc) const;

    /**
     Gets occupation numbers for 3s elements.

     @param nocc number of 3s orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc3s(const double nocc) const;

    /**
     Gets occupation numbers for 3p elements.

     @param nocc number of 3s3p orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc3s3p(const double nocc) const;

    /**
     Gets occupation numbers for 4s elements.

     @param nocc number of 4s orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc4s(const double nocc) const;

    /**
     Gets occupation numbers for 3d elements.

     @param nocc number of 3d orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc3d(const double nocc) const;

    /**
     Gets occupation numbers for 4p elements.

     @param nocc number of 4s4p orbitals.
     @return vector of occupation numbers.
     */
    std::vector<double> _getOcc4s4p(const double nocc) const;

    /**
     Computes SAD initial guess.

     @param molecule the molecule.
     @param basis_1 the minimal (smaller) basis set.
     @param basis_2 the molecular (larger) basis set.
     @param S12 the crossing overlap matrix between basis_1 and basis_2.
     @param S22 the overlap matrix computed from basis_2.
     @param restricted the flag for generating restricted initial guess
     @return the density matrix of SAD guess.
     */
    CAODensityMatrix _compSADGuess(const CMolecule&       molecule,
                                   const CMolecularBasis& basis_1,
                                   const CMolecularBasis& basis_2,
                                   const COverlapMatrix&  S12,
                                   const COverlapMatrix&  S22,
                                   const bool             restricted) const;

   public:
    /**
     Creates a SAD guess driver object using MPI info.

     @param comm the MPI communicator.
     */
    CSADGuessDriver(MPI_Comm comm);

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
     @param restricted the flag for generating restricted initial guess
     @return the density matrix of SAD guess.
     */
    CAODensityMatrix compute(const CMolecule&       molecule,
                             const CMolecularBasis& basis_1,
                             const CMolecularBasis& basis_2,
                             const COverlapMatrix&  S12,
                             const COverlapMatrix&  S22,
                             const bool             restricted) const;

    /**
     Computes indicies of atomic orbitals that are located on each atom.

     @param molecule the molecule.
     @param basis the molecular basis set.
     @return a vector of vector containing the atomic orbital indicies for each atom in the molecule.
     */
    std::vector<std::vector<int32_t>> getAOIndicesOfAtoms(const CMolecule& molecule, const CMolecularBasis& basis) const;

    /**
     Computes occupation numbers for a given element.

     @param elem_id the element id (nuclear charge).
     @param nelec the number of excessive alpha or beta electrons.
     @return a vector of occupation numbers.
     */
    std::vector<double> getOccupationNumbersForElement(const int32_t elem_id, const double nelec) const;

    /**
     Computes occupation numbers for a given molecule.

     @param molecule the molecule.
     @param nelec the number of excessive alpha or beta electrons.
     @return a vector of vector containing the occupation numbers for each atom in the molecule.
     */
    std::vector<std::vector<double>> getOccupationNumbersForMolecule(const CMolecule& molecule, const double nelec) const;
};

#endif /* SADGuessDriver_hpp */
