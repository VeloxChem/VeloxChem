//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularOrbitals_hpp
#define MolecularOrbitals_hpp

#include <cstdint>
#include <vector>

#include "MolecularOrbitalsType.hpp"
#include "DenseMatrix.hpp"
#include "MemBlock.hpp"
#include "AODensityMatrix.hpp"
#include "MolecularBasis.hpp"
#include "SpinBlock.hpp"

/**
 Class CMolecularOrbitals stores data about molecular orbitals and provides set
 of methods for handling of molecular orbitals data.
 
 @author Z. Rinkevicius
 */
class CMolecularOrbitals
{
    /**
     The type of molecular orbitals.
     */
    molorb _orbitalsType;
    
    /**
     The vector of dense matrices for storing molecular orbitals.
     */
    std::vector<CDenseMatrix> _orbitals;
    
    /**
     The vector of vectors for storing molecular orbital energies.
     */
    std::vector<CMemBlock<double>> _energies;
    
public:
    
    /**
     Creates an empty molecular orbitals object.
     */
    CMolecularOrbitals();
    
    /**
     Creates a molecular orbitals object.
     
     @param orbitals the vector of dense matrices with molecular orbitals.
     @param energies the vector of orbital energies vectors.
     @param orbitalsType the type of molecular orbitals.
     */
    CMolecularOrbitals(const std::vector<CDenseMatrix>&      orbitals,
                       const std::vector<CMemBlock<double>>& energies,
                       const molorb                          orbitalsType);
    
    /**
     Creates a molecular orbitals object by copying other molecular orbitals object.
     
     @param source the molecular orbitals object.
     */
    CMolecularOrbitals(const CMolecularOrbitals& source);
    
    /**
     Creates a molecular orbitals object by moving other molecular orbitals object.
     
     @param source the molecular basis object.
     */
    CMolecularOrbitals(CMolecularOrbitals&& source) noexcept;
    
    /**
     Destroys a molecular orbitals object.
     */
    ~CMolecularOrbitals();
    
    /**
     Assigns a molecular orbitals object by copying other molecular orbitals object.
     
     @param source the molecular orbitals object.
     */
    CMolecularOrbitals& operator=(const CMolecularOrbitals& source);
    
    /**
     Assigns a molecular orbitals object by moving other molecular orbitals object.
     
     @param source the molecular orbitals object.
     */
    CMolecularOrbitals& operator=(CMolecularOrbitals&& source) noexcept;
    
    /**
     Compares molecular orbitals object with other molecular orbitals object.
     
     @param other the molecular orbitals object.
     @return true if molecular orbitals objects are equal, false otherwise.
     */
    bool operator==(const CMolecularOrbitals& other) const;
    
    /**
     Compares molecular orbitals object with other molecular orbitals object.
     
     @param other the molecular orbitals object.
     @return true if molecular orbitals objects are not equal, false otherwise.
     */
    bool operator!=(const CMolecularOrbitals& other) const;
    
    /**
     Creates molecular orbitals objects from this molecular orbitals object
     according to given basis sets pair.

     @param molecule the molecule.
     @param aoBasis the molecular basis for created orbitals.
     @param minBasis the molecular basis for this orbitals.
     @return the molecular orbitals.
     */
    CMolecularOrbitals insert(const CMolecule&       molecule,
                              const CMolecularBasis& aoBasis,
                              const CMolecularBasis& minBasis) const;
    
    /**
     Gets type of molecular orbital matrix.

     @return the type of molecular orbital matrix.
     */
    molorb getOrbitalsType() const;
    
    /**
     Gets number of rows in specific molecular orbital matrix.
     
     @return the number of rows.
     */
    int32_t getNumberOfRows() const;
    
    /**
     Gets number of columns in specific molecular orbital matrix.
     
     @return the number of columns.
     */
    int32_t getNumberOfColumns() const;
    
    /**
     Gets constant pointer to first element of specific alpha molecular orbital
     matrix.
     
     @return the constant pointer to first element of alpha molecular orbital
     matrix.
     */
    const double* alphaOrbitals() const;
    
    /**
     Gets constant pointer to first element of specific beta molecular orbital
     matrix.
     
     @return the constant pointer to first element of beta molecular orbital
     matrix.
     */
    const double* betaOrbitals() const;
    
    /**
     Gets constant pointer to first element of specific alpha energy eigenvalues.
     
     @return the constant pointer to first element of alpha energy eigenvalues.
     */
    const double* alphaEnergies() const;
    
    /**
     Gets constant pointer to first element of specific beta energy eigenvalues.
     
     @return the constant pointer to first element of beta energy eigenvalues.
     */
    const double* betaEnergies() const;
    
    /**
     Gets alpha orbitals within specific range.
     
     @param iMolecularOrbital the index of starting orbital.
     @param nMolecularOrbitals the number of molecular orbitals.
     @return the alpha orbitals.
     */
    CDenseMatrix alphaOrbitals(const int32_t iMolecularOrbital,
                               const int32_t nMolecularOrbitals) const;
    
    /**
     Gets beta orbitals within specific range.
     
     @param iMolecularOrbital the index of starting orbital.
     @param nMolecularOrbitals the number of molecular orbitals.
     @return the beta orbitals.
     */
    CDenseMatrix betaOrbitals(const int32_t iMolecularOrbital,
                              const int32_t nMolecularOrbitals) const;
    
    /**
     Gets alpha orbitals within given vector of molecular orbital indexes.

     @param iMolecularOrbitals the vector of molecular orbital indexes.
     @return the alpha orbitals.
     */
    CDenseMatrix alphaOrbitals(const std::vector<int32_t>& iMolecularOrbitals) const;
    
    /**
     Gets beta orbitals within given vector of molecular orbital indexes.
     
     @param iMolecularOrbitals the vector of molecular orbital indexes.
     @return the beta orbitals.
     */
    CDenseMatrix betaOrbitals(const std::vector<int32_t>& iMolecularOrbitals) const;
    
    /**
     Gets string representation of density matrix object.

     @return the string representation.
     */
    std::string getString() const;
    
    /**
     Computes spin restricted electron density matrix in AO basis for specific
     number of electrons.

     @param nElectrons the total number of electrons.
     @return the AO density matrix.
     */
    CAODensityMatrix getAODensity(const int32_t nElectrons) const;
    
    /**
     Computes spin unrestricted electron density matrix in AO basis for specific
     number of alpha and beta electrons.
     
     @param nAlphaElectrons the number of alpha electrons.
     @param nBetaElectrons the number of beta electrons.
     @return the AO density matrix.
     */
    CAODensityMatrix getAODensity(const int32_t nAlphaElectrons,
                                  const int32_t nBetaElectrons) const;
    
    /**
     Computes restricted pair C_i C_j^T density matrix in AO basis.
     
     @param iMolecularOrbital the index of i-th molecular orbital C_i.
     @param jMolecularOrbital the index of j-th molecular orbital C_j.
     @return the AO pair density matrix.
     */
    CAODensityMatrix getRestrictedPairDensity(const int32_t iMolecularOrbital,
                                              const int32_t jMolecularOrbital) const;
    
    /**
     Computes set of restricted pair C_i C_j^T density matrices in AO basis.
     
     @param iMolecularOrbitals the vector of index of i-th molecular orbital C_i.
     @param jMolecularOrbitals the vector of index of j-th molecular orbital C_j.
     @return the AO pair density matrix.
     */
    CAODensityMatrix getRestrictedPairDensity(const std::vector<int32_t>& iMolecularOrbitals,
                                              const std::vector<int32_t>& jMolecularOrbitals) const;
    
    
    /**
     Transforms matrix in AO basis to matrix in MO basis using selected
     molecular orbitals.

     @param aoMatrix the matrix in AO basis.
     @param spinPair the spin block of molecular orbitals.
     @return the matrix in MO basis.
     */
    CDenseMatrix transform(const CDenseMatrix& aoMatrix,
                           const szblock       spinPair) const;

    /**
     Broadcasts molecular orbitals object within domain of MPI communicator.
     
     @param rank the rank of MPI process.
     @param comm the MPI communicator.
     */
    void broadcast(int32_t rank, MPI_Comm comm);
    
    /**
     Converts molecular orbitals object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the molecular orbitals object.
     */
    friend std::ostream& operator<<(      std::ostream&       output,
                                    const CMolecularOrbitals& source);
};

#endif /* MolecularOrbitals_hpp */
