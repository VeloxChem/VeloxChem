//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef MolecularOrbitals_hpp
#define MolecularOrbitals_hpp

#include <cstdint>
#include <vector>

#include "MolecularOrbitalsType.hpp"
#include "DenseMatrix.hpp"
#include "MemBlock.hpp"
#include "AODensityMatrix.hpp"

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
     Gets constant pointer to first element of specific unrestricted alpha
     molecular orbital matrix.
     
     @return the constant pointer to first element of unrestricted alpha
     molecular orbital matrix.
     */
    const double* alphaOrbitals() const;
    
    /**
     Gets constant pointer to first element of specific unrestricted beta
     molecular orbital matrix.
     
     @return the constant pointer to first element of unrestricted beta
     molecular orbital matrix.
     */
    const double* betaOrbitals() const;
    
    /**
     Gets constant pointer to first element of specific unrestricted alpha
     energy eigenvalues.
     
     @return the constant pointer to first element of unrestricted alpha
     energy eigenvalues.
     */
    const double* alphaEnergies() const;
    
    /**
     Gets constant pointer to first element of specific unrestricted beta
     energy eigenvalues.
     
     @return the constant pointer to first element of unrestricted beta
     energy eigenvalues.
     */
    const double* betaEnergies() const;
    
    /**
     Gets string representation of density matrix object.

     @return the string representation.
     */
    std::string getString() const;

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
     Converts molecular orbitals object to text output and insert it into output
     text stream.
     
     @param output the output text stream.
     @param source the molecular orbitals object.
     */
    friend std::ostream& operator<<(      std::ostream&       output,
                                    const CMolecularOrbitals& source);
};

#endif /* MolecularOrbitals_hpp */
