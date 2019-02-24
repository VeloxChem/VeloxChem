//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef TDASigmaVectorDriver_hpp
#define TDASigmaVectorDriver_hpp

#include <cstdint>
#include <vector>

#include "mpi.h"

#include "MpiFunc.hpp"
#include "MemBlock.hpp"
#include "DenseMatrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MolecularOrbitals.hpp"
#include "ScreeningContainer.hpp"
#include "ExcitationVector.hpp"

/**
 Class CTDASigmaVectorDriver class computes sigma vector i.e. sigma = A * Z.
 
 @author Z. Rinkevicius
 */
class CTDASigmaVectorDriver
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
     Allocates and initializes set of sigma vectors according to given set of Z
     vectors.

     @param zVectors the set of Z vectors.
     @return the set of sigma vectors.
     */
    std::vector<CDenseMatrix> _allocSigmaVectors(const std::vector<CExcitationVector>& zVectors) const;
    
    /**
     Adds canonical Fock matrix contribution to A * Z vectors.

     @param sigmaVectors the set of sigma vectors.
     @param zVectors the set of Z vectors.
     @param molecularOrbitals the molecular orbitals.
     */
    void _addCanonicalFockContribution(      std::vector<CDenseMatrix>&      sigmaVectors,
                                       const std::vector<CExcitationVector>& zVectors,
                                       const CMolecularOrbitals&             molecularOrbitals) const;
    
    /**
     Adds first order Fock matrix contribution to A * Z vectors.

     @param sigmaVectors the set of sigma vectors.
     @param zVectors the set of Z vectors.
     @param screeningContainer the electron repulsion integrals screeners
            container.
     @param molecularOrbitals the molecular orbitals.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param comm the MPI communicator.
     */
    void _addFirstOrderFockContribution(      std::vector<CDenseMatrix>&      sigmaVectors,
                                        const std::vector<CExcitationVector>& zVectors,
                                        const CScreeningContainer&            screeningContainer,
                                        const CMolecularOrbitals&             molecularOrbitals,
                                        const CMolecule&                      molecule,
                                        const CMolecularBasis&                basis,
                                              MPI_Comm                        comm) const;
    
public:
    
    /**
     Creates a TDA sigma vector driver object.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CTDASigmaVectorDriver(const int32_t  globRank,
                          const int32_t  globNodes,
                                MPI_Comm comm);
    
    /**
     Destroys a TDA sigma vector driver object.
     */
    ~CTDASigmaVectorDriver();
    
    /**
     Computes sigma = A * Z vectors.

     @param zVectors the vector of excitation vectors.
     @param screeningContainer the electron repulsion integrals screeners
            container.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param comm the MPI communication.
     @return the vector of sigma vectors.
     */
    std::vector<CDenseMatrix> compute(const std::vector<CExcitationVector>& zVectors,
                                      const CScreeningContainer&            screeningContainer,
                                      const CMolecularOrbitals&             molecularOrbitals,
                                      const CMolecule&                      molecule,
                                      const CMolecularBasis&                basis,
                                            MPI_Comm                        comm) const;
};

#endif /* TDASigmaVectorDriver_hpp */
