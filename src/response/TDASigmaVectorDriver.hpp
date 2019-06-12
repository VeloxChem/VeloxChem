//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef TDASigmaVectorDriver_hpp
#define TDASigmaVectorDriver_hpp

#include <cstdint>
#include <vector>

#include "mpi.h"

#include "DenseMatrix.hpp"
#include "ExcitationVector.hpp"
#include "MemBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularOrbitals.hpp"
#include "Molecule.hpp"
#include "MpiFunc.hpp"
#include "ScreeningContainer.hpp"

/**
 Class CTDASigmaVectorDriver class computes sigma vector i.e. sigma = A * Z.

 @author Z. Rinkevicius
 */
class CTDASigmaVectorDriver
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
     The MPI communicator.
     */
    MPI_Comm _locComm;

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
    void _addCanonicalFockContribution(std::vector<CDenseMatrix>&            sigmaVectors,
                                       const std::vector<CExcitationVector>& zVectors,
                                       const CMolecularOrbitals&             molecularOrbitals) const;

    /**
     Adds first order Fock matrix contribution to A * Z vectors.

     @param sigmaVectors the set of sigma vectors.
     @param zVectors the set of Z vectors.
     @param isTripletStates the flag indicating A matrix construction for
            triplet excited states, instead of default singlet excited states.
     @param screeningContainer the electron repulsion integrals screeners
            container.
     @param molecularOrbitals the molecular orbitals.
     @param molecule the molecule.
     @param basis the molecular basis.
     */
    void _addFirstOrderFockContribution(std::vector<CDenseMatrix>&            sigmaVectors,
                                        const std::vector<CExcitationVector>& zVectors,
                                        const bool                            isTripletStates,
                                        const CScreeningContainer&            screeningContainer,
                                        const CMolecularOrbitals&             molecularOrbitals,
                                        const CMolecule&                      molecule,
                                        const CMolecularBasis&                basis) const;

   public:
    /**
     Creates a TDA sigma vector driver object.

     @param comm the MPI communicator.
     */
    CTDASigmaVectorDriver(MPI_Comm comm);

    /**
     Destroys a TDA sigma vector driver object.
     */
    ~CTDASigmaVectorDriver();

    /**
     Computes sigma = A * Z vectors.

     @param zVectors the vector of excitation vectors.
     @param isTripletStates the flag indicating A matrix construction for
            triplet excited states, instead of default singlet excited states.
     @param screeningContainer the electron repulsion integrals screeners
     container.
     @param molecule the molecule.
     @param basis the molecular basis.
     @return the vector of sigma vectors.
     */
    std::vector<CDenseMatrix> compute(const std::vector<CExcitationVector>& zVectors,
                                      const bool                            isTripletStates,
                                      const CScreeningContainer&            screeningContainer,
                                      const CMolecularOrbitals&             molecularOrbitals,
                                      const CMolecule&                      molecule,
                                      const CMolecularBasis&                basis) const;
};

#endif /* TDASigmaVectorDriver_hpp */
