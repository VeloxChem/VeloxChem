//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef LinearMomemtumIntegralsDriver_hpp
#define LinearMomemtumIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "ExecMode.hpp"
#include "LinearMomentumMatrix.hpp"
#include "GtoContainer.hpp"
#include "OneIntsDistributor.hpp"
#include "RecursionMap.hpp"

/**
 Class CLinearMomentumIntegralsDriver computes one-electron linear momentum
 integrals.
 
 @author Z. Rinkevicius
 */
class CLinearMomentumIntegralsDriver
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
     Comutes linear momentum integrals for pair of GTOs containers.
     
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the linear momentum matrix object.
     */
    CLinearMomentumMatrix _compLinearMomentumIntegrals(const CGtoContainer* braGtoContainer,
                                                       const CGtoContainer* ketGtoContainer) const;
    
    /**
     Computes linear momentum integrals for specific pair of GTOs blocks and
     stores integrals in given distribution buffer.
     
     @param distPatternX the pointer to X component of integrals distribution
     pattern.
     @param distPatternY the pointer to Y component of integrals distribution
     pattern.
     @param distPatternZ the pointer to Z component of integrals distribution
     pattern.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compLinearMomentumForGtoBlocks(      COneIntsDistribution* distPatternX,
                                               COneIntsDistribution* distPatternY,
                                               COneIntsDistribution* distPatternZ,
                                         const CGtoBlock&            braGtoBlock,
                                         const CGtoBlock&            ketGtoBlock) const;
    
    /**
     Computes batch of primitive linear momentum integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTO on bra side) x (all contracted GTOs on ket
     side).
     
     @param primBuffer the primitives buffer.
     @param recursionMap the recursion map for Obara-Saika recursion.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimLinearMomentumInts(      CMemBlock2D<double>&  primBuffer,
                                     const CRecursionMap&        recursionMap,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  abDistances,
                                     const CMemBlock2D<double>&  paDistances,
                                     const CMemBlock2D<double>&  pbDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoBlock&            ketGtoBlock,
                                     const int32_t               iContrGto) const;
    
    /**
     Sets recursion map object for linear momentum integrals of specified angular
     momentum.
     
     @param braAngularMomentum the angular momentum of bra side.
     @param ketAngularMomentum the angular momentum of ket side.
     @return the recursion map for linear momentum integrals.
     */
    CRecursionMap _setRecursionMap(const int32_t braAngularMomentum,
                                   const int32_t ketAngularMomentum,
                                   const int32_t maxNumberOfPrimitives) const;
    
public:
    
    /**
     Creates a linear momentum integrals driver object using MPI info.
     
     @param comm the MPI communicator.
     */
    CLinearMomentumIntegralsDriver(MPI_Comm comm);
    
    /**
     Destroys a linear momentum integrals driver object.
     */
    ~CLinearMomentumIntegralsDriver();
    
    /**
     Computes linear momentum integrals for molecule in specific basis set and
     stores results in electric dipole matrix object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @return the linear momentum matrix object.
     */
    CLinearMomentumMatrix compute(const CMolecule&       molecule,
                                  const CMolecularBasis& basis) const;
    
    /**
     Computes linear momentum integrals for molecule in two basis sets and stores
     results in electric dipole matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the linear momentum matrix object.
     */
    CLinearMomentumMatrix compute(const CMolecule&       molecule,
                                  const CMolecularBasis& braBasis,
                                  const CMolecularBasis& ketBasis) const;
    
    /**
     Computes linear momentum integrals for two molecules in basis set and
     stores results in electric dipole matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param basis the molecular basis.
     @return the linear momentum matrix object.
     */
    CLinearMomentumMatrix compute(const CMolecule&       braMolecule,
                                  const CMolecule&       ketMolecule,
                                  const CMolecularBasis& basis) const;
    
    /**
     Computes linear momentum integrals for two molecules in different basis
     sets and stores results in electric dipole matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the linear momentum matrix object.
     */
    CLinearMomentumMatrix compute(const CMolecule&       braMolecule,
                                  const CMolecule&       ketMolecule,
                                  const CMolecularBasis& braBasis,
                                  const CMolecularBasis& ketBasis) const;
    
    /**
     Computes linear momentum integrals blocks for pair of GTOs blocks and
     stores them into integrals batch.
     
     @param intsBatchX the pointer to batch buffer of X component of electric
            dipole integrals.
     @param intsBatchY the pointer to batch buffer of Y component of electric
            dipole integrals.
     @param intsBatchZ the pointer to batch buffer of Z component of electric
            dipole integrals.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(      double*    intsBatchX,
                       double*    intsBatchY,
                       double*    intsBatchZ,
                 const CGtoBlock& braGtoBlock,
                 const CGtoBlock& ketGtoBlock) const;
};

#endif /* LinearMomemtumIntegralsDriver_hpp */
