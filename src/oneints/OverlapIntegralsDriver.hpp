//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef OverlapIntegralsDriver_hpp
#define OverlapIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "ExecMode.hpp"
#include "OverlapMatrix.hpp"
#include "GtoContainer.hpp"
#include "OneIntsDistributor.hpp"

/**
 Class COverlapIntegralsDriver computes one-electron overlap integrals.
 
 @author Z. Rinkevicius
 */
class COverlapIntegralsDriver
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
     Comutes overlap integrals for pair of GTOs containers.

     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the overlap matrix object.
     */
    COverlapMatrix _compOverlapIntegrals(const CGtoContainer* braGtoContainer,
                                         const CGtoContainer* ketGtoContainer) const;
    
    /**
     Computes overlap integrals for specific pair of GTOs blocks and stores
     integrals in given distribution buffer.
     
     @param distPattern the pointer to integrals distribution pattern.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compOverlapForGtoBlocks(      COneIntsDistribution* distPattern,
                                  const CGtoBlock&            braGtoBlock,
                                  const CGtoBlock&            ketGtoBlock) const;
    
    /**
     Computes batch of primitive overlap integrals using Obara-Saika recursion
     and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTO on bra side) x (all contracted GTOs on ket
     side).
     
     @param primBuffer the primitive integrals buffer.
     @param auxBuffer the auxilary integrals buffer.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param pa2bDistances the vector of distances products R(PA) x R(PB).
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimOverlapInts(      CMemBlock2D<double>&  primBuffer,
                                    CMemBlock2D<double>&  auxBuffer,
                              const CMemBlock2D<double>&  osFactors,
                              const CMemBlock2D<double>&  abDistances,
                              const CMemBlock2D<double>&  paDistances,
                              const CMemBlock2D<double>&  pbDistances,
                              const CMemBlock2D<double>&  pa2bDistances,
                              const CGtoBlock&            braGtoBlock,
                              const CGtoBlock&            ketGtoBlock,
                              const int32_t               iContrGto) const;
    
public:
    
    /**
     Creates a overlap integrals driver object using MPI info.
     
     @param comm the MPI communicator.
     */
    COverlapIntegralsDriver(MPI_Comm comm);
    
    /**
     Destroys a overlap integrals driver object.
     */
    ~COverlapIntegralsDriver();
    
    /**
     Computes overlap integrals for molecule in specific basis set and stores
     results in overlap matrix object.

     @param molecule the molecule.
     @param basis the molecular basis.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       molecule,
                           const CMolecularBasis& basis) const;
    
    /**
     Computes overlap integrals for molecule in two basis sets and stores
     results in overlap matrix object.
     
     @param molecule the molecule.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       molecule,
                           const CMolecularBasis& braBasis,
                           const CMolecularBasis& ketBasis) const;
    
    /**
     Computes overlap integrals for two molecules in basis set and stores
     results in overlap matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param basis the molecular basis.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       braMolecule,
                           const CMolecule&       ketMolecule,
                           const CMolecularBasis& basis) const;
    
    /**
     Computes overlap integrals for two molecules in different basis sets and
     stores results in overlap matrix object.
     
     @param braMolecule the molecule for bra side of overlap matrix.
     @param ketMolecule the molecule for ket side of overlap matrix.
     @param braBasis the molecular basis for bra side of overlap matrix.
     @param ketBasis the molecular basis for ket side of overlap matrix.
     @return the overlap matrix object.
     */
    COverlapMatrix compute(const CMolecule&       braMolecule,
                           const CMolecule&       ketMolecule,
                           const CMolecularBasis& braBasis,
                           const CMolecularBasis& ketBasis) const;
    
    /**
     Computes overlap integrals blocks for pair of GTOs blocks and stores them
     into integrals batch.
 
     @param intsBatch the pointer to integrals batch buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(      double*    intsBatch,
                 const CGtoBlock& braGtoBlock,
                 const CGtoBlock& ketGtoBlock) const;
};

#endif /* OverlapIntegralsDriver_hpp */
