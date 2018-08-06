//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef NuclearPotentialIntegralsDriver_hpp
#define NuclearPotentialIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "NuclearPotentialMatrix.hpp"
#include "GtoContainer.hpp"
#include "SparseMatrix.hpp"
#include "ThreeIndexes.hpp"
#include "VecIndexes.hpp"

/**
 Class CNuclearPotentialIntegralsDriver computes one-electron nuclear potential
 integrals.
 
 @author Z. Rinkevicius
 */
class CNuclearPotentialIntegralsDriver
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
     Comutes nuclear potential integrals for pair of GTOs containers.
     
     @param charges the vector of point charges.
     @param coordinates the vector of point charges coordines.
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix _compNuclearPotentialIntegrals(const CMemBlock<double>*   charges,
                                                           const CMemBlock2D<double>* coordinates,
                                                           const CGtoContainer*       braGtoContainer,
                                                           const CGtoContainer*       ketGtoContainer) const;
    
    /**
     Computes nuclear potential integrals for specific pair of GTOs blocks.
     
     @param sparseBuffer the buffer of sparce matrices.
     @param charges the vector of point charges.
     @param coordinates the vector of point charges coordines.
     @param braGtoContainer the GTOs container for bra side.
     @param iBraGtoBlock the index of GTOs block object in GTOs container for
     bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @param iKetGtoBlock the index of GTOs block object in GTOs container for
     ket side.
     */
    void _compNuclearPotentialForGtoBlocks(      CSparseMatrix*       sparseBuffer,
                                           const CMemBlock<double>*   charges,
                                           const CMemBlock2D<double>* coordinates,
                                           const CGtoContainer*       braGtoContainer,
                                           const int32_t              iBraGtoBlock,
                                           const CGtoContainer*       ketGtoContainer,
                                           const int32_t              iKetGtoBlock) const;
    
    /**
     Gets Obara-Saika recursion pattern for specific combination of GTOs blocks
     on bra and ket sides.
     
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @return the vector of three indexes object with recursion pattern.
     */
    CVecThreeIndexes _getRecursionPattern(const CGtoBlock& braGtoBlock,
                                          const CGtoBlock& ketGtoBlock) const;
    
    /**
     Gets vector of unified indexes of primitive GTOs buffer for specific
     Obara-Saika recursion pattern.
     
     @param recIndexes the vector of starting indexes of data blocks in recursion
     pattern.
     @param recPattern the recursion pattern.
     @param maxPrimGtos the maximum number of primitive GTOs in contracted
     GTO on bra side.
     @return the total number of blocks in recursion pattern.
     */
    int32_t _getIndexesForRecursionPattern(      std::vector<int32_t>& recIndexes,
                                           const CVecThreeIndexes&     recPattern,
                                           const int32_t               maxPrimGtos) const;
    
    /**
     Creates buffer of sparse matrices for all allowed combinations of GTOs
     blocks from GTOs containers on bra and ket sides.
     
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the vector of sparse matrices.
     */
    CSparseMatrix* _createSparseBuffer(const CGtoContainer* braGtoContainer,
                                       const CGtoContainer* ketGtoContainer) const;
    
    /**
     Adds single point charge contribution from primitives recursion buffer to
     primitives accumulation buffer, which contains primitive nuclear potential
     integrals.

     @param accBuffer the primitive integrals accumulation buffer.
     @param primBuffer the primitives recursion buffer.
     @param primIndex the index of specific integrals in primitives recursion
            buffer.
     @param charges the vector of point charges.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     @param iPointCharge the index of point charge in vector point charges.
     */
    void _addPointChargeContribution(      CMemBlock2D<double>& accBuffer,
                                     const CMemBlock2D<double>& primBuffer,
                                     const int32_t              primIndex,
                                     const CMemBlock<double>&   charges,
                                     const CGtoBlock&           braGtoBlock,
                                     const CGtoBlock&           ketGtoBlock,
                                     const int32_t              iContrGto,
                                     const int32_t              iPointCharge) const;
    
public:
    
    /**
     Creates a nuclear potential integrals driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CNuclearPotentialIntegralsDriver(const int32_t  globRank,
                                     const int32_t  globNodes,
                                     MPI_Comm comm);
    
    /**
     Destroys a nuclear potential integrals driver object.
     */
    ~CNuclearPotentialIntegralsDriver();
    
    /**
     Computes nuclear potential integrals for molecule in specific basis set and
     stores results in nuclear potential matrix object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param comm the MPI communicator.
     @return the nuclear potential matrix object.
     */
    CNuclearPotentialMatrix compute(const CMolecule&       molecule,
                                    const CMolecularBasis& basis,
                                          MPI_Comm         comm) const;
};

#endif /* NuclearPotentialIntegralsDriver_hpp */
