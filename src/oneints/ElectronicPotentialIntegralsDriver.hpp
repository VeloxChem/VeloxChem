//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ElectronicPotentialIntegralsDriver_hpp
#define ElectronicPotentialIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "OutputStream.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "OutputStream.hpp"
#include "ExecMode.hpp"
#include "ElectronicPotentialMatrix.hpp"
#include "GtoContainer.hpp"
#include "VecIndexes.hpp"
#include "SparseMatrix.hpp"
#include "BoysFunction.hpp"
#include "OutputStream.hpp"
#include "SystemClock.hpp"
#include "OneIntsDistributor.hpp"

/**
 Class CElectronicPotentialIntegralsDriver computes electronic potential
 <f(r)| 1/|r-r'||g(r')> integrals.
 
 @author Z. Rinkevicius
 */
class CElectronicPotentialIntegralsDriver
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
     Comutes electronic potential integrals for pair of GTOs containers.
     
     @param braGtoContainer the GTOs container for bra side.
     @param ketGtoContainer the GTOs container for ket side.
     @return the kinetic energy matrix object.
     */
    CElectronicPotentialMatrix _compElectronicPotentialIntegrals(const CGtoContainer* braGtoContainer,
                                                                 const CGtoContainer* ketGtoContainer) const;
    
    /**
     Computes electronic potential integrals for specific pair of GTOs blocks and
     stores integrals in matrix of predefined size.
     
     @param distPattern the pointer to integrals distribution pattern.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void _compElectronicPotentialForGtoBlocks(      COneIntsDistribution* distPattern,
                                              const CGtoBlock&            braGtoBlock,
                                              const CGtoBlock&            ketGtoBlock) const;
    
    /**
     Computes batch of primitive electronic potential integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTO on bra side) x (all contracted GTOs on ket
     side).
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void _compPrimElectronicPotentialInts(      CMemBlock2D<double>&  primBuffer,
                                          const CVecThreeIndexes&     recPattern,
                                          const std::vector<int32_t>& recIndexes,
                                          const CBoysFunction&        bfTable,
                                                CMemBlock<double>&    bfArguments,
                                                CMemBlock2D<double>&  bfValues,
                                          const int32_t               bfOrder,
                                          const CMemBlock2D<double>&  osFactors,
                                          const CMemBlock2D<double>&  abDistances,
                                          const CMemBlock2D<double>&  paDistances,
                                          const CMemBlock2D<double>&  pbDistances,
                                          const CGtoBlock&            braGtoBlock,
                                          const CGtoBlock&            ketGtoBlock,
                                          const int32_t               iContrGto) const;
    
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
     Prints electronic potential integrals computation time to output stream.
     
     @param timer the system clock timer.
     @param oStream the output stream.
     */
    void _printComputationTime(const CSystemClock&  timer,
                                     COutputStream& oStream) const;
    
public:
    
    /**
     Creates a electronic potential integrals driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CElectronicPotentialIntegralsDriver(const int32_t  globRank,
                                        const int32_t  globNodes,
                                              MPI_Comm comm);
    
    /**
     Destroys a electronic potential integrals driver object.
     */
    ~CElectronicPotentialIntegralsDriver();
    
    /**
     Computes electronic potential integrals for molecule in specific basis set
     and stores results in electronic potential matrix object.
     
     @param molecule the molecule.
     @param basis the molecular basis.
     @param oStream the output stream.
     @param comm the MPI communicator.
     @return the electronic potential matrix object.
     */
    CElectronicPotentialMatrix compute(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                             COutputStream&   oStream, 
                                             MPI_Comm         comm) const;
    
    /**
     Computes electronic potential integrals blocks for pair of GTOs blocks and
     stores them into integrals batch.
     
     @param intsBatch the pointer to integrals batch buffer.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     */
    void compute(      double*    intsBatch,
                 const CGtoBlock& braGtoBlock,
                 const CGtoBlock& ketGtoBlock) const;
};


#endif /* ElectronicPotentialIntegralsDriver_hpp */
