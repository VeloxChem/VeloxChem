//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ElectronRepulsionIntegralsDriver_hpp
#define ElectronRepulsionIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "OutputStream.hpp"
#include "GtoPairsBlock.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"

/**
 Class CElectronicRepulsionIntegralsDriver computes electron repulsion
 <f(r)g(r)| 1/|r-r'||h(r') i(r')> integrals.
 
 @author Z. Rinkevicius
 */
class CElectronRepulsionIntegralsDriver
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
     Gets Obara-Saika bra side horizontal recursion pattern for specific
     combination of GTOs pairs blocks on bra and ket sides.
     
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @return the vector of three indexes object with recursion pattern.
     */
    CVecFourIndexes _getBraHorizontalRecursionPattern(const CGtoPairsBlock& braGtoPairsBlock,
                                                      const CGtoPairsBlock& ketGtoPairsBlock) const;
    
    /**
     Gets Obara-Saika ket side horizontal recursion pattern for given set of leading terms.
     
     @param leadTerms the vector of leading terms in recursion pattern.
     @return the vector of three indexes object with recursion pattern.
     */
    CVecThreeIndexes _getKetHorizontalRecursionPattern(const CVecThreeIndexes& leadTerms) const;
    
    /**
     Gets Obara-Saika vertical recursion pattern for given set of leading terms.
     
     @param leadTerms the vector of leading terms in recursion pattern.
     @return the vector of three indexes object with recursion pattern.
     */
    CVecThreeIndexes _getVerticalRecursionPattern(const CVecThreeIndexes& leadTerms) const;
    
    
    /**
     Gets vector of unified indexes of primitive GTOs pairs buffer for vertical
     Obara-Saika recursion pattern.
     
     @param recIndexes the vector of starting indexes of data blocks in recursion
            pattern.
     @param recPattern the recursion pattern.
     @param maxPrimPairs the maximum number of primitive GTOs pairs in contracted
            GTOs pair on bra side.
     @return the total number of blocks in recursion pattern.
     */
    int32_t _getIndexesForVerticalRecursionPattern(      std::vector<int32_t>& recIndexes,
                                                   const CVecThreeIndexes&     recPattern,
                                                   const int32_t               maxPrimPairs) const;
    
    /**
     Gets vector of unified indexes of contracted GTOs pairs buffer.
     
     @param contrIndexes the vector of starting indexes in contracted GTOs buffer.
     @param contrListing the contracted integrals listing.
     @return the total number of blocks in contracted integrals listing.
     */
    int32_t _getIndexesForContractedIntegrals(      std::vector<int32_t>& contrIndexes,
                                              const CVecThreeIndexes&     contrListing) const;
    
    /**
     Computes batch of primitive electron repulsion integrals using Obara-Saika
     recursion and stores results in primitives buffer.
     Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
     
     Batch size: (one contracted GTOs pair on bra side) x (all contracted GTOs
                  pairs on ket side).
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param pqDistances the vector of distances R(PQ) = P - Q.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
            blocks.
     @param iContrPair the index of contracted GTO pair on bra side.
     */
    void _compPrimElectronRepulsionInts(      CMemBlock2D<double>&  primBuffer,
                                        const CVecThreeIndexes&     recPattern,
                                        const std::vector<int32_t>& recIndexes,
                                        const CBoysFunction&        bfTable,
                                              CMemBlock<double>&    bfArguments,
                                              CMemBlock2D<double>&  bfValues,
                                        const int32_t               bfOrder,
                                        const CMemBlock2D<double>&  osFactors,
                                        const CMemBlock2D<double>&  pqDistances,
                                        const CMemBlock2D<double>&  wpDistances,
                                        const CMemBlock2D<double>&  wqDistances,
                                        const CGtoPairsBlock&       braGtoPairsBlock,
                                        const CGtoPairsBlock&       ketGtoPairsBlock,
                                        const bool                  isBraEqualKet,
                                        const int32_t               iContrPair) const;
    
public:
    
    /**
     Creates an electron repulsion integrals driver object using MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CElectronRepulsionIntegralsDriver(const int32_t  globRank,
                                      const int32_t  globNodes,
                                            MPI_Comm comm);
    
    /**
     Destroys an electron repulsion integrals driver object.
     */
    ~CElectronRepulsionIntegralsDriver();
    
    /**
     Computes electron repulsion integrals for molecule with specific AO basis
     set and process results according to provided distribution function.
     
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @param threshold the integrals cut-off threshold.
     @param oStream the output stream.
     @param comm the MPI communicator.
     */
    void  compute(const CMolecule&       molecule,
                  const CMolecularBasis& aoBasis,
                  const double           threshold,
                        COutputStream&   oStream,
                        MPI_Comm         comm) const;
    
    /**
     Computes electronic repulsion integrals for combination of GTOs pairs
     blocks.
     
     @param braGtoPairsBlock the GTOs pairsblock on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForGtoPairsBlocks(const CGtoPairsBlock& braGtoPairsBlock,
                                                const CGtoPairsBlock& ketGtoPairsBlock) const;
};


#endif /* ElectronRepulsionIntegralsDriver_hpp */
