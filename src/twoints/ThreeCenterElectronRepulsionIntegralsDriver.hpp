//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ThreeCenterElectronRepulsionIntegralsDriver_hpp
#define ThreeCenterElectronRepulsionIntegralsDriver_hpp

#include <cstdint>

#include "mpi.h"

#include "Molecule.hpp"
#include "MolecularBasis.hpp"
#include "MemBlock2D.hpp"
#include "GtoPairsContainer.hpp"
#include "GtoContainer.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"

/**
 Class CThreeCenterElectronicRepulsionIntegralsDriver computes electron repulsion
 <f(r)| 1/|r-r'||g(r') h(r')> integrals.
 
 @author Z. Rinkevicius
 */
class CThreeCenterElectronRepulsionIntegralsDriver
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
     Comutes electronic repulsion integrals for pair of GTOs and GTOs
     pairs containers.
     
     @param braGtoContainer the GTOs container on bra side.
     @param ketGtoPairsContainer the GTOs pairs container on ket side.
     */
    void _compElectronRepulsionIntegrals(const CGtoContainer*      braGtoContainer,
                                         const CGtoPairsContainer* ketGtoPairsContainer) const;
    
   /**
    Computes batch of primitive electron repulsion integrals using Obara-Saika
    recursion and stores results in primitives buffer.
    Reference: S. Obara, A. Saika, J. Chem. Phys. 84, 3963 (1986).
    
    Batch size: (one contracted GTO on bra side) x (all contracted GTOs pairs on
                                                    ket side).
    
    @param primBuffer the primitives buffer.
    @param recPattern the recursion pattern.
    @param recIndexes the indexes of data blocks in recursion pattern.
    @param bfTable the Boys function evaluator.
    @param bfArguments the vector of Boys function arguments.
    @param bfValues the vector of Boys function values.
    @param bfOrder the order of Boys function.
    @param osFactors the Obara-Saika recursion factors.
    @param aqDistances the vector of distances R(AQ) = A - Q.
    @param waDistances the vector of distances R(WA) = W - A.
    @param wqDistances the vector of distances R(WQ) = W - Q.
    @param braGtoBlock the GTOs block on bra side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param iContrGto the index of contracted GTO on bra side.
    */
    void _compPrimElectronRepulsionInts(      CMemBlock2D<double>&  primBuffer,
                                        const CVecThreeIndexes&     recPattern,
                                        const std::vector<int32_t>& recIndexes,
                                        const CBoysFunction&        bfTable,
                                              CMemBlock<double>&    bfArguments,
                                              CMemBlock2D<double>&  bfValues,
                                        const int32_t               bfOrder,
                                        const CMemBlock2D<double>&  osFactors,
                                        const CMemBlock2D<double>&  aqDistances,
                                        const CMemBlock2D<double>&  waDistances,
                                        const CMemBlock2D<double>&  wqDistances,
                                        const CGtoBlock&            braGtoBlock,
                                        const CGtoPairsBlock&       ketGtoPairsBlock,
                                        const int32_t               iContrGto) const;
    
    /**
     Computes batch of contracted, half-transformed electron repulsion integrals
     using horizontal recursion and stores results in half-transformed integrals
     buffer.

     @param contrBuffer the contracted, half-transformed integrals buffer.
     @param recPattern the horizontal recursion pattern.
     @param recIndexes the indexes of data blocks in horizontal recursion pattern.
     @param cdDistances the vector of distances R(CD) = C - D.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void _compContrElectronRepulsionInts(      CMemBlock2D<double>&  contrBuffer,
                                         const CVecThreeIndexes&     recPattern,
                                         const std::vector<int32_t>& recIndexes,
                                         const CMemBlock2D<double>&  cdDistances,
                                         const CGtoBlock&            braGtoBlock,
                                         const CGtoPairsBlock&       ketGtoPairsBlock) const;
    
    /**
     Creates atoms list splitting pattern for generation of GTO blocks on each
     MPI process. The splitting pattern enables effective task based computation
     of integrals on systems with variable number of CPU cores.

     @param molecule the molecule.
     @param riBasis the molecular RI basis.
     @param gtoPairsContainer the GTOs pairs container.
     @return the splitting pattern for atoms list.
     */
    CMemBlock2D<int32_t> _getBatchesOfGtoBlocks(const CMolecule&          molecule,
                                                const CMolecularBasis&    riBasis,
                                                const CGtoPairsContainer& gtoPairsContainer) const;
    
    /**
     Gets Obara-Saika horizontal recursion pattern for specific combination of
     GTOs blocks on bra and ket sides.
     
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @return the vector of three indexes object with recursion pattern.
     */
    CVecThreeIndexes _getHorizontalRecursionPattern(const CGtoBlock&      braGtoBlock,
                                                    const CGtoPairsBlock& ketGtoPairsBlock) const;
    
    
    /**
     Gets Obara-Saika vertical recursion pattern for given set of leading terms.
     
     @param leadTerms the vector of leading terms in recursion pattern.
     @return the vector of three indexes object with recursion pattern.
     */
    CVecThreeIndexes _getVerticalRecursionPattern(const CVecThreeIndexes& leadTerms) const;
    
    /**
     Gets vector of unified indexes of primitive GTOs buffer for vertical
     Obara-Saika recursion pattern.
     
     @param recIndexes the vector of starting indexes of data blocks in recursion
            pattern.
     @param recPattern the recursion pattern.
     @param maxPrimGtos the maximum number of primitive GTOs in contracted
     GTO on bra side.
     @return the total number of blocks in recursion pattern.
     */
    int32_t _getIndexesForVerticalRecursionPattern(      std::vector<int32_t>& recIndexes,
                                                   const CVecThreeIndexes&     recPattern,
                                                   const int32_t               maxPrimGtos) const;
    
    /**
     Gets vector of unified indexes of contracted GTOs buffer.

     @param contrIndexes the vector of starting indexes in contracted GTOs buffer.
     @param contrListing the contracted integrals listing.
     @return the total number of blocks in contracted integrals listing.
     */
    int32_t _getIndexesForContractedIntegrals(      std::vector<int32_t>& contrIndexes,
                                              const CVecThreeIndexes&     contrListing) const;
    
    /**
     Gets vector of unified indexes of half transformed integrals buffer.
     
     @param intsIndexes the vector of starting indexes in half transformed
            integrals buffer.
     @param intsListing the half transformed integrals listing.
     @return the total number of blocks in half transformed integrals listing.
     */
    int32_t _getIndexesForHalfTransformedIntegrals(      std::vector<int32_t>& intsIndexes,
                                                   const CVecThreeIndexes&     intsListing) const;

public:
    
    /**
     Creates a three center electron repulsion integrals driver object using
     MPI info.
     
     @param globRank the the rank of MPI process.
     @param globNodes the total number of MPI processes.
     @param comm the MPI communicator.
     */
    CThreeCenterElectronRepulsionIntegralsDriver(const int32_t  globRank,
                                                 const int32_t  globNodes,
                                                       MPI_Comm comm);
    
    /**
     Destroys a three center electron repulsion integrals driver object.
     */
    ~CThreeCenterElectronRepulsionIntegralsDriver();
    
    /**
     Computes three center electron repulsion integrals for molecule with
     specific combination of AO and RI basis sets and stores results in
     compressed format on each MPI process within domain of MPI communicator.
     
     @param molecule the molecule.
     @param aoBasis the molecular AO basis.
     @param riBasis the molecular RI basis.
     @param threshold the integrals cut-off threshold. 
     @param comm the MPI communicator.
     */
    void  compute(const CMolecule&       molecule,
                  const CMolecularBasis& aoBasis,
                  const CMolecularBasis& riBasis,
                  const double           threshold,
                        MPI_Comm         comm) const;
    
    /**
     Computes electronic repulsion integrals for specific pair of GTOs and
     GTOs pairs blocks.
     
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     */
    void compElectronRepulsionForGtoBlocks(const CGtoBlock&      braGtoBlock,
                                           const CGtoPairsBlock& ketGtoPairsBlock) const;
};

#endif /* ThreeCenterElectronRepulsionIntegralsDriver_hpp */
