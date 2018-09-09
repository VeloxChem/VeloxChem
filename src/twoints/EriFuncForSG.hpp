//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef EriFuncForSG_hpp
#define EriFuncForSG_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"
#include "GtoPairsBlock.hpp"

namespace erifunc { // erifunc namespace
    
    /**
     Computes batch of primitive (SS|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param pqDistances the vector of distances R(PQ) = P - Q.
     @param braGtoPairsBlock the GTOs psirs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
            blocks.
     @param iContrPair the index of contracted GTOs apir on bra side.
     */
    void compElectronRepulsionForSSSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CBoysFunction&        bfTable,
                                            CMemBlock<double>&    bfArguments,
                                            CMemBlock2D<double>&  bfValues,
                                      const int32_t               bfOrder,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  pqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SS|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
            blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSSSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
            blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSPSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSPSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SS|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSSSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSDSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSPSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSDSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSDSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SS|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSSSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSFSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSPSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSFSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSDSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSFSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSFSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SS|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSSSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSGSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSPSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSGSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSDSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSGSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSFSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSGSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param isBraEqualKet the flag for equality for bra and ket GTOs pairs
     blocks.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSGSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
} // erifunc namespace

#endif /* EriFuncForSG_hpp */
