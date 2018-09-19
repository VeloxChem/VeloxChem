//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef EriFuncForH_hpp
#define EriFuncForH_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"
#include "GtoPairsBlock.hpp"

namespace erifunc { // erifunc namespace

    /**
     Computes batch of primitive (SS|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSSSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SS)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSPSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SP)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSDSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SD)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSFSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SF)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSGSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SG)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
} // erifunc namespace

#endif /* EriFuncForH_hpp */
