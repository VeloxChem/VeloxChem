//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef EriFuncForL_hpp
#define EriFuncForL_hpp

#include "MemBlock2D.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"
#include "GtoPairsBlock.hpp"

namespace erifunc { // erifunc namespace
    
    /**
     Computes batch of primitive (SS|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSSSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SS)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSPSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SP)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSDSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SD)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSFSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SF)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSGSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SG)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSKSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SK)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SL|g(r,r')|SL)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSLSL(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const bool                  isBraEqualKet,
                                      const int32_t               iContrPair);
    
} // erifunc namespace
    

#endif /* EriFuncForL_hpp */
