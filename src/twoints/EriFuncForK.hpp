//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef EriFuncForK_hpp
#define EriFuncForK_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"
#include "GtoPairsBlock.hpp"

namespace erifunc { // erifunc namespace

    /**
     Computes batch of primitive (SS|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSSSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSPSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSDSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSFSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSGSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSHSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SH)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSISK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SI)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SK|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wpDistances the vector of distances R(WP) = W - P.
     @param braGtoPairsBlock the GTOs pairs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param nKetPrimPairs the number of primitive GTOs pairs on ket side.
     @param iContrPair the index of contracted GTOs pair on bra side.
     */
    void compElectronRepulsionForSKSK(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    
} // erifunc namespace

#endif /* EriFuncForK_hpp */
