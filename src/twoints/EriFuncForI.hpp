//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef EriFuncForI_hpp
#define EriFuncForI_hpp

#include <cstdint>
#include <vector>

#include "MemBlock2D.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"
#include "GtoPairsBlock.hpp"

namespace erifunc { // erifunc namespace

    /**
     Computes batch of primitive (SS|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSSSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wqDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SS)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SP|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSPSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SP)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SD|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSDSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SD)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SF|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSFSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SF)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SG|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSGSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SG)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SH|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSHSI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SH)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISH(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
    
    /**
     Computes batch of primitive (SI|g(r,r')|SI)^(m) electron repulsion integrals
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
    void compElectronRepulsionForSISI(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  wpDistances,
                                      const CGtoPairsBlock&       braGtoPairsBlock,
                                      const CGtoPairsBlock&       ketGtoPairsBlock,
                                      const int32_t               nKetPrimPairs,
                                      const int32_t               iContrPair);
} // erifunc namespace

#endif /* EriFuncForI_hpp */
