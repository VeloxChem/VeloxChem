//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ThreeCenterEriFunc_hpp
#define ThreeCenterEriFunc_hpp

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "BoysFunction.hpp"
#include "VecIndexes.hpp"
#include "GtoBlock.hpp"
#include "GtoPairsBlock.hpp"

namespace t3erifunc { // t3erifunc namespace
    
    /**
     Computes batch of primitive (S|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param aqDistances the vector of distances R(AQ) = A - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSS(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CBoysFunction&        bfTable,
                                           CMemBlock<double>&    bfArguments,
                                           CMemBlock2D<double>&  bfValues,
                                     const int32_t               bfOrder,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  aqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSP(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSS(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSP(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSD(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSS(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSD(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSP(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSD(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSF(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSS(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSF(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSP(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSF(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSD(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSF(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSG(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SS)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSS(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSG(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SP)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSP(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSG(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SD)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSD(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSG(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SF)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSF(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SG)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSG(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SH)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSH(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SH)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSH(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SH)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSH(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SH)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSH(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SH)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSH(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SI)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSI(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SI)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSI(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SI)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSI(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SI)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSI(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SI)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSI(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSK(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSK(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSK(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSK(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SK)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSK(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|SL)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param wqDistances the vector of distances R(WQ) = W - Q.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForSSL(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  wqDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|SL)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForPSL(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|SL)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForDSL(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|SL)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForFSL(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|SL)^(m) electron repulsion integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param waDistances the vector of distances R(WA) = W - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoPairsBlock the GTOs pairs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronRepulsionForGSL(      CMemBlock2D<double>&  primBuffer,
                                     const CVecThreeIndexes&     recPattern,
                                     const std::vector<int32_t>& recIndexes,
                                     const CMemBlock2D<double>&  osFactors,
                                     const CMemBlock2D<double>&  waDistances,
                                     const CGtoBlock&            braGtoBlock,
                                     const CGtoPairsBlock&       ketGtoPairsBlock,
                                     const int32_t               iContrGto);

}  // t3erifunc namespace

#endif /* ThreeCenterEriFunc_hpp */
