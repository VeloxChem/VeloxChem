//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef ElectronicPotentialRecFunc_hpp
#define ElectronicPotentialRecFunc_hpp

#include "MemBlock.hpp"
#include "MemBlock2D.hpp"
#include "BoysFunction.hpp"
#include "VecIndexes.hpp"
#include "GtoBlock.hpp"

namespace epotrecfunc { // epotrecfunc namespace
    
    /**
     Computes batch of primitive (S|g(r,r')|S)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForSS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CBoysFunction&        bfTable,
                                            CMemBlock<double>&    bfArguments,
                                            CMemBlock2D<double>&  bfValues,
                                      const int32_t               bfOrder,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  abDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|P)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForSP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  pbDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|S)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForPS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|P)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForPP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|D)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForSD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  pbDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|S)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForDS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|D)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForPD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|P)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForDP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|D)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForDD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|F)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForSF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  pbDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|S)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForFS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|F)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForPF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|P)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForFP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|F)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForDF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|D)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForFD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|F)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForFF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|g(r,r')|G)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForSG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  pbDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|S)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForGS(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|G)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForPG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|P)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForGP(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|g(r,r')|G)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForDG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|D)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForGD(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (F|g(r,r')|G)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForFG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|F)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForGF(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|g(r,r')|G)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForPG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (G|g(r,r')|G)^(m) electronic potential integrals
     and stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compElectronicPotentialForGG(      CMemBlock2D<double>&  primBuffer,
                                      const CVecThreeIndexes&     recPattern,
                                      const std::vector<int32_t>& recIndexes,
                                      const CMemBlock2D<double>&  osFactors,
                                      const CMemBlock2D<double>&  paDistances,
                                      const CGtoBlock&            braGtoBlock,
                                      const CGtoBlock&            ketGtoBlock,
                                      const int32_t               iContrGto);

} // epotrecfunc namespace

#endif /* ElectronicPotentialRecFunc_hpp */
