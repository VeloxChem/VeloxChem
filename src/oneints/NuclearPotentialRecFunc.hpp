//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

#ifndef NuclearPotentialRecFunc_hpp
#define NuclearPotentialRecFunc_hpp

#include "MemBlock2D.hpp"
#include "MemBlock.hpp"
#include "GtoBlock.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"

namespace npotrecfunc { // npotrecfunc namespace
  
    /**
     Computes batch of primitive (S|A(0)|S)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param bfTable the Boys function evaluator.
     @param bfArguments the vector of Boys function arguments.
     @param bfValues the vector of Boys function values.
     @param bfOrder the order of Boys function.
     @param osFactors the Obara-Saika recursion factors.
     @param abDistances the vector of distances R(AB) = A - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForSS(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CBoysFunction&        bfTable,
                                         CMemBlock<double>&    bfArguments,
                                         CMemBlock2D<double>&  bfValues,
                                   const int32_t               bfOrder,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  abDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|A(0)|P)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForSP(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  pbDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|A(0)|S)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForPS(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  paDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|A(0)|P)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForPP(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  paDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (S|A(0)|D)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param pbDistances the vector of distances R(PB) = P - B.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForSD(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  pbDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|A(0)|S)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForDS(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  paDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (P|A(0)|D)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForPD(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  paDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|A(0)|P)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForDP(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  paDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
    
    /**
     Computes batch of primitive (D|A(0)|D)^(m) nuclear potential integrals and
     stores results in primitives buffer.
     
     @param primBuffer the primitives buffer.
     @param recPattern the recursion pattern.
     @param recIndexes the indexes of data blocks in recursion pattern.
     @param osFactors the Obara-Saika recursion factors.
     @param paDistances the vector of distances R(PA) = P - A.
     @param pcDistances the vector of distances R(PC) = P - C.
     @param braGtoBlock the GTOs block on bra side.
     @param ketGtoBlock the GTOs block on ket side.
     @param iContrGto the index of contracted GTO on bra side.
     */
    void compNuclearPotentialForDD(      CMemBlock2D<double>&  primBuffer,
                                   const CVecThreeIndexes&     recPattern,
                                   const std::vector<int32_t>& recIndexes,
                                   const CMemBlock2D<double>&  osFactors,
                                   const CMemBlock2D<double>&  paDistances,
                                   const CMemBlock2D<double>&  pcDistances,
                                   const CGtoBlock&            braGtoBlock,
                                   const CGtoBlock&            ketGtoBlock,
                                   const int32_t               iContrGto);
} // npotrecfunc namespace

#endif /* NuclearPotentialRecFunc_hpp */
