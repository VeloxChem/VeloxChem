//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Created by Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.

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
    void compElectronicPotentialForSSS(      CMemBlock2D<double>&  primBuffer,
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
    
}  // t3erifunc namespace

#endif /* ThreeCenterEriFunc_hpp */
