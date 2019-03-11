//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#ifndef ElectricFieldRecFunc_hpp
#define ElectricFieldRecFunc_hpp

#include "MemBlock2D.hpp"
#include "MemBlock.hpp"
#include "GtoBlock.hpp"
#include "VecIndexes.hpp"
#include "BoysFunction.hpp"

namespace efieldrecfunc { // npotrecfunc namespace
    
    /**
     Computes batch of primitive (S|A(1)|S)^(m) electric field integrals and
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
    void compElectricFieldForSS(      CMemBlock2D<double>&  primBuffer,
                                const CVecFourIndexes&      recPattern,
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
    
} // efieldrecfunc namespace

#endif /* ElectricFieldRecFunc_hpp */
