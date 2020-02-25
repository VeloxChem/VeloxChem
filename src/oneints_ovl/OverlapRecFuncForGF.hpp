//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace ovlrecfunc {  // ovlrecfunc namespace

/**
Computes batch of primitive (G|S|F) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForGF(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes sub-batch (0,50) of primitive (G|S|F) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForGF_0_50(CMemBlock2D<double>&       primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const int32_t              nOSFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto);

/**
Computes sub-batch (50,100) of primitive (G|S|F) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForGF_50_100(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const int32_t              nOSFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes sub-batch (100,150) of primitive (G|S|F) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForGF_100_150(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const int32_t              nOSFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto);

}  // namespace ovlrecfunc
