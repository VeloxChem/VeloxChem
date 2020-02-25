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

namespace lmomrecfunc {  // lmomrecfunc namespace

/**
Computes batch of primitive (S|P|S) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForSS(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (S|P|P) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForSP(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (P|P|S) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForPS(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (S|P|D) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForSD(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (D|P|S) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForDS(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (S|P|F) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForSF(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (F|P|S) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForFS(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (S|P|G) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForSG(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& pbDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

/**
Computes batch of primitive (G|P|S) linear momentum integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compLinearMomentumForGS(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto);

}  // namespace lmomrecfunc
