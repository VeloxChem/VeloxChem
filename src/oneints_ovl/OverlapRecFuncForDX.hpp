//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace ovlrecfunc {  // ovlrecfunc namespace

/**
Computes batch of primitive (D|S|D) overlap integrals and stores
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
void compOverlapForDD(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (D|S|F) overlap integrals and stores
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
void compOverlapForDF(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes sub-batch (0,30) of primitive (D|S|F) overlap integrals and stores
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
void compOverlapForDF_0_30(CMemBlock2D<double>&       primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const int32_t              nOSFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto);

/**
Computes sub-batch (30,60) of primitive (D|S|F) overlap integrals and stores
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
void compOverlapForDF_30_60(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (F|S|D) overlap integrals and stores
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
void compOverlapForFD(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes sub-batch (0,30) of primitive (F|S|D) overlap integrals and stores
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
void compOverlapForFD_0_30(CMemBlock2D<double>&       primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const int32_t              nOSFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto);

/**
Computes sub-batch (30,60) of primitive (F|S|D) overlap integrals and stores
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
void compOverlapForFD_30_60(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (D|S|G) overlap integrals and stores
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
void compOverlapForDG(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (D|S|G) overlap integrals and stores
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
void compOverlapForDG_0_45(CMemBlock2D<double>&       primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const int32_t              nOSFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (D|S|G) overlap integrals and stores
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
void compOverlapForDG_45_90(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes batch of primitive (G|S|D) overlap integrals and stores
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
void compOverlapForGD(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes sub-batch (0,45) of primitive (G|S|D) overlap integrals and stores
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
void compOverlapForGD_0_45(CMemBlock2D<double>&       primBuffer,
                           const CRecursionMap&       recursionMap,
                           const CMemBlock2D<double>& osFactors,
                           const int32_t              nOSFactors,
                           const CMemBlock2D<double>& paDistances,
                           const CGtoBlock&           braGtoBlock,
                           const CGtoBlock&           ketGtoBlock,
                           const int32_t              iContrGto);

/**
Computes sub-batch (45,90) of primitive (G|S|D) overlap integrals and stores
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
void compOverlapForGD_45_90(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const int32_t              nOSFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

}  // namespace ovlrecfunc
