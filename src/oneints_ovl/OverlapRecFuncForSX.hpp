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
Computes batch of primitive (S|S|S) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param abDistances the vector of distances R(AB) = A - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForSS(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& abDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (S|S|P) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForSP(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& pbDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (P|S|S) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param paDistances the vector of distances R(PA) = P - A.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForPS(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (S|S|D) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForSD(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& pbDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (D|S|S) overlap integrals and stores
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
void compOverlapForDS(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (S|S|F) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForSF(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& pbDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (F|S|S) overlap integrals and stores
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
void compOverlapForFS(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (S|S|G) overlap integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compOverlapForSG(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& pbDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

/**
Computes batch of primitive (G|S|S) overlap integrals and stores
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
void compOverlapForGS(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto);

}  // namespace ovlrecfunc
