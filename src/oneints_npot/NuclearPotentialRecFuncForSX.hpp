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

#include "BoysFunction.hpp"
#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

/**
Computes batch of primitive (S|A|S) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param abDistances the vector of distances R(AB) = A - B.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForSS(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CBoysFunction&       bfTable,
                               CMemBlock<double>&         bfArguments,
                               CMemBlock2D<double>&       bfValues,
                               const int32_t              bfOrder,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& abDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (S|A|P) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param pbDistances the vector of distances R(PB) = P - B.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForSP(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& pbDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (P|A|S) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param paDistances the vector of distances R(PA) = P - A.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForPS(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (S|A|D) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForSD(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& pbDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (D|A|S) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForDS(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (S|A|F) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForSF(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& pbDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (F|A|S) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForFS(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (S|A|G) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param pbDistances the vector of distances R(PB) = P - B.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForSG(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& pbDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

/**
Computes batch of primitive (G|A|S) nuclear potential integrals and stores
results in primitives buffer.

@param primBuffer the primitives buffer.
@param recursionMap the recursion map object.
@param osFactors the Obara-Saika recursion factors.
@param paDistances the vector of distances R(PA) = P - A.
@param pcDistances the vector of distances R(PC) = P - C.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
*/
void compNuclearPotentialForGS(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto);

}  // namespace npotrecfunc
