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

namespace efieldrecfunc {  // efieldrecfunc namespace

/**
Computes batch of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CMemBlock2D<double>& pcDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto);

/**
Computes sub-batch (0,50) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_0_50(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
                                 const CGtoBlock&           braGtoBlock,
                                 const CGtoBlock&           ketGtoBlock,
                                 const int32_t              iContrGto);

/**
Computes sub-batch (50,100) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_50_100(CMemBlock2D<double>&       primBuffer,
                                   const CRecursionMap&       recursionMap,
                                   const CMemBlock2D<double>& osFactors,
                                   const CMemBlock2D<double>& paDistances,
                                   const CMemBlock2D<double>& pcDistances,
                                   const CGtoBlock&           braGtoBlock,
                                   const CGtoBlock&           ketGtoBlock,
                                   const int32_t              iContrGto);

/**
Computes sub-batch (100,150) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_100_150(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (150,200) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_150_200(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (200,250) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_200_250(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (250,300) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_250_300(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (300,350) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_300_350(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (350,400) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_350_400(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

/**
Computes sub-batch (400,450) of primitive (G|E|F) electric field integrals and stores
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
void compElectricFieldForGF_400_450(CMemBlock2D<double>&       primBuffer,
                                    const CRecursionMap&       recursionMap,
                                    const CMemBlock2D<double>& osFactors,
                                    const CMemBlock2D<double>& paDistances,
                                    const CMemBlock2D<double>& pcDistances,
                                    const CGtoBlock&           braGtoBlock,
                                    const CGtoBlock&           ketGtoBlock,
                                    const int32_t              iContrGto);

}  // namespace efieldrecfunc
