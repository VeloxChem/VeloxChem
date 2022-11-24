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

#ifndef GeomRecFunc_hpp
#define GeomRecFunc_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "MemBlock2D.hpp"
#include "RecursionMap.hpp"

namespace geomrecfunc {  // geomrecfunc namespace

/**
Computes batch of primitive geometrical derivatives of generic (S|Z|X)  integrals
and stores results in primitives buffer.

@param primBuffer the primitives buffer.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param iPrimBuffSX the index of primitive SX integrals buffer.
@param iPrimBuffPX the index of primitive PX integrals buffer.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
@param axis the Cartesian axis to compute geometrical gradient.
*/
void compGeomForSX(CMemBlock2D<double>&       primBuffer,
                   const CMemBlock2D<double>& osFactors,
                   const int32_t              nOSFactors,
                   const int32_t              iPrimBuffSX,
                   const int32_t              iPrimBuffPX,
                   const CGtoBlock&           braGtoBlock,
                   const CGtoBlock&           ketGtoBlock,
                   const int32_t              iContrGto,
                   const char                 axis);

/**
Computes batch of primitive geometrical derivatives of generic (P|Z|X)  integrals
and stores results in primitives buffer.

@param primBuffer the primitives buffer.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param iPrimBuffPX the index of primitive PX integrals buffer.
@param iPrimBuffDX the index of primitive DX integrals buffer.
@param iPrimBuffSX the index of primitive SX integrals buffer.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
@param axis the Cartesian axis to compute geometrical gradient.
*/
void compGeomForPX(CMemBlock2D<double>&       primBuffer,
                   const CMemBlock2D<double>& osFactors,
                   const int32_t              nOSFactors,
                   const int32_t              iPrimBuffPX,
                   const int32_t              iPrimBuffDX,
                   const int32_t              iPrimBuffSX,
                   const CGtoBlock&           braGtoBlock,
                   const CGtoBlock&           ketGtoBlock,
                   const int32_t              iContrGto,
                   const char                 axis);

/**
Computes batch of primitive geometrical derivatives of generic (D|Z|X)  integrals
and stores results in primitives buffer.

@param primBuffer the primitives buffer.
@param osFactors the Obara-Saika recursion factors.
@param nOSFactors the number of Obara-Saika recursion factors.
@param iPrimBuffDX the index of primitive DX integrals buffer.
@param iPrimBuffFX the index of primitive FX integrals buffer.
@param iPrimBuffPX the index of primitive PX integrals buffer.
@param braGtoBlock the GTOs block on bra side.
@param ketGtoBlock the GTOs block on ket side.
@param iContrGto the index of contracted GTO on bra side.
@param axis the Cartesian axis to compute geometrical gradient.
*/
void compGeomForDX(CMemBlock2D<double>&       primBuffer,
                   const CMemBlock2D<double>& osFactors,
                   const int32_t              nOSFactors,
                   const int32_t              iPrimBuffDX,
                   const int32_t              iPrimBuffFX,
                   const int32_t              iPrimBuffPX,
                   const CGtoBlock&           braGtoBlock,
                   const CGtoBlock&           ketGtoBlock,
                   const int32_t              iContrGto,
                   const char                 axis);

}  // namespace geomrecfunc

#endif /* GeomRecFunc_hpp */
