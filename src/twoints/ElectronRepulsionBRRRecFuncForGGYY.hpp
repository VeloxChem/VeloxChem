//
//                           VELOXCHEM 1.0-RC
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

#include<cstdint>

#include "MemBlock2D.hpp"
#include "GtoPairsBlock.hpp"
#include "RecursionMap.hpp"

namespace eribrrfunc { // eribrrfunc namespace

    /**
    Computes batch of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY(      CMemBlock2D<double>& braBuffer,
                                      const CRecursionMap&       recursionMap,
                                      const CMemBlock2D<double>& abDistances,
                                      const CGtoPairsBlock&      braGtoPairsBlock,
                                      const CGtoPairsBlock&      ketGtoPairsBlock,
                                      const int32_t              nKetContrPairs,
                                      const int32_t              iContrPair);

    /**
    Computes sub-batch (0,75) of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY_0_75(      CMemBlock2D<double>& braBuffer,
                                           const CRecursionMap&       recursionMap,
                                           const CMemBlock2D<double>& abDistances,
                                           const CGtoPairsBlock&      braGtoPairsBlock,
                                           const CGtoPairsBlock&      ketGtoPairsBlock,
                                           const int32_t              nKetContrPairs,
                                           const int32_t              iContrPair);

    /**
    Computes sub-batch (75,150) of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY_75_150(      CMemBlock2D<double>& braBuffer,
                                             const CRecursionMap&       recursionMap,
                                             const CMemBlock2D<double>& abDistances,
                                             const CGtoPairsBlock&      braGtoPairsBlock,
                                             const CGtoPairsBlock&      ketGtoPairsBlock,
                                             const int32_t              nKetContrPairs,
                                             const int32_t              iContrPair);

    /**
    Computes sub-batch (150,225) of contracted (GG|G|XY electron repulsion integrals and stores
    results in integrals buffer.

    @param braBuffer the horizontal recursion buffer for bra side.
    @param recursionMap the recursion map object.
    @param abDistances the vector of distances R(AB) = A - B.
    @param braGtoPairsBlock the GTOs pairs block on ket side.
    @param ketGtoPairsBlock the GTOs pairs block on ket side.
    @param nKetContrPairs the number of contractes GTOs pairs on ket side.
    @param iContrPair the index of contracted GTO pair on bra side.
    */
    void compElectronRepulsionForGGXY_150_225(      CMemBlock2D<double>& braBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& abDistances,
                                              const CGtoPairsBlock&      braGtoPairsBlock,
                                              const CGtoPairsBlock&      ketGtoPairsBlock,
                                              const int32_t              nKetContrPairs,
                                              const int32_t              iContrPair);


} // eribrrfunc namespace

