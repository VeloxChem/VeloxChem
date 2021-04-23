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

#include "LinearMomentumRecFuncForGG.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForGG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForGG_0_49(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_49_98(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_98_147(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_147_195(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_195_243(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_243_291(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_291_339(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_339_387(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_387_435(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_435_483(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_483_531(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_531_579(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_579_627(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGG_627_675(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForGG_0_49(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
{
    // Batch of Integrals (0,49)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_xxx_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx);

        auto tpy_xxx_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx);

        auto tpz_xxx_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx);

        auto tpx_xxx_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 1);

        auto tpy_xxx_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 1);

        auto tpz_xxx_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 1);

        auto tpx_xxx_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 2);

        auto tpy_xxx_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 2);

        auto tpz_xxx_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 2);

        auto tpx_xxx_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 3);

        auto tpy_xxx_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 3);

        auto tpz_xxx_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 3);

        auto tpx_xxx_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 4);

        auto tpy_xxx_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 4);

        auto tpz_xxx_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 4);

        auto tpx_xxx_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 5);

        auto tpy_xxx_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 5);

        auto tpz_xxx_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 5);

        auto tpx_xxx_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 6);

        auto tpy_xxx_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 6);

        auto tpz_xxx_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 6);

        auto tpx_xxx_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 7);

        auto tpy_xxx_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 7);

        auto tpz_xxx_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 7);

        auto tpx_xxx_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 8);

        auto tpy_xxx_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 8);

        auto tpz_xxx_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 8);

        auto tpx_xxx_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 9);

        auto tpy_xxx_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 9);

        auto tpz_xxx_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 9);

        auto tpx_xxx_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 10);

        auto tpy_xxx_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 10);

        auto tpz_xxx_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 10);

        auto tpx_xxx_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 11);

        auto tpy_xxx_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 11);

        auto tpz_xxx_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 11);

        auto tpx_xxx_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 12);

        auto tpy_xxx_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 12);

        auto tpz_xxx_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 12);

        auto tpx_xxx_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 13);

        auto tpy_xxx_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 13);

        auto tpz_xxx_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 13);

        auto tpx_xxx_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 14);

        auto tpy_xxx_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 14);

        auto tpz_xxx_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 14);

        auto tpx_xxy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 15);

        auto tpy_xxy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 15);

        auto tpz_xxy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 15);

        auto tpx_xxy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 16);

        auto tpx_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx);

        auto tpy_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx);

        auto tpz_xx_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx);

        auto tpx_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 1);

        auto tpy_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tpz_xx_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tpx_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 2);

        auto tpy_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tpz_xx_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tpx_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 3);

        auto tpy_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tpz_xx_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tpx_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 4);

        auto tpy_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tpz_xx_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tpx_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 5);

        auto tpy_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tpz_xx_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tpx_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 6);

        auto tpy_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tpz_xx_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tpx_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 7);

        auto tpy_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tpz_xx_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tpx_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 8);

        auto tpy_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tpz_xx_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tpx_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 9);

        auto tpy_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tpz_xx_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tpx_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 10);

        auto tpy_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tpz_xx_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tpx_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 11);

        auto tpy_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tpz_xx_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tpx_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 12);

        auto tpy_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tpz_xx_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tpx_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 13);

        auto tpy_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tpz_xx_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tpx_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 14);

        auto tpy_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tpz_xx_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 14);

        auto tpx_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 15);

        auto tpy_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tpz_xy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tpx_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 16);

        auto tpx_xxx_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx);

        auto tpy_xxx_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx);

        auto tpz_xxx_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx);

        auto tpx_xxx_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 1);

        auto tpy_xxx_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 1);

        auto tpz_xxx_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 1);

        auto tpx_xxx_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 2);

        auto tpy_xxx_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 2);

        auto tpz_xxx_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 2);

        auto tpx_xxx_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 3);

        auto tpy_xxx_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 3);

        auto tpz_xxx_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 3);

        auto tpx_xxx_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 4);

        auto tpy_xxx_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 4);

        auto tpz_xxx_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 4);

        auto tpx_xxx_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 5);

        auto tpy_xxx_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 5);

        auto tpz_xxx_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 5);

        auto tpx_xxx_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 6);

        auto tpy_xxx_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 6);

        auto tpz_xxx_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 6);

        auto tpx_xxx_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 7);

        auto tpy_xxx_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 7);

        auto tpz_xxx_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 7);

        auto tpx_xxx_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 8);

        auto tpy_xxx_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 8);

        auto tpz_xxx_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 8);

        auto tpx_xxx_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 9);

        auto tpy_xxx_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 9);

        auto tpz_xxx_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 9);

        auto tpx_xxy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 10);

        auto tpy_xxy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 10);

        auto tpz_xxy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 10);

        auto tpx_xxy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 11);

        auto ts_xxx_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx);

        auto ts_xxx_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 1);

        auto ts_xxx_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 2);

        auto ts_xxx_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 3);

        auto ts_xxx_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 4);

        auto ts_xxx_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 5);

        auto ts_xxx_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 6);

        auto ts_xxx_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 7);

        auto ts_xxx_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 8);

        auto ts_xxx_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 9);

        auto ts_xxx_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 10);

        auto ts_xxx_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 11);

        auto ts_xxx_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 12);

        auto ts_xxx_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 13);

        auto ts_xxx_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 14);

        auto ts_xxy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 15);

        auto ts_xxy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 16);

        // set up pointers to integrals

        auto tpx_xxxx_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx);

        auto tpy_xxxx_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx);

        auto tpz_xxxx_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx);

        auto tpx_xxxx_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 1);

        auto tpy_xxxx_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 1);

        auto tpz_xxxx_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 1);

        auto tpx_xxxx_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 2);

        auto tpy_xxxx_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 2);

        auto tpz_xxxx_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 2);

        auto tpx_xxxx_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 3);

        auto tpy_xxxx_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 3);

        auto tpz_xxxx_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 3);

        auto tpx_xxxx_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 4);

        auto tpy_xxxx_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 4);

        auto tpz_xxxx_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 4);

        auto tpx_xxxx_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 5);

        auto tpy_xxxx_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 5);

        auto tpz_xxxx_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 5);

        auto tpx_xxxx_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 6);

        auto tpy_xxxx_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 6);

        auto tpz_xxxx_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 6);

        auto tpx_xxxx_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 7);

        auto tpy_xxxx_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 7);

        auto tpz_xxxx_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 7);

        auto tpx_xxxx_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 8);

        auto tpy_xxxx_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 8);

        auto tpz_xxxx_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 8);

        auto tpx_xxxx_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 9);

        auto tpy_xxxx_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 9);

        auto tpz_xxxx_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 9);

        auto tpx_xxxx_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 10);

        auto tpy_xxxx_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 10);

        auto tpz_xxxx_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 10);

        auto tpx_xxxx_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 11);

        auto tpy_xxxx_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 11);

        auto tpz_xxxx_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 11);

        auto tpx_xxxx_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 12);

        auto tpy_xxxx_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 12);

        auto tpz_xxxx_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 12);

        auto tpx_xxxx_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 13);

        auto tpy_xxxx_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 13);

        auto tpz_xxxx_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 13);

        auto tpx_xxxx_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 14);

        auto tpy_xxxx_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 14);

        auto tpz_xxxx_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 14);

        auto tpx_xxxy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 15);

        auto tpy_xxxy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 15);

        auto tpz_xxxy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 15);

        auto tpx_xxxy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 16);

        // Batch of Integrals (0,49)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xx_xxxx_0, tpx_xx_xxxy_0, tpx_xx_xxxz_0, \
                                     tpx_xx_xxyy_0, tpx_xx_xxyz_0, tpx_xx_xxzz_0, tpx_xx_xyyy_0, tpx_xx_xyyz_0, \
                                     tpx_xx_xyzz_0, tpx_xx_xzzz_0, tpx_xx_yyyy_0, tpx_xx_yyyz_0, tpx_xx_yyzz_0, \
                                     tpx_xx_yzzz_0, tpx_xx_zzzz_0, tpx_xxx_xxx_0, tpx_xxx_xxxx_0, tpx_xxx_xxxy_0, \
                                     tpx_xxx_xxxz_0, tpx_xxx_xxy_0, tpx_xxx_xxyy_0, tpx_xxx_xxyz_0, tpx_xxx_xxz_0, \
                                     tpx_xxx_xxzz_0, tpx_xxx_xyy_0, tpx_xxx_xyyy_0, tpx_xxx_xyyz_0, tpx_xxx_xyz_0, \
                                     tpx_xxx_xyzz_0, tpx_xxx_xzz_0, tpx_xxx_xzzz_0, tpx_xxx_yyy_0, tpx_xxx_yyyy_0, \
                                     tpx_xxx_yyyz_0, tpx_xxx_yyz_0, tpx_xxx_yyzz_0, tpx_xxx_yzz_0, tpx_xxx_yzzz_0, \
                                     tpx_xxx_zzz_0, tpx_xxx_zzzz_0, tpx_xxxx_xxxx_0, tpx_xxxx_xxxy_0, tpx_xxxx_xxxz_0, \
                                     tpx_xxxx_xxyy_0, tpx_xxxx_xxyz_0, tpx_xxxx_xxzz_0, tpx_xxxx_xyyy_0, tpx_xxxx_xyyz_0, \
                                     tpx_xxxx_xyzz_0, tpx_xxxx_xzzz_0, tpx_xxxx_yyyy_0, tpx_xxxx_yyyz_0, tpx_xxxx_yyzz_0, \
                                     tpx_xxxx_yzzz_0, tpx_xxxx_zzzz_0, tpx_xxxy_xxxx_0, tpx_xxxy_xxxy_0, tpx_xxy_xxx_0, \
                                     tpx_xxy_xxxx_0, tpx_xxy_xxxy_0, tpx_xxy_xxy_0, tpx_xy_xxxx_0, tpx_xy_xxxy_0, \
                                     tpy_xx_xxxx_0, tpy_xx_xxxy_0, tpy_xx_xxxz_0, tpy_xx_xxyy_0, tpy_xx_xxyz_0, \
                                     tpy_xx_xxzz_0, tpy_xx_xyyy_0, tpy_xx_xyyz_0, tpy_xx_xyzz_0, tpy_xx_xzzz_0, \
                                     tpy_xx_yyyy_0, tpy_xx_yyyz_0, tpy_xx_yyzz_0, tpy_xx_yzzz_0, tpy_xx_zzzz_0, \
                                     tpy_xxx_xxx_0, tpy_xxx_xxxx_0, tpy_xxx_xxxy_0, tpy_xxx_xxxz_0, tpy_xxx_xxy_0, \
                                     tpy_xxx_xxyy_0, tpy_xxx_xxyz_0, tpy_xxx_xxz_0, tpy_xxx_xxzz_0, tpy_xxx_xyy_0, \
                                     tpy_xxx_xyyy_0, tpy_xxx_xyyz_0, tpy_xxx_xyz_0, tpy_xxx_xyzz_0, tpy_xxx_xzz_0, \
                                     tpy_xxx_xzzz_0, tpy_xxx_yyy_0, tpy_xxx_yyyy_0, tpy_xxx_yyyz_0, tpy_xxx_yyz_0, \
                                     tpy_xxx_yyzz_0, tpy_xxx_yzz_0, tpy_xxx_yzzz_0, tpy_xxx_zzz_0, tpy_xxx_zzzz_0, \
                                     tpy_xxxx_xxxx_0, tpy_xxxx_xxxy_0, tpy_xxxx_xxxz_0, tpy_xxxx_xxyy_0, tpy_xxxx_xxyz_0, \
                                     tpy_xxxx_xxzz_0, tpy_xxxx_xyyy_0, tpy_xxxx_xyyz_0, tpy_xxxx_xyzz_0, tpy_xxxx_xzzz_0, \
                                     tpy_xxxx_yyyy_0, tpy_xxxx_yyyz_0, tpy_xxxx_yyzz_0, tpy_xxxx_yzzz_0, tpy_xxxx_zzzz_0, \
                                     tpy_xxxy_xxxx_0, tpy_xxy_xxx_0, tpy_xxy_xxxx_0, tpy_xy_xxxx_0, tpz_xx_xxxx_0, \
                                     tpz_xx_xxxy_0, tpz_xx_xxxz_0, tpz_xx_xxyy_0, tpz_xx_xxyz_0, tpz_xx_xxzz_0, \
                                     tpz_xx_xyyy_0, tpz_xx_xyyz_0, tpz_xx_xyzz_0, tpz_xx_xzzz_0, tpz_xx_yyyy_0, \
                                     tpz_xx_yyyz_0, tpz_xx_yyzz_0, tpz_xx_yzzz_0, tpz_xx_zzzz_0, tpz_xxx_xxx_0, \
                                     tpz_xxx_xxxx_0, tpz_xxx_xxxy_0, tpz_xxx_xxxz_0, tpz_xxx_xxy_0, tpz_xxx_xxyy_0, \
                                     tpz_xxx_xxyz_0, tpz_xxx_xxz_0, tpz_xxx_xxzz_0, tpz_xxx_xyy_0, tpz_xxx_xyyy_0, \
                                     tpz_xxx_xyyz_0, tpz_xxx_xyz_0, tpz_xxx_xyzz_0, tpz_xxx_xzz_0, tpz_xxx_xzzz_0, \
                                     tpz_xxx_yyy_0, tpz_xxx_yyyy_0, tpz_xxx_yyyz_0, tpz_xxx_yyz_0, tpz_xxx_yyzz_0, \
                                     tpz_xxx_yzz_0, tpz_xxx_yzzz_0, tpz_xxx_zzz_0, tpz_xxx_zzzz_0, tpz_xxxx_xxxx_0, \
                                     tpz_xxxx_xxxy_0, tpz_xxxx_xxxz_0, tpz_xxxx_xxyy_0, tpz_xxxx_xxyz_0, tpz_xxxx_xxzz_0, \
                                     tpz_xxxx_xyyy_0, tpz_xxxx_xyyz_0, tpz_xxxx_xyzz_0, tpz_xxxx_xzzz_0, tpz_xxxx_yyyy_0, \
                                     tpz_xxxx_yyyz_0, tpz_xxxx_yyzz_0, tpz_xxxx_yzzz_0, tpz_xxxx_zzzz_0, tpz_xxxy_xxxx_0, \
                                     tpz_xxy_xxx_0, tpz_xxy_xxxx_0, tpz_xy_xxxx_0, ts_xxx_xxxx_0, ts_xxx_xxxy_0, \
                                     ts_xxx_xxxz_0, ts_xxx_xxyy_0, ts_xxx_xxyz_0, ts_xxx_xxzz_0, ts_xxx_xyyy_0, \
                                     ts_xxx_xyyz_0, ts_xxx_xyzz_0, ts_xxx_xzzz_0, ts_xxx_yyyy_0, ts_xxx_yyyz_0, \
                                     ts_xxx_yyzz_0, ts_xxx_yzzz_0, ts_xxx_zzzz_0, ts_xxy_xxxx_0, ts_xxy_xxxy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxxx_xxxx_0[j] =
                pa_x[j] * tpx_xxx_xxxx_0[j] + 1.5 * fl1_fx * tpx_xx_xxxx_0[j] + 2.0 * fl1_fx * tpx_xxx_xxx_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxxx_0[j];

            tpy_xxxx_xxxx_0[j] = pa_x[j] * tpy_xxx_xxxx_0[j] + 1.5 * fl1_fx * tpy_xx_xxxx_0[j] + 2.0 * fl1_fx * tpy_xxx_xxx_0[j];

            tpz_xxxx_xxxx_0[j] = pa_x[j] * tpz_xxx_xxxx_0[j] + 1.5 * fl1_fx * tpz_xx_xxxx_0[j] + 2.0 * fl1_fx * tpz_xxx_xxx_0[j];

            tpx_xxxx_xxxy_0[j] =
                pa_x[j] * tpx_xxx_xxxy_0[j] + 1.5 * fl1_fx * tpx_xx_xxxy_0[j] + 1.5 * fl1_fx * tpx_xxx_xxy_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxxy_0[j];

            tpy_xxxx_xxxy_0[j] = pa_x[j] * tpy_xxx_xxxy_0[j] + 1.5 * fl1_fx * tpy_xx_xxxy_0[j] + 1.5 * fl1_fx * tpy_xxx_xxy_0[j];

            tpz_xxxx_xxxy_0[j] = pa_x[j] * tpz_xxx_xxxy_0[j] + 1.5 * fl1_fx * tpz_xx_xxxy_0[j] + 1.5 * fl1_fx * tpz_xxx_xxy_0[j];

            tpx_xxxx_xxxz_0[j] =
                pa_x[j] * tpx_xxx_xxxz_0[j] + 1.5 * fl1_fx * tpx_xx_xxxz_0[j] + 1.5 * fl1_fx * tpx_xxx_xxz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxxz_0[j];

            tpy_xxxx_xxxz_0[j] = pa_x[j] * tpy_xxx_xxxz_0[j] + 1.5 * fl1_fx * tpy_xx_xxxz_0[j] + 1.5 * fl1_fx * tpy_xxx_xxz_0[j];

            tpz_xxxx_xxxz_0[j] = pa_x[j] * tpz_xxx_xxxz_0[j] + 1.5 * fl1_fx * tpz_xx_xxxz_0[j] + 1.5 * fl1_fx * tpz_xxx_xxz_0[j];

            tpx_xxxx_xxyy_0[j] =
                pa_x[j] * tpx_xxx_xxyy_0[j] + 1.5 * fl1_fx * tpx_xx_xxyy_0[j] + fl1_fx * tpx_xxx_xyy_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxyy_0[j];

            tpy_xxxx_xxyy_0[j] = pa_x[j] * tpy_xxx_xxyy_0[j] + 1.5 * fl1_fx * tpy_xx_xxyy_0[j] + fl1_fx * tpy_xxx_xyy_0[j];

            tpz_xxxx_xxyy_0[j] = pa_x[j] * tpz_xxx_xxyy_0[j] + 1.5 * fl1_fx * tpz_xx_xxyy_0[j] + fl1_fx * tpz_xxx_xyy_0[j];

            tpx_xxxx_xxyz_0[j] =
                pa_x[j] * tpx_xxx_xxyz_0[j] + 1.5 * fl1_fx * tpx_xx_xxyz_0[j] + fl1_fx * tpx_xxx_xyz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxyz_0[j];

            tpy_xxxx_xxyz_0[j] = pa_x[j] * tpy_xxx_xxyz_0[j] + 1.5 * fl1_fx * tpy_xx_xxyz_0[j] + fl1_fx * tpy_xxx_xyz_0[j];

            tpz_xxxx_xxyz_0[j] = pa_x[j] * tpz_xxx_xxyz_0[j] + 1.5 * fl1_fx * tpz_xx_xxyz_0[j] + fl1_fx * tpz_xxx_xyz_0[j];

            tpx_xxxx_xxzz_0[j] =
                pa_x[j] * tpx_xxx_xxzz_0[j] + 1.5 * fl1_fx * tpx_xx_xxzz_0[j] + fl1_fx * tpx_xxx_xzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxzz_0[j];

            tpy_xxxx_xxzz_0[j] = pa_x[j] * tpy_xxx_xxzz_0[j] + 1.5 * fl1_fx * tpy_xx_xxzz_0[j] + fl1_fx * tpy_xxx_xzz_0[j];

            tpz_xxxx_xxzz_0[j] = pa_x[j] * tpz_xxx_xxzz_0[j] + 1.5 * fl1_fx * tpz_xx_xxzz_0[j] + fl1_fx * tpz_xxx_xzz_0[j];

            tpx_xxxx_xyyy_0[j] =
                pa_x[j] * tpx_xxx_xyyy_0[j] + 1.5 * fl1_fx * tpx_xx_xyyy_0[j] + 0.5 * fl1_fx * tpx_xxx_yyy_0[j] - fl1_fgb * fl1_fx * ts_xxx_xyyy_0[j];

            tpy_xxxx_xyyy_0[j] = pa_x[j] * tpy_xxx_xyyy_0[j] + 1.5 * fl1_fx * tpy_xx_xyyy_0[j] + 0.5 * fl1_fx * tpy_xxx_yyy_0[j];

            tpz_xxxx_xyyy_0[j] = pa_x[j] * tpz_xxx_xyyy_0[j] + 1.5 * fl1_fx * tpz_xx_xyyy_0[j] + 0.5 * fl1_fx * tpz_xxx_yyy_0[j];

            tpx_xxxx_xyyz_0[j] =
                pa_x[j] * tpx_xxx_xyyz_0[j] + 1.5 * fl1_fx * tpx_xx_xyyz_0[j] + 0.5 * fl1_fx * tpx_xxx_yyz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xyyz_0[j];

            tpy_xxxx_xyyz_0[j] = pa_x[j] * tpy_xxx_xyyz_0[j] + 1.5 * fl1_fx * tpy_xx_xyyz_0[j] + 0.5 * fl1_fx * tpy_xxx_yyz_0[j];

            tpz_xxxx_xyyz_0[j] = pa_x[j] * tpz_xxx_xyyz_0[j] + 1.5 * fl1_fx * tpz_xx_xyyz_0[j] + 0.5 * fl1_fx * tpz_xxx_yyz_0[j];

            tpx_xxxx_xyzz_0[j] =
                pa_x[j] * tpx_xxx_xyzz_0[j] + 1.5 * fl1_fx * tpx_xx_xyzz_0[j] + 0.5 * fl1_fx * tpx_xxx_yzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xyzz_0[j];

            tpy_xxxx_xyzz_0[j] = pa_x[j] * tpy_xxx_xyzz_0[j] + 1.5 * fl1_fx * tpy_xx_xyzz_0[j] + 0.5 * fl1_fx * tpy_xxx_yzz_0[j];

            tpz_xxxx_xyzz_0[j] = pa_x[j] * tpz_xxx_xyzz_0[j] + 1.5 * fl1_fx * tpz_xx_xyzz_0[j] + 0.5 * fl1_fx * tpz_xxx_yzz_0[j];

            tpx_xxxx_xzzz_0[j] =
                pa_x[j] * tpx_xxx_xzzz_0[j] + 1.5 * fl1_fx * tpx_xx_xzzz_0[j] + 0.5 * fl1_fx * tpx_xxx_zzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xzzz_0[j];

            tpy_xxxx_xzzz_0[j] = pa_x[j] * tpy_xxx_xzzz_0[j] + 1.5 * fl1_fx * tpy_xx_xzzz_0[j] + 0.5 * fl1_fx * tpy_xxx_zzz_0[j];

            tpz_xxxx_xzzz_0[j] = pa_x[j] * tpz_xxx_xzzz_0[j] + 1.5 * fl1_fx * tpz_xx_xzzz_0[j] + 0.5 * fl1_fx * tpz_xxx_zzz_0[j];

            tpx_xxxx_yyyy_0[j] = pa_x[j] * tpx_xxx_yyyy_0[j] + 1.5 * fl1_fx * tpx_xx_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xxx_yyyy_0[j];

            tpy_xxxx_yyyy_0[j] = pa_x[j] * tpy_xxx_yyyy_0[j] + 1.5 * fl1_fx * tpy_xx_yyyy_0[j];

            tpz_xxxx_yyyy_0[j] = pa_x[j] * tpz_xxx_yyyy_0[j] + 1.5 * fl1_fx * tpz_xx_yyyy_0[j];

            tpx_xxxx_yyyz_0[j] = pa_x[j] * tpx_xxx_yyyz_0[j] + 1.5 * fl1_fx * tpx_xx_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xxx_yyyz_0[j];

            tpy_xxxx_yyyz_0[j] = pa_x[j] * tpy_xxx_yyyz_0[j] + 1.5 * fl1_fx * tpy_xx_yyyz_0[j];

            tpz_xxxx_yyyz_0[j] = pa_x[j] * tpz_xxx_yyyz_0[j] + 1.5 * fl1_fx * tpz_xx_yyyz_0[j];

            tpx_xxxx_yyzz_0[j] = pa_x[j] * tpx_xxx_yyzz_0[j] + 1.5 * fl1_fx * tpx_xx_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_yyzz_0[j];

            tpy_xxxx_yyzz_0[j] = pa_x[j] * tpy_xxx_yyzz_0[j] + 1.5 * fl1_fx * tpy_xx_yyzz_0[j];

            tpz_xxxx_yyzz_0[j] = pa_x[j] * tpz_xxx_yyzz_0[j] + 1.5 * fl1_fx * tpz_xx_yyzz_0[j];

            tpx_xxxx_yzzz_0[j] = pa_x[j] * tpx_xxx_yzzz_0[j] + 1.5 * fl1_fx * tpx_xx_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_yzzz_0[j];

            tpy_xxxx_yzzz_0[j] = pa_x[j] * tpy_xxx_yzzz_0[j] + 1.5 * fl1_fx * tpy_xx_yzzz_0[j];

            tpz_xxxx_yzzz_0[j] = pa_x[j] * tpz_xxx_yzzz_0[j] + 1.5 * fl1_fx * tpz_xx_yzzz_0[j];

            tpx_xxxx_zzzz_0[j] = pa_x[j] * tpx_xxx_zzzz_0[j] + 1.5 * fl1_fx * tpx_xx_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_zzzz_0[j];

            tpy_xxxx_zzzz_0[j] = pa_x[j] * tpy_xxx_zzzz_0[j] + 1.5 * fl1_fx * tpy_xx_zzzz_0[j];

            tpz_xxxx_zzzz_0[j] = pa_x[j] * tpz_xxx_zzzz_0[j] + 1.5 * fl1_fx * tpz_xx_zzzz_0[j];

            tpx_xxxy_xxxx_0[j] =
                pa_x[j] * tpx_xxy_xxxx_0[j] + fl1_fx * tpx_xy_xxxx_0[j] + 2.0 * fl1_fx * tpx_xxy_xxx_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxxx_0[j];

            tpy_xxxy_xxxx_0[j] = pa_x[j] * tpy_xxy_xxxx_0[j] + fl1_fx * tpy_xy_xxxx_0[j] + 2.0 * fl1_fx * tpy_xxy_xxx_0[j];

            tpz_xxxy_xxxx_0[j] = pa_x[j] * tpz_xxy_xxxx_0[j] + fl1_fx * tpz_xy_xxxx_0[j] + 2.0 * fl1_fx * tpz_xxy_xxx_0[j];

            tpx_xxxy_xxxy_0[j] =
                pa_x[j] * tpx_xxy_xxxy_0[j] + fl1_fx * tpx_xy_xxxy_0[j] + 1.5 * fl1_fx * tpx_xxy_xxy_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxxy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_49_98(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (49,98)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpy_xxy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 16);

        auto tpz_xxy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 16);

        auto tpx_xxy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 17);

        auto tpy_xxy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 17);

        auto tpz_xxy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 17);

        auto tpx_xxy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 18);

        auto tpy_xxy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 18);

        auto tpz_xxy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 18);

        auto tpx_xxy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 19);

        auto tpy_xxy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 19);

        auto tpz_xxy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 19);

        auto tpx_xxy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 20);

        auto tpy_xxy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 20);

        auto tpz_xxy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 20);

        auto tpx_xxy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 21);

        auto tpy_xxy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 21);

        auto tpz_xxy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 21);

        auto tpx_xxy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 22);

        auto tpy_xxy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 22);

        auto tpz_xxy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 22);

        auto tpx_xxy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 23);

        auto tpy_xxy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 23);

        auto tpz_xxy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 23);

        auto tpx_xxy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 24);

        auto tpy_xxy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 24);

        auto tpz_xxy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 24);

        auto tpx_xxy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 25);

        auto tpy_xxy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 25);

        auto tpz_xxy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 25);

        auto tpx_xxy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 26);

        auto tpy_xxy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 26);

        auto tpz_xxy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 26);

        auto tpx_xxy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 27);

        auto tpy_xxy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 27);

        auto tpz_xxy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 27);

        auto tpx_xxy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 28);

        auto tpy_xxy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 28);

        auto tpz_xxy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 28);

        auto tpx_xxy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 29);

        auto tpy_xxy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 29);

        auto tpz_xxy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 29);

        auto tpx_xxz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 30);

        auto tpy_xxz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 30);

        auto tpz_xxz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 30);

        auto tpx_xxz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 31);

        auto tpy_xxz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 31);

        auto tpz_xxz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 31);

        auto tpx_xxz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 32);

        auto tpy_xxz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 32);

        auto tpy_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tpz_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 16);

        auto tpx_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 17);

        auto tpy_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tpz_xy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tpx_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 18);

        auto tpy_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tpz_xy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tpx_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 19);

        auto tpy_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tpz_xy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tpx_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 20);

        auto tpy_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tpz_xy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tpx_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 21);

        auto tpy_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tpz_xy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tpx_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 22);

        auto tpy_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tpz_xy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tpx_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 23);

        auto tpy_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tpz_xy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tpx_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 24);

        auto tpy_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tpz_xy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tpx_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 25);

        auto tpy_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tpz_xy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tpx_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 26);

        auto tpy_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tpz_xy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tpx_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 27);

        auto tpy_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tpz_xy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tpx_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 28);

        auto tpy_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tpz_xy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tpx_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 29);

        auto tpy_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tpz_xy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 29);

        auto tpx_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 30);

        auto tpy_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tpz_xz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tpx_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 31);

        auto tpy_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tpz_xz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tpx_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 32);

        auto tpy_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tpy_xxy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 11);

        auto tpz_xxy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 11);

        auto tpx_xxy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 12);

        auto tpy_xxy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 12);

        auto tpz_xxy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 12);

        auto tpx_xxy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 13);

        auto tpy_xxy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 13);

        auto tpz_xxy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 13);

        auto tpx_xxy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 14);

        auto tpy_xxy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 14);

        auto tpz_xxy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 14);

        auto tpx_xxy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 15);

        auto tpy_xxy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 15);

        auto tpz_xxy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 15);

        auto tpx_xxy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 16);

        auto tpy_xxy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 16);

        auto tpz_xxy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 16);

        auto tpx_xxy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 17);

        auto tpy_xxy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 17);

        auto tpz_xxy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 17);

        auto tpx_xxy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 18);

        auto tpy_xxy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 18);

        auto tpz_xxy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 18);

        auto tpx_xxy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 19);

        auto tpy_xxy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 19);

        auto tpz_xxy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 19);

        auto tpx_xxz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 20);

        auto tpy_xxz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 20);

        auto tpz_xxz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 20);

        auto tpx_xxz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 21);

        auto tpy_xxz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 21);

        auto tpz_xxz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 21);

        auto tpx_xxz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 22);

        auto tpy_xxz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 22);

        auto ts_xxy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 17);

        auto ts_xxy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 18);

        auto ts_xxy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 19);

        auto ts_xxy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 20);

        auto ts_xxy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 21);

        auto ts_xxy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 22);

        auto ts_xxy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 23);

        auto ts_xxy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 24);

        auto ts_xxy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 25);

        auto ts_xxy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 26);

        auto ts_xxy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 27);

        auto ts_xxy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 28);

        auto ts_xxy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 29);

        auto ts_xxz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 30);

        auto ts_xxz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 31);

        auto ts_xxz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 32);

        // set up pointers to integrals

        auto tpy_xxxy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 16);

        auto tpz_xxxy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 16);

        auto tpx_xxxy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 17);

        auto tpy_xxxy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 17);

        auto tpz_xxxy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 17);

        auto tpx_xxxy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 18);

        auto tpy_xxxy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 18);

        auto tpz_xxxy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 18);

        auto tpx_xxxy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 19);

        auto tpy_xxxy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 19);

        auto tpz_xxxy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 19);

        auto tpx_xxxy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 20);

        auto tpy_xxxy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 20);

        auto tpz_xxxy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 20);

        auto tpx_xxxy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 21);

        auto tpy_xxxy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 21);

        auto tpz_xxxy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 21);

        auto tpx_xxxy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 22);

        auto tpy_xxxy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 22);

        auto tpz_xxxy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 22);

        auto tpx_xxxy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 23);

        auto tpy_xxxy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 23);

        auto tpz_xxxy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 23);

        auto tpx_xxxy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 24);

        auto tpy_xxxy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 24);

        auto tpz_xxxy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 24);

        auto tpx_xxxy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 25);

        auto tpy_xxxy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 25);

        auto tpz_xxxy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 25);

        auto tpx_xxxy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 26);

        auto tpy_xxxy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 26);

        auto tpz_xxxy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 26);

        auto tpx_xxxy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 27);

        auto tpy_xxxy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 27);

        auto tpz_xxxy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 27);

        auto tpx_xxxy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 28);

        auto tpy_xxxy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 28);

        auto tpz_xxxy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 28);

        auto tpx_xxxy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 29);

        auto tpy_xxxy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 29);

        auto tpz_xxxy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 29);

        auto tpx_xxxz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 30);

        auto tpy_xxxz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 30);

        auto tpz_xxxz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 30);

        auto tpx_xxxz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 31);

        auto tpy_xxxz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 31);

        auto tpz_xxxz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 31);

        auto tpx_xxxz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 32);

        auto tpy_xxxz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 32);

        // Batch of Integrals (49,98)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxxy_xxxz_0, tpx_xxxy_xxyy_0, tpx_xxxy_xxyz_0, \
                                     tpx_xxxy_xxzz_0, tpx_xxxy_xyyy_0, tpx_xxxy_xyyz_0, tpx_xxxy_xyzz_0, tpx_xxxy_xzzz_0, \
                                     tpx_xxxy_yyyy_0, tpx_xxxy_yyyz_0, tpx_xxxy_yyzz_0, tpx_xxxy_yzzz_0, tpx_xxxy_zzzz_0, \
                                     tpx_xxxz_xxxx_0, tpx_xxxz_xxxy_0, tpx_xxxz_xxxz_0, tpx_xxy_xxxz_0, tpx_xxy_xxyy_0, \
                                     tpx_xxy_xxyz_0, tpx_xxy_xxz_0, tpx_xxy_xxzz_0, tpx_xxy_xyy_0, tpx_xxy_xyyy_0, \
                                     tpx_xxy_xyyz_0, tpx_xxy_xyz_0, tpx_xxy_xyzz_0, tpx_xxy_xzz_0, tpx_xxy_xzzz_0, \
                                     tpx_xxy_yyy_0, tpx_xxy_yyyy_0, tpx_xxy_yyyz_0, tpx_xxy_yyz_0, tpx_xxy_yyzz_0, \
                                     tpx_xxy_yzz_0, tpx_xxy_yzzz_0, tpx_xxy_zzz_0, tpx_xxy_zzzz_0, tpx_xxz_xxx_0, \
                                     tpx_xxz_xxxx_0, tpx_xxz_xxxy_0, tpx_xxz_xxxz_0, tpx_xxz_xxy_0, tpx_xxz_xxz_0, \
                                     tpx_xy_xxxz_0, tpx_xy_xxyy_0, tpx_xy_xxyz_0, tpx_xy_xxzz_0, tpx_xy_xyyy_0, \
                                     tpx_xy_xyyz_0, tpx_xy_xyzz_0, tpx_xy_xzzz_0, tpx_xy_yyyy_0, tpx_xy_yyyz_0, \
                                     tpx_xy_yyzz_0, tpx_xy_yzzz_0, tpx_xy_zzzz_0, tpx_xz_xxxx_0, tpx_xz_xxxy_0, \
                                     tpx_xz_xxxz_0, tpy_xxxy_xxxy_0, tpy_xxxy_xxxz_0, tpy_xxxy_xxyy_0, tpy_xxxy_xxyz_0, \
                                     tpy_xxxy_xxzz_0, tpy_xxxy_xyyy_0, tpy_xxxy_xyyz_0, tpy_xxxy_xyzz_0, tpy_xxxy_xzzz_0, \
                                     tpy_xxxy_yyyy_0, tpy_xxxy_yyyz_0, tpy_xxxy_yyzz_0, tpy_xxxy_yzzz_0, tpy_xxxy_zzzz_0, \
                                     tpy_xxxz_xxxx_0, tpy_xxxz_xxxy_0, tpy_xxxz_xxxz_0, tpy_xxy_xxxy_0, tpy_xxy_xxxz_0, \
                                     tpy_xxy_xxy_0, tpy_xxy_xxyy_0, tpy_xxy_xxyz_0, tpy_xxy_xxz_0, tpy_xxy_xxzz_0, \
                                     tpy_xxy_xyy_0, tpy_xxy_xyyy_0, tpy_xxy_xyyz_0, tpy_xxy_xyz_0, tpy_xxy_xyzz_0, \
                                     tpy_xxy_xzz_0, tpy_xxy_xzzz_0, tpy_xxy_yyy_0, tpy_xxy_yyyy_0, tpy_xxy_yyyz_0, \
                                     tpy_xxy_yyz_0, tpy_xxy_yyzz_0, tpy_xxy_yzz_0, tpy_xxy_yzzz_0, tpy_xxy_zzz_0, \
                                     tpy_xxy_zzzz_0, tpy_xxz_xxx_0, tpy_xxz_xxxx_0, tpy_xxz_xxxy_0, tpy_xxz_xxxz_0, \
                                     tpy_xxz_xxy_0, tpy_xxz_xxz_0, tpy_xy_xxxy_0, tpy_xy_xxxz_0, tpy_xy_xxyy_0, \
                                     tpy_xy_xxyz_0, tpy_xy_xxzz_0, tpy_xy_xyyy_0, tpy_xy_xyyz_0, tpy_xy_xyzz_0, \
                                     tpy_xy_xzzz_0, tpy_xy_yyyy_0, tpy_xy_yyyz_0, tpy_xy_yyzz_0, tpy_xy_yzzz_0, \
                                     tpy_xy_zzzz_0, tpy_xz_xxxx_0, tpy_xz_xxxy_0, tpy_xz_xxxz_0, tpz_xxxy_xxxy_0, \
                                     tpz_xxxy_xxxz_0, tpz_xxxy_xxyy_0, tpz_xxxy_xxyz_0, tpz_xxxy_xxzz_0, tpz_xxxy_xyyy_0, \
                                     tpz_xxxy_xyyz_0, tpz_xxxy_xyzz_0, tpz_xxxy_xzzz_0, tpz_xxxy_yyyy_0, tpz_xxxy_yyyz_0, \
                                     tpz_xxxy_yyzz_0, tpz_xxxy_yzzz_0, tpz_xxxy_zzzz_0, tpz_xxxz_xxxx_0, tpz_xxxz_xxxy_0, \
                                     tpz_xxy_xxxy_0, tpz_xxy_xxxz_0, tpz_xxy_xxy_0, tpz_xxy_xxyy_0, tpz_xxy_xxyz_0, \
                                     tpz_xxy_xxz_0, tpz_xxy_xxzz_0, tpz_xxy_xyy_0, tpz_xxy_xyyy_0, tpz_xxy_xyyz_0, \
                                     tpz_xxy_xyz_0, tpz_xxy_xyzz_0, tpz_xxy_xzz_0, tpz_xxy_xzzz_0, tpz_xxy_yyy_0, \
                                     tpz_xxy_yyyy_0, tpz_xxy_yyyz_0, tpz_xxy_yyz_0, tpz_xxy_yyzz_0, tpz_xxy_yzz_0, \
                                     tpz_xxy_yzzz_0, tpz_xxy_zzz_0, tpz_xxy_zzzz_0, tpz_xxz_xxx_0, tpz_xxz_xxxx_0, \
                                     tpz_xxz_xxxy_0, tpz_xxz_xxy_0, tpz_xy_xxxy_0, tpz_xy_xxxz_0, tpz_xy_xxyy_0, \
                                     tpz_xy_xxyz_0, tpz_xy_xxzz_0, tpz_xy_xyyy_0, tpz_xy_xyyz_0, tpz_xy_xyzz_0, \
                                     tpz_xy_xzzz_0, tpz_xy_yyyy_0, tpz_xy_yyyz_0, tpz_xy_yyzz_0, tpz_xy_yzzz_0, \
                                     tpz_xy_zzzz_0, tpz_xz_xxxx_0, tpz_xz_xxxy_0, ts_xxy_xxxz_0, ts_xxy_xxyy_0, \
                                     ts_xxy_xxyz_0, ts_xxy_xxzz_0, ts_xxy_xyyy_0, ts_xxy_xyyz_0, ts_xxy_xyzz_0, \
                                     ts_xxy_xzzz_0, ts_xxy_yyyy_0, ts_xxy_yyyz_0, ts_xxy_yyzz_0, ts_xxy_yzzz_0, \
                                     ts_xxy_zzzz_0, ts_xxz_xxxx_0, ts_xxz_xxxy_0, ts_xxz_xxxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_xxxy_xxxy_0[j] = pa_x[j] * tpy_xxy_xxxy_0[j] + fl1_fx * tpy_xy_xxxy_0[j] + 1.5 * fl1_fx * tpy_xxy_xxy_0[j];

            tpz_xxxy_xxxy_0[j] = pa_x[j] * tpz_xxy_xxxy_0[j] + fl1_fx * tpz_xy_xxxy_0[j] + 1.5 * fl1_fx * tpz_xxy_xxy_0[j];

            tpx_xxxy_xxxz_0[j] =
                pa_x[j] * tpx_xxy_xxxz_0[j] + fl1_fx * tpx_xy_xxxz_0[j] + 1.5 * fl1_fx * tpx_xxy_xxz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxxz_0[j];

            tpy_xxxy_xxxz_0[j] = pa_x[j] * tpy_xxy_xxxz_0[j] + fl1_fx * tpy_xy_xxxz_0[j] + 1.5 * fl1_fx * tpy_xxy_xxz_0[j];

            tpz_xxxy_xxxz_0[j] = pa_x[j] * tpz_xxy_xxxz_0[j] + fl1_fx * tpz_xy_xxxz_0[j] + 1.5 * fl1_fx * tpz_xxy_xxz_0[j];

            tpx_xxxy_xxyy_0[j] =
                pa_x[j] * tpx_xxy_xxyy_0[j] + fl1_fx * tpx_xy_xxyy_0[j] + fl1_fx * tpx_xxy_xyy_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxyy_0[j];

            tpy_xxxy_xxyy_0[j] = pa_x[j] * tpy_xxy_xxyy_0[j] + fl1_fx * tpy_xy_xxyy_0[j] + fl1_fx * tpy_xxy_xyy_0[j];

            tpz_xxxy_xxyy_0[j] = pa_x[j] * tpz_xxy_xxyy_0[j] + fl1_fx * tpz_xy_xxyy_0[j] + fl1_fx * tpz_xxy_xyy_0[j];

            tpx_xxxy_xxyz_0[j] =
                pa_x[j] * tpx_xxy_xxyz_0[j] + fl1_fx * tpx_xy_xxyz_0[j] + fl1_fx * tpx_xxy_xyz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxyz_0[j];

            tpy_xxxy_xxyz_0[j] = pa_x[j] * tpy_xxy_xxyz_0[j] + fl1_fx * tpy_xy_xxyz_0[j] + fl1_fx * tpy_xxy_xyz_0[j];

            tpz_xxxy_xxyz_0[j] = pa_x[j] * tpz_xxy_xxyz_0[j] + fl1_fx * tpz_xy_xxyz_0[j] + fl1_fx * tpz_xxy_xyz_0[j];

            tpx_xxxy_xxzz_0[j] =
                pa_x[j] * tpx_xxy_xxzz_0[j] + fl1_fx * tpx_xy_xxzz_0[j] + fl1_fx * tpx_xxy_xzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxzz_0[j];

            tpy_xxxy_xxzz_0[j] = pa_x[j] * tpy_xxy_xxzz_0[j] + fl1_fx * tpy_xy_xxzz_0[j] + fl1_fx * tpy_xxy_xzz_0[j];

            tpz_xxxy_xxzz_0[j] = pa_x[j] * tpz_xxy_xxzz_0[j] + fl1_fx * tpz_xy_xxzz_0[j] + fl1_fx * tpz_xxy_xzz_0[j];

            tpx_xxxy_xyyy_0[j] =
                pa_x[j] * tpx_xxy_xyyy_0[j] + fl1_fx * tpx_xy_xyyy_0[j] + 0.5 * fl1_fx * tpx_xxy_yyy_0[j] - fl1_fgb * fl1_fx * ts_xxy_xyyy_0[j];

            tpy_xxxy_xyyy_0[j] = pa_x[j] * tpy_xxy_xyyy_0[j] + fl1_fx * tpy_xy_xyyy_0[j] + 0.5 * fl1_fx * tpy_xxy_yyy_0[j];

            tpz_xxxy_xyyy_0[j] = pa_x[j] * tpz_xxy_xyyy_0[j] + fl1_fx * tpz_xy_xyyy_0[j] + 0.5 * fl1_fx * tpz_xxy_yyy_0[j];

            tpx_xxxy_xyyz_0[j] =
                pa_x[j] * tpx_xxy_xyyz_0[j] + fl1_fx * tpx_xy_xyyz_0[j] + 0.5 * fl1_fx * tpx_xxy_yyz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xyyz_0[j];

            tpy_xxxy_xyyz_0[j] = pa_x[j] * tpy_xxy_xyyz_0[j] + fl1_fx * tpy_xy_xyyz_0[j] + 0.5 * fl1_fx * tpy_xxy_yyz_0[j];

            tpz_xxxy_xyyz_0[j] = pa_x[j] * tpz_xxy_xyyz_0[j] + fl1_fx * tpz_xy_xyyz_0[j] + 0.5 * fl1_fx * tpz_xxy_yyz_0[j];

            tpx_xxxy_xyzz_0[j] =
                pa_x[j] * tpx_xxy_xyzz_0[j] + fl1_fx * tpx_xy_xyzz_0[j] + 0.5 * fl1_fx * tpx_xxy_yzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xyzz_0[j];

            tpy_xxxy_xyzz_0[j] = pa_x[j] * tpy_xxy_xyzz_0[j] + fl1_fx * tpy_xy_xyzz_0[j] + 0.5 * fl1_fx * tpy_xxy_yzz_0[j];

            tpz_xxxy_xyzz_0[j] = pa_x[j] * tpz_xxy_xyzz_0[j] + fl1_fx * tpz_xy_xyzz_0[j] + 0.5 * fl1_fx * tpz_xxy_yzz_0[j];

            tpx_xxxy_xzzz_0[j] =
                pa_x[j] * tpx_xxy_xzzz_0[j] + fl1_fx * tpx_xy_xzzz_0[j] + 0.5 * fl1_fx * tpx_xxy_zzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xzzz_0[j];

            tpy_xxxy_xzzz_0[j] = pa_x[j] * tpy_xxy_xzzz_0[j] + fl1_fx * tpy_xy_xzzz_0[j] + 0.5 * fl1_fx * tpy_xxy_zzz_0[j];

            tpz_xxxy_xzzz_0[j] = pa_x[j] * tpz_xxy_xzzz_0[j] + fl1_fx * tpz_xy_xzzz_0[j] + 0.5 * fl1_fx * tpz_xxy_zzz_0[j];

            tpx_xxxy_yyyy_0[j] = pa_x[j] * tpx_xxy_yyyy_0[j] + fl1_fx * tpx_xy_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xxy_yyyy_0[j];

            tpy_xxxy_yyyy_0[j] = pa_x[j] * tpy_xxy_yyyy_0[j] + fl1_fx * tpy_xy_yyyy_0[j];

            tpz_xxxy_yyyy_0[j] = pa_x[j] * tpz_xxy_yyyy_0[j] + fl1_fx * tpz_xy_yyyy_0[j];

            tpx_xxxy_yyyz_0[j] = pa_x[j] * tpx_xxy_yyyz_0[j] + fl1_fx * tpx_xy_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xxy_yyyz_0[j];

            tpy_xxxy_yyyz_0[j] = pa_x[j] * tpy_xxy_yyyz_0[j] + fl1_fx * tpy_xy_yyyz_0[j];

            tpz_xxxy_yyyz_0[j] = pa_x[j] * tpz_xxy_yyyz_0[j] + fl1_fx * tpz_xy_yyyz_0[j];

            tpx_xxxy_yyzz_0[j] = pa_x[j] * tpx_xxy_yyzz_0[j] + fl1_fx * tpx_xy_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_yyzz_0[j];

            tpy_xxxy_yyzz_0[j] = pa_x[j] * tpy_xxy_yyzz_0[j] + fl1_fx * tpy_xy_yyzz_0[j];

            tpz_xxxy_yyzz_0[j] = pa_x[j] * tpz_xxy_yyzz_0[j] + fl1_fx * tpz_xy_yyzz_0[j];

            tpx_xxxy_yzzz_0[j] = pa_x[j] * tpx_xxy_yzzz_0[j] + fl1_fx * tpx_xy_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_yzzz_0[j];

            tpy_xxxy_yzzz_0[j] = pa_x[j] * tpy_xxy_yzzz_0[j] + fl1_fx * tpy_xy_yzzz_0[j];

            tpz_xxxy_yzzz_0[j] = pa_x[j] * tpz_xxy_yzzz_0[j] + fl1_fx * tpz_xy_yzzz_0[j];

            tpx_xxxy_zzzz_0[j] = pa_x[j] * tpx_xxy_zzzz_0[j] + fl1_fx * tpx_xy_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_zzzz_0[j];

            tpy_xxxy_zzzz_0[j] = pa_x[j] * tpy_xxy_zzzz_0[j] + fl1_fx * tpy_xy_zzzz_0[j];

            tpz_xxxy_zzzz_0[j] = pa_x[j] * tpz_xxy_zzzz_0[j] + fl1_fx * tpz_xy_zzzz_0[j];

            tpx_xxxz_xxxx_0[j] =
                pa_x[j] * tpx_xxz_xxxx_0[j] + fl1_fx * tpx_xz_xxxx_0[j] + 2.0 * fl1_fx * tpx_xxz_xxx_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxxx_0[j];

            tpy_xxxz_xxxx_0[j] = pa_x[j] * tpy_xxz_xxxx_0[j] + fl1_fx * tpy_xz_xxxx_0[j] + 2.0 * fl1_fx * tpy_xxz_xxx_0[j];

            tpz_xxxz_xxxx_0[j] = pa_x[j] * tpz_xxz_xxxx_0[j] + fl1_fx * tpz_xz_xxxx_0[j] + 2.0 * fl1_fx * tpz_xxz_xxx_0[j];

            tpx_xxxz_xxxy_0[j] =
                pa_x[j] * tpx_xxz_xxxy_0[j] + fl1_fx * tpx_xz_xxxy_0[j] + 1.5 * fl1_fx * tpx_xxz_xxy_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxxy_0[j];

            tpy_xxxz_xxxy_0[j] = pa_x[j] * tpy_xxz_xxxy_0[j] + fl1_fx * tpy_xz_xxxy_0[j] + 1.5 * fl1_fx * tpy_xxz_xxy_0[j];

            tpz_xxxz_xxxy_0[j] = pa_x[j] * tpz_xxz_xxxy_0[j] + fl1_fx * tpz_xz_xxxy_0[j] + 1.5 * fl1_fx * tpz_xxz_xxy_0[j];

            tpx_xxxz_xxxz_0[j] =
                pa_x[j] * tpx_xxz_xxxz_0[j] + fl1_fx * tpx_xz_xxxz_0[j] + 1.5 * fl1_fx * tpx_xxz_xxz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxxz_0[j];

            tpy_xxxz_xxxz_0[j] = pa_x[j] * tpy_xxz_xxxz_0[j] + fl1_fx * tpy_xz_xxxz_0[j] + 1.5 * fl1_fx * tpy_xxz_xxz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_98_147(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (98,147)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpz_xxz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tpx_xxz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 33);

        auto tpy_xxz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 33);

        auto tpz_xxz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 33);

        auto tpx_xxz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 34);

        auto tpy_xxz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 34);

        auto tpz_xxz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 34);

        auto tpx_xxz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 35);

        auto tpy_xxz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 35);

        auto tpz_xxz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 35);

        auto tpx_xxz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 36);

        auto tpy_xxz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 36);

        auto tpz_xxz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 36);

        auto tpx_xxz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 37);

        auto tpy_xxz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 37);

        auto tpz_xxz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 37);

        auto tpx_xxz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 38);

        auto tpy_xxz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 38);

        auto tpz_xxz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 38);

        auto tpx_xxz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 39);

        auto tpy_xxz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 39);

        auto tpz_xxz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 39);

        auto tpx_xxz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 40);

        auto tpy_xxz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 40);

        auto tpz_xxz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 40);

        auto tpx_xxz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 41);

        auto tpy_xxz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 41);

        auto tpz_xxz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 41);

        auto tpx_xxz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 42);

        auto tpy_xxz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 42);

        auto tpz_xxz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 42);

        auto tpx_xxz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 43);

        auto tpy_xxz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 43);

        auto tpz_xxz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 43);

        auto tpx_xxz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 44);

        auto tpy_xxz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 44);

        auto tpz_xxz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 44);

        auto tpx_xyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 45);

        auto tpy_xyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 45);

        auto tpz_xyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 45);

        auto tpx_xyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 46);

        auto tpy_xyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 46);

        auto tpz_xyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 46);

        auto tpx_xyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 47);

        auto tpy_xyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 47);

        auto tpz_xyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 47);

        auto tpx_xyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 48);

        auto tpy_xyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 48);

        auto tpz_xyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 48);

        auto tpz_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tpx_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 33);

        auto tpy_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tpz_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tpx_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 34);

        auto tpy_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tpz_xz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tpx_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 35);

        auto tpy_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tpz_xz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tpx_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 36);

        auto tpy_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tpz_xz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tpx_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 37);

        auto tpy_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tpz_xz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tpx_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 38);

        auto tpy_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tpz_xz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tpx_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 39);

        auto tpy_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tpz_xz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tpx_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 40);

        auto tpy_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tpz_xz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tpx_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 41);

        auto tpy_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tpz_xz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tpx_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 42);

        auto tpy_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tpz_xz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tpx_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 43);

        auto tpy_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tpz_xz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tpx_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 44);

        auto tpy_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tpz_xz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 44);

        auto tpx_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 45);

        auto tpy_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tpz_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tpx_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 46);

        auto tpy_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tpz_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tpx_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 47);

        auto tpy_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tpz_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tpx_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 48);

        auto tpy_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tpz_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tpz_xxz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 22);

        auto tpx_xxz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 23);

        auto tpy_xxz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 23);

        auto tpz_xxz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 23);

        auto tpx_xxz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 24);

        auto tpy_xxz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 24);

        auto tpz_xxz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 24);

        auto tpx_xxz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 25);

        auto tpy_xxz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 25);

        auto tpz_xxz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 25);

        auto tpx_xxz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 26);

        auto tpy_xxz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 26);

        auto tpz_xxz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 26);

        auto tpx_xxz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 27);

        auto tpy_xxz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 27);

        auto tpz_xxz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 27);

        auto tpx_xxz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 28);

        auto tpy_xxz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 28);

        auto tpz_xxz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 28);

        auto tpx_xxz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 29);

        auto tpy_xxz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 29);

        auto tpz_xxz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 29);

        auto tpx_xyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 30);

        auto tpy_xyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 30);

        auto tpz_xyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 30);

        auto tpx_xyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 31);

        auto tpy_xyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 31);

        auto tpz_xyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 31);

        auto tpx_xyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 32);

        auto tpy_xyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 32);

        auto tpz_xyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 32);

        auto tpx_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 33);

        auto tpy_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tpz_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 33);

        auto ts_xxz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 33);

        auto ts_xxz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 34);

        auto ts_xxz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 35);

        auto ts_xxz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 36);

        auto ts_xxz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 37);

        auto ts_xxz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 38);

        auto ts_xxz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 39);

        auto ts_xxz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 40);

        auto ts_xxz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 41);

        auto ts_xxz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 42);

        auto ts_xxz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 43);

        auto ts_xxz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 44);

        auto ts_xyy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 45);

        auto ts_xyy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 46);

        auto ts_xyy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 47);

        auto ts_xyy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 48);

        // set up pointers to integrals

        auto tpz_xxxz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 32);

        auto tpx_xxxz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 33);

        auto tpy_xxxz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 33);

        auto tpz_xxxz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 33);

        auto tpx_xxxz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 34);

        auto tpy_xxxz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 34);

        auto tpz_xxxz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 34);

        auto tpx_xxxz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 35);

        auto tpy_xxxz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 35);

        auto tpz_xxxz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 35);

        auto tpx_xxxz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 36);

        auto tpy_xxxz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 36);

        auto tpz_xxxz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 36);

        auto tpx_xxxz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 37);

        auto tpy_xxxz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 37);

        auto tpz_xxxz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 37);

        auto tpx_xxxz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 38);

        auto tpy_xxxz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 38);

        auto tpz_xxxz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 38);

        auto tpx_xxxz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 39);

        auto tpy_xxxz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 39);

        auto tpz_xxxz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 39);

        auto tpx_xxxz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 40);

        auto tpy_xxxz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 40);

        auto tpz_xxxz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 40);

        auto tpx_xxxz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 41);

        auto tpy_xxxz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 41);

        auto tpz_xxxz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 41);

        auto tpx_xxxz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 42);

        auto tpy_xxxz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 42);

        auto tpz_xxxz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 42);

        auto tpx_xxxz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 43);

        auto tpy_xxxz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 43);

        auto tpz_xxxz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 43);

        auto tpx_xxxz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 44);

        auto tpy_xxxz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 44);

        auto tpz_xxxz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 44);

        auto tpx_xxyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 45);

        auto tpy_xxyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 45);

        auto tpz_xxyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 45);

        auto tpx_xxyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 46);

        auto tpy_xxyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 46);

        auto tpz_xxyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 46);

        auto tpx_xxyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 47);

        auto tpy_xxyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 47);

        auto tpz_xxyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 47);

        auto tpx_xxyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 48);

        auto tpy_xxyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 48);

        auto tpz_xxyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 48);

        // Batch of Integrals (98,147)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxxz_xxyy_0, tpx_xxxz_xxyz_0, tpx_xxxz_xxzz_0, \
                                     tpx_xxxz_xyyy_0, tpx_xxxz_xyyz_0, tpx_xxxz_xyzz_0, tpx_xxxz_xzzz_0, tpx_xxxz_yyyy_0, \
                                     tpx_xxxz_yyyz_0, tpx_xxxz_yyzz_0, tpx_xxxz_yzzz_0, tpx_xxxz_zzzz_0, tpx_xxyy_xxxx_0, \
                                     tpx_xxyy_xxxy_0, tpx_xxyy_xxxz_0, tpx_xxyy_xxyy_0, tpx_xxz_xxyy_0, tpx_xxz_xxyz_0, \
                                     tpx_xxz_xxzz_0, tpx_xxz_xyy_0, tpx_xxz_xyyy_0, tpx_xxz_xyyz_0, tpx_xxz_xyz_0, \
                                     tpx_xxz_xyzz_0, tpx_xxz_xzz_0, tpx_xxz_xzzz_0, tpx_xxz_yyy_0, tpx_xxz_yyyy_0, \
                                     tpx_xxz_yyyz_0, tpx_xxz_yyz_0, tpx_xxz_yyzz_0, tpx_xxz_yzz_0, tpx_xxz_yzzz_0, \
                                     tpx_xxz_zzz_0, tpx_xxz_zzzz_0, tpx_xyy_xxx_0, tpx_xyy_xxxx_0, tpx_xyy_xxxy_0, \
                                     tpx_xyy_xxxz_0, tpx_xyy_xxy_0, tpx_xyy_xxyy_0, tpx_xyy_xxz_0, tpx_xyy_xyy_0, \
                                     tpx_xz_xxyy_0, tpx_xz_xxyz_0, tpx_xz_xxzz_0, tpx_xz_xyyy_0, tpx_xz_xyyz_0, \
                                     tpx_xz_xyzz_0, tpx_xz_xzzz_0, tpx_xz_yyyy_0, tpx_xz_yyyz_0, tpx_xz_yyzz_0, \
                                     tpx_xz_yzzz_0, tpx_xz_zzzz_0, tpx_yy_xxxx_0, tpx_yy_xxxy_0, tpx_yy_xxxz_0, \
                                     tpx_yy_xxyy_0, tpy_xxxz_xxyy_0, tpy_xxxz_xxyz_0, tpy_xxxz_xxzz_0, tpy_xxxz_xyyy_0, \
                                     tpy_xxxz_xyyz_0, tpy_xxxz_xyzz_0, tpy_xxxz_xzzz_0, tpy_xxxz_yyyy_0, tpy_xxxz_yyyz_0, \
                                     tpy_xxxz_yyzz_0, tpy_xxxz_yzzz_0, tpy_xxxz_zzzz_0, tpy_xxyy_xxxx_0, tpy_xxyy_xxxy_0, \
                                     tpy_xxyy_xxxz_0, tpy_xxyy_xxyy_0, tpy_xxz_xxyy_0, tpy_xxz_xxyz_0, tpy_xxz_xxzz_0, \
                                     tpy_xxz_xyy_0, tpy_xxz_xyyy_0, tpy_xxz_xyyz_0, tpy_xxz_xyz_0, tpy_xxz_xyzz_0, \
                                     tpy_xxz_xzz_0, tpy_xxz_xzzz_0, tpy_xxz_yyy_0, tpy_xxz_yyyy_0, tpy_xxz_yyyz_0, \
                                     tpy_xxz_yyz_0, tpy_xxz_yyzz_0, tpy_xxz_yzz_0, tpy_xxz_yzzz_0, tpy_xxz_zzz_0, \
                                     tpy_xxz_zzzz_0, tpy_xyy_xxx_0, tpy_xyy_xxxx_0, tpy_xyy_xxxy_0, tpy_xyy_xxxz_0, \
                                     tpy_xyy_xxy_0, tpy_xyy_xxyy_0, tpy_xyy_xxz_0, tpy_xyy_xyy_0, tpy_xz_xxyy_0, \
                                     tpy_xz_xxyz_0, tpy_xz_xxzz_0, tpy_xz_xyyy_0, tpy_xz_xyyz_0, tpy_xz_xyzz_0, \
                                     tpy_xz_xzzz_0, tpy_xz_yyyy_0, tpy_xz_yyyz_0, tpy_xz_yyzz_0, tpy_xz_yzzz_0, \
                                     tpy_xz_zzzz_0, tpy_yy_xxxx_0, tpy_yy_xxxy_0, tpy_yy_xxxz_0, tpy_yy_xxyy_0, \
                                     tpz_xxxz_xxxz_0, tpz_xxxz_xxyy_0, tpz_xxxz_xxyz_0, tpz_xxxz_xxzz_0, tpz_xxxz_xyyy_0, \
                                     tpz_xxxz_xyyz_0, tpz_xxxz_xyzz_0, tpz_xxxz_xzzz_0, tpz_xxxz_yyyy_0, tpz_xxxz_yyyz_0, \
                                     tpz_xxxz_yyzz_0, tpz_xxxz_yzzz_0, tpz_xxxz_zzzz_0, tpz_xxyy_xxxx_0, tpz_xxyy_xxxy_0, \
                                     tpz_xxyy_xxxz_0, tpz_xxyy_xxyy_0, tpz_xxz_xxxz_0, tpz_xxz_xxyy_0, tpz_xxz_xxyz_0, \
                                     tpz_xxz_xxz_0, tpz_xxz_xxzz_0, tpz_xxz_xyy_0, tpz_xxz_xyyy_0, tpz_xxz_xyyz_0, \
                                     tpz_xxz_xyz_0, tpz_xxz_xyzz_0, tpz_xxz_xzz_0, tpz_xxz_xzzz_0, tpz_xxz_yyy_0, \
                                     tpz_xxz_yyyy_0, tpz_xxz_yyyz_0, tpz_xxz_yyz_0, tpz_xxz_yyzz_0, tpz_xxz_yzz_0, \
                                     tpz_xxz_yzzz_0, tpz_xxz_zzz_0, tpz_xxz_zzzz_0, tpz_xyy_xxx_0, tpz_xyy_xxxx_0, \
                                     tpz_xyy_xxxy_0, tpz_xyy_xxxz_0, tpz_xyy_xxy_0, tpz_xyy_xxyy_0, tpz_xyy_xxz_0, \
                                     tpz_xyy_xyy_0, tpz_xz_xxxz_0, tpz_xz_xxyy_0, tpz_xz_xxyz_0, tpz_xz_xxzz_0, \
                                     tpz_xz_xyyy_0, tpz_xz_xyyz_0, tpz_xz_xyzz_0, tpz_xz_xzzz_0, tpz_xz_yyyy_0, \
                                     tpz_xz_yyyz_0, tpz_xz_yyzz_0, tpz_xz_yzzz_0, tpz_xz_zzzz_0, tpz_yy_xxxx_0, \
                                     tpz_yy_xxxy_0, tpz_yy_xxxz_0, tpz_yy_xxyy_0, ts_xxz_xxyy_0, ts_xxz_xxyz_0, \
                                     ts_xxz_xxzz_0, ts_xxz_xyyy_0, ts_xxz_xyyz_0, ts_xxz_xyzz_0, ts_xxz_xzzz_0, \
                                     ts_xxz_yyyy_0, ts_xxz_yyyz_0, ts_xxz_yyzz_0, ts_xxz_yzzz_0, ts_xxz_zzzz_0, \
                                     ts_xyy_xxxx_0, ts_xyy_xxxy_0, ts_xyy_xxxz_0, ts_xyy_xxyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_xxxz_xxxz_0[j] = pa_x[j] * tpz_xxz_xxxz_0[j] + fl1_fx * tpz_xz_xxxz_0[j] + 1.5 * fl1_fx * tpz_xxz_xxz_0[j];

            tpx_xxxz_xxyy_0[j] =
                pa_x[j] * tpx_xxz_xxyy_0[j] + fl1_fx * tpx_xz_xxyy_0[j] + fl1_fx * tpx_xxz_xyy_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxyy_0[j];

            tpy_xxxz_xxyy_0[j] = pa_x[j] * tpy_xxz_xxyy_0[j] + fl1_fx * tpy_xz_xxyy_0[j] + fl1_fx * tpy_xxz_xyy_0[j];

            tpz_xxxz_xxyy_0[j] = pa_x[j] * tpz_xxz_xxyy_0[j] + fl1_fx * tpz_xz_xxyy_0[j] + fl1_fx * tpz_xxz_xyy_0[j];

            tpx_xxxz_xxyz_0[j] =
                pa_x[j] * tpx_xxz_xxyz_0[j] + fl1_fx * tpx_xz_xxyz_0[j] + fl1_fx * tpx_xxz_xyz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxyz_0[j];

            tpy_xxxz_xxyz_0[j] = pa_x[j] * tpy_xxz_xxyz_0[j] + fl1_fx * tpy_xz_xxyz_0[j] + fl1_fx * tpy_xxz_xyz_0[j];

            tpz_xxxz_xxyz_0[j] = pa_x[j] * tpz_xxz_xxyz_0[j] + fl1_fx * tpz_xz_xxyz_0[j] + fl1_fx * tpz_xxz_xyz_0[j];

            tpx_xxxz_xxzz_0[j] =
                pa_x[j] * tpx_xxz_xxzz_0[j] + fl1_fx * tpx_xz_xxzz_0[j] + fl1_fx * tpx_xxz_xzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxzz_0[j];

            tpy_xxxz_xxzz_0[j] = pa_x[j] * tpy_xxz_xxzz_0[j] + fl1_fx * tpy_xz_xxzz_0[j] + fl1_fx * tpy_xxz_xzz_0[j];

            tpz_xxxz_xxzz_0[j] = pa_x[j] * tpz_xxz_xxzz_0[j] + fl1_fx * tpz_xz_xxzz_0[j] + fl1_fx * tpz_xxz_xzz_0[j];

            tpx_xxxz_xyyy_0[j] =
                pa_x[j] * tpx_xxz_xyyy_0[j] + fl1_fx * tpx_xz_xyyy_0[j] + 0.5 * fl1_fx * tpx_xxz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xxz_xyyy_0[j];

            tpy_xxxz_xyyy_0[j] = pa_x[j] * tpy_xxz_xyyy_0[j] + fl1_fx * tpy_xz_xyyy_0[j] + 0.5 * fl1_fx * tpy_xxz_yyy_0[j];

            tpz_xxxz_xyyy_0[j] = pa_x[j] * tpz_xxz_xyyy_0[j] + fl1_fx * tpz_xz_xyyy_0[j] + 0.5 * fl1_fx * tpz_xxz_yyy_0[j];

            tpx_xxxz_xyyz_0[j] =
                pa_x[j] * tpx_xxz_xyyz_0[j] + fl1_fx * tpx_xz_xyyz_0[j] + 0.5 * fl1_fx * tpx_xxz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xyyz_0[j];

            tpy_xxxz_xyyz_0[j] = pa_x[j] * tpy_xxz_xyyz_0[j] + fl1_fx * tpy_xz_xyyz_0[j] + 0.5 * fl1_fx * tpy_xxz_yyz_0[j];

            tpz_xxxz_xyyz_0[j] = pa_x[j] * tpz_xxz_xyyz_0[j] + fl1_fx * tpz_xz_xyyz_0[j] + 0.5 * fl1_fx * tpz_xxz_yyz_0[j];

            tpx_xxxz_xyzz_0[j] =
                pa_x[j] * tpx_xxz_xyzz_0[j] + fl1_fx * tpx_xz_xyzz_0[j] + 0.5 * fl1_fx * tpx_xxz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xyzz_0[j];

            tpy_xxxz_xyzz_0[j] = pa_x[j] * tpy_xxz_xyzz_0[j] + fl1_fx * tpy_xz_xyzz_0[j] + 0.5 * fl1_fx * tpy_xxz_yzz_0[j];

            tpz_xxxz_xyzz_0[j] = pa_x[j] * tpz_xxz_xyzz_0[j] + fl1_fx * tpz_xz_xyzz_0[j] + 0.5 * fl1_fx * tpz_xxz_yzz_0[j];

            tpx_xxxz_xzzz_0[j] =
                pa_x[j] * tpx_xxz_xzzz_0[j] + fl1_fx * tpx_xz_xzzz_0[j] + 0.5 * fl1_fx * tpx_xxz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xzzz_0[j];

            tpy_xxxz_xzzz_0[j] = pa_x[j] * tpy_xxz_xzzz_0[j] + fl1_fx * tpy_xz_xzzz_0[j] + 0.5 * fl1_fx * tpy_xxz_zzz_0[j];

            tpz_xxxz_xzzz_0[j] = pa_x[j] * tpz_xxz_xzzz_0[j] + fl1_fx * tpz_xz_xzzz_0[j] + 0.5 * fl1_fx * tpz_xxz_zzz_0[j];

            tpx_xxxz_yyyy_0[j] = pa_x[j] * tpx_xxz_yyyy_0[j] + fl1_fx * tpx_xz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xxz_yyyy_0[j];

            tpy_xxxz_yyyy_0[j] = pa_x[j] * tpy_xxz_yyyy_0[j] + fl1_fx * tpy_xz_yyyy_0[j];

            tpz_xxxz_yyyy_0[j] = pa_x[j] * tpz_xxz_yyyy_0[j] + fl1_fx * tpz_xz_yyyy_0[j];

            tpx_xxxz_yyyz_0[j] = pa_x[j] * tpx_xxz_yyyz_0[j] + fl1_fx * tpx_xz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xxz_yyyz_0[j];

            tpy_xxxz_yyyz_0[j] = pa_x[j] * tpy_xxz_yyyz_0[j] + fl1_fx * tpy_xz_yyyz_0[j];

            tpz_xxxz_yyyz_0[j] = pa_x[j] * tpz_xxz_yyyz_0[j] + fl1_fx * tpz_xz_yyyz_0[j];

            tpx_xxxz_yyzz_0[j] = pa_x[j] * tpx_xxz_yyzz_0[j] + fl1_fx * tpx_xz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_yyzz_0[j];

            tpy_xxxz_yyzz_0[j] = pa_x[j] * tpy_xxz_yyzz_0[j] + fl1_fx * tpy_xz_yyzz_0[j];

            tpz_xxxz_yyzz_0[j] = pa_x[j] * tpz_xxz_yyzz_0[j] + fl1_fx * tpz_xz_yyzz_0[j];

            tpx_xxxz_yzzz_0[j] = pa_x[j] * tpx_xxz_yzzz_0[j] + fl1_fx * tpx_xz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_yzzz_0[j];

            tpy_xxxz_yzzz_0[j] = pa_x[j] * tpy_xxz_yzzz_0[j] + fl1_fx * tpy_xz_yzzz_0[j];

            tpz_xxxz_yzzz_0[j] = pa_x[j] * tpz_xxz_yzzz_0[j] + fl1_fx * tpz_xz_yzzz_0[j];

            tpx_xxxz_zzzz_0[j] = pa_x[j] * tpx_xxz_zzzz_0[j] + fl1_fx * tpx_xz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_zzzz_0[j];

            tpy_xxxz_zzzz_0[j] = pa_x[j] * tpy_xxz_zzzz_0[j] + fl1_fx * tpy_xz_zzzz_0[j];

            tpz_xxxz_zzzz_0[j] = pa_x[j] * tpz_xxz_zzzz_0[j] + fl1_fx * tpz_xz_zzzz_0[j];

            tpx_xxyy_xxxx_0[j] =
                pa_x[j] * tpx_xyy_xxxx_0[j] + 0.5 * fl1_fx * tpx_yy_xxxx_0[j] + 2.0 * fl1_fx * tpx_xyy_xxx_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxxx_0[j];

            tpy_xxyy_xxxx_0[j] = pa_x[j] * tpy_xyy_xxxx_0[j] + 0.5 * fl1_fx * tpy_yy_xxxx_0[j] + 2.0 * fl1_fx * tpy_xyy_xxx_0[j];

            tpz_xxyy_xxxx_0[j] = pa_x[j] * tpz_xyy_xxxx_0[j] + 0.5 * fl1_fx * tpz_yy_xxxx_0[j] + 2.0 * fl1_fx * tpz_xyy_xxx_0[j];

            tpx_xxyy_xxxy_0[j] =
                pa_x[j] * tpx_xyy_xxxy_0[j] + 0.5 * fl1_fx * tpx_yy_xxxy_0[j] + 1.5 * fl1_fx * tpx_xyy_xxy_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxxy_0[j];

            tpy_xxyy_xxxy_0[j] = pa_x[j] * tpy_xyy_xxxy_0[j] + 0.5 * fl1_fx * tpy_yy_xxxy_0[j] + 1.5 * fl1_fx * tpy_xyy_xxy_0[j];

            tpz_xxyy_xxxy_0[j] = pa_x[j] * tpz_xyy_xxxy_0[j] + 0.5 * fl1_fx * tpz_yy_xxxy_0[j] + 1.5 * fl1_fx * tpz_xyy_xxy_0[j];

            tpx_xxyy_xxxz_0[j] =
                pa_x[j] * tpx_xyy_xxxz_0[j] + 0.5 * fl1_fx * tpx_yy_xxxz_0[j] + 1.5 * fl1_fx * tpx_xyy_xxz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxxz_0[j];

            tpy_xxyy_xxxz_0[j] = pa_x[j] * tpy_xyy_xxxz_0[j] + 0.5 * fl1_fx * tpy_yy_xxxz_0[j] + 1.5 * fl1_fx * tpy_xyy_xxz_0[j];

            tpz_xxyy_xxxz_0[j] = pa_x[j] * tpz_xyy_xxxz_0[j] + 0.5 * fl1_fx * tpz_yy_xxxz_0[j] + 1.5 * fl1_fx * tpz_xyy_xxz_0[j];

            tpx_xxyy_xxyy_0[j] =
                pa_x[j] * tpx_xyy_xxyy_0[j] + 0.5 * fl1_fx * tpx_yy_xxyy_0[j] + fl1_fx * tpx_xyy_xyy_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxyy_0[j];

            tpy_xxyy_xxyy_0[j] = pa_x[j] * tpy_xyy_xxyy_0[j] + 0.5 * fl1_fx * tpy_yy_xxyy_0[j] + fl1_fx * tpy_xyy_xyy_0[j];

            tpz_xxyy_xxyy_0[j] = pa_x[j] * tpz_xyy_xxyy_0[j] + 0.5 * fl1_fx * tpz_yy_xxyy_0[j] + fl1_fx * tpz_xyy_xyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_147_195(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (147,195)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 49);

        auto tpy_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tpz_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 49);

        auto tpx_xyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 50);

        auto tpy_xyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 50);

        auto tpz_xyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 50);

        auto tpx_xyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 51);

        auto tpy_xyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 51);

        auto tpz_xyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 51);

        auto tpx_xyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 52);

        auto tpy_xyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 52);

        auto tpz_xyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 52);

        auto tpx_xyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 53);

        auto tpy_xyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 53);

        auto tpz_xyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 53);

        auto tpx_xyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 54);

        auto tpy_xyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 54);

        auto tpz_xyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 54);

        auto tpx_xyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 55);

        auto tpy_xyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 55);

        auto tpz_xyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 55);

        auto tpx_xyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 56);

        auto tpy_xyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 56);

        auto tpz_xyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 56);

        auto tpx_xyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 57);

        auto tpy_xyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 57);

        auto tpz_xyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 57);

        auto tpx_xyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 58);

        auto tpy_xyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 58);

        auto tpz_xyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 58);

        auto tpx_xyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 59);

        auto tpy_xyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 59);

        auto tpz_xyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 59);

        auto tpx_xyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 60);

        auto tpy_xyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 60);

        auto tpz_xyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 60);

        auto tpx_xyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 61);

        auto tpy_xyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 61);

        auto tpz_xyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 61);

        auto tpx_xyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 62);

        auto tpy_xyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 62);

        auto tpz_xyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 62);

        auto tpx_xyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 63);

        auto tpy_xyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 63);

        auto tpz_xyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 63);

        auto tpx_xyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 64);

        auto tpy_xyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 64);

        auto tpz_xyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 64);

        auto tpx_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 49);

        auto tpy_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tpz_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tpx_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 50);

        auto tpy_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tpz_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tpx_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 51);

        auto tpy_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tpz_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tpx_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 52);

        auto tpy_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tpz_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tpx_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 53);

        auto tpy_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tpz_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tpx_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 54);

        auto tpy_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tpz_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tpx_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 55);

        auto tpy_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tpz_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tpx_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 56);

        auto tpy_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tpz_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tpx_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 57);

        auto tpy_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tpz_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tpx_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 58);

        auto tpy_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tpz_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tpx_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 59);

        auto tpy_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tpz_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tpx_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 60);

        auto tpy_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tpz_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tpx_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 61);

        auto tpy_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tpz_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tpx_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 62);

        auto tpy_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tpz_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tpx_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 63);

        auto tpy_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tpz_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tpx_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 64);

        auto tpy_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tpz_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tpx_xyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 34);

        auto tpy_xyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 34);

        auto tpz_xyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 34);

        auto tpx_xyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 35);

        auto tpy_xyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 35);

        auto tpz_xyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 35);

        auto tpx_xyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 36);

        auto tpy_xyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 36);

        auto tpz_xyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 36);

        auto tpx_xyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 37);

        auto tpy_xyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 37);

        auto tpz_xyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 37);

        auto tpx_xyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 38);

        auto tpy_xyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 38);

        auto tpz_xyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 38);

        auto tpx_xyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 39);

        auto tpy_xyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 39);

        auto tpz_xyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 39);

        auto tpx_xyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 40);

        auto tpy_xyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 40);

        auto tpz_xyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 40);

        auto tpx_xyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 41);

        auto tpy_xyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 41);

        auto tpz_xyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 41);

        auto tpx_xyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 42);

        auto tpy_xyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 42);

        auto tpz_xyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 42);

        auto tpx_xyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 43);

        auto tpy_xyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 43);

        auto tpz_xyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 43);

        auto tpx_xyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 44);

        auto tpy_xyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 44);

        auto tpz_xyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 44);

        auto ts_xyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 49);

        auto ts_xyy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 50);

        auto ts_xyy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 51);

        auto ts_xyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 52);

        auto ts_xyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 53);

        auto ts_xyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 54);

        auto ts_xyy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 55);

        auto ts_xyy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 56);

        auto ts_xyy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 57);

        auto ts_xyy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 58);

        auto ts_xyy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 59);

        auto ts_xyz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 60);

        auto ts_xyz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 61);

        auto ts_xyz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 62);

        auto ts_xyz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 63);

        auto ts_xyz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 64);

        // set up pointers to integrals

        auto tpx_xxyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 49);

        auto tpy_xxyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 49);

        auto tpz_xxyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 49);

        auto tpx_xxyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 50);

        auto tpy_xxyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 50);

        auto tpz_xxyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 50);

        auto tpx_xxyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 51);

        auto tpy_xxyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 51);

        auto tpz_xxyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 51);

        auto tpx_xxyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 52);

        auto tpy_xxyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 52);

        auto tpz_xxyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 52);

        auto tpx_xxyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 53);

        auto tpy_xxyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 53);

        auto tpz_xxyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 53);

        auto tpx_xxyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 54);

        auto tpy_xxyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 54);

        auto tpz_xxyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 54);

        auto tpx_xxyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 55);

        auto tpy_xxyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 55);

        auto tpz_xxyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 55);

        auto tpx_xxyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 56);

        auto tpy_xxyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 56);

        auto tpz_xxyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 56);

        auto tpx_xxyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 57);

        auto tpy_xxyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 57);

        auto tpz_xxyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 57);

        auto tpx_xxyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 58);

        auto tpy_xxyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 58);

        auto tpz_xxyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 58);

        auto tpx_xxyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 59);

        auto tpy_xxyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 59);

        auto tpz_xxyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 59);

        auto tpx_xxyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 60);

        auto tpy_xxyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 60);

        auto tpz_xxyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 60);

        auto tpx_xxyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 61);

        auto tpy_xxyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 61);

        auto tpz_xxyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 61);

        auto tpx_xxyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 62);

        auto tpy_xxyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 62);

        auto tpz_xxyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 62);

        auto tpx_xxyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 63);

        auto tpy_xxyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 63);

        auto tpz_xxyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 63);

        auto tpx_xxyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 64);

        auto tpy_xxyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 64);

        auto tpz_xxyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 64);

        // Batch of Integrals (147,195)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxyy_xxyz_0, tpx_xxyy_xxzz_0, tpx_xxyy_xyyy_0, \
                                     tpx_xxyy_xyyz_0, tpx_xxyy_xyzz_0, tpx_xxyy_xzzz_0, tpx_xxyy_yyyy_0, tpx_xxyy_yyyz_0, \
                                     tpx_xxyy_yyzz_0, tpx_xxyy_yzzz_0, tpx_xxyy_zzzz_0, tpx_xxyz_xxxx_0, tpx_xxyz_xxxy_0, \
                                     tpx_xxyz_xxxz_0, tpx_xxyz_xxyy_0, tpx_xxyz_xxyz_0, tpx_xyy_xxyz_0, tpx_xyy_xxzz_0, \
                                     tpx_xyy_xyyy_0, tpx_xyy_xyyz_0, tpx_xyy_xyz_0, tpx_xyy_xyzz_0, tpx_xyy_xzz_0, \
                                     tpx_xyy_xzzz_0, tpx_xyy_yyy_0, tpx_xyy_yyyy_0, tpx_xyy_yyyz_0, tpx_xyy_yyz_0, \
                                     tpx_xyy_yyzz_0, tpx_xyy_yzz_0, tpx_xyy_yzzz_0, tpx_xyy_zzz_0, tpx_xyy_zzzz_0, \
                                     tpx_xyz_xxx_0, tpx_xyz_xxxx_0, tpx_xyz_xxxy_0, tpx_xyz_xxxz_0, tpx_xyz_xxy_0, \
                                     tpx_xyz_xxyy_0, tpx_xyz_xxyz_0, tpx_xyz_xxz_0, tpx_xyz_xyy_0, tpx_xyz_xyz_0, \
                                     tpx_yy_xxyz_0, tpx_yy_xxzz_0, tpx_yy_xyyy_0, tpx_yy_xyyz_0, tpx_yy_xyzz_0, \
                                     tpx_yy_xzzz_0, tpx_yy_yyyy_0, tpx_yy_yyyz_0, tpx_yy_yyzz_0, tpx_yy_yzzz_0, \
                                     tpx_yy_zzzz_0, tpx_yz_xxxx_0, tpx_yz_xxxy_0, tpx_yz_xxxz_0, tpx_yz_xxyy_0, \
                                     tpx_yz_xxyz_0, tpy_xxyy_xxyz_0, tpy_xxyy_xxzz_0, tpy_xxyy_xyyy_0, tpy_xxyy_xyyz_0, \
                                     tpy_xxyy_xyzz_0, tpy_xxyy_xzzz_0, tpy_xxyy_yyyy_0, tpy_xxyy_yyyz_0, tpy_xxyy_yyzz_0, \
                                     tpy_xxyy_yzzz_0, tpy_xxyy_zzzz_0, tpy_xxyz_xxxx_0, tpy_xxyz_xxxy_0, tpy_xxyz_xxxz_0, \
                                     tpy_xxyz_xxyy_0, tpy_xxyz_xxyz_0, tpy_xyy_xxyz_0, tpy_xyy_xxzz_0, tpy_xyy_xyyy_0, \
                                     tpy_xyy_xyyz_0, tpy_xyy_xyz_0, tpy_xyy_xyzz_0, tpy_xyy_xzz_0, tpy_xyy_xzzz_0, \
                                     tpy_xyy_yyy_0, tpy_xyy_yyyy_0, tpy_xyy_yyyz_0, tpy_xyy_yyz_0, tpy_xyy_yyzz_0, \
                                     tpy_xyy_yzz_0, tpy_xyy_yzzz_0, tpy_xyy_zzz_0, tpy_xyy_zzzz_0, tpy_xyz_xxx_0, \
                                     tpy_xyz_xxxx_0, tpy_xyz_xxxy_0, tpy_xyz_xxxz_0, tpy_xyz_xxy_0, tpy_xyz_xxyy_0, \
                                     tpy_xyz_xxyz_0, tpy_xyz_xxz_0, tpy_xyz_xyy_0, tpy_xyz_xyz_0, tpy_yy_xxyz_0, \
                                     tpy_yy_xxzz_0, tpy_yy_xyyy_0, tpy_yy_xyyz_0, tpy_yy_xyzz_0, tpy_yy_xzzz_0, \
                                     tpy_yy_yyyy_0, tpy_yy_yyyz_0, tpy_yy_yyzz_0, tpy_yy_yzzz_0, tpy_yy_zzzz_0, \
                                     tpy_yz_xxxx_0, tpy_yz_xxxy_0, tpy_yz_xxxz_0, tpy_yz_xxyy_0, tpy_yz_xxyz_0, \
                                     tpz_xxyy_xxyz_0, tpz_xxyy_xxzz_0, tpz_xxyy_xyyy_0, tpz_xxyy_xyyz_0, tpz_xxyy_xyzz_0, \
                                     tpz_xxyy_xzzz_0, tpz_xxyy_yyyy_0, tpz_xxyy_yyyz_0, tpz_xxyy_yyzz_0, tpz_xxyy_yzzz_0, \
                                     tpz_xxyy_zzzz_0, tpz_xxyz_xxxx_0, tpz_xxyz_xxxy_0, tpz_xxyz_xxxz_0, tpz_xxyz_xxyy_0, \
                                     tpz_xxyz_xxyz_0, tpz_xyy_xxyz_0, tpz_xyy_xxzz_0, tpz_xyy_xyyy_0, tpz_xyy_xyyz_0, \
                                     tpz_xyy_xyz_0, tpz_xyy_xyzz_0, tpz_xyy_xzz_0, tpz_xyy_xzzz_0, tpz_xyy_yyy_0, \
                                     tpz_xyy_yyyy_0, tpz_xyy_yyyz_0, tpz_xyy_yyz_0, tpz_xyy_yyzz_0, tpz_xyy_yzz_0, \
                                     tpz_xyy_yzzz_0, tpz_xyy_zzz_0, tpz_xyy_zzzz_0, tpz_xyz_xxx_0, tpz_xyz_xxxx_0, \
                                     tpz_xyz_xxxy_0, tpz_xyz_xxxz_0, tpz_xyz_xxy_0, tpz_xyz_xxyy_0, tpz_xyz_xxyz_0, \
                                     tpz_xyz_xxz_0, tpz_xyz_xyy_0, tpz_xyz_xyz_0, tpz_yy_xxyz_0, tpz_yy_xxzz_0, \
                                     tpz_yy_xyyy_0, tpz_yy_xyyz_0, tpz_yy_xyzz_0, tpz_yy_xzzz_0, tpz_yy_yyyy_0, \
                                     tpz_yy_yyyz_0, tpz_yy_yyzz_0, tpz_yy_yzzz_0, tpz_yy_zzzz_0, tpz_yz_xxxx_0, \
                                     tpz_yz_xxxy_0, tpz_yz_xxxz_0, tpz_yz_xxyy_0, tpz_yz_xxyz_0, ts_xyy_xxyz_0, \
                                     ts_xyy_xxzz_0, ts_xyy_xyyy_0, ts_xyy_xyyz_0, ts_xyy_xyzz_0, ts_xyy_xzzz_0, \
                                     ts_xyy_yyyy_0, ts_xyy_yyyz_0, ts_xyy_yyzz_0, ts_xyy_yzzz_0, ts_xyy_zzzz_0, \
                                     ts_xyz_xxxx_0, ts_xyz_xxxy_0, ts_xyz_xxxz_0, ts_xyz_xxyy_0, ts_xyz_xxyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxyy_xxyz_0[j] =
                pa_x[j] * tpx_xyy_xxyz_0[j] + 0.5 * fl1_fx * tpx_yy_xxyz_0[j] + fl1_fx * tpx_xyy_xyz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxyz_0[j];

            tpy_xxyy_xxyz_0[j] = pa_x[j] * tpy_xyy_xxyz_0[j] + 0.5 * fl1_fx * tpy_yy_xxyz_0[j] + fl1_fx * tpy_xyy_xyz_0[j];

            tpz_xxyy_xxyz_0[j] = pa_x[j] * tpz_xyy_xxyz_0[j] + 0.5 * fl1_fx * tpz_yy_xxyz_0[j] + fl1_fx * tpz_xyy_xyz_0[j];

            tpx_xxyy_xxzz_0[j] =
                pa_x[j] * tpx_xyy_xxzz_0[j] + 0.5 * fl1_fx * tpx_yy_xxzz_0[j] + fl1_fx * tpx_xyy_xzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxzz_0[j];

            tpy_xxyy_xxzz_0[j] = pa_x[j] * tpy_xyy_xxzz_0[j] + 0.5 * fl1_fx * tpy_yy_xxzz_0[j] + fl1_fx * tpy_xyy_xzz_0[j];

            tpz_xxyy_xxzz_0[j] = pa_x[j] * tpz_xyy_xxzz_0[j] + 0.5 * fl1_fx * tpz_yy_xxzz_0[j] + fl1_fx * tpz_xyy_xzz_0[j];

            tpx_xxyy_xyyy_0[j] =
                pa_x[j] * tpx_xyy_xyyy_0[j] + 0.5 * fl1_fx * tpx_yy_xyyy_0[j] + 0.5 * fl1_fx * tpx_xyy_yyy_0[j] - fl1_fgb * fl1_fx * ts_xyy_xyyy_0[j];

            tpy_xxyy_xyyy_0[j] = pa_x[j] * tpy_xyy_xyyy_0[j] + 0.5 * fl1_fx * tpy_yy_xyyy_0[j] + 0.5 * fl1_fx * tpy_xyy_yyy_0[j];

            tpz_xxyy_xyyy_0[j] = pa_x[j] * tpz_xyy_xyyy_0[j] + 0.5 * fl1_fx * tpz_yy_xyyy_0[j] + 0.5 * fl1_fx * tpz_xyy_yyy_0[j];

            tpx_xxyy_xyyz_0[j] =
                pa_x[j] * tpx_xyy_xyyz_0[j] + 0.5 * fl1_fx * tpx_yy_xyyz_0[j] + 0.5 * fl1_fx * tpx_xyy_yyz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xyyz_0[j];

            tpy_xxyy_xyyz_0[j] = pa_x[j] * tpy_xyy_xyyz_0[j] + 0.5 * fl1_fx * tpy_yy_xyyz_0[j] + 0.5 * fl1_fx * tpy_xyy_yyz_0[j];

            tpz_xxyy_xyyz_0[j] = pa_x[j] * tpz_xyy_xyyz_0[j] + 0.5 * fl1_fx * tpz_yy_xyyz_0[j] + 0.5 * fl1_fx * tpz_xyy_yyz_0[j];

            tpx_xxyy_xyzz_0[j] =
                pa_x[j] * tpx_xyy_xyzz_0[j] + 0.5 * fl1_fx * tpx_yy_xyzz_0[j] + 0.5 * fl1_fx * tpx_xyy_yzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xyzz_0[j];

            tpy_xxyy_xyzz_0[j] = pa_x[j] * tpy_xyy_xyzz_0[j] + 0.5 * fl1_fx * tpy_yy_xyzz_0[j] + 0.5 * fl1_fx * tpy_xyy_yzz_0[j];

            tpz_xxyy_xyzz_0[j] = pa_x[j] * tpz_xyy_xyzz_0[j] + 0.5 * fl1_fx * tpz_yy_xyzz_0[j] + 0.5 * fl1_fx * tpz_xyy_yzz_0[j];

            tpx_xxyy_xzzz_0[j] =
                pa_x[j] * tpx_xyy_xzzz_0[j] + 0.5 * fl1_fx * tpx_yy_xzzz_0[j] + 0.5 * fl1_fx * tpx_xyy_zzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xzzz_0[j];

            tpy_xxyy_xzzz_0[j] = pa_x[j] * tpy_xyy_xzzz_0[j] + 0.5 * fl1_fx * tpy_yy_xzzz_0[j] + 0.5 * fl1_fx * tpy_xyy_zzz_0[j];

            tpz_xxyy_xzzz_0[j] = pa_x[j] * tpz_xyy_xzzz_0[j] + 0.5 * fl1_fx * tpz_yy_xzzz_0[j] + 0.5 * fl1_fx * tpz_xyy_zzz_0[j];

            tpx_xxyy_yyyy_0[j] = pa_x[j] * tpx_xyy_yyyy_0[j] + 0.5 * fl1_fx * tpx_yy_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xyy_yyyy_0[j];

            tpy_xxyy_yyyy_0[j] = pa_x[j] * tpy_xyy_yyyy_0[j] + 0.5 * fl1_fx * tpy_yy_yyyy_0[j];

            tpz_xxyy_yyyy_0[j] = pa_x[j] * tpz_xyy_yyyy_0[j] + 0.5 * fl1_fx * tpz_yy_yyyy_0[j];

            tpx_xxyy_yyyz_0[j] = pa_x[j] * tpx_xyy_yyyz_0[j] + 0.5 * fl1_fx * tpx_yy_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xyy_yyyz_0[j];

            tpy_xxyy_yyyz_0[j] = pa_x[j] * tpy_xyy_yyyz_0[j] + 0.5 * fl1_fx * tpy_yy_yyyz_0[j];

            tpz_xxyy_yyyz_0[j] = pa_x[j] * tpz_xyy_yyyz_0[j] + 0.5 * fl1_fx * tpz_yy_yyyz_0[j];

            tpx_xxyy_yyzz_0[j] = pa_x[j] * tpx_xyy_yyzz_0[j] + 0.5 * fl1_fx * tpx_yy_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_yyzz_0[j];

            tpy_xxyy_yyzz_0[j] = pa_x[j] * tpy_xyy_yyzz_0[j] + 0.5 * fl1_fx * tpy_yy_yyzz_0[j];

            tpz_xxyy_yyzz_0[j] = pa_x[j] * tpz_xyy_yyzz_0[j] + 0.5 * fl1_fx * tpz_yy_yyzz_0[j];

            tpx_xxyy_yzzz_0[j] = pa_x[j] * tpx_xyy_yzzz_0[j] + 0.5 * fl1_fx * tpx_yy_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_yzzz_0[j];

            tpy_xxyy_yzzz_0[j] = pa_x[j] * tpy_xyy_yzzz_0[j] + 0.5 * fl1_fx * tpy_yy_yzzz_0[j];

            tpz_xxyy_yzzz_0[j] = pa_x[j] * tpz_xyy_yzzz_0[j] + 0.5 * fl1_fx * tpz_yy_yzzz_0[j];

            tpx_xxyy_zzzz_0[j] = pa_x[j] * tpx_xyy_zzzz_0[j] + 0.5 * fl1_fx * tpx_yy_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_zzzz_0[j];

            tpy_xxyy_zzzz_0[j] = pa_x[j] * tpy_xyy_zzzz_0[j] + 0.5 * fl1_fx * tpy_yy_zzzz_0[j];

            tpz_xxyy_zzzz_0[j] = pa_x[j] * tpz_xyy_zzzz_0[j] + 0.5 * fl1_fx * tpz_yy_zzzz_0[j];

            tpx_xxyz_xxxx_0[j] =
                pa_x[j] * tpx_xyz_xxxx_0[j] + 0.5 * fl1_fx * tpx_yz_xxxx_0[j] + 2.0 * fl1_fx * tpx_xyz_xxx_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxxx_0[j];

            tpy_xxyz_xxxx_0[j] = pa_x[j] * tpy_xyz_xxxx_0[j] + 0.5 * fl1_fx * tpy_yz_xxxx_0[j] + 2.0 * fl1_fx * tpy_xyz_xxx_0[j];

            tpz_xxyz_xxxx_0[j] = pa_x[j] * tpz_xyz_xxxx_0[j] + 0.5 * fl1_fx * tpz_yz_xxxx_0[j] + 2.0 * fl1_fx * tpz_xyz_xxx_0[j];

            tpx_xxyz_xxxy_0[j] =
                pa_x[j] * tpx_xyz_xxxy_0[j] + 0.5 * fl1_fx * tpx_yz_xxxy_0[j] + 1.5 * fl1_fx * tpx_xyz_xxy_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxxy_0[j];

            tpy_xxyz_xxxy_0[j] = pa_x[j] * tpy_xyz_xxxy_0[j] + 0.5 * fl1_fx * tpy_yz_xxxy_0[j] + 1.5 * fl1_fx * tpy_xyz_xxy_0[j];

            tpz_xxyz_xxxy_0[j] = pa_x[j] * tpz_xyz_xxxy_0[j] + 0.5 * fl1_fx * tpz_yz_xxxy_0[j] + 1.5 * fl1_fx * tpz_xyz_xxy_0[j];

            tpx_xxyz_xxxz_0[j] =
                pa_x[j] * tpx_xyz_xxxz_0[j] + 0.5 * fl1_fx * tpx_yz_xxxz_0[j] + 1.5 * fl1_fx * tpx_xyz_xxz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxxz_0[j];

            tpy_xxyz_xxxz_0[j] = pa_x[j] * tpy_xyz_xxxz_0[j] + 0.5 * fl1_fx * tpy_yz_xxxz_0[j] + 1.5 * fl1_fx * tpy_xyz_xxz_0[j];

            tpz_xxyz_xxxz_0[j] = pa_x[j] * tpz_xyz_xxxz_0[j] + 0.5 * fl1_fx * tpz_yz_xxxz_0[j] + 1.5 * fl1_fx * tpz_xyz_xxz_0[j];

            tpx_xxyz_xxyy_0[j] =
                pa_x[j] * tpx_xyz_xxyy_0[j] + 0.5 * fl1_fx * tpx_yz_xxyy_0[j] + fl1_fx * tpx_xyz_xyy_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxyy_0[j];

            tpy_xxyz_xxyy_0[j] = pa_x[j] * tpy_xyz_xxyy_0[j] + 0.5 * fl1_fx * tpy_yz_xxyy_0[j] + fl1_fx * tpy_xyz_xyy_0[j];

            tpz_xxyz_xxyy_0[j] = pa_x[j] * tpz_xyz_xxyy_0[j] + 0.5 * fl1_fx * tpz_yz_xxyy_0[j] + fl1_fx * tpz_xyz_xyy_0[j];

            tpx_xxyz_xxyz_0[j] =
                pa_x[j] * tpx_xyz_xxyz_0[j] + 0.5 * fl1_fx * tpx_yz_xxyz_0[j] + fl1_fx * tpx_xyz_xyz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxyz_0[j];

            tpy_xxyz_xxyz_0[j] = pa_x[j] * tpy_xyz_xxyz_0[j] + 0.5 * fl1_fx * tpy_yz_xxyz_0[j] + fl1_fx * tpy_xyz_xyz_0[j];

            tpz_xxyz_xxyz_0[j] = pa_x[j] * tpz_xyz_xxyz_0[j] + 0.5 * fl1_fx * tpz_yz_xxyz_0[j] + fl1_fx * tpz_xyz_xyz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_195_243(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (195,243)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 65);

        auto tpy_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tpz_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tpx_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 66);

        auto tpy_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 66);

        auto tpz_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 66);

        auto tpx_xyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 67);

        auto tpy_xyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 67);

        auto tpz_xyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 67);

        auto tpx_xyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 68);

        auto tpy_xyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 68);

        auto tpz_xyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 68);

        auto tpx_xyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 69);

        auto tpy_xyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 69);

        auto tpz_xyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 69);

        auto tpx_xyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 70);

        auto tpy_xyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 70);

        auto tpz_xyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 70);

        auto tpx_xyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 71);

        auto tpy_xyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 71);

        auto tpz_xyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 71);

        auto tpx_xyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 72);

        auto tpy_xyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 72);

        auto tpz_xyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 72);

        auto tpx_xyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 73);

        auto tpy_xyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 73);

        auto tpz_xyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 73);

        auto tpx_xyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 74);

        auto tpy_xyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 74);

        auto tpz_xyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 74);

        auto tpx_xzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 75);

        auto tpy_xzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 75);

        auto tpz_xzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 75);

        auto tpx_xzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 76);

        auto tpy_xzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 76);

        auto tpz_xzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 76);

        auto tpx_xzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 77);

        auto tpy_xzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 77);

        auto tpz_xzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 77);

        auto tpx_xzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 78);

        auto tpy_xzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 78);

        auto tpz_xzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 78);

        auto tpx_xzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 79);

        auto tpy_xzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 79);

        auto tpz_xzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 79);

        auto tpx_xzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 80);

        auto tpy_xzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 80);

        auto tpz_xzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 80);

        auto tpx_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 65);

        auto tpy_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tpz_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tpx_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 66);

        auto tpy_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tpz_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tpx_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 67);

        auto tpy_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tpz_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tpx_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 68);

        auto tpy_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tpz_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tpx_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 69);

        auto tpy_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tpz_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tpx_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 70);

        auto tpy_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tpz_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tpx_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 71);

        auto tpy_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tpz_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tpx_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 72);

        auto tpy_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tpz_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tpx_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 73);

        auto tpy_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tpz_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tpx_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 74);

        auto tpy_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tpz_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tpx_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 75);

        auto tpy_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tpz_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tpx_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 76);

        auto tpy_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tpz_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tpx_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 77);

        auto tpy_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tpz_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tpx_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 78);

        auto tpy_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tpz_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tpx_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 79);

        auto tpy_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tpz_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tpx_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 80);

        auto tpy_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tpz_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tpx_xyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 45);

        auto tpy_xyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 45);

        auto tpz_xyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 45);

        auto tpx_xyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 46);

        auto tpy_xyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 46);

        auto tpz_xyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 46);

        auto tpx_xyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 47);

        auto tpy_xyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 47);

        auto tpz_xyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 47);

        auto tpx_xyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 48);

        auto tpy_xyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 48);

        auto tpz_xyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 48);

        auto tpx_xyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 49);

        auto tpy_xyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 49);

        auto tpz_xyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 49);

        auto tpx_xzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 50);

        auto tpy_xzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 50);

        auto tpz_xzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 50);

        auto tpx_xzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 51);

        auto tpy_xzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 51);

        auto tpz_xzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 51);

        auto tpx_xzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 52);

        auto tpy_xzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 52);

        auto tpz_xzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 52);

        auto tpx_xzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 53);

        auto tpy_xzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 53);

        auto tpz_xzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 53);

        auto tpx_xzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 54);

        auto tpy_xzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 54);

        auto tpz_xzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 54);

        auto tpx_xzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 55);

        auto tpy_xzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 55);

        auto tpz_xzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 55);

        auto ts_xyz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 65);

        auto ts_xyz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 66);

        auto ts_xyz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 67);

        auto ts_xyz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 68);

        auto ts_xyz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 69);

        auto ts_xyz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 70);

        auto ts_xyz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 71);

        auto ts_xyz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 72);

        auto ts_xyz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 73);

        auto ts_xyz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 74);

        auto ts_xzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 75);

        auto ts_xzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 76);

        auto ts_xzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 77);

        auto ts_xzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 78);

        auto ts_xzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 79);

        auto ts_xzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 80);

        // set up pointers to integrals

        auto tpx_xxyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 65);

        auto tpy_xxyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 65);

        auto tpz_xxyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 65);

        auto tpx_xxyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 66);

        auto tpy_xxyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 66);

        auto tpz_xxyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 66);

        auto tpx_xxyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 67);

        auto tpy_xxyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 67);

        auto tpz_xxyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 67);

        auto tpx_xxyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 68);

        auto tpy_xxyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 68);

        auto tpz_xxyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 68);

        auto tpx_xxyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 69);

        auto tpy_xxyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 69);

        auto tpz_xxyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 69);

        auto tpx_xxyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 70);

        auto tpy_xxyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 70);

        auto tpz_xxyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 70);

        auto tpx_xxyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 71);

        auto tpy_xxyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 71);

        auto tpz_xxyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 71);

        auto tpx_xxyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 72);

        auto tpy_xxyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 72);

        auto tpz_xxyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 72);

        auto tpx_xxyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 73);

        auto tpy_xxyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 73);

        auto tpz_xxyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 73);

        auto tpx_xxyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 74);

        auto tpy_xxyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 74);

        auto tpz_xxyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 74);

        auto tpx_xxzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 75);

        auto tpy_xxzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 75);

        auto tpz_xxzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 75);

        auto tpx_xxzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 76);

        auto tpy_xxzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 76);

        auto tpz_xxzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 76);

        auto tpx_xxzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 77);

        auto tpy_xxzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 77);

        auto tpz_xxzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 77);

        auto tpx_xxzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 78);

        auto tpy_xxzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 78);

        auto tpz_xxzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 78);

        auto tpx_xxzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 79);

        auto tpy_xxzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 79);

        auto tpz_xxzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 79);

        auto tpx_xxzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 80);

        auto tpy_xxzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 80);

        auto tpz_xxzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 80);

        // Batch of Integrals (195,243)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxyz_xxzz_0, tpx_xxyz_xyyy_0, tpx_xxyz_xyyz_0, \
                                     tpx_xxyz_xyzz_0, tpx_xxyz_xzzz_0, tpx_xxyz_yyyy_0, tpx_xxyz_yyyz_0, tpx_xxyz_yyzz_0, \
                                     tpx_xxyz_yzzz_0, tpx_xxyz_zzzz_0, tpx_xxzz_xxxx_0, tpx_xxzz_xxxy_0, tpx_xxzz_xxxz_0, \
                                     tpx_xxzz_xxyy_0, tpx_xxzz_xxyz_0, tpx_xxzz_xxzz_0, tpx_xyz_xxzz_0, tpx_xyz_xyyy_0, \
                                     tpx_xyz_xyyz_0, tpx_xyz_xyzz_0, tpx_xyz_xzz_0, tpx_xyz_xzzz_0, tpx_xyz_yyy_0, \
                                     tpx_xyz_yyyy_0, tpx_xyz_yyyz_0, tpx_xyz_yyz_0, tpx_xyz_yyzz_0, tpx_xyz_yzz_0, \
                                     tpx_xyz_yzzz_0, tpx_xyz_zzz_0, tpx_xyz_zzzz_0, tpx_xzz_xxx_0, tpx_xzz_xxxx_0, \
                                     tpx_xzz_xxxy_0, tpx_xzz_xxxz_0, tpx_xzz_xxy_0, tpx_xzz_xxyy_0, tpx_xzz_xxyz_0, \
                                     tpx_xzz_xxz_0, tpx_xzz_xxzz_0, tpx_xzz_xyy_0, tpx_xzz_xyz_0, tpx_xzz_xzz_0, \
                                     tpx_yz_xxzz_0, tpx_yz_xyyy_0, tpx_yz_xyyz_0, tpx_yz_xyzz_0, tpx_yz_xzzz_0, \
                                     tpx_yz_yyyy_0, tpx_yz_yyyz_0, tpx_yz_yyzz_0, tpx_yz_yzzz_0, tpx_yz_zzzz_0, \
                                     tpx_zz_xxxx_0, tpx_zz_xxxy_0, tpx_zz_xxxz_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, \
                                     tpx_zz_xxzz_0, tpy_xxyz_xxzz_0, tpy_xxyz_xyyy_0, tpy_xxyz_xyyz_0, tpy_xxyz_xyzz_0, \
                                     tpy_xxyz_xzzz_0, tpy_xxyz_yyyy_0, tpy_xxyz_yyyz_0, tpy_xxyz_yyzz_0, tpy_xxyz_yzzz_0, \
                                     tpy_xxyz_zzzz_0, tpy_xxzz_xxxx_0, tpy_xxzz_xxxy_0, tpy_xxzz_xxxz_0, tpy_xxzz_xxyy_0, \
                                     tpy_xxzz_xxyz_0, tpy_xxzz_xxzz_0, tpy_xyz_xxzz_0, tpy_xyz_xyyy_0, tpy_xyz_xyyz_0, \
                                     tpy_xyz_xyzz_0, tpy_xyz_xzz_0, tpy_xyz_xzzz_0, tpy_xyz_yyy_0, tpy_xyz_yyyy_0, \
                                     tpy_xyz_yyyz_0, tpy_xyz_yyz_0, tpy_xyz_yyzz_0, tpy_xyz_yzz_0, tpy_xyz_yzzz_0, \
                                     tpy_xyz_zzz_0, tpy_xyz_zzzz_0, tpy_xzz_xxx_0, tpy_xzz_xxxx_0, tpy_xzz_xxxy_0, \
                                     tpy_xzz_xxxz_0, tpy_xzz_xxy_0, tpy_xzz_xxyy_0, tpy_xzz_xxyz_0, tpy_xzz_xxz_0, \
                                     tpy_xzz_xxzz_0, tpy_xzz_xyy_0, tpy_xzz_xyz_0, tpy_xzz_xzz_0, tpy_yz_xxzz_0, \
                                     tpy_yz_xyyy_0, tpy_yz_xyyz_0, tpy_yz_xyzz_0, tpy_yz_xzzz_0, tpy_yz_yyyy_0, \
                                     tpy_yz_yyyz_0, tpy_yz_yyzz_0, tpy_yz_yzzz_0, tpy_yz_zzzz_0, tpy_zz_xxxx_0, \
                                     tpy_zz_xxxy_0, tpy_zz_xxxz_0, tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxzz_0, \
                                     tpz_xxyz_xxzz_0, tpz_xxyz_xyyy_0, tpz_xxyz_xyyz_0, tpz_xxyz_xyzz_0, tpz_xxyz_xzzz_0, \
                                     tpz_xxyz_yyyy_0, tpz_xxyz_yyyz_0, tpz_xxyz_yyzz_0, tpz_xxyz_yzzz_0, tpz_xxyz_zzzz_0, \
                                     tpz_xxzz_xxxx_0, tpz_xxzz_xxxy_0, tpz_xxzz_xxxz_0, tpz_xxzz_xxyy_0, tpz_xxzz_xxyz_0, \
                                     tpz_xxzz_xxzz_0, tpz_xyz_xxzz_0, tpz_xyz_xyyy_0, tpz_xyz_xyyz_0, tpz_xyz_xyzz_0, \
                                     tpz_xyz_xzz_0, tpz_xyz_xzzz_0, tpz_xyz_yyy_0, tpz_xyz_yyyy_0, tpz_xyz_yyyz_0, \
                                     tpz_xyz_yyz_0, tpz_xyz_yyzz_0, tpz_xyz_yzz_0, tpz_xyz_yzzz_0, tpz_xyz_zzz_0, \
                                     tpz_xyz_zzzz_0, tpz_xzz_xxx_0, tpz_xzz_xxxx_0, tpz_xzz_xxxy_0, tpz_xzz_xxxz_0, \
                                     tpz_xzz_xxy_0, tpz_xzz_xxyy_0, tpz_xzz_xxyz_0, tpz_xzz_xxz_0, tpz_xzz_xxzz_0, \
                                     tpz_xzz_xyy_0, tpz_xzz_xyz_0, tpz_xzz_xzz_0, tpz_yz_xxzz_0, tpz_yz_xyyy_0, \
                                     tpz_yz_xyyz_0, tpz_yz_xyzz_0, tpz_yz_xzzz_0, tpz_yz_yyyy_0, tpz_yz_yyyz_0, \
                                     tpz_yz_yyzz_0, tpz_yz_yzzz_0, tpz_yz_zzzz_0, tpz_zz_xxxx_0, tpz_zz_xxxy_0, \
                                     tpz_zz_xxxz_0, tpz_zz_xxyy_0, tpz_zz_xxyz_0, tpz_zz_xxzz_0, ts_xyz_xxzz_0, \
                                     ts_xyz_xyyy_0, ts_xyz_xyyz_0, ts_xyz_xyzz_0, ts_xyz_xzzz_0, ts_xyz_yyyy_0, \
                                     ts_xyz_yyyz_0, ts_xyz_yyzz_0, ts_xyz_yzzz_0, ts_xyz_zzzz_0, ts_xzz_xxxx_0, \
                                     ts_xzz_xxxy_0, ts_xzz_xxxz_0, ts_xzz_xxyy_0, ts_xzz_xxyz_0, ts_xzz_xxzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxyz_xxzz_0[j] =
                pa_x[j] * tpx_xyz_xxzz_0[j] + 0.5 * fl1_fx * tpx_yz_xxzz_0[j] + fl1_fx * tpx_xyz_xzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxzz_0[j];

            tpy_xxyz_xxzz_0[j] = pa_x[j] * tpy_xyz_xxzz_0[j] + 0.5 * fl1_fx * tpy_yz_xxzz_0[j] + fl1_fx * tpy_xyz_xzz_0[j];

            tpz_xxyz_xxzz_0[j] = pa_x[j] * tpz_xyz_xxzz_0[j] + 0.5 * fl1_fx * tpz_yz_xxzz_0[j] + fl1_fx * tpz_xyz_xzz_0[j];

            tpx_xxyz_xyyy_0[j] =
                pa_x[j] * tpx_xyz_xyyy_0[j] + 0.5 * fl1_fx * tpx_yz_xyyy_0[j] + 0.5 * fl1_fx * tpx_xyz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xyz_xyyy_0[j];

            tpy_xxyz_xyyy_0[j] = pa_x[j] * tpy_xyz_xyyy_0[j] + 0.5 * fl1_fx * tpy_yz_xyyy_0[j] + 0.5 * fl1_fx * tpy_xyz_yyy_0[j];

            tpz_xxyz_xyyy_0[j] = pa_x[j] * tpz_xyz_xyyy_0[j] + 0.5 * fl1_fx * tpz_yz_xyyy_0[j] + 0.5 * fl1_fx * tpz_xyz_yyy_0[j];

            tpx_xxyz_xyyz_0[j] =
                pa_x[j] * tpx_xyz_xyyz_0[j] + 0.5 * fl1_fx * tpx_yz_xyyz_0[j] + 0.5 * fl1_fx * tpx_xyz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xyyz_0[j];

            tpy_xxyz_xyyz_0[j] = pa_x[j] * tpy_xyz_xyyz_0[j] + 0.5 * fl1_fx * tpy_yz_xyyz_0[j] + 0.5 * fl1_fx * tpy_xyz_yyz_0[j];

            tpz_xxyz_xyyz_0[j] = pa_x[j] * tpz_xyz_xyyz_0[j] + 0.5 * fl1_fx * tpz_yz_xyyz_0[j] + 0.5 * fl1_fx * tpz_xyz_yyz_0[j];

            tpx_xxyz_xyzz_0[j] =
                pa_x[j] * tpx_xyz_xyzz_0[j] + 0.5 * fl1_fx * tpx_yz_xyzz_0[j] + 0.5 * fl1_fx * tpx_xyz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xyzz_0[j];

            tpy_xxyz_xyzz_0[j] = pa_x[j] * tpy_xyz_xyzz_0[j] + 0.5 * fl1_fx * tpy_yz_xyzz_0[j] + 0.5 * fl1_fx * tpy_xyz_yzz_0[j];

            tpz_xxyz_xyzz_0[j] = pa_x[j] * tpz_xyz_xyzz_0[j] + 0.5 * fl1_fx * tpz_yz_xyzz_0[j] + 0.5 * fl1_fx * tpz_xyz_yzz_0[j];

            tpx_xxyz_xzzz_0[j] =
                pa_x[j] * tpx_xyz_xzzz_0[j] + 0.5 * fl1_fx * tpx_yz_xzzz_0[j] + 0.5 * fl1_fx * tpx_xyz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xzzz_0[j];

            tpy_xxyz_xzzz_0[j] = pa_x[j] * tpy_xyz_xzzz_0[j] + 0.5 * fl1_fx * tpy_yz_xzzz_0[j] + 0.5 * fl1_fx * tpy_xyz_zzz_0[j];

            tpz_xxyz_xzzz_0[j] = pa_x[j] * tpz_xyz_xzzz_0[j] + 0.5 * fl1_fx * tpz_yz_xzzz_0[j] + 0.5 * fl1_fx * tpz_xyz_zzz_0[j];

            tpx_xxyz_yyyy_0[j] = pa_x[j] * tpx_xyz_yyyy_0[j] + 0.5 * fl1_fx * tpx_yz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xyz_yyyy_0[j];

            tpy_xxyz_yyyy_0[j] = pa_x[j] * tpy_xyz_yyyy_0[j] + 0.5 * fl1_fx * tpy_yz_yyyy_0[j];

            tpz_xxyz_yyyy_0[j] = pa_x[j] * tpz_xyz_yyyy_0[j] + 0.5 * fl1_fx * tpz_yz_yyyy_0[j];

            tpx_xxyz_yyyz_0[j] = pa_x[j] * tpx_xyz_yyyz_0[j] + 0.5 * fl1_fx * tpx_yz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xyz_yyyz_0[j];

            tpy_xxyz_yyyz_0[j] = pa_x[j] * tpy_xyz_yyyz_0[j] + 0.5 * fl1_fx * tpy_yz_yyyz_0[j];

            tpz_xxyz_yyyz_0[j] = pa_x[j] * tpz_xyz_yyyz_0[j] + 0.5 * fl1_fx * tpz_yz_yyyz_0[j];

            tpx_xxyz_yyzz_0[j] = pa_x[j] * tpx_xyz_yyzz_0[j] + 0.5 * fl1_fx * tpx_yz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_yyzz_0[j];

            tpy_xxyz_yyzz_0[j] = pa_x[j] * tpy_xyz_yyzz_0[j] + 0.5 * fl1_fx * tpy_yz_yyzz_0[j];

            tpz_xxyz_yyzz_0[j] = pa_x[j] * tpz_xyz_yyzz_0[j] + 0.5 * fl1_fx * tpz_yz_yyzz_0[j];

            tpx_xxyz_yzzz_0[j] = pa_x[j] * tpx_xyz_yzzz_0[j] + 0.5 * fl1_fx * tpx_yz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_yzzz_0[j];

            tpy_xxyz_yzzz_0[j] = pa_x[j] * tpy_xyz_yzzz_0[j] + 0.5 * fl1_fx * tpy_yz_yzzz_0[j];

            tpz_xxyz_yzzz_0[j] = pa_x[j] * tpz_xyz_yzzz_0[j] + 0.5 * fl1_fx * tpz_yz_yzzz_0[j];

            tpx_xxyz_zzzz_0[j] = pa_x[j] * tpx_xyz_zzzz_0[j] + 0.5 * fl1_fx * tpx_yz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_zzzz_0[j];

            tpy_xxyz_zzzz_0[j] = pa_x[j] * tpy_xyz_zzzz_0[j] + 0.5 * fl1_fx * tpy_yz_zzzz_0[j];

            tpz_xxyz_zzzz_0[j] = pa_x[j] * tpz_xyz_zzzz_0[j] + 0.5 * fl1_fx * tpz_yz_zzzz_0[j];

            tpx_xxzz_xxxx_0[j] =
                pa_x[j] * tpx_xzz_xxxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxxx_0[j] + 2.0 * fl1_fx * tpx_xzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxxx_0[j];

            tpy_xxzz_xxxx_0[j] = pa_x[j] * tpy_xzz_xxxx_0[j] + 0.5 * fl1_fx * tpy_zz_xxxx_0[j] + 2.0 * fl1_fx * tpy_xzz_xxx_0[j];

            tpz_xxzz_xxxx_0[j] = pa_x[j] * tpz_xzz_xxxx_0[j] + 0.5 * fl1_fx * tpz_zz_xxxx_0[j] + 2.0 * fl1_fx * tpz_xzz_xxx_0[j];

            tpx_xxzz_xxxy_0[j] =
                pa_x[j] * tpx_xzz_xxxy_0[j] + 0.5 * fl1_fx * tpx_zz_xxxy_0[j] + 1.5 * fl1_fx * tpx_xzz_xxy_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxxy_0[j];

            tpy_xxzz_xxxy_0[j] = pa_x[j] * tpy_xzz_xxxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxxy_0[j] + 1.5 * fl1_fx * tpy_xzz_xxy_0[j];

            tpz_xxzz_xxxy_0[j] = pa_x[j] * tpz_xzz_xxxy_0[j] + 0.5 * fl1_fx * tpz_zz_xxxy_0[j] + 1.5 * fl1_fx * tpz_xzz_xxy_0[j];

            tpx_xxzz_xxxz_0[j] =
                pa_x[j] * tpx_xzz_xxxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxxz_0[j] + 1.5 * fl1_fx * tpx_xzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxxz_0[j];

            tpy_xxzz_xxxz_0[j] = pa_x[j] * tpy_xzz_xxxz_0[j] + 0.5 * fl1_fx * tpy_zz_xxxz_0[j] + 1.5 * fl1_fx * tpy_xzz_xxz_0[j];

            tpz_xxzz_xxxz_0[j] = pa_x[j] * tpz_xzz_xxxz_0[j] + 0.5 * fl1_fx * tpz_zz_xxxz_0[j] + 1.5 * fl1_fx * tpz_xzz_xxz_0[j];

            tpx_xxzz_xxyy_0[j] =
                pa_x[j] * tpx_xzz_xxyy_0[j] + 0.5 * fl1_fx * tpx_zz_xxyy_0[j] + fl1_fx * tpx_xzz_xyy_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxyy_0[j];

            tpy_xxzz_xxyy_0[j] = pa_x[j] * tpy_xzz_xxyy_0[j] + 0.5 * fl1_fx * tpy_zz_xxyy_0[j] + fl1_fx * tpy_xzz_xyy_0[j];

            tpz_xxzz_xxyy_0[j] = pa_x[j] * tpz_xzz_xxyy_0[j] + 0.5 * fl1_fx * tpz_zz_xxyy_0[j] + fl1_fx * tpz_xzz_xyy_0[j];

            tpx_xxzz_xxyz_0[j] =
                pa_x[j] * tpx_xzz_xxyz_0[j] + 0.5 * fl1_fx * tpx_zz_xxyz_0[j] + fl1_fx * tpx_xzz_xyz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxyz_0[j];

            tpy_xxzz_xxyz_0[j] = pa_x[j] * tpy_xzz_xxyz_0[j] + 0.5 * fl1_fx * tpy_zz_xxyz_0[j] + fl1_fx * tpy_xzz_xyz_0[j];

            tpz_xxzz_xxyz_0[j] = pa_x[j] * tpz_xzz_xxyz_0[j] + 0.5 * fl1_fx * tpz_zz_xxyz_0[j] + fl1_fx * tpz_xzz_xyz_0[j];

            tpx_xxzz_xxzz_0[j] =
                pa_x[j] * tpx_xzz_xxzz_0[j] + 0.5 * fl1_fx * tpx_zz_xxzz_0[j] + fl1_fx * tpx_xzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxzz_0[j];

            tpy_xxzz_xxzz_0[j] = pa_x[j] * tpy_xzz_xxzz_0[j] + 0.5 * fl1_fx * tpy_zz_xxzz_0[j] + fl1_fx * tpy_xzz_xzz_0[j];

            tpz_xxzz_xxzz_0[j] = pa_x[j] * tpz_xzz_xxzz_0[j] + 0.5 * fl1_fx * tpz_zz_xxzz_0[j] + fl1_fx * tpz_xzz_xzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_243_291(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (243,291)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 81);

        auto tpy_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tpz_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tpx_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 82);

        auto tpy_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tpz_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tpx_xzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 83);

        auto tpy_xzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 83);

        auto tpz_xzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 83);

        auto tpx_xzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 84);

        auto tpy_xzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 84);

        auto tpz_xzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 84);

        auto tpx_xzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 85);

        auto tpy_xzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 85);

        auto tpz_xzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 85);

        auto tpx_xzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 86);

        auto tpy_xzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 86);

        auto tpz_xzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 86);

        auto tpx_xzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 87);

        auto tpy_xzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 87);

        auto tpz_xzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 87);

        auto tpx_xzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 88);

        auto tpy_xzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 88);

        auto tpz_xzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 88);

        auto tpx_xzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 89);

        auto tpy_xzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 89);

        auto tpz_xzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 89);

        auto tpx_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 90);

        auto tpy_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tpz_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tpx_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 91);

        auto tpy_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tpz_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tpx_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 92);

        auto tpy_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tpz_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tpx_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 93);

        auto tpy_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tpz_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tpx_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 94);

        auto tpy_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tpz_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tpx_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 95);

        auto tpy_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tpz_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tpx_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 96);

        auto tpy_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tpz_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpy_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tpz_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tpx_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 84);

        auto tpy_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tpz_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tpx_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 85);

        auto tpy_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tpz_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tpx_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 86);

        auto tpy_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tpz_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tpx_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 87);

        auto tpy_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tpz_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tpx_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 88);

        auto tpy_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tpz_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tpx_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 89);

        auto tpy_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tpx_xzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 56);

        auto tpy_xzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 56);

        auto tpz_xzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 56);

        auto tpx_xzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 57);

        auto tpy_xzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 57);

        auto tpz_xzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 57);

        auto tpx_xzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 58);

        auto tpy_xzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 58);

        auto tpz_xzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 58);

        auto tpx_xzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 59);

        auto tpy_xzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 59);

        auto tpz_xzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 59);

        auto tpx_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 60);

        auto tpy_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tpz_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tpx_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 61);

        auto tpy_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tpz_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tpx_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 62);

        auto tpy_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tpz_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tpx_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 63);

        auto tpy_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tpz_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tpx_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 64);

        auto tpy_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tpz_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tpx_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 65);

        auto tpy_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tpz_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tpx_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 66);

        auto tpy_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tpz_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto ts_xzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 81);

        auto ts_xzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 82);

        auto ts_xzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 83);

        auto ts_xzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 84);

        auto ts_xzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 85);

        auto ts_xzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 86);

        auto ts_xzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 87);

        auto ts_xzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 88);

        auto ts_xzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 89);

        auto ts_yyy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 90);

        auto ts_yyy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 91);

        auto ts_yyy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 92);

        auto ts_yyy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 93);

        auto ts_yyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 94);

        auto ts_yyy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 95);

        auto ts_yyy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 96);

        // set up pointers to integrals

        auto tpx_xxzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 81);

        auto tpy_xxzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 81);

        auto tpz_xxzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 81);

        auto tpx_xxzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 82);

        auto tpy_xxzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 82);

        auto tpz_xxzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 82);

        auto tpx_xxzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 83);

        auto tpy_xxzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 83);

        auto tpz_xxzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 83);

        auto tpx_xxzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 84);

        auto tpy_xxzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 84);

        auto tpz_xxzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 84);

        auto tpx_xxzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 85);

        auto tpy_xxzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 85);

        auto tpz_xxzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 85);

        auto tpx_xxzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 86);

        auto tpy_xxzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 86);

        auto tpz_xxzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 86);

        auto tpx_xxzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 87);

        auto tpy_xxzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 87);

        auto tpz_xxzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 87);

        auto tpx_xxzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 88);

        auto tpy_xxzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 88);

        auto tpz_xxzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 88);

        auto tpx_xxzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 89);

        auto tpy_xxzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 89);

        auto tpz_xxzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 89);

        auto tpx_xyyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 90);

        auto tpy_xyyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 90);

        auto tpz_xyyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 90);

        auto tpx_xyyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 91);

        auto tpy_xyyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 91);

        auto tpz_xyyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 91);

        auto tpx_xyyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 92);

        auto tpy_xyyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 92);

        auto tpz_xyyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 92);

        auto tpx_xyyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 93);

        auto tpy_xyyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 93);

        auto tpz_xyyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 93);

        auto tpx_xyyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 94);

        auto tpy_xyyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 94);

        auto tpz_xyyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 94);

        auto tpx_xyyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 95);

        auto tpy_xyyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 95);

        auto tpz_xyyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 95);

        auto tpx_xyyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 96);

        auto tpy_xyyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 96);

        auto tpz_xyyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 96);

        // Batch of Integrals (243,291)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxzz_xyyy_0, tpx_xxzz_xyyz_0, tpx_xxzz_xyzz_0, \
                                     tpx_xxzz_xzzz_0, tpx_xxzz_yyyy_0, tpx_xxzz_yyyz_0, tpx_xxzz_yyzz_0, tpx_xxzz_yzzz_0, \
                                     tpx_xxzz_zzzz_0, tpx_xyyy_xxxx_0, tpx_xyyy_xxxy_0, tpx_xyyy_xxxz_0, tpx_xyyy_xxyy_0, \
                                     tpx_xyyy_xxyz_0, tpx_xyyy_xxzz_0, tpx_xyyy_xyyy_0, tpx_xzz_xyyy_0, tpx_xzz_xyyz_0, \
                                     tpx_xzz_xyzz_0, tpx_xzz_xzzz_0, tpx_xzz_yyy_0, tpx_xzz_yyyy_0, tpx_xzz_yyyz_0, \
                                     tpx_xzz_yyz_0, tpx_xzz_yyzz_0, tpx_xzz_yzz_0, tpx_xzz_yzzz_0, tpx_xzz_zzz_0, \
                                     tpx_xzz_zzzz_0, tpx_yyy_xxx_0, tpx_yyy_xxxx_0, tpx_yyy_xxxy_0, tpx_yyy_xxxz_0, \
                                     tpx_yyy_xxy_0, tpx_yyy_xxyy_0, tpx_yyy_xxyz_0, tpx_yyy_xxz_0, tpx_yyy_xxzz_0, \
                                     tpx_yyy_xyy_0, tpx_yyy_xyyy_0, tpx_yyy_xyz_0, tpx_yyy_xzz_0, tpx_yyy_yyy_0, \
                                     tpx_zz_xyyy_0, tpx_zz_xyyz_0, tpx_zz_xyzz_0, tpx_zz_xzzz_0, tpx_zz_yyyy_0, \
                                     tpx_zz_yyyz_0, tpx_zz_yyzz_0, tpx_zz_yzzz_0, tpx_zz_zzzz_0, tpy_xxzz_xyyy_0, \
                                     tpy_xxzz_xyyz_0, tpy_xxzz_xyzz_0, tpy_xxzz_xzzz_0, tpy_xxzz_yyyy_0, tpy_xxzz_yyyz_0, \
                                     tpy_xxzz_yyzz_0, tpy_xxzz_yzzz_0, tpy_xxzz_zzzz_0, tpy_xyyy_xxxx_0, tpy_xyyy_xxxy_0, \
                                     tpy_xyyy_xxxz_0, tpy_xyyy_xxyy_0, tpy_xyyy_xxyz_0, tpy_xyyy_xxzz_0, tpy_xyyy_xyyy_0, \
                                     tpy_xzz_xyyy_0, tpy_xzz_xyyz_0, tpy_xzz_xyzz_0, tpy_xzz_xzzz_0, tpy_xzz_yyy_0, \
                                     tpy_xzz_yyyy_0, tpy_xzz_yyyz_0, tpy_xzz_yyz_0, tpy_xzz_yyzz_0, tpy_xzz_yzz_0, \
                                     tpy_xzz_yzzz_0, tpy_xzz_zzz_0, tpy_xzz_zzzz_0, tpy_yyy_xxx_0, tpy_yyy_xxxx_0, \
                                     tpy_yyy_xxxy_0, tpy_yyy_xxxz_0, tpy_yyy_xxy_0, tpy_yyy_xxyy_0, tpy_yyy_xxyz_0, \
                                     tpy_yyy_xxz_0, tpy_yyy_xxzz_0, tpy_yyy_xyy_0, tpy_yyy_xyyy_0, tpy_yyy_xyz_0, \
                                     tpy_yyy_xzz_0, tpy_yyy_yyy_0, tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpy_zz_xyzz_0, \
                                     tpy_zz_xzzz_0, tpy_zz_yyyy_0, tpy_zz_yyyz_0, tpy_zz_yyzz_0, tpy_zz_yzzz_0, \
                                     tpy_zz_zzzz_0, tpz_xxzz_xyyy_0, tpz_xxzz_xyyz_0, tpz_xxzz_xyzz_0, tpz_xxzz_xzzz_0, \
                                     tpz_xxzz_yyyy_0, tpz_xxzz_yyyz_0, tpz_xxzz_yyzz_0, tpz_xxzz_yzzz_0, tpz_xxzz_zzzz_0, \
                                     tpz_xyyy_xxxx_0, tpz_xyyy_xxxy_0, tpz_xyyy_xxxz_0, tpz_xyyy_xxyy_0, tpz_xyyy_xxyz_0, \
                                     tpz_xyyy_xxzz_0, tpz_xyyy_xyyy_0, tpz_xzz_xyyy_0, tpz_xzz_xyyz_0, tpz_xzz_xyzz_0, \
                                     tpz_xzz_xzzz_0, tpz_xzz_yyy_0, tpz_xzz_yyyy_0, tpz_xzz_yyyz_0, tpz_xzz_yyz_0, \
                                     tpz_xzz_yyzz_0, tpz_xzz_yzz_0, tpz_xzz_yzzz_0, tpz_xzz_zzz_0, tpz_xzz_zzzz_0, \
                                     tpz_yyy_xxx_0, tpz_yyy_xxxx_0, tpz_yyy_xxxy_0, tpz_yyy_xxxz_0, tpz_yyy_xxy_0, \
                                     tpz_yyy_xxyy_0, tpz_yyy_xxyz_0, tpz_yyy_xxz_0, tpz_yyy_xxzz_0, tpz_yyy_xyy_0, \
                                     tpz_yyy_xyyy_0, tpz_yyy_xyz_0, tpz_yyy_xzz_0, tpz_yyy_yyy_0, tpz_zz_xyyy_0, \
                                     tpz_zz_xyyz_0, tpz_zz_xyzz_0, tpz_zz_xzzz_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, \
                                     tpz_zz_yyzz_0, tpz_zz_yzzz_0, tpz_zz_zzzz_0, ts_xzz_xyyy_0, ts_xzz_xyyz_0, \
                                     ts_xzz_xyzz_0, ts_xzz_xzzz_0, ts_xzz_yyyy_0, ts_xzz_yyyz_0, ts_xzz_yyzz_0, \
                                     ts_xzz_yzzz_0, ts_xzz_zzzz_0, ts_yyy_xxxx_0, ts_yyy_xxxy_0, ts_yyy_xxxz_0, \
                                     ts_yyy_xxyy_0, ts_yyy_xxyz_0, ts_yyy_xxzz_0, ts_yyy_xyyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxzz_xyyy_0[j] =
                pa_x[j] * tpx_xzz_xyyy_0[j] + 0.5 * fl1_fx * tpx_zz_xyyy_0[j] + 0.5 * fl1_fx * tpx_xzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xzz_xyyy_0[j];

            tpy_xxzz_xyyy_0[j] = pa_x[j] * tpy_xzz_xyyy_0[j] + 0.5 * fl1_fx * tpy_zz_xyyy_0[j] + 0.5 * fl1_fx * tpy_xzz_yyy_0[j];

            tpz_xxzz_xyyy_0[j] = pa_x[j] * tpz_xzz_xyyy_0[j] + 0.5 * fl1_fx * tpz_zz_xyyy_0[j] + 0.5 * fl1_fx * tpz_xzz_yyy_0[j];

            tpx_xxzz_xyyz_0[j] =
                pa_x[j] * tpx_xzz_xyyz_0[j] + 0.5 * fl1_fx * tpx_zz_xyyz_0[j] + 0.5 * fl1_fx * tpx_xzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xyyz_0[j];

            tpy_xxzz_xyyz_0[j] = pa_x[j] * tpy_xzz_xyyz_0[j] + 0.5 * fl1_fx * tpy_zz_xyyz_0[j] + 0.5 * fl1_fx * tpy_xzz_yyz_0[j];

            tpz_xxzz_xyyz_0[j] = pa_x[j] * tpz_xzz_xyyz_0[j] + 0.5 * fl1_fx * tpz_zz_xyyz_0[j] + 0.5 * fl1_fx * tpz_xzz_yyz_0[j];

            tpx_xxzz_xyzz_0[j] =
                pa_x[j] * tpx_xzz_xyzz_0[j] + 0.5 * fl1_fx * tpx_zz_xyzz_0[j] + 0.5 * fl1_fx * tpx_xzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xyzz_0[j];

            tpy_xxzz_xyzz_0[j] = pa_x[j] * tpy_xzz_xyzz_0[j] + 0.5 * fl1_fx * tpy_zz_xyzz_0[j] + 0.5 * fl1_fx * tpy_xzz_yzz_0[j];

            tpz_xxzz_xyzz_0[j] = pa_x[j] * tpz_xzz_xyzz_0[j] + 0.5 * fl1_fx * tpz_zz_xyzz_0[j] + 0.5 * fl1_fx * tpz_xzz_yzz_0[j];

            tpx_xxzz_xzzz_0[j] =
                pa_x[j] * tpx_xzz_xzzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzzz_0[j] + 0.5 * fl1_fx * tpx_xzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xzzz_0[j];

            tpy_xxzz_xzzz_0[j] = pa_x[j] * tpy_xzz_xzzz_0[j] + 0.5 * fl1_fx * tpy_zz_xzzz_0[j] + 0.5 * fl1_fx * tpy_xzz_zzz_0[j];

            tpz_xxzz_xzzz_0[j] = pa_x[j] * tpz_xzz_xzzz_0[j] + 0.5 * fl1_fx * tpz_zz_xzzz_0[j] + 0.5 * fl1_fx * tpz_xzz_zzz_0[j];

            tpx_xxzz_yyyy_0[j] = pa_x[j] * tpx_xzz_yyyy_0[j] + 0.5 * fl1_fx * tpx_zz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xzz_yyyy_0[j];

            tpy_xxzz_yyyy_0[j] = pa_x[j] * tpy_xzz_yyyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyyy_0[j];

            tpz_xxzz_yyyy_0[j] = pa_x[j] * tpz_xzz_yyyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyyy_0[j];

            tpx_xxzz_yyyz_0[j] = pa_x[j] * tpx_xzz_yyyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xzz_yyyz_0[j];

            tpy_xxzz_yyyz_0[j] = pa_x[j] * tpy_xzz_yyyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyyz_0[j];

            tpz_xxzz_yyyz_0[j] = pa_x[j] * tpz_xzz_yyyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyyz_0[j];

            tpx_xxzz_yyzz_0[j] = pa_x[j] * tpx_xzz_yyzz_0[j] + 0.5 * fl1_fx * tpx_zz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_yyzz_0[j];

            tpy_xxzz_yyzz_0[j] = pa_x[j] * tpy_xzz_yyzz_0[j] + 0.5 * fl1_fx * tpy_zz_yyzz_0[j];

            tpz_xxzz_yyzz_0[j] = pa_x[j] * tpz_xzz_yyzz_0[j] + 0.5 * fl1_fx * tpz_zz_yyzz_0[j];

            tpx_xxzz_yzzz_0[j] = pa_x[j] * tpx_xzz_yzzz_0[j] + 0.5 * fl1_fx * tpx_zz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_yzzz_0[j];

            tpy_xxzz_yzzz_0[j] = pa_x[j] * tpy_xzz_yzzz_0[j] + 0.5 * fl1_fx * tpy_zz_yzzz_0[j];

            tpz_xxzz_yzzz_0[j] = pa_x[j] * tpz_xzz_yzzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzzz_0[j];

            tpx_xxzz_zzzz_0[j] = pa_x[j] * tpx_xzz_zzzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_zzzz_0[j];

            tpy_xxzz_zzzz_0[j] = pa_x[j] * tpy_xzz_zzzz_0[j] + 0.5 * fl1_fx * tpy_zz_zzzz_0[j];

            tpz_xxzz_zzzz_0[j] = pa_x[j] * tpz_xzz_zzzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzzz_0[j];

            tpx_xyyy_xxxx_0[j] = pa_x[j] * tpx_yyy_xxxx_0[j] + 2.0 * fl1_fx * tpx_yyy_xxx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxxx_0[j];

            tpy_xyyy_xxxx_0[j] = pa_x[j] * tpy_yyy_xxxx_0[j] + 2.0 * fl1_fx * tpy_yyy_xxx_0[j];

            tpz_xyyy_xxxx_0[j] = pa_x[j] * tpz_yyy_xxxx_0[j] + 2.0 * fl1_fx * tpz_yyy_xxx_0[j];

            tpx_xyyy_xxxy_0[j] = pa_x[j] * tpx_yyy_xxxy_0[j] + 1.5 * fl1_fx * tpx_yyy_xxy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxxy_0[j];

            tpy_xyyy_xxxy_0[j] = pa_x[j] * tpy_yyy_xxxy_0[j] + 1.5 * fl1_fx * tpy_yyy_xxy_0[j];

            tpz_xyyy_xxxy_0[j] = pa_x[j] * tpz_yyy_xxxy_0[j] + 1.5 * fl1_fx * tpz_yyy_xxy_0[j];

            tpx_xyyy_xxxz_0[j] = pa_x[j] * tpx_yyy_xxxz_0[j] + 1.5 * fl1_fx * tpx_yyy_xxz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxxz_0[j];

            tpy_xyyy_xxxz_0[j] = pa_x[j] * tpy_yyy_xxxz_0[j] + 1.5 * fl1_fx * tpy_yyy_xxz_0[j];

            tpz_xyyy_xxxz_0[j] = pa_x[j] * tpz_yyy_xxxz_0[j] + 1.5 * fl1_fx * tpz_yyy_xxz_0[j];

            tpx_xyyy_xxyy_0[j] = pa_x[j] * tpx_yyy_xxyy_0[j] + fl1_fx * tpx_yyy_xyy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxyy_0[j];

            tpy_xyyy_xxyy_0[j] = pa_x[j] * tpy_yyy_xxyy_0[j] + fl1_fx * tpy_yyy_xyy_0[j];

            tpz_xyyy_xxyy_0[j] = pa_x[j] * tpz_yyy_xxyy_0[j] + fl1_fx * tpz_yyy_xyy_0[j];

            tpx_xyyy_xxyz_0[j] = pa_x[j] * tpx_yyy_xxyz_0[j] + fl1_fx * tpx_yyy_xyz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxyz_0[j];

            tpy_xyyy_xxyz_0[j] = pa_x[j] * tpy_yyy_xxyz_0[j] + fl1_fx * tpy_yyy_xyz_0[j];

            tpz_xyyy_xxyz_0[j] = pa_x[j] * tpz_yyy_xxyz_0[j] + fl1_fx * tpz_yyy_xyz_0[j];

            tpx_xyyy_xxzz_0[j] = pa_x[j] * tpx_yyy_xxzz_0[j] + fl1_fx * tpx_yyy_xzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxzz_0[j];

            tpy_xyyy_xxzz_0[j] = pa_x[j] * tpy_yyy_xxzz_0[j] + fl1_fx * tpy_yyy_xzz_0[j];

            tpz_xyyy_xxzz_0[j] = pa_x[j] * tpz_yyy_xxzz_0[j] + fl1_fx * tpz_yyy_xzz_0[j];

            tpx_xyyy_xyyy_0[j] = pa_x[j] * tpx_yyy_xyyy_0[j] + 0.5 * fl1_fx * tpx_yyy_yyy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyyy_0[j];

            tpy_xyyy_xyyy_0[j] = pa_x[j] * tpy_yyy_xyyy_0[j] + 0.5 * fl1_fx * tpy_yyy_yyy_0[j];

            tpz_xyyy_xyyy_0[j] = pa_x[j] * tpz_yyy_xyyy_0[j] + 0.5 * fl1_fx * tpz_yyy_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_291_339(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (291,339)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 97);

        auto tpy_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tpz_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tpx_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 98);

        auto tpy_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tpz_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tpx_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 99);

        auto tpy_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tpz_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tpx_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 100);

        auto tpy_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tpz_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tpx_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 101);

        auto tpy_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tpz_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tpx_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 102);

        auto tpy_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tpz_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tpx_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 103);

        auto tpy_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tpz_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tpx_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 104);

        auto tpy_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tpz_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tpx_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 105);

        auto tpy_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tpz_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tpx_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 106);

        auto tpy_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tpz_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tpx_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 107);

        auto tpy_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tpz_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tpx_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 108);

        auto tpy_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tpz_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tpx_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 109);

        auto tpy_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tpz_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tpx_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 110);

        auto tpy_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tpz_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tpx_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 111);

        auto tpy_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tpz_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tpx_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 112);

        auto tpy_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tpz_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tpx_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 67);

        auto tpy_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tpz_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tpx_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 68);

        auto tpy_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tpz_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tpx_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 69);

        auto tpy_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tpz_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tpx_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 70);

        auto tpy_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tpz_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tpx_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 71);

        auto tpy_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tpz_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tpx_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 72);

        auto tpy_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tpz_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tpx_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 73);

        auto tpy_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tpz_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tpx_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 74);

        auto tpy_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tpz_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tpx_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 75);

        auto tpy_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tpz_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tpx_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 76);

        auto tpy_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tpz_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tpx_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 77);

        auto tpy_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tpz_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto ts_yyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 97);

        auto ts_yyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 98);

        auto ts_yyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 99);

        auto ts_yyy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 100);

        auto ts_yyy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 101);

        auto ts_yyy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 102);

        auto ts_yyy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 103);

        auto ts_yyy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 104);

        auto ts_yyz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 105);

        auto ts_yyz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 106);

        auto ts_yyz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 107);

        auto ts_yyz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 108);

        auto ts_yyz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 109);

        auto ts_yyz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 110);

        auto ts_yyz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 111);

        auto ts_yyz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 112);

        // set up pointers to integrals

        auto tpx_xyyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 97);

        auto tpy_xyyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 97);

        auto tpz_xyyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 97);

        auto tpx_xyyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 98);

        auto tpy_xyyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 98);

        auto tpz_xyyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 98);

        auto tpx_xyyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 99);

        auto tpy_xyyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 99);

        auto tpz_xyyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 99);

        auto tpx_xyyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 100);

        auto tpy_xyyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 100);

        auto tpz_xyyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 100);

        auto tpx_xyyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 101);

        auto tpy_xyyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 101);

        auto tpz_xyyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 101);

        auto tpx_xyyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 102);

        auto tpy_xyyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 102);

        auto tpz_xyyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 102);

        auto tpx_xyyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 103);

        auto tpy_xyyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 103);

        auto tpz_xyyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 103);

        auto tpx_xyyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 104);

        auto tpy_xyyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 104);

        auto tpz_xyyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 104);

        auto tpx_xyyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 105);

        auto tpy_xyyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 105);

        auto tpz_xyyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 105);

        auto tpx_xyyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 106);

        auto tpy_xyyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 106);

        auto tpz_xyyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 106);

        auto tpx_xyyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 107);

        auto tpy_xyyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 107);

        auto tpz_xyyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 107);

        auto tpx_xyyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 108);

        auto tpy_xyyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 108);

        auto tpz_xyyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 108);

        auto tpx_xyyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 109);

        auto tpy_xyyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 109);

        auto tpz_xyyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 109);

        auto tpx_xyyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 110);

        auto tpy_xyyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 110);

        auto tpz_xyyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 110);

        auto tpx_xyyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 111);

        auto tpy_xyyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 111);

        auto tpz_xyyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 111);

        auto tpx_xyyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 112);

        auto tpy_xyyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 112);

        auto tpz_xyyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 112);

        // Batch of Integrals (291,339)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyyy_xyyz_0, tpx_xyyy_xyzz_0, tpx_xyyy_xzzz_0, \
                                     tpx_xyyy_yyyy_0, tpx_xyyy_yyyz_0, tpx_xyyy_yyzz_0, tpx_xyyy_yzzz_0, tpx_xyyy_zzzz_0, \
                                     tpx_xyyz_xxxx_0, tpx_xyyz_xxxy_0, tpx_xyyz_xxxz_0, tpx_xyyz_xxyy_0, tpx_xyyz_xxyz_0, \
                                     tpx_xyyz_xxzz_0, tpx_xyyz_xyyy_0, tpx_xyyz_xyyz_0, tpx_yyy_xyyz_0, tpx_yyy_xyzz_0, \
                                     tpx_yyy_xzzz_0, tpx_yyy_yyyy_0, tpx_yyy_yyyz_0, tpx_yyy_yyz_0, tpx_yyy_yyzz_0, \
                                     tpx_yyy_yzz_0, tpx_yyy_yzzz_0, tpx_yyy_zzz_0, tpx_yyy_zzzz_0, tpx_yyz_xxx_0, \
                                     tpx_yyz_xxxx_0, tpx_yyz_xxxy_0, tpx_yyz_xxxz_0, tpx_yyz_xxy_0, tpx_yyz_xxyy_0, \
                                     tpx_yyz_xxyz_0, tpx_yyz_xxz_0, tpx_yyz_xxzz_0, tpx_yyz_xyy_0, tpx_yyz_xyyy_0, \
                                     tpx_yyz_xyyz_0, tpx_yyz_xyz_0, tpx_yyz_xzz_0, tpx_yyz_yyy_0, tpx_yyz_yyz_0, \
                                     tpy_xyyy_xyyz_0, tpy_xyyy_xyzz_0, tpy_xyyy_xzzz_0, tpy_xyyy_yyyy_0, tpy_xyyy_yyyz_0, \
                                     tpy_xyyy_yyzz_0, tpy_xyyy_yzzz_0, tpy_xyyy_zzzz_0, tpy_xyyz_xxxx_0, tpy_xyyz_xxxy_0, \
                                     tpy_xyyz_xxxz_0, tpy_xyyz_xxyy_0, tpy_xyyz_xxyz_0, tpy_xyyz_xxzz_0, tpy_xyyz_xyyy_0, \
                                     tpy_xyyz_xyyz_0, tpy_yyy_xyyz_0, tpy_yyy_xyzz_0, tpy_yyy_xzzz_0, tpy_yyy_yyyy_0, \
                                     tpy_yyy_yyyz_0, tpy_yyy_yyz_0, tpy_yyy_yyzz_0, tpy_yyy_yzz_0, tpy_yyy_yzzz_0, \
                                     tpy_yyy_zzz_0, tpy_yyy_zzzz_0, tpy_yyz_xxx_0, tpy_yyz_xxxx_0, tpy_yyz_xxxy_0, \
                                     tpy_yyz_xxxz_0, tpy_yyz_xxy_0, tpy_yyz_xxyy_0, tpy_yyz_xxyz_0, tpy_yyz_xxz_0, \
                                     tpy_yyz_xxzz_0, tpy_yyz_xyy_0, tpy_yyz_xyyy_0, tpy_yyz_xyyz_0, tpy_yyz_xyz_0, \
                                     tpy_yyz_xzz_0, tpy_yyz_yyy_0, tpy_yyz_yyz_0, tpz_xyyy_xyyz_0, tpz_xyyy_xyzz_0, \
                                     tpz_xyyy_xzzz_0, tpz_xyyy_yyyy_0, tpz_xyyy_yyyz_0, tpz_xyyy_yyzz_0, tpz_xyyy_yzzz_0, \
                                     tpz_xyyy_zzzz_0, tpz_xyyz_xxxx_0, tpz_xyyz_xxxy_0, tpz_xyyz_xxxz_0, tpz_xyyz_xxyy_0, \
                                     tpz_xyyz_xxyz_0, tpz_xyyz_xxzz_0, tpz_xyyz_xyyy_0, tpz_xyyz_xyyz_0, tpz_yyy_xyyz_0, \
                                     tpz_yyy_xyzz_0, tpz_yyy_xzzz_0, tpz_yyy_yyyy_0, tpz_yyy_yyyz_0, tpz_yyy_yyz_0, \
                                     tpz_yyy_yyzz_0, tpz_yyy_yzz_0, tpz_yyy_yzzz_0, tpz_yyy_zzz_0, tpz_yyy_zzzz_0, \
                                     tpz_yyz_xxx_0, tpz_yyz_xxxx_0, tpz_yyz_xxxy_0, tpz_yyz_xxxz_0, tpz_yyz_xxy_0, \
                                     tpz_yyz_xxyy_0, tpz_yyz_xxyz_0, tpz_yyz_xxz_0, tpz_yyz_xxzz_0, tpz_yyz_xyy_0, \
                                     tpz_yyz_xyyy_0, tpz_yyz_xyyz_0, tpz_yyz_xyz_0, tpz_yyz_xzz_0, tpz_yyz_yyy_0, \
                                     tpz_yyz_yyz_0, ts_yyy_xyyz_0, ts_yyy_xyzz_0, ts_yyy_xzzz_0, ts_yyy_yyyy_0, \
                                     ts_yyy_yyyz_0, ts_yyy_yyzz_0, ts_yyy_yzzz_0, ts_yyy_zzzz_0, ts_yyz_xxxx_0, \
                                     ts_yyz_xxxy_0, ts_yyz_xxxz_0, ts_yyz_xxyy_0, ts_yyz_xxyz_0, ts_yyz_xxzz_0, \
                                     ts_yyz_xyyy_0, ts_yyz_xyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xyyy_xyyz_0[j] = pa_x[j] * tpx_yyy_xyyz_0[j] + 0.5 * fl1_fx * tpx_yyy_yyz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyyz_0[j];

            tpy_xyyy_xyyz_0[j] = pa_x[j] * tpy_yyy_xyyz_0[j] + 0.5 * fl1_fx * tpy_yyy_yyz_0[j];

            tpz_xyyy_xyyz_0[j] = pa_x[j] * tpz_yyy_xyyz_0[j] + 0.5 * fl1_fx * tpz_yyy_yyz_0[j];

            tpx_xyyy_xyzz_0[j] = pa_x[j] * tpx_yyy_xyzz_0[j] + 0.5 * fl1_fx * tpx_yyy_yzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyzz_0[j];

            tpy_xyyy_xyzz_0[j] = pa_x[j] * tpy_yyy_xyzz_0[j] + 0.5 * fl1_fx * tpy_yyy_yzz_0[j];

            tpz_xyyy_xyzz_0[j] = pa_x[j] * tpz_yyy_xyzz_0[j] + 0.5 * fl1_fx * tpz_yyy_yzz_0[j];

            tpx_xyyy_xzzz_0[j] = pa_x[j] * tpx_yyy_xzzz_0[j] + 0.5 * fl1_fx * tpx_yyy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xzzz_0[j];

            tpy_xyyy_xzzz_0[j] = pa_x[j] * tpy_yyy_xzzz_0[j] + 0.5 * fl1_fx * tpy_yyy_zzz_0[j];

            tpz_xyyy_xzzz_0[j] = pa_x[j] * tpz_yyy_xzzz_0[j] + 0.5 * fl1_fx * tpz_yyy_zzz_0[j];

            tpx_xyyy_yyyy_0[j] = pa_x[j] * tpx_yyy_yyyy_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyyy_0[j];

            tpy_xyyy_yyyy_0[j] = pa_x[j] * tpy_yyy_yyyy_0[j];

            tpz_xyyy_yyyy_0[j] = pa_x[j] * tpz_yyy_yyyy_0[j];

            tpx_xyyy_yyyz_0[j] = pa_x[j] * tpx_yyy_yyyz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyyz_0[j];

            tpy_xyyy_yyyz_0[j] = pa_x[j] * tpy_yyy_yyyz_0[j];

            tpz_xyyy_yyyz_0[j] = pa_x[j] * tpz_yyy_yyyz_0[j];

            tpx_xyyy_yyzz_0[j] = pa_x[j] * tpx_yyy_yyzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyzz_0[j];

            tpy_xyyy_yyzz_0[j] = pa_x[j] * tpy_yyy_yyzz_0[j];

            tpz_xyyy_yyzz_0[j] = pa_x[j] * tpz_yyy_yyzz_0[j];

            tpx_xyyy_yzzz_0[j] = pa_x[j] * tpx_yyy_yzzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yzzz_0[j];

            tpy_xyyy_yzzz_0[j] = pa_x[j] * tpy_yyy_yzzz_0[j];

            tpz_xyyy_yzzz_0[j] = pa_x[j] * tpz_yyy_yzzz_0[j];

            tpx_xyyy_zzzz_0[j] = pa_x[j] * tpx_yyy_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_zzzz_0[j];

            tpy_xyyy_zzzz_0[j] = pa_x[j] * tpy_yyy_zzzz_0[j];

            tpz_xyyy_zzzz_0[j] = pa_x[j] * tpz_yyy_zzzz_0[j];

            tpx_xyyz_xxxx_0[j] = pa_x[j] * tpx_yyz_xxxx_0[j] + 2.0 * fl1_fx * tpx_yyz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxxx_0[j];

            tpy_xyyz_xxxx_0[j] = pa_x[j] * tpy_yyz_xxxx_0[j] + 2.0 * fl1_fx * tpy_yyz_xxx_0[j];

            tpz_xyyz_xxxx_0[j] = pa_x[j] * tpz_yyz_xxxx_0[j] + 2.0 * fl1_fx * tpz_yyz_xxx_0[j];

            tpx_xyyz_xxxy_0[j] = pa_x[j] * tpx_yyz_xxxy_0[j] + 1.5 * fl1_fx * tpx_yyz_xxy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxxy_0[j];

            tpy_xyyz_xxxy_0[j] = pa_x[j] * tpy_yyz_xxxy_0[j] + 1.5 * fl1_fx * tpy_yyz_xxy_0[j];

            tpz_xyyz_xxxy_0[j] = pa_x[j] * tpz_yyz_xxxy_0[j] + 1.5 * fl1_fx * tpz_yyz_xxy_0[j];

            tpx_xyyz_xxxz_0[j] = pa_x[j] * tpx_yyz_xxxz_0[j] + 1.5 * fl1_fx * tpx_yyz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxxz_0[j];

            tpy_xyyz_xxxz_0[j] = pa_x[j] * tpy_yyz_xxxz_0[j] + 1.5 * fl1_fx * tpy_yyz_xxz_0[j];

            tpz_xyyz_xxxz_0[j] = pa_x[j] * tpz_yyz_xxxz_0[j] + 1.5 * fl1_fx * tpz_yyz_xxz_0[j];

            tpx_xyyz_xxyy_0[j] = pa_x[j] * tpx_yyz_xxyy_0[j] + fl1_fx * tpx_yyz_xyy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxyy_0[j];

            tpy_xyyz_xxyy_0[j] = pa_x[j] * tpy_yyz_xxyy_0[j] + fl1_fx * tpy_yyz_xyy_0[j];

            tpz_xyyz_xxyy_0[j] = pa_x[j] * tpz_yyz_xxyy_0[j] + fl1_fx * tpz_yyz_xyy_0[j];

            tpx_xyyz_xxyz_0[j] = pa_x[j] * tpx_yyz_xxyz_0[j] + fl1_fx * tpx_yyz_xyz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxyz_0[j];

            tpy_xyyz_xxyz_0[j] = pa_x[j] * tpy_yyz_xxyz_0[j] + fl1_fx * tpy_yyz_xyz_0[j];

            tpz_xyyz_xxyz_0[j] = pa_x[j] * tpz_yyz_xxyz_0[j] + fl1_fx * tpz_yyz_xyz_0[j];

            tpx_xyyz_xxzz_0[j] = pa_x[j] * tpx_yyz_xxzz_0[j] + fl1_fx * tpx_yyz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxzz_0[j];

            tpy_xyyz_xxzz_0[j] = pa_x[j] * tpy_yyz_xxzz_0[j] + fl1_fx * tpy_yyz_xzz_0[j];

            tpz_xyyz_xxzz_0[j] = pa_x[j] * tpz_yyz_xxzz_0[j] + fl1_fx * tpz_yyz_xzz_0[j];

            tpx_xyyz_xyyy_0[j] = pa_x[j] * tpx_yyz_xyyy_0[j] + 0.5 * fl1_fx * tpx_yyz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyyy_0[j];

            tpy_xyyz_xyyy_0[j] = pa_x[j] * tpy_yyz_xyyy_0[j] + 0.5 * fl1_fx * tpy_yyz_yyy_0[j];

            tpz_xyyz_xyyy_0[j] = pa_x[j] * tpz_yyz_xyyy_0[j] + 0.5 * fl1_fx * tpz_yyz_yyy_0[j];

            tpx_xyyz_xyyz_0[j] = pa_x[j] * tpx_yyz_xyyz_0[j] + 0.5 * fl1_fx * tpx_yyz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyyz_0[j];

            tpy_xyyz_xyyz_0[j] = pa_x[j] * tpy_yyz_xyyz_0[j] + 0.5 * fl1_fx * tpy_yyz_yyz_0[j];

            tpz_xyyz_xyyz_0[j] = pa_x[j] * tpz_yyz_xyyz_0[j] + 0.5 * fl1_fx * tpz_yyz_yyz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_339_387(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (339,387)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 113);

        auto tpy_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tpz_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tpx_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 114);

        auto tpy_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tpz_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tpx_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 115);

        auto tpy_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tpz_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tpx_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 116);

        auto tpy_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 116);

        auto tpz_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tpx_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 117);

        auto tpy_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tpz_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tpx_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 118);

        auto tpy_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tpz_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tpx_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 119);

        auto tpy_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tpz_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tpx_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 120);

        auto tpy_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tpz_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tpx_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 121);

        auto tpy_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tpz_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tpx_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 122);

        auto tpy_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tpz_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tpx_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 123);

        auto tpy_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tpz_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tpx_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 124);

        auto tpy_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tpz_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tpx_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 125);

        auto tpy_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tpz_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tpx_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 126);

        auto tpy_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tpz_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tpx_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 127);

        auto tpy_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tpz_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tpx_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 128);

        auto tpy_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tpz_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tpx_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 78);

        auto tpy_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tpz_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tpx_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 79);

        auto tpy_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tpz_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tpx_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 80);

        auto tpy_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tpz_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tpx_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 81);

        auto tpy_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tpz_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tpx_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 82);

        auto tpy_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tpz_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tpx_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 83);

        auto tpy_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tpz_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tpx_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 84);

        auto tpy_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tpz_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tpx_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 85);

        auto tpy_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tpz_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tpx_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 86);

        auto tpy_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tpz_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tpx_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 87);

        auto tpy_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tpz_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tpx_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 88);

        auto tpy_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tpz_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto ts_yyz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 113);

        auto ts_yyz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 114);

        auto ts_yyz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 115);

        auto ts_yyz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 116);

        auto ts_yyz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 117);

        auto ts_yyz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 118);

        auto ts_yyz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 119);

        auto ts_yzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 120);

        auto ts_yzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 121);

        auto ts_yzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 122);

        auto ts_yzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 123);

        auto ts_yzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 124);

        auto ts_yzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 125);

        auto ts_yzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 126);

        auto ts_yzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 127);

        auto ts_yzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 128);

        // set up pointers to integrals

        auto tpx_xyyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 113);

        auto tpy_xyyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 113);

        auto tpz_xyyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 113);

        auto tpx_xyyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 114);

        auto tpy_xyyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 114);

        auto tpz_xyyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 114);

        auto tpx_xyyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 115);

        auto tpy_xyyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 115);

        auto tpz_xyyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 115);

        auto tpx_xyyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 116);

        auto tpy_xyyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 116);

        auto tpz_xyyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 116);

        auto tpx_xyyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 117);

        auto tpy_xyyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 117);

        auto tpz_xyyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 117);

        auto tpx_xyyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 118);

        auto tpy_xyyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 118);

        auto tpz_xyyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 118);

        auto tpx_xyyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 119);

        auto tpy_xyyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 119);

        auto tpz_xyyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 119);

        auto tpx_xyzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 120);

        auto tpy_xyzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 120);

        auto tpz_xyzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 120);

        auto tpx_xyzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 121);

        auto tpy_xyzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 121);

        auto tpz_xyzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 121);

        auto tpx_xyzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 122);

        auto tpy_xyzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 122);

        auto tpz_xyzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 122);

        auto tpx_xyzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 123);

        auto tpy_xyzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 123);

        auto tpz_xyzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 123);

        auto tpx_xyzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 124);

        auto tpy_xyzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 124);

        auto tpz_xyzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 124);

        auto tpx_xyzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 125);

        auto tpy_xyzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 125);

        auto tpz_xyzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 125);

        auto tpx_xyzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 126);

        auto tpy_xyzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 126);

        auto tpz_xyzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 126);

        auto tpx_xyzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 127);

        auto tpy_xyzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 127);

        auto tpz_xyzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 127);

        auto tpx_xyzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 128);

        auto tpy_xyzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 128);

        auto tpz_xyzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 128);

        // Batch of Integrals (339,387)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyyz_xyzz_0, tpx_xyyz_xzzz_0, tpx_xyyz_yyyy_0, \
                                     tpx_xyyz_yyyz_0, tpx_xyyz_yyzz_0, tpx_xyyz_yzzz_0, tpx_xyyz_zzzz_0, tpx_xyzz_xxxx_0, \
                                     tpx_xyzz_xxxy_0, tpx_xyzz_xxxz_0, tpx_xyzz_xxyy_0, tpx_xyzz_xxyz_0, tpx_xyzz_xxzz_0, \
                                     tpx_xyzz_xyyy_0, tpx_xyzz_xyyz_0, tpx_xyzz_xyzz_0, tpx_yyz_xyzz_0, tpx_yyz_xzzz_0, \
                                     tpx_yyz_yyyy_0, tpx_yyz_yyyz_0, tpx_yyz_yyzz_0, tpx_yyz_yzz_0, tpx_yyz_yzzz_0, \
                                     tpx_yyz_zzz_0, tpx_yyz_zzzz_0, tpx_yzz_xxx_0, tpx_yzz_xxxx_0, tpx_yzz_xxxy_0, \
                                     tpx_yzz_xxxz_0, tpx_yzz_xxy_0, tpx_yzz_xxyy_0, tpx_yzz_xxyz_0, tpx_yzz_xxz_0, \
                                     tpx_yzz_xxzz_0, tpx_yzz_xyy_0, tpx_yzz_xyyy_0, tpx_yzz_xyyz_0, tpx_yzz_xyz_0, \
                                     tpx_yzz_xyzz_0, tpx_yzz_xzz_0, tpx_yzz_yyy_0, tpx_yzz_yyz_0, tpx_yzz_yzz_0, \
                                     tpy_xyyz_xyzz_0, tpy_xyyz_xzzz_0, tpy_xyyz_yyyy_0, tpy_xyyz_yyyz_0, tpy_xyyz_yyzz_0, \
                                     tpy_xyyz_yzzz_0, tpy_xyyz_zzzz_0, tpy_xyzz_xxxx_0, tpy_xyzz_xxxy_0, tpy_xyzz_xxxz_0, \
                                     tpy_xyzz_xxyy_0, tpy_xyzz_xxyz_0, tpy_xyzz_xxzz_0, tpy_xyzz_xyyy_0, tpy_xyzz_xyyz_0, \
                                     tpy_xyzz_xyzz_0, tpy_yyz_xyzz_0, tpy_yyz_xzzz_0, tpy_yyz_yyyy_0, tpy_yyz_yyyz_0, \
                                     tpy_yyz_yyzz_0, tpy_yyz_yzz_0, tpy_yyz_yzzz_0, tpy_yyz_zzz_0, tpy_yyz_zzzz_0, \
                                     tpy_yzz_xxx_0, tpy_yzz_xxxx_0, tpy_yzz_xxxy_0, tpy_yzz_xxxz_0, tpy_yzz_xxy_0, \
                                     tpy_yzz_xxyy_0, tpy_yzz_xxyz_0, tpy_yzz_xxz_0, tpy_yzz_xxzz_0, tpy_yzz_xyy_0, \
                                     tpy_yzz_xyyy_0, tpy_yzz_xyyz_0, tpy_yzz_xyz_0, tpy_yzz_xyzz_0, tpy_yzz_xzz_0, \
                                     tpy_yzz_yyy_0, tpy_yzz_yyz_0, tpy_yzz_yzz_0, tpz_xyyz_xyzz_0, tpz_xyyz_xzzz_0, \
                                     tpz_xyyz_yyyy_0, tpz_xyyz_yyyz_0, tpz_xyyz_yyzz_0, tpz_xyyz_yzzz_0, tpz_xyyz_zzzz_0, \
                                     tpz_xyzz_xxxx_0, tpz_xyzz_xxxy_0, tpz_xyzz_xxxz_0, tpz_xyzz_xxyy_0, tpz_xyzz_xxyz_0, \
                                     tpz_xyzz_xxzz_0, tpz_xyzz_xyyy_0, tpz_xyzz_xyyz_0, tpz_xyzz_xyzz_0, tpz_yyz_xyzz_0, \
                                     tpz_yyz_xzzz_0, tpz_yyz_yyyy_0, tpz_yyz_yyyz_0, tpz_yyz_yyzz_0, tpz_yyz_yzz_0, \
                                     tpz_yyz_yzzz_0, tpz_yyz_zzz_0, tpz_yyz_zzzz_0, tpz_yzz_xxx_0, tpz_yzz_xxxx_0, \
                                     tpz_yzz_xxxy_0, tpz_yzz_xxxz_0, tpz_yzz_xxy_0, tpz_yzz_xxyy_0, tpz_yzz_xxyz_0, \
                                     tpz_yzz_xxz_0, tpz_yzz_xxzz_0, tpz_yzz_xyy_0, tpz_yzz_xyyy_0, tpz_yzz_xyyz_0, \
                                     tpz_yzz_xyz_0, tpz_yzz_xyzz_0, tpz_yzz_xzz_0, tpz_yzz_yyy_0, tpz_yzz_yyz_0, \
                                     tpz_yzz_yzz_0, ts_yyz_xyzz_0, ts_yyz_xzzz_0, ts_yyz_yyyy_0, ts_yyz_yyyz_0, \
                                     ts_yyz_yyzz_0, ts_yyz_yzzz_0, ts_yyz_zzzz_0, ts_yzz_xxxx_0, ts_yzz_xxxy_0, \
                                     ts_yzz_xxxz_0, ts_yzz_xxyy_0, ts_yzz_xxyz_0, ts_yzz_xxzz_0, ts_yzz_xyyy_0, \
                                     ts_yzz_xyyz_0, ts_yzz_xyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xyyz_xyzz_0[j] = pa_x[j] * tpx_yyz_xyzz_0[j] + 0.5 * fl1_fx * tpx_yyz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyzz_0[j];

            tpy_xyyz_xyzz_0[j] = pa_x[j] * tpy_yyz_xyzz_0[j] + 0.5 * fl1_fx * tpy_yyz_yzz_0[j];

            tpz_xyyz_xyzz_0[j] = pa_x[j] * tpz_yyz_xyzz_0[j] + 0.5 * fl1_fx * tpz_yyz_yzz_0[j];

            tpx_xyyz_xzzz_0[j] = pa_x[j] * tpx_yyz_xzzz_0[j] + 0.5 * fl1_fx * tpx_yyz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xzzz_0[j];

            tpy_xyyz_xzzz_0[j] = pa_x[j] * tpy_yyz_xzzz_0[j] + 0.5 * fl1_fx * tpy_yyz_zzz_0[j];

            tpz_xyyz_xzzz_0[j] = pa_x[j] * tpz_yyz_xzzz_0[j] + 0.5 * fl1_fx * tpz_yyz_zzz_0[j];

            tpx_xyyz_yyyy_0[j] = pa_x[j] * tpx_yyz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyyy_0[j];

            tpy_xyyz_yyyy_0[j] = pa_x[j] * tpy_yyz_yyyy_0[j];

            tpz_xyyz_yyyy_0[j] = pa_x[j] * tpz_yyz_yyyy_0[j];

            tpx_xyyz_yyyz_0[j] = pa_x[j] * tpx_yyz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyyz_0[j];

            tpy_xyyz_yyyz_0[j] = pa_x[j] * tpy_yyz_yyyz_0[j];

            tpz_xyyz_yyyz_0[j] = pa_x[j] * tpz_yyz_yyyz_0[j];

            tpx_xyyz_yyzz_0[j] = pa_x[j] * tpx_yyz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyzz_0[j];

            tpy_xyyz_yyzz_0[j] = pa_x[j] * tpy_yyz_yyzz_0[j];

            tpz_xyyz_yyzz_0[j] = pa_x[j] * tpz_yyz_yyzz_0[j];

            tpx_xyyz_yzzz_0[j] = pa_x[j] * tpx_yyz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yzzz_0[j];

            tpy_xyyz_yzzz_0[j] = pa_x[j] * tpy_yyz_yzzz_0[j];

            tpz_xyyz_yzzz_0[j] = pa_x[j] * tpz_yyz_yzzz_0[j];

            tpx_xyyz_zzzz_0[j] = pa_x[j] * tpx_yyz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_zzzz_0[j];

            tpy_xyyz_zzzz_0[j] = pa_x[j] * tpy_yyz_zzzz_0[j];

            tpz_xyyz_zzzz_0[j] = pa_x[j] * tpz_yyz_zzzz_0[j];

            tpx_xyzz_xxxx_0[j] = pa_x[j] * tpx_yzz_xxxx_0[j] + 2.0 * fl1_fx * tpx_yzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxxx_0[j];

            tpy_xyzz_xxxx_0[j] = pa_x[j] * tpy_yzz_xxxx_0[j] + 2.0 * fl1_fx * tpy_yzz_xxx_0[j];

            tpz_xyzz_xxxx_0[j] = pa_x[j] * tpz_yzz_xxxx_0[j] + 2.0 * fl1_fx * tpz_yzz_xxx_0[j];

            tpx_xyzz_xxxy_0[j] = pa_x[j] * tpx_yzz_xxxy_0[j] + 1.5 * fl1_fx * tpx_yzz_xxy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxxy_0[j];

            tpy_xyzz_xxxy_0[j] = pa_x[j] * tpy_yzz_xxxy_0[j] + 1.5 * fl1_fx * tpy_yzz_xxy_0[j];

            tpz_xyzz_xxxy_0[j] = pa_x[j] * tpz_yzz_xxxy_0[j] + 1.5 * fl1_fx * tpz_yzz_xxy_0[j];

            tpx_xyzz_xxxz_0[j] = pa_x[j] * tpx_yzz_xxxz_0[j] + 1.5 * fl1_fx * tpx_yzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxxz_0[j];

            tpy_xyzz_xxxz_0[j] = pa_x[j] * tpy_yzz_xxxz_0[j] + 1.5 * fl1_fx * tpy_yzz_xxz_0[j];

            tpz_xyzz_xxxz_0[j] = pa_x[j] * tpz_yzz_xxxz_0[j] + 1.5 * fl1_fx * tpz_yzz_xxz_0[j];

            tpx_xyzz_xxyy_0[j] = pa_x[j] * tpx_yzz_xxyy_0[j] + fl1_fx * tpx_yzz_xyy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxyy_0[j];

            tpy_xyzz_xxyy_0[j] = pa_x[j] * tpy_yzz_xxyy_0[j] + fl1_fx * tpy_yzz_xyy_0[j];

            tpz_xyzz_xxyy_0[j] = pa_x[j] * tpz_yzz_xxyy_0[j] + fl1_fx * tpz_yzz_xyy_0[j];

            tpx_xyzz_xxyz_0[j] = pa_x[j] * tpx_yzz_xxyz_0[j] + fl1_fx * tpx_yzz_xyz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxyz_0[j];

            tpy_xyzz_xxyz_0[j] = pa_x[j] * tpy_yzz_xxyz_0[j] + fl1_fx * tpy_yzz_xyz_0[j];

            tpz_xyzz_xxyz_0[j] = pa_x[j] * tpz_yzz_xxyz_0[j] + fl1_fx * tpz_yzz_xyz_0[j];

            tpx_xyzz_xxzz_0[j] = pa_x[j] * tpx_yzz_xxzz_0[j] + fl1_fx * tpx_yzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxzz_0[j];

            tpy_xyzz_xxzz_0[j] = pa_x[j] * tpy_yzz_xxzz_0[j] + fl1_fx * tpy_yzz_xzz_0[j];

            tpz_xyzz_xxzz_0[j] = pa_x[j] * tpz_yzz_xxzz_0[j] + fl1_fx * tpz_yzz_xzz_0[j];

            tpx_xyzz_xyyy_0[j] = pa_x[j] * tpx_yzz_xyyy_0[j] + 0.5 * fl1_fx * tpx_yzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyyy_0[j];

            tpy_xyzz_xyyy_0[j] = pa_x[j] * tpy_yzz_xyyy_0[j] + 0.5 * fl1_fx * tpy_yzz_yyy_0[j];

            tpz_xyzz_xyyy_0[j] = pa_x[j] * tpz_yzz_xyyy_0[j] + 0.5 * fl1_fx * tpz_yzz_yyy_0[j];

            tpx_xyzz_xyyz_0[j] = pa_x[j] * tpx_yzz_xyyz_0[j] + 0.5 * fl1_fx * tpx_yzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyyz_0[j];

            tpy_xyzz_xyyz_0[j] = pa_x[j] * tpy_yzz_xyyz_0[j] + 0.5 * fl1_fx * tpy_yzz_yyz_0[j];

            tpz_xyzz_xyyz_0[j] = pa_x[j] * tpz_yzz_xyyz_0[j] + 0.5 * fl1_fx * tpz_yzz_yyz_0[j];

            tpx_xyzz_xyzz_0[j] = pa_x[j] * tpx_yzz_xyzz_0[j] + 0.5 * fl1_fx * tpx_yzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyzz_0[j];

            tpy_xyzz_xyzz_0[j] = pa_x[j] * tpy_yzz_xyzz_0[j] + 0.5 * fl1_fx * tpy_yzz_yzz_0[j];

            tpz_xyzz_xyzz_0[j] = pa_x[j] * tpz_yzz_xyzz_0[j] + 0.5 * fl1_fx * tpz_yzz_yzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_387_435(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (387,435)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tpx_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 129);

        auto tpy_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tpz_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tpx_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 130);

        auto tpy_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tpz_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tpx_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 131);

        auto tpy_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tpz_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tpx_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 132);

        auto tpy_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tpz_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tpx_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 133);

        auto tpy_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tpz_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tpx_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 134);

        auto tpy_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tpz_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tpx_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 135);

        auto tpy_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tpz_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tpx_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 136);

        auto tpy_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tpz_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tpx_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 137);

        auto tpy_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tpz_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tpx_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 138);

        auto tpy_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tpz_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tpx_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 139);

        auto tpy_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tpz_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tpx_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 140);

        auto tpy_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tpz_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tpx_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 141);

        auto tpy_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tpz_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tpx_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 142);

        auto tpy_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tpz_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tpx_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 143);

        auto tpy_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tpz_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tpx_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 144);

        auto tpy_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tpz_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tpx_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 89);

        auto tpy_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tpz_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tpx_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 90);

        auto tpy_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tpz_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tpx_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 91);

        auto tpy_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tpz_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tpx_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 92);

        auto tpy_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tpz_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tpx_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 93);

        auto tpy_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tpz_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tpx_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 94);

        auto tpy_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tpz_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tpx_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 95);

        auto tpy_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tpz_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tpx_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 96);

        auto tpy_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tpz_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tpx_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 97);

        auto tpy_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tpz_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tpx_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 98);

        auto tpy_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tpz_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tpx_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 99);

        auto tpy_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tpz_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto ts_yzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 129);

        auto ts_yzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 130);

        auto ts_yzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 131);

        auto ts_yzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 132);

        auto ts_yzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 133);

        auto ts_yzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 134);

        auto ts_zzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 135);

        auto ts_zzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 136);

        auto ts_zzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 137);

        auto ts_zzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 138);

        auto ts_zzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 139);

        auto ts_zzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 140);

        auto ts_zzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 141);

        auto ts_zzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 142);

        auto ts_zzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 143);

        auto ts_zzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 144);

        // set up pointers to integrals

        auto tpx_xyzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 129);

        auto tpy_xyzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 129);

        auto tpz_xyzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 129);

        auto tpx_xyzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 130);

        auto tpy_xyzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 130);

        auto tpz_xyzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 130);

        auto tpx_xyzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 131);

        auto tpy_xyzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 131);

        auto tpz_xyzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 131);

        auto tpx_xyzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 132);

        auto tpy_xyzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 132);

        auto tpz_xyzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 132);

        auto tpx_xyzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 133);

        auto tpy_xyzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 133);

        auto tpz_xyzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 133);

        auto tpx_xyzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 134);

        auto tpy_xyzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 134);

        auto tpz_xyzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 134);

        auto tpx_xzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 135);

        auto tpy_xzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 135);

        auto tpz_xzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 135);

        auto tpx_xzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 136);

        auto tpy_xzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 136);

        auto tpz_xzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 136);

        auto tpx_xzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 137);

        auto tpy_xzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 137);

        auto tpz_xzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 137);

        auto tpx_xzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 138);

        auto tpy_xzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 138);

        auto tpz_xzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 138);

        auto tpx_xzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 139);

        auto tpy_xzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 139);

        auto tpz_xzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 139);

        auto tpx_xzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 140);

        auto tpy_xzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 140);

        auto tpz_xzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 140);

        auto tpx_xzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 141);

        auto tpy_xzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 141);

        auto tpz_xzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 141);

        auto tpx_xzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 142);

        auto tpy_xzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 142);

        auto tpz_xzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 142);

        auto tpx_xzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 143);

        auto tpy_xzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 143);

        auto tpz_xzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 143);

        auto tpx_xzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 144);

        auto tpy_xzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 144);

        auto tpz_xzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 144);

        // Batch of Integrals (387,435)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyzz_xzzz_0, tpx_xyzz_yyyy_0, tpx_xyzz_yyyz_0, \
                                     tpx_xyzz_yyzz_0, tpx_xyzz_yzzz_0, tpx_xyzz_zzzz_0, tpx_xzzz_xxxx_0, tpx_xzzz_xxxy_0, \
                                     tpx_xzzz_xxxz_0, tpx_xzzz_xxyy_0, tpx_xzzz_xxyz_0, tpx_xzzz_xxzz_0, tpx_xzzz_xyyy_0, \
                                     tpx_xzzz_xyyz_0, tpx_xzzz_xyzz_0, tpx_xzzz_xzzz_0, tpx_yzz_xzzz_0, tpx_yzz_yyyy_0, \
                                     tpx_yzz_yyyz_0, tpx_yzz_yyzz_0, tpx_yzz_yzzz_0, tpx_yzz_zzz_0, tpx_yzz_zzzz_0, \
                                     tpx_zzz_xxx_0, tpx_zzz_xxxx_0, tpx_zzz_xxxy_0, tpx_zzz_xxxz_0, tpx_zzz_xxy_0, \
                                     tpx_zzz_xxyy_0, tpx_zzz_xxyz_0, tpx_zzz_xxz_0, tpx_zzz_xxzz_0, tpx_zzz_xyy_0, \
                                     tpx_zzz_xyyy_0, tpx_zzz_xyyz_0, tpx_zzz_xyz_0, tpx_zzz_xyzz_0, tpx_zzz_xzz_0, \
                                     tpx_zzz_xzzz_0, tpx_zzz_yyy_0, tpx_zzz_yyz_0, tpx_zzz_yzz_0, tpx_zzz_zzz_0, \
                                     tpy_xyzz_xzzz_0, tpy_xyzz_yyyy_0, tpy_xyzz_yyyz_0, tpy_xyzz_yyzz_0, tpy_xyzz_yzzz_0, \
                                     tpy_xyzz_zzzz_0, tpy_xzzz_xxxx_0, tpy_xzzz_xxxy_0, tpy_xzzz_xxxz_0, tpy_xzzz_xxyy_0, \
                                     tpy_xzzz_xxyz_0, tpy_xzzz_xxzz_0, tpy_xzzz_xyyy_0, tpy_xzzz_xyyz_0, tpy_xzzz_xyzz_0, \
                                     tpy_xzzz_xzzz_0, tpy_yzz_xzzz_0, tpy_yzz_yyyy_0, tpy_yzz_yyyz_0, tpy_yzz_yyzz_0, \
                                     tpy_yzz_yzzz_0, tpy_yzz_zzz_0, tpy_yzz_zzzz_0, tpy_zzz_xxx_0, tpy_zzz_xxxx_0, \
                                     tpy_zzz_xxxy_0, tpy_zzz_xxxz_0, tpy_zzz_xxy_0, tpy_zzz_xxyy_0, tpy_zzz_xxyz_0, \
                                     tpy_zzz_xxz_0, tpy_zzz_xxzz_0, tpy_zzz_xyy_0, tpy_zzz_xyyy_0, tpy_zzz_xyyz_0, \
                                     tpy_zzz_xyz_0, tpy_zzz_xyzz_0, tpy_zzz_xzz_0, tpy_zzz_xzzz_0, tpy_zzz_yyy_0, \
                                     tpy_zzz_yyz_0, tpy_zzz_yzz_0, tpy_zzz_zzz_0, tpz_xyzz_xzzz_0, tpz_xyzz_yyyy_0, \
                                     tpz_xyzz_yyyz_0, tpz_xyzz_yyzz_0, tpz_xyzz_yzzz_0, tpz_xyzz_zzzz_0, tpz_xzzz_xxxx_0, \
                                     tpz_xzzz_xxxy_0, tpz_xzzz_xxxz_0, tpz_xzzz_xxyy_0, tpz_xzzz_xxyz_0, tpz_xzzz_xxzz_0, \
                                     tpz_xzzz_xyyy_0, tpz_xzzz_xyyz_0, tpz_xzzz_xyzz_0, tpz_xzzz_xzzz_0, tpz_yzz_xzzz_0, \
                                     tpz_yzz_yyyy_0, tpz_yzz_yyyz_0, tpz_yzz_yyzz_0, tpz_yzz_yzzz_0, tpz_yzz_zzz_0, \
                                     tpz_yzz_zzzz_0, tpz_zzz_xxx_0, tpz_zzz_xxxx_0, tpz_zzz_xxxy_0, tpz_zzz_xxxz_0, \
                                     tpz_zzz_xxy_0, tpz_zzz_xxyy_0, tpz_zzz_xxyz_0, tpz_zzz_xxz_0, tpz_zzz_xxzz_0, \
                                     tpz_zzz_xyy_0, tpz_zzz_xyyy_0, tpz_zzz_xyyz_0, tpz_zzz_xyz_0, tpz_zzz_xyzz_0, \
                                     tpz_zzz_xzz_0, tpz_zzz_xzzz_0, tpz_zzz_yyy_0, tpz_zzz_yyz_0, tpz_zzz_yzz_0, \
                                     tpz_zzz_zzz_0, ts_yzz_xzzz_0, ts_yzz_yyyy_0, ts_yzz_yyyz_0, ts_yzz_yyzz_0, \
                                     ts_yzz_yzzz_0, ts_yzz_zzzz_0, ts_zzz_xxxx_0, ts_zzz_xxxy_0, ts_zzz_xxxz_0, \
                                     ts_zzz_xxyy_0, ts_zzz_xxyz_0, ts_zzz_xxzz_0, ts_zzz_xyyy_0, ts_zzz_xyyz_0, \
                                     ts_zzz_xyzz_0, ts_zzz_xzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xyzz_xzzz_0[j] = pa_x[j] * tpx_yzz_xzzz_0[j] + 0.5 * fl1_fx * tpx_yzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xzzz_0[j];

            tpy_xyzz_xzzz_0[j] = pa_x[j] * tpy_yzz_xzzz_0[j] + 0.5 * fl1_fx * tpy_yzz_zzz_0[j];

            tpz_xyzz_xzzz_0[j] = pa_x[j] * tpz_yzz_xzzz_0[j] + 0.5 * fl1_fx * tpz_yzz_zzz_0[j];

            tpx_xyzz_yyyy_0[j] = pa_x[j] * tpx_yzz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyyy_0[j];

            tpy_xyzz_yyyy_0[j] = pa_x[j] * tpy_yzz_yyyy_0[j];

            tpz_xyzz_yyyy_0[j] = pa_x[j] * tpz_yzz_yyyy_0[j];

            tpx_xyzz_yyyz_0[j] = pa_x[j] * tpx_yzz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyyz_0[j];

            tpy_xyzz_yyyz_0[j] = pa_x[j] * tpy_yzz_yyyz_0[j];

            tpz_xyzz_yyyz_0[j] = pa_x[j] * tpz_yzz_yyyz_0[j];

            tpx_xyzz_yyzz_0[j] = pa_x[j] * tpx_yzz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyzz_0[j];

            tpy_xyzz_yyzz_0[j] = pa_x[j] * tpy_yzz_yyzz_0[j];

            tpz_xyzz_yyzz_0[j] = pa_x[j] * tpz_yzz_yyzz_0[j];

            tpx_xyzz_yzzz_0[j] = pa_x[j] * tpx_yzz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yzzz_0[j];

            tpy_xyzz_yzzz_0[j] = pa_x[j] * tpy_yzz_yzzz_0[j];

            tpz_xyzz_yzzz_0[j] = pa_x[j] * tpz_yzz_yzzz_0[j];

            tpx_xyzz_zzzz_0[j] = pa_x[j] * tpx_yzz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_zzzz_0[j];

            tpy_xyzz_zzzz_0[j] = pa_x[j] * tpy_yzz_zzzz_0[j];

            tpz_xyzz_zzzz_0[j] = pa_x[j] * tpz_yzz_zzzz_0[j];

            tpx_xzzz_xxxx_0[j] = pa_x[j] * tpx_zzz_xxxx_0[j] + 2.0 * fl1_fx * tpx_zzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxx_0[j];

            tpy_xzzz_xxxx_0[j] = pa_x[j] * tpy_zzz_xxxx_0[j] + 2.0 * fl1_fx * tpy_zzz_xxx_0[j];

            tpz_xzzz_xxxx_0[j] = pa_x[j] * tpz_zzz_xxxx_0[j] + 2.0 * fl1_fx * tpz_zzz_xxx_0[j];

            tpx_xzzz_xxxy_0[j] = pa_x[j] * tpx_zzz_xxxy_0[j] + 1.5 * fl1_fx * tpx_zzz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxy_0[j];

            tpy_xzzz_xxxy_0[j] = pa_x[j] * tpy_zzz_xxxy_0[j] + 1.5 * fl1_fx * tpy_zzz_xxy_0[j];

            tpz_xzzz_xxxy_0[j] = pa_x[j] * tpz_zzz_xxxy_0[j] + 1.5 * fl1_fx * tpz_zzz_xxy_0[j];

            tpx_xzzz_xxxz_0[j] = pa_x[j] * tpx_zzz_xxxz_0[j] + 1.5 * fl1_fx * tpx_zzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxz_0[j];

            tpy_xzzz_xxxz_0[j] = pa_x[j] * tpy_zzz_xxxz_0[j] + 1.5 * fl1_fx * tpy_zzz_xxz_0[j];

            tpz_xzzz_xxxz_0[j] = pa_x[j] * tpz_zzz_xxxz_0[j] + 1.5 * fl1_fx * tpz_zzz_xxz_0[j];

            tpx_xzzz_xxyy_0[j] = pa_x[j] * tpx_zzz_xxyy_0[j] + fl1_fx * tpx_zzz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxyy_0[j];

            tpy_xzzz_xxyy_0[j] = pa_x[j] * tpy_zzz_xxyy_0[j] + fl1_fx * tpy_zzz_xyy_0[j];

            tpz_xzzz_xxyy_0[j] = pa_x[j] * tpz_zzz_xxyy_0[j] + fl1_fx * tpz_zzz_xyy_0[j];

            tpx_xzzz_xxyz_0[j] = pa_x[j] * tpx_zzz_xxyz_0[j] + fl1_fx * tpx_zzz_xyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxyz_0[j];

            tpy_xzzz_xxyz_0[j] = pa_x[j] * tpy_zzz_xxyz_0[j] + fl1_fx * tpy_zzz_xyz_0[j];

            tpz_xzzz_xxyz_0[j] = pa_x[j] * tpz_zzz_xxyz_0[j] + fl1_fx * tpz_zzz_xyz_0[j];

            tpx_xzzz_xxzz_0[j] = pa_x[j] * tpx_zzz_xxzz_0[j] + fl1_fx * tpx_zzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxzz_0[j];

            tpy_xzzz_xxzz_0[j] = pa_x[j] * tpy_zzz_xxzz_0[j] + fl1_fx * tpy_zzz_xzz_0[j];

            tpz_xzzz_xxzz_0[j] = pa_x[j] * tpz_zzz_xxzz_0[j] + fl1_fx * tpz_zzz_xzz_0[j];

            tpx_xzzz_xyyy_0[j] = pa_x[j] * tpx_zzz_xyyy_0[j] + 0.5 * fl1_fx * tpx_zzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyyy_0[j];

            tpy_xzzz_xyyy_0[j] = pa_x[j] * tpy_zzz_xyyy_0[j] + 0.5 * fl1_fx * tpy_zzz_yyy_0[j];

            tpz_xzzz_xyyy_0[j] = pa_x[j] * tpz_zzz_xyyy_0[j] + 0.5 * fl1_fx * tpz_zzz_yyy_0[j];

            tpx_xzzz_xyyz_0[j] = pa_x[j] * tpx_zzz_xyyz_0[j] + 0.5 * fl1_fx * tpx_zzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyyz_0[j];

            tpy_xzzz_xyyz_0[j] = pa_x[j] * tpy_zzz_xyyz_0[j] + 0.5 * fl1_fx * tpy_zzz_yyz_0[j];

            tpz_xzzz_xyyz_0[j] = pa_x[j] * tpz_zzz_xyyz_0[j] + 0.5 * fl1_fx * tpz_zzz_yyz_0[j];

            tpx_xzzz_xyzz_0[j] = pa_x[j] * tpx_zzz_xyzz_0[j] + 0.5 * fl1_fx * tpx_zzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyzz_0[j];

            tpy_xzzz_xyzz_0[j] = pa_x[j] * tpy_zzz_xyzz_0[j] + 0.5 * fl1_fx * tpy_zzz_yzz_0[j];

            tpz_xzzz_xyzz_0[j] = pa_x[j] * tpz_zzz_xyzz_0[j] + 0.5 * fl1_fx * tpz_zzz_yzz_0[j];

            tpx_xzzz_xzzz_0[j] = pa_x[j] * tpx_zzz_xzzz_0[j] + 0.5 * fl1_fx * tpx_zzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xzzz_0[j];

            tpy_xzzz_xzzz_0[j] = pa_x[j] * tpy_zzz_xzzz_0[j] + 0.5 * fl1_fx * tpy_zzz_zzz_0[j];

            tpz_xzzz_xzzz_0[j] = pa_x[j] * tpz_zzz_xzzz_0[j] + 0.5 * fl1_fx * tpz_zzz_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_435_483(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (435,483)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tpx_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 90);

        auto tpy_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tpz_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tpx_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 91);

        auto tpy_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tpz_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tpx_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 92);

        auto tpy_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tpz_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tpx_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 93);

        auto tpy_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tpz_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tpx_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 94);

        auto tpy_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tpz_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tpx_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 95);

        auto tpy_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tpz_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tpx_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 96);

        auto tpy_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tpz_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tpx_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 97);

        auto tpy_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tpz_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tpx_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 98);

        auto tpy_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tpz_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tpx_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 99);

        auto tpy_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tpz_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tpx_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 100);

        auto tpy_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tpz_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tpx_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 145);

        auto tpy_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tpz_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tpx_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 146);

        auto tpy_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tpz_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tpx_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 147);

        auto tpy_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tpz_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tpx_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 148);

        auto tpy_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tpz_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tpx_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 149);

        auto tpy_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tpz_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 149);

        auto tpx_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 45);

        auto tpy_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tpz_yy_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tpx_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 46);

        auto tpy_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tpz_yy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tpx_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 47);

        auto tpy_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tpz_yy_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tpx_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 48);

        auto tpy_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tpz_yy_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tpx_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 49);

        auto tpy_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tpz_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tpx_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 50);

        auto tpy_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tpz_yy_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tpx_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 51);

        auto tpy_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tpz_yy_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tpx_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 52);

        auto tpy_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tpz_yy_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tpx_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 53);

        auto tpy_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tpz_yy_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tpx_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 54);

        auto tpy_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tpz_yy_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tpx_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 55);

        auto tpy_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tpz_yy_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tpx_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 60);

        auto tpy_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tpz_yyy_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tpx_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 61);

        auto tpy_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tpz_yyy_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tpx_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 62);

        auto tpy_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tpz_yyy_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tpx_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 63);

        auto tpy_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tpz_yyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tpx_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 64);

        auto tpy_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tpz_yyy_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tpx_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 65);

        auto tpy_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tpz_yyy_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tpx_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 66);

        auto tpy_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tpz_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto ts_yyy_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 90);

        auto ts_yyy_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 91);

        auto ts_yyy_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 92);

        auto ts_yyy_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 93);

        auto ts_yyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 94);

        auto ts_yyy_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 95);

        auto ts_yyy_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 96);

        auto ts_yyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 97);

        auto ts_yyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 98);

        auto ts_yyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 99);

        auto ts_yyy_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 100);

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        auto ts_zzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 149);

        // set up pointers to integrals

        auto tpx_xzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 145);

        auto tpy_xzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 145);

        auto tpz_xzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 145);

        auto tpx_xzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 146);

        auto tpy_xzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 146);

        auto tpz_xzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 146);

        auto tpx_xzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 147);

        auto tpy_xzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 147);

        auto tpz_xzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 147);

        auto tpx_xzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 148);

        auto tpy_xzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 148);

        auto tpz_xzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 148);

        auto tpx_xzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 149);

        auto tpy_xzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 149);

        auto tpz_xzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 149);

        auto tpx_yyyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 150);

        auto tpy_yyyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 150);

        auto tpz_yyyy_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 150);

        auto tpx_yyyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 151);

        auto tpy_yyyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 151);

        auto tpz_yyyy_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 151);

        auto tpx_yyyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 152);

        auto tpy_yyyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 152);

        auto tpz_yyyy_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 152);

        auto tpx_yyyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 153);

        auto tpy_yyyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 153);

        auto tpz_yyyy_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 153);

        auto tpx_yyyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 154);

        auto tpy_yyyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 154);

        auto tpz_yyyy_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 154);

        auto tpx_yyyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 155);

        auto tpy_yyyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 155);

        auto tpz_yyyy_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 155);

        auto tpx_yyyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 156);

        auto tpy_yyyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 156);

        auto tpz_yyyy_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 156);

        auto tpx_yyyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 157);

        auto tpy_yyyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 157);

        auto tpz_yyyy_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 157);

        auto tpx_yyyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 158);

        auto tpy_yyyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 158);

        auto tpz_yyyy_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 158);

        auto tpx_yyyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 159);

        auto tpy_yyyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 159);

        auto tpz_yyyy_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 159);

        auto tpx_yyyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 160);

        auto tpy_yyyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 160);

        auto tpz_yyyy_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 160);

        // Batch of Integrals (435,483)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_xzzz_yyyy_0, tpx_xzzz_yyyz_0, tpx_xzzz_yyzz_0, \
                                     tpx_xzzz_yzzz_0, tpx_xzzz_zzzz_0, tpx_yy_xxxx_0, tpx_yy_xxxy_0, tpx_yy_xxxz_0, \
                                     tpx_yy_xxyy_0, tpx_yy_xxyz_0, tpx_yy_xxzz_0, tpx_yy_xyyy_0, tpx_yy_xyyz_0, \
                                     tpx_yy_xyzz_0, tpx_yy_xzzz_0, tpx_yy_yyyy_0, tpx_yyy_xxx_0, tpx_yyy_xxxx_0, \
                                     tpx_yyy_xxxy_0, tpx_yyy_xxxz_0, tpx_yyy_xxy_0, tpx_yyy_xxyy_0, tpx_yyy_xxyz_0, \
                                     tpx_yyy_xxz_0, tpx_yyy_xxzz_0, tpx_yyy_xyy_0, tpx_yyy_xyyy_0, tpx_yyy_xyyz_0, \
                                     tpx_yyy_xyz_0, tpx_yyy_xyzz_0, tpx_yyy_xzz_0, tpx_yyy_xzzz_0, tpx_yyy_yyy_0, \
                                     tpx_yyy_yyyy_0, tpx_yyyy_xxxx_0, tpx_yyyy_xxxy_0, tpx_yyyy_xxxz_0, tpx_yyyy_xxyy_0, \
                                     tpx_yyyy_xxyz_0, tpx_yyyy_xxzz_0, tpx_yyyy_xyyy_0, tpx_yyyy_xyyz_0, tpx_yyyy_xyzz_0, \
                                     tpx_yyyy_xzzz_0, tpx_yyyy_yyyy_0, tpx_zzz_yyyy_0, tpx_zzz_yyyz_0, tpx_zzz_yyzz_0, \
                                     tpx_zzz_yzzz_0, tpx_zzz_zzzz_0, tpy_xzzz_yyyy_0, tpy_xzzz_yyyz_0, tpy_xzzz_yyzz_0, \
                                     tpy_xzzz_yzzz_0, tpy_xzzz_zzzz_0, tpy_yy_xxxx_0, tpy_yy_xxxy_0, tpy_yy_xxxz_0, \
                                     tpy_yy_xxyy_0, tpy_yy_xxyz_0, tpy_yy_xxzz_0, tpy_yy_xyyy_0, tpy_yy_xyyz_0, \
                                     tpy_yy_xyzz_0, tpy_yy_xzzz_0, tpy_yy_yyyy_0, tpy_yyy_xxx_0, tpy_yyy_xxxx_0, \
                                     tpy_yyy_xxxy_0, tpy_yyy_xxxz_0, tpy_yyy_xxy_0, tpy_yyy_xxyy_0, tpy_yyy_xxyz_0, \
                                     tpy_yyy_xxz_0, tpy_yyy_xxzz_0, tpy_yyy_xyy_0, tpy_yyy_xyyy_0, tpy_yyy_xyyz_0, \
                                     tpy_yyy_xyz_0, tpy_yyy_xyzz_0, tpy_yyy_xzz_0, tpy_yyy_xzzz_0, tpy_yyy_yyy_0, \
                                     tpy_yyy_yyyy_0, tpy_yyyy_xxxx_0, tpy_yyyy_xxxy_0, tpy_yyyy_xxxz_0, tpy_yyyy_xxyy_0, \
                                     tpy_yyyy_xxyz_0, tpy_yyyy_xxzz_0, tpy_yyyy_xyyy_0, tpy_yyyy_xyyz_0, tpy_yyyy_xyzz_0, \
                                     tpy_yyyy_xzzz_0, tpy_yyyy_yyyy_0, tpy_zzz_yyyy_0, tpy_zzz_yyyz_0, tpy_zzz_yyzz_0, \
                                     tpy_zzz_yzzz_0, tpy_zzz_zzzz_0, tpz_xzzz_yyyy_0, tpz_xzzz_yyyz_0, tpz_xzzz_yyzz_0, \
                                     tpz_xzzz_yzzz_0, tpz_xzzz_zzzz_0, tpz_yy_xxxx_0, tpz_yy_xxxy_0, tpz_yy_xxxz_0, \
                                     tpz_yy_xxyy_0, tpz_yy_xxyz_0, tpz_yy_xxzz_0, tpz_yy_xyyy_0, tpz_yy_xyyz_0, \
                                     tpz_yy_xyzz_0, tpz_yy_xzzz_0, tpz_yy_yyyy_0, tpz_yyy_xxx_0, tpz_yyy_xxxx_0, \
                                     tpz_yyy_xxxy_0, tpz_yyy_xxxz_0, tpz_yyy_xxy_0, tpz_yyy_xxyy_0, tpz_yyy_xxyz_0, \
                                     tpz_yyy_xxz_0, tpz_yyy_xxzz_0, tpz_yyy_xyy_0, tpz_yyy_xyyy_0, tpz_yyy_xyyz_0, \
                                     tpz_yyy_xyz_0, tpz_yyy_xyzz_0, tpz_yyy_xzz_0, tpz_yyy_xzzz_0, tpz_yyy_yyy_0, \
                                     tpz_yyy_yyyy_0, tpz_yyyy_xxxx_0, tpz_yyyy_xxxy_0, tpz_yyyy_xxxz_0, tpz_yyyy_xxyy_0, \
                                     tpz_yyyy_xxyz_0, tpz_yyyy_xxzz_0, tpz_yyyy_xyyy_0, tpz_yyyy_xyyz_0, tpz_yyyy_xyzz_0, \
                                     tpz_yyyy_xzzz_0, tpz_yyyy_yyyy_0, tpz_zzz_yyyy_0, tpz_zzz_yyyz_0, tpz_zzz_yyzz_0, \
                                     tpz_zzz_yzzz_0, tpz_zzz_zzzz_0, ts_yyy_xxxx_0, ts_yyy_xxxy_0, ts_yyy_xxxz_0, \
                                     ts_yyy_xxyy_0, ts_yyy_xxyz_0, ts_yyy_xxzz_0, ts_yyy_xyyy_0, ts_yyy_xyyz_0, \
                                     ts_yyy_xyzz_0, ts_yyy_xzzz_0, ts_yyy_yyyy_0, ts_zzz_yyyy_0, ts_zzz_yyyz_0, \
                                     ts_zzz_yyzz_0, ts_zzz_yzzz_0, ts_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xzzz_yyyy_0[j] = pa_x[j] * tpx_zzz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyyy_0[j];

            tpy_xzzz_yyyy_0[j] = pa_x[j] * tpy_zzz_yyyy_0[j];

            tpz_xzzz_yyyy_0[j] = pa_x[j] * tpz_zzz_yyyy_0[j];

            tpx_xzzz_yyyz_0[j] = pa_x[j] * tpx_zzz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyyz_0[j];

            tpy_xzzz_yyyz_0[j] = pa_x[j] * tpy_zzz_yyyz_0[j];

            tpz_xzzz_yyyz_0[j] = pa_x[j] * tpz_zzz_yyyz_0[j];

            tpx_xzzz_yyzz_0[j] = pa_x[j] * tpx_zzz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyzz_0[j];

            tpy_xzzz_yyzz_0[j] = pa_x[j] * tpy_zzz_yyzz_0[j];

            tpz_xzzz_yyzz_0[j] = pa_x[j] * tpz_zzz_yyzz_0[j];

            tpx_xzzz_yzzz_0[j] = pa_x[j] * tpx_zzz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yzzz_0[j];

            tpy_xzzz_yzzz_0[j] = pa_x[j] * tpy_zzz_yzzz_0[j];

            tpz_xzzz_yzzz_0[j] = pa_x[j] * tpz_zzz_yzzz_0[j];

            tpx_xzzz_zzzz_0[j] = pa_x[j] * tpx_zzz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zzzz_0[j];

            tpy_xzzz_zzzz_0[j] = pa_x[j] * tpy_zzz_zzzz_0[j];

            tpz_xzzz_zzzz_0[j] = pa_x[j] * tpz_zzz_zzzz_0[j];

            tpx_yyyy_xxxx_0[j] = pa_y[j] * tpx_yyy_xxxx_0[j] + 1.5 * fl1_fx * tpx_yy_xxxx_0[j];

            tpy_yyyy_xxxx_0[j] = pa_y[j] * tpy_yyy_xxxx_0[j] + 1.5 * fl1_fx * tpy_yy_xxxx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxxx_0[j];

            tpz_yyyy_xxxx_0[j] = pa_y[j] * tpz_yyy_xxxx_0[j] + 1.5 * fl1_fx * tpz_yy_xxxx_0[j];

            tpx_yyyy_xxxy_0[j] = pa_y[j] * tpx_yyy_xxxy_0[j] + 1.5 * fl1_fx * tpx_yy_xxxy_0[j] + 0.5 * fl1_fx * tpx_yyy_xxx_0[j];

            tpy_yyyy_xxxy_0[j] =
                pa_y[j] * tpy_yyy_xxxy_0[j] + 1.5 * fl1_fx * tpy_yy_xxxy_0[j] + 0.5 * fl1_fx * tpy_yyy_xxx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxxy_0[j];

            tpz_yyyy_xxxy_0[j] = pa_y[j] * tpz_yyy_xxxy_0[j] + 1.5 * fl1_fx * tpz_yy_xxxy_0[j] + 0.5 * fl1_fx * tpz_yyy_xxx_0[j];

            tpx_yyyy_xxxz_0[j] = pa_y[j] * tpx_yyy_xxxz_0[j] + 1.5 * fl1_fx * tpx_yy_xxxz_0[j];

            tpy_yyyy_xxxz_0[j] = pa_y[j] * tpy_yyy_xxxz_0[j] + 1.5 * fl1_fx * tpy_yy_xxxz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxxz_0[j];

            tpz_yyyy_xxxz_0[j] = pa_y[j] * tpz_yyy_xxxz_0[j] + 1.5 * fl1_fx * tpz_yy_xxxz_0[j];

            tpx_yyyy_xxyy_0[j] = pa_y[j] * tpx_yyy_xxyy_0[j] + 1.5 * fl1_fx * tpx_yy_xxyy_0[j] + fl1_fx * tpx_yyy_xxy_0[j];

            tpy_yyyy_xxyy_0[j] =
                pa_y[j] * tpy_yyy_xxyy_0[j] + 1.5 * fl1_fx * tpy_yy_xxyy_0[j] + fl1_fx * tpy_yyy_xxy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxyy_0[j];

            tpz_yyyy_xxyy_0[j] = pa_y[j] * tpz_yyy_xxyy_0[j] + 1.5 * fl1_fx * tpz_yy_xxyy_0[j] + fl1_fx * tpz_yyy_xxy_0[j];

            tpx_yyyy_xxyz_0[j] = pa_y[j] * tpx_yyy_xxyz_0[j] + 1.5 * fl1_fx * tpx_yy_xxyz_0[j] + 0.5 * fl1_fx * tpx_yyy_xxz_0[j];

            tpy_yyyy_xxyz_0[j] =
                pa_y[j] * tpy_yyy_xxyz_0[j] + 1.5 * fl1_fx * tpy_yy_xxyz_0[j] + 0.5 * fl1_fx * tpy_yyy_xxz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxyz_0[j];

            tpz_yyyy_xxyz_0[j] = pa_y[j] * tpz_yyy_xxyz_0[j] + 1.5 * fl1_fx * tpz_yy_xxyz_0[j] + 0.5 * fl1_fx * tpz_yyy_xxz_0[j];

            tpx_yyyy_xxzz_0[j] = pa_y[j] * tpx_yyy_xxzz_0[j] + 1.5 * fl1_fx * tpx_yy_xxzz_0[j];

            tpy_yyyy_xxzz_0[j] = pa_y[j] * tpy_yyy_xxzz_0[j] + 1.5 * fl1_fx * tpy_yy_xxzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxzz_0[j];

            tpz_yyyy_xxzz_0[j] = pa_y[j] * tpz_yyy_xxzz_0[j] + 1.5 * fl1_fx * tpz_yy_xxzz_0[j];

            tpx_yyyy_xyyy_0[j] = pa_y[j] * tpx_yyy_xyyy_0[j] + 1.5 * fl1_fx * tpx_yy_xyyy_0[j] + 1.5 * fl1_fx * tpx_yyy_xyy_0[j];

            tpy_yyyy_xyyy_0[j] =
                pa_y[j] * tpy_yyy_xyyy_0[j] + 1.5 * fl1_fx * tpy_yy_xyyy_0[j] + 1.5 * fl1_fx * tpy_yyy_xyy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyyy_0[j];

            tpz_yyyy_xyyy_0[j] = pa_y[j] * tpz_yyy_xyyy_0[j] + 1.5 * fl1_fx * tpz_yy_xyyy_0[j] + 1.5 * fl1_fx * tpz_yyy_xyy_0[j];

            tpx_yyyy_xyyz_0[j] = pa_y[j] * tpx_yyy_xyyz_0[j] + 1.5 * fl1_fx * tpx_yy_xyyz_0[j] + fl1_fx * tpx_yyy_xyz_0[j];

            tpy_yyyy_xyyz_0[j] =
                pa_y[j] * tpy_yyy_xyyz_0[j] + 1.5 * fl1_fx * tpy_yy_xyyz_0[j] + fl1_fx * tpy_yyy_xyz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyyz_0[j];

            tpz_yyyy_xyyz_0[j] = pa_y[j] * tpz_yyy_xyyz_0[j] + 1.5 * fl1_fx * tpz_yy_xyyz_0[j] + fl1_fx * tpz_yyy_xyz_0[j];

            tpx_yyyy_xyzz_0[j] = pa_y[j] * tpx_yyy_xyzz_0[j] + 1.5 * fl1_fx * tpx_yy_xyzz_0[j] + 0.5 * fl1_fx * tpx_yyy_xzz_0[j];

            tpy_yyyy_xyzz_0[j] =
                pa_y[j] * tpy_yyy_xyzz_0[j] + 1.5 * fl1_fx * tpy_yy_xyzz_0[j] + 0.5 * fl1_fx * tpy_yyy_xzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyzz_0[j];

            tpz_yyyy_xyzz_0[j] = pa_y[j] * tpz_yyy_xyzz_0[j] + 1.5 * fl1_fx * tpz_yy_xyzz_0[j] + 0.5 * fl1_fx * tpz_yyy_xzz_0[j];

            tpx_yyyy_xzzz_0[j] = pa_y[j] * tpx_yyy_xzzz_0[j] + 1.5 * fl1_fx * tpx_yy_xzzz_0[j];

            tpy_yyyy_xzzz_0[j] = pa_y[j] * tpy_yyy_xzzz_0[j] + 1.5 * fl1_fx * tpy_yy_xzzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xzzz_0[j];

            tpz_yyyy_xzzz_0[j] = pa_y[j] * tpz_yyy_xzzz_0[j] + 1.5 * fl1_fx * tpz_yy_xzzz_0[j];

            tpx_yyyy_yyyy_0[j] = pa_y[j] * tpx_yyy_yyyy_0[j] + 1.5 * fl1_fx * tpx_yy_yyyy_0[j] + 2.0 * fl1_fx * tpx_yyy_yyy_0[j];

            tpy_yyyy_yyyy_0[j] =
                pa_y[j] * tpy_yyy_yyyy_0[j] + 1.5 * fl1_fx * tpy_yy_yyyy_0[j] + 2.0 * fl1_fx * tpy_yyy_yyy_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyyy_0[j];

            tpz_yyyy_yyyy_0[j] = pa_y[j] * tpz_yyy_yyyy_0[j] + 1.5 * fl1_fx * tpz_yy_yyyy_0[j] + 2.0 * fl1_fx * tpz_yyy_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_483_531(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (483,531)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tpx_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 101);

        auto tpy_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tpz_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tpx_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 102);

        auto tpy_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tpz_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tpx_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 103);

        auto tpy_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tpz_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tpx_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 104);

        auto tpy_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tpz_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tpx_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 105);

        auto tpy_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tpz_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tpx_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 106);

        auto tpy_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tpz_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tpx_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 107);

        auto tpy_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tpz_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tpx_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 108);

        auto tpy_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tpz_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tpx_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 109);

        auto tpy_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tpz_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tpx_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 110);

        auto tpy_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tpz_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tpx_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 111);

        auto tpy_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tpz_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tpx_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 112);

        auto tpy_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tpz_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tpx_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 113);

        auto tpy_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tpz_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tpx_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 114);

        auto tpy_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tpz_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tpx_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 115);

        auto tpy_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tpz_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tpx_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 116);

        auto tpy_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 116);

        auto tpz_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tpx_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 56);

        auto tpy_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tpz_yy_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tpx_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 57);

        auto tpy_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tpz_yy_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tpx_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 58);

        auto tpy_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tpz_yy_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tpx_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 59);

        auto tpy_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tpz_yy_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tpx_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 60);

        auto tpy_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tpz_yz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tpx_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 61);

        auto tpy_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tpz_yz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tpx_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 62);

        auto tpy_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tpz_yz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tpx_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 63);

        auto tpy_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tpz_yz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tpx_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 64);

        auto tpy_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tpz_yz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tpx_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 65);

        auto tpy_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tpz_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tpx_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 66);

        auto tpy_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tpz_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tpx_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 67);

        auto tpy_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tpz_yz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tpx_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 68);

        auto tpy_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tpz_yz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tpx_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 69);

        auto tpy_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tpz_yz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tpx_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 70);

        auto tpy_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tpz_yz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tpx_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 71);

        auto tpy_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tpz_yz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tpx_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 67);

        auto tpy_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tpz_yyy_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tpx_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 68);

        auto tpy_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tpz_yyy_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tpx_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 69);

        auto tpy_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tpz_yyy_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tpx_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 70);

        auto tpy_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tpz_yyz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tpx_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 71);

        auto tpy_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tpz_yyz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tpx_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 72);

        auto tpy_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tpz_yyz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tpx_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 73);

        auto tpy_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tpz_yyz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tpx_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 74);

        auto tpy_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tpz_yyz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tpx_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 75);

        auto tpy_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tpz_yyz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tpx_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 76);

        auto tpy_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tpz_yyz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tpx_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 77);

        auto tpy_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tpz_yyz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto ts_yyy_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 101);

        auto ts_yyy_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 102);

        auto ts_yyy_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 103);

        auto ts_yyy_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 104);

        auto ts_yyz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 105);

        auto ts_yyz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 106);

        auto ts_yyz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 107);

        auto ts_yyz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 108);

        auto ts_yyz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 109);

        auto ts_yyz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 110);

        auto ts_yyz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 111);

        auto ts_yyz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 112);

        auto ts_yyz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 113);

        auto ts_yyz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 114);

        auto ts_yyz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 115);

        auto ts_yyz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 116);

        // set up pointers to integrals

        auto tpx_yyyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 161);

        auto tpy_yyyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 161);

        auto tpz_yyyy_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 161);

        auto tpx_yyyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 162);

        auto tpy_yyyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 162);

        auto tpz_yyyy_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 162);

        auto tpx_yyyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 163);

        auto tpy_yyyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 163);

        auto tpz_yyyy_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 163);

        auto tpx_yyyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 164);

        auto tpy_yyyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 164);

        auto tpz_yyyy_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 164);

        auto tpx_yyyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 165);

        auto tpy_yyyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 165);

        auto tpz_yyyz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 165);

        auto tpx_yyyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 166);

        auto tpy_yyyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 166);

        auto tpz_yyyz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 166);

        auto tpx_yyyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 167);

        auto tpy_yyyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 167);

        auto tpz_yyyz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 167);

        auto tpx_yyyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 168);

        auto tpy_yyyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 168);

        auto tpz_yyyz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 168);

        auto tpx_yyyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 169);

        auto tpy_yyyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 169);

        auto tpz_yyyz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 169);

        auto tpx_yyyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 170);

        auto tpy_yyyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 170);

        auto tpz_yyyz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 170);

        auto tpx_yyyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 171);

        auto tpy_yyyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 171);

        auto tpz_yyyz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 171);

        auto tpx_yyyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 172);

        auto tpy_yyyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 172);

        auto tpz_yyyz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 172);

        auto tpx_yyyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 173);

        auto tpy_yyyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 173);

        auto tpz_yyyz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 173);

        auto tpx_yyyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 174);

        auto tpy_yyyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 174);

        auto tpz_yyyz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 174);

        auto tpx_yyyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 175);

        auto tpy_yyyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 175);

        auto tpz_yyyz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 175);

        auto tpx_yyyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 176);

        auto tpy_yyyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 176);

        auto tpz_yyyz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 176);

        // Batch of Integrals (483,531)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yy_yyyz_0, tpx_yy_yyzz_0, tpx_yy_yzzz_0, \
                                     tpx_yy_zzzz_0, tpx_yyy_yyyz_0, tpx_yyy_yyz_0, tpx_yyy_yyzz_0, tpx_yyy_yzz_0, \
                                     tpx_yyy_yzzz_0, tpx_yyy_zzz_0, tpx_yyy_zzzz_0, tpx_yyyy_yyyz_0, tpx_yyyy_yyzz_0, \
                                     tpx_yyyy_yzzz_0, tpx_yyyy_zzzz_0, tpx_yyyz_xxxx_0, tpx_yyyz_xxxy_0, tpx_yyyz_xxxz_0, \
                                     tpx_yyyz_xxyy_0, tpx_yyyz_xxyz_0, tpx_yyyz_xxzz_0, tpx_yyyz_xyyy_0, tpx_yyyz_xyyz_0, \
                                     tpx_yyyz_xyzz_0, tpx_yyyz_xzzz_0, tpx_yyyz_yyyy_0, tpx_yyyz_yyyz_0, tpx_yyz_xxx_0, \
                                     tpx_yyz_xxxx_0, tpx_yyz_xxxy_0, tpx_yyz_xxxz_0, tpx_yyz_xxy_0, tpx_yyz_xxyy_0, \
                                     tpx_yyz_xxyz_0, tpx_yyz_xxz_0, tpx_yyz_xxzz_0, tpx_yyz_xyy_0, tpx_yyz_xyyy_0, \
                                     tpx_yyz_xyyz_0, tpx_yyz_xyz_0, tpx_yyz_xyzz_0, tpx_yyz_xzz_0, tpx_yyz_xzzz_0, \
                                     tpx_yyz_yyy_0, tpx_yyz_yyyy_0, tpx_yyz_yyyz_0, tpx_yyz_yyz_0, tpx_yz_xxxx_0, \
                                     tpx_yz_xxxy_0, tpx_yz_xxxz_0, tpx_yz_xxyy_0, tpx_yz_xxyz_0, tpx_yz_xxzz_0, \
                                     tpx_yz_xyyy_0, tpx_yz_xyyz_0, tpx_yz_xyzz_0, tpx_yz_xzzz_0, tpx_yz_yyyy_0, \
                                     tpx_yz_yyyz_0, tpy_yy_yyyz_0, tpy_yy_yyzz_0, tpy_yy_yzzz_0, tpy_yy_zzzz_0, \
                                     tpy_yyy_yyyz_0, tpy_yyy_yyz_0, tpy_yyy_yyzz_0, tpy_yyy_yzz_0, tpy_yyy_yzzz_0, \
                                     tpy_yyy_zzz_0, tpy_yyy_zzzz_0, tpy_yyyy_yyyz_0, tpy_yyyy_yyzz_0, tpy_yyyy_yzzz_0, \
                                     tpy_yyyy_zzzz_0, tpy_yyyz_xxxx_0, tpy_yyyz_xxxy_0, tpy_yyyz_xxxz_0, tpy_yyyz_xxyy_0, \
                                     tpy_yyyz_xxyz_0, tpy_yyyz_xxzz_0, tpy_yyyz_xyyy_0, tpy_yyyz_xyyz_0, tpy_yyyz_xyzz_0, \
                                     tpy_yyyz_xzzz_0, tpy_yyyz_yyyy_0, tpy_yyyz_yyyz_0, tpy_yyz_xxx_0, tpy_yyz_xxxx_0, \
                                     tpy_yyz_xxxy_0, tpy_yyz_xxxz_0, tpy_yyz_xxy_0, tpy_yyz_xxyy_0, tpy_yyz_xxyz_0, \
                                     tpy_yyz_xxz_0, tpy_yyz_xxzz_0, tpy_yyz_xyy_0, tpy_yyz_xyyy_0, tpy_yyz_xyyz_0, \
                                     tpy_yyz_xyz_0, tpy_yyz_xyzz_0, tpy_yyz_xzz_0, tpy_yyz_xzzz_0, tpy_yyz_yyy_0, \
                                     tpy_yyz_yyyy_0, tpy_yyz_yyyz_0, tpy_yyz_yyz_0, tpy_yz_xxxx_0, tpy_yz_xxxy_0, \
                                     tpy_yz_xxxz_0, tpy_yz_xxyy_0, tpy_yz_xxyz_0, tpy_yz_xxzz_0, tpy_yz_xyyy_0, \
                                     tpy_yz_xyyz_0, tpy_yz_xyzz_0, tpy_yz_xzzz_0, tpy_yz_yyyy_0, tpy_yz_yyyz_0, \
                                     tpz_yy_yyyz_0, tpz_yy_yyzz_0, tpz_yy_yzzz_0, tpz_yy_zzzz_0, tpz_yyy_yyyz_0, \
                                     tpz_yyy_yyz_0, tpz_yyy_yyzz_0, tpz_yyy_yzz_0, tpz_yyy_yzzz_0, tpz_yyy_zzz_0, \
                                     tpz_yyy_zzzz_0, tpz_yyyy_yyyz_0, tpz_yyyy_yyzz_0, tpz_yyyy_yzzz_0, tpz_yyyy_zzzz_0, \
                                     tpz_yyyz_xxxx_0, tpz_yyyz_xxxy_0, tpz_yyyz_xxxz_0, tpz_yyyz_xxyy_0, tpz_yyyz_xxyz_0, \
                                     tpz_yyyz_xxzz_0, tpz_yyyz_xyyy_0, tpz_yyyz_xyyz_0, tpz_yyyz_xyzz_0, tpz_yyyz_xzzz_0, \
                                     tpz_yyyz_yyyy_0, tpz_yyyz_yyyz_0, tpz_yyz_xxx_0, tpz_yyz_xxxx_0, tpz_yyz_xxxy_0, \
                                     tpz_yyz_xxxz_0, tpz_yyz_xxy_0, tpz_yyz_xxyy_0, tpz_yyz_xxyz_0, tpz_yyz_xxz_0, \
                                     tpz_yyz_xxzz_0, tpz_yyz_xyy_0, tpz_yyz_xyyy_0, tpz_yyz_xyyz_0, tpz_yyz_xyz_0, \
                                     tpz_yyz_xyzz_0, tpz_yyz_xzz_0, tpz_yyz_xzzz_0, tpz_yyz_yyy_0, tpz_yyz_yyyy_0, \
                                     tpz_yyz_yyyz_0, tpz_yyz_yyz_0, tpz_yz_xxxx_0, tpz_yz_xxxy_0, tpz_yz_xxxz_0, \
                                     tpz_yz_xxyy_0, tpz_yz_xxyz_0, tpz_yz_xxzz_0, tpz_yz_xyyy_0, tpz_yz_xyyz_0, \
                                     tpz_yz_xyzz_0, tpz_yz_xzzz_0, tpz_yz_yyyy_0, tpz_yz_yyyz_0, ts_yyy_yyyz_0, \
                                     ts_yyy_yyzz_0, ts_yyy_yzzz_0, ts_yyy_zzzz_0, ts_yyz_xxxx_0, ts_yyz_xxxy_0, \
                                     ts_yyz_xxxz_0, ts_yyz_xxyy_0, ts_yyz_xxyz_0, ts_yyz_xxzz_0, ts_yyz_xyyy_0, \
                                     ts_yyz_xyyz_0, ts_yyz_xyzz_0, ts_yyz_xzzz_0, ts_yyz_yyyy_0, ts_yyz_yyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyyy_yyyz_0[j] = pa_y[j] * tpx_yyy_yyyz_0[j] + 1.5 * fl1_fx * tpx_yy_yyyz_0[j] + 1.5 * fl1_fx * tpx_yyy_yyz_0[j];

            tpy_yyyy_yyyz_0[j] =
                pa_y[j] * tpy_yyy_yyyz_0[j] + 1.5 * fl1_fx * tpy_yy_yyyz_0[j] + 1.5 * fl1_fx * tpy_yyy_yyz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyyz_0[j];

            tpz_yyyy_yyyz_0[j] = pa_y[j] * tpz_yyy_yyyz_0[j] + 1.5 * fl1_fx * tpz_yy_yyyz_0[j] + 1.5 * fl1_fx * tpz_yyy_yyz_0[j];

            tpx_yyyy_yyzz_0[j] = pa_y[j] * tpx_yyy_yyzz_0[j] + 1.5 * fl1_fx * tpx_yy_yyzz_0[j] + fl1_fx * tpx_yyy_yzz_0[j];

            tpy_yyyy_yyzz_0[j] =
                pa_y[j] * tpy_yyy_yyzz_0[j] + 1.5 * fl1_fx * tpy_yy_yyzz_0[j] + fl1_fx * tpy_yyy_yzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyzz_0[j];

            tpz_yyyy_yyzz_0[j] = pa_y[j] * tpz_yyy_yyzz_0[j] + 1.5 * fl1_fx * tpz_yy_yyzz_0[j] + fl1_fx * tpz_yyy_yzz_0[j];

            tpx_yyyy_yzzz_0[j] = pa_y[j] * tpx_yyy_yzzz_0[j] + 1.5 * fl1_fx * tpx_yy_yzzz_0[j] + 0.5 * fl1_fx * tpx_yyy_zzz_0[j];

            tpy_yyyy_yzzz_0[j] =
                pa_y[j] * tpy_yyy_yzzz_0[j] + 1.5 * fl1_fx * tpy_yy_yzzz_0[j] + 0.5 * fl1_fx * tpy_yyy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yzzz_0[j];

            tpz_yyyy_yzzz_0[j] = pa_y[j] * tpz_yyy_yzzz_0[j] + 1.5 * fl1_fx * tpz_yy_yzzz_0[j] + 0.5 * fl1_fx * tpz_yyy_zzz_0[j];

            tpx_yyyy_zzzz_0[j] = pa_y[j] * tpx_yyy_zzzz_0[j] + 1.5 * fl1_fx * tpx_yy_zzzz_0[j];

            tpy_yyyy_zzzz_0[j] = pa_y[j] * tpy_yyy_zzzz_0[j] + 1.5 * fl1_fx * tpy_yy_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_zzzz_0[j];

            tpz_yyyy_zzzz_0[j] = pa_y[j] * tpz_yyy_zzzz_0[j] + 1.5 * fl1_fx * tpz_yy_zzzz_0[j];

            tpx_yyyz_xxxx_0[j] = pa_y[j] * tpx_yyz_xxxx_0[j] + fl1_fx * tpx_yz_xxxx_0[j];

            tpy_yyyz_xxxx_0[j] = pa_y[j] * tpy_yyz_xxxx_0[j] + fl1_fx * tpy_yz_xxxx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxxx_0[j];

            tpz_yyyz_xxxx_0[j] = pa_y[j] * tpz_yyz_xxxx_0[j] + fl1_fx * tpz_yz_xxxx_0[j];

            tpx_yyyz_xxxy_0[j] = pa_y[j] * tpx_yyz_xxxy_0[j] + fl1_fx * tpx_yz_xxxy_0[j] + 0.5 * fl1_fx * tpx_yyz_xxx_0[j];

            tpy_yyyz_xxxy_0[j] =
                pa_y[j] * tpy_yyz_xxxy_0[j] + fl1_fx * tpy_yz_xxxy_0[j] + 0.5 * fl1_fx * tpy_yyz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxxy_0[j];

            tpz_yyyz_xxxy_0[j] = pa_y[j] * tpz_yyz_xxxy_0[j] + fl1_fx * tpz_yz_xxxy_0[j] + 0.5 * fl1_fx * tpz_yyz_xxx_0[j];

            tpx_yyyz_xxxz_0[j] = pa_y[j] * tpx_yyz_xxxz_0[j] + fl1_fx * tpx_yz_xxxz_0[j];

            tpy_yyyz_xxxz_0[j] = pa_y[j] * tpy_yyz_xxxz_0[j] + fl1_fx * tpy_yz_xxxz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxxz_0[j];

            tpz_yyyz_xxxz_0[j] = pa_y[j] * tpz_yyz_xxxz_0[j] + fl1_fx * tpz_yz_xxxz_0[j];

            tpx_yyyz_xxyy_0[j] = pa_y[j] * tpx_yyz_xxyy_0[j] + fl1_fx * tpx_yz_xxyy_0[j] + fl1_fx * tpx_yyz_xxy_0[j];

            tpy_yyyz_xxyy_0[j] =
                pa_y[j] * tpy_yyz_xxyy_0[j] + fl1_fx * tpy_yz_xxyy_0[j] + fl1_fx * tpy_yyz_xxy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxyy_0[j];

            tpz_yyyz_xxyy_0[j] = pa_y[j] * tpz_yyz_xxyy_0[j] + fl1_fx * tpz_yz_xxyy_0[j] + fl1_fx * tpz_yyz_xxy_0[j];

            tpx_yyyz_xxyz_0[j] = pa_y[j] * tpx_yyz_xxyz_0[j] + fl1_fx * tpx_yz_xxyz_0[j] + 0.5 * fl1_fx * tpx_yyz_xxz_0[j];

            tpy_yyyz_xxyz_0[j] =
                pa_y[j] * tpy_yyz_xxyz_0[j] + fl1_fx * tpy_yz_xxyz_0[j] + 0.5 * fl1_fx * tpy_yyz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxyz_0[j];

            tpz_yyyz_xxyz_0[j] = pa_y[j] * tpz_yyz_xxyz_0[j] + fl1_fx * tpz_yz_xxyz_0[j] + 0.5 * fl1_fx * tpz_yyz_xxz_0[j];

            tpx_yyyz_xxzz_0[j] = pa_y[j] * tpx_yyz_xxzz_0[j] + fl1_fx * tpx_yz_xxzz_0[j];

            tpy_yyyz_xxzz_0[j] = pa_y[j] * tpy_yyz_xxzz_0[j] + fl1_fx * tpy_yz_xxzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxzz_0[j];

            tpz_yyyz_xxzz_0[j] = pa_y[j] * tpz_yyz_xxzz_0[j] + fl1_fx * tpz_yz_xxzz_0[j];

            tpx_yyyz_xyyy_0[j] = pa_y[j] * tpx_yyz_xyyy_0[j] + fl1_fx * tpx_yz_xyyy_0[j] + 1.5 * fl1_fx * tpx_yyz_xyy_0[j];

            tpy_yyyz_xyyy_0[j] =
                pa_y[j] * tpy_yyz_xyyy_0[j] + fl1_fx * tpy_yz_xyyy_0[j] + 1.5 * fl1_fx * tpy_yyz_xyy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyyy_0[j];

            tpz_yyyz_xyyy_0[j] = pa_y[j] * tpz_yyz_xyyy_0[j] + fl1_fx * tpz_yz_xyyy_0[j] + 1.5 * fl1_fx * tpz_yyz_xyy_0[j];

            tpx_yyyz_xyyz_0[j] = pa_y[j] * tpx_yyz_xyyz_0[j] + fl1_fx * tpx_yz_xyyz_0[j] + fl1_fx * tpx_yyz_xyz_0[j];

            tpy_yyyz_xyyz_0[j] =
                pa_y[j] * tpy_yyz_xyyz_0[j] + fl1_fx * tpy_yz_xyyz_0[j] + fl1_fx * tpy_yyz_xyz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyyz_0[j];

            tpz_yyyz_xyyz_0[j] = pa_y[j] * tpz_yyz_xyyz_0[j] + fl1_fx * tpz_yz_xyyz_0[j] + fl1_fx * tpz_yyz_xyz_0[j];

            tpx_yyyz_xyzz_0[j] = pa_y[j] * tpx_yyz_xyzz_0[j] + fl1_fx * tpx_yz_xyzz_0[j] + 0.5 * fl1_fx * tpx_yyz_xzz_0[j];

            tpy_yyyz_xyzz_0[j] =
                pa_y[j] * tpy_yyz_xyzz_0[j] + fl1_fx * tpy_yz_xyzz_0[j] + 0.5 * fl1_fx * tpy_yyz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyzz_0[j];

            tpz_yyyz_xyzz_0[j] = pa_y[j] * tpz_yyz_xyzz_0[j] + fl1_fx * tpz_yz_xyzz_0[j] + 0.5 * fl1_fx * tpz_yyz_xzz_0[j];

            tpx_yyyz_xzzz_0[j] = pa_y[j] * tpx_yyz_xzzz_0[j] + fl1_fx * tpx_yz_xzzz_0[j];

            tpy_yyyz_xzzz_0[j] = pa_y[j] * tpy_yyz_xzzz_0[j] + fl1_fx * tpy_yz_xzzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xzzz_0[j];

            tpz_yyyz_xzzz_0[j] = pa_y[j] * tpz_yyz_xzzz_0[j] + fl1_fx * tpz_yz_xzzz_0[j];

            tpx_yyyz_yyyy_0[j] = pa_y[j] * tpx_yyz_yyyy_0[j] + fl1_fx * tpx_yz_yyyy_0[j] + 2.0 * fl1_fx * tpx_yyz_yyy_0[j];

            tpy_yyyz_yyyy_0[j] =
                pa_y[j] * tpy_yyz_yyyy_0[j] + fl1_fx * tpy_yz_yyyy_0[j] + 2.0 * fl1_fx * tpy_yyz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyyy_0[j];

            tpz_yyyz_yyyy_0[j] = pa_y[j] * tpz_yyz_yyyy_0[j] + fl1_fx * tpz_yz_yyyy_0[j] + 2.0 * fl1_fx * tpz_yyz_yyy_0[j];

            tpx_yyyz_yyyz_0[j] = pa_y[j] * tpx_yyz_yyyz_0[j] + fl1_fx * tpx_yz_yyyz_0[j] + 1.5 * fl1_fx * tpx_yyz_yyz_0[j];

            tpy_yyyz_yyyz_0[j] =
                pa_y[j] * tpy_yyz_yyyz_0[j] + fl1_fx * tpy_yz_yyyz_0[j] + 1.5 * fl1_fx * tpy_yyz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyyz_0[j];

            tpz_yyyz_yyyz_0[j] = pa_y[j] * tpz_yyz_yyyz_0[j] + fl1_fx * tpz_yz_yyyz_0[j] + 1.5 * fl1_fx * tpz_yyz_yyz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_531_579(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (531,579)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tpx_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 117);

        auto tpy_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tpz_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tpx_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 118);

        auto tpy_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tpz_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tpx_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 119);

        auto tpy_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tpz_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tpx_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 120);

        auto tpy_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tpz_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tpx_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 121);

        auto tpy_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tpz_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tpx_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 122);

        auto tpy_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tpz_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tpx_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 123);

        auto tpy_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tpz_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tpx_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 124);

        auto tpy_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tpz_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tpx_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 125);

        auto tpy_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tpz_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tpx_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 126);

        auto tpy_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tpz_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tpx_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 127);

        auto tpy_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tpz_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tpx_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 128);

        auto tpy_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tpz_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tpx_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 129);

        auto tpy_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tpz_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tpx_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 130);

        auto tpy_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tpz_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tpx_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 131);

        auto tpy_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tpz_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tpx_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 132);

        auto tpy_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tpz_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tpx_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 72);

        auto tpy_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tpz_yz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tpx_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 73);

        auto tpy_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tpz_yz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tpx_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 74);

        auto tpy_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tpz_yz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tpx_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 75);

        auto tpy_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tpz_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tpx_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 76);

        auto tpy_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tpz_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tpx_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 77);

        auto tpy_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tpz_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tpx_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 78);

        auto tpy_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tpz_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tpx_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 79);

        auto tpy_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tpz_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tpx_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 80);

        auto tpy_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tpz_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpy_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tpz_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tpx_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 84);

        auto tpy_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tpz_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tpx_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 85);

        auto tpy_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tpz_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tpx_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 86);

        auto tpy_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tpz_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tpx_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 87);

        auto tpy_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tpz_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tpx_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 78);

        auto tpy_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tpz_yyz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tpx_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 79);

        auto tpy_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tpz_yyz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tpx_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 80);

        auto tpy_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tpz_yzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tpx_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 81);

        auto tpy_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tpz_yzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tpx_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 82);

        auto tpy_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tpz_yzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tpx_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 83);

        auto tpy_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tpz_yzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tpx_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 84);

        auto tpy_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tpz_yzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tpx_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 85);

        auto tpy_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tpz_yzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tpx_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 86);

        auto tpy_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tpz_yzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tpx_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 87);

        auto tpy_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tpz_yzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tpx_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 88);

        auto tpy_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tpz_yzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto ts_yyz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 117);

        auto ts_yyz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 118);

        auto ts_yyz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 119);

        auto ts_yzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 120);

        auto ts_yzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 121);

        auto ts_yzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 122);

        auto ts_yzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 123);

        auto ts_yzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 124);

        auto ts_yzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 125);

        auto ts_yzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 126);

        auto ts_yzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 127);

        auto ts_yzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 128);

        auto ts_yzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 129);

        auto ts_yzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 130);

        auto ts_yzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 131);

        auto ts_yzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 132);

        // set up pointers to integrals

        auto tpx_yyyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 177);

        auto tpy_yyyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 177);

        auto tpz_yyyz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 177);

        auto tpx_yyyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 178);

        auto tpy_yyyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 178);

        auto tpz_yyyz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 178);

        auto tpx_yyyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 179);

        auto tpy_yyyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 179);

        auto tpz_yyyz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 179);

        auto tpx_yyzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 180);

        auto tpy_yyzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 180);

        auto tpz_yyzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 180);

        auto tpx_yyzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 181);

        auto tpy_yyzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 181);

        auto tpz_yyzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 181);

        auto tpx_yyzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 182);

        auto tpy_yyzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 182);

        auto tpz_yyzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 182);

        auto tpx_yyzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 183);

        auto tpy_yyzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 183);

        auto tpz_yyzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 183);

        auto tpx_yyzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 184);

        auto tpy_yyzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 184);

        auto tpz_yyzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 184);

        auto tpx_yyzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 185);

        auto tpy_yyzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 185);

        auto tpz_yyzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 185);

        auto tpx_yyzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 186);

        auto tpy_yyzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 186);

        auto tpz_yyzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 186);

        auto tpx_yyzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 187);

        auto tpy_yyzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 187);

        auto tpz_yyzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 187);

        auto tpx_yyzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 188);

        auto tpy_yyzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 188);

        auto tpz_yyzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 188);

        auto tpx_yyzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 189);

        auto tpy_yyzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 189);

        auto tpz_yyzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 189);

        auto tpx_yyzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 190);

        auto tpy_yyzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 190);

        auto tpz_yyzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 190);

        auto tpx_yyzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 191);

        auto tpy_yyzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 191);

        auto tpz_yyzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 191);

        auto tpx_yyzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 192);

        auto tpy_yyzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 192);

        auto tpz_yyzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 192);

        // Batch of Integrals (531,579)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yyyz_yyzz_0, tpx_yyyz_yzzz_0, tpx_yyyz_zzzz_0, \
                                     tpx_yyz_yyzz_0, tpx_yyz_yzz_0, tpx_yyz_yzzz_0, tpx_yyz_zzz_0, tpx_yyz_zzzz_0, \
                                     tpx_yyzz_xxxx_0, tpx_yyzz_xxxy_0, tpx_yyzz_xxxz_0, tpx_yyzz_xxyy_0, tpx_yyzz_xxyz_0, \
                                     tpx_yyzz_xxzz_0, tpx_yyzz_xyyy_0, tpx_yyzz_xyyz_0, tpx_yyzz_xyzz_0, tpx_yyzz_xzzz_0, \
                                     tpx_yyzz_yyyy_0, tpx_yyzz_yyyz_0, tpx_yyzz_yyzz_0, tpx_yz_yyzz_0, tpx_yz_yzzz_0, \
                                     tpx_yz_zzzz_0, tpx_yzz_xxx_0, tpx_yzz_xxxx_0, tpx_yzz_xxxy_0, tpx_yzz_xxxz_0, \
                                     tpx_yzz_xxy_0, tpx_yzz_xxyy_0, tpx_yzz_xxyz_0, tpx_yzz_xxz_0, tpx_yzz_xxzz_0, \
                                     tpx_yzz_xyy_0, tpx_yzz_xyyy_0, tpx_yzz_xyyz_0, tpx_yzz_xyz_0, tpx_yzz_xyzz_0, \
                                     tpx_yzz_xzz_0, tpx_yzz_xzzz_0, tpx_yzz_yyy_0, tpx_yzz_yyyy_0, tpx_yzz_yyyz_0, \
                                     tpx_yzz_yyz_0, tpx_yzz_yyzz_0, tpx_yzz_yzz_0, tpx_zz_xxxx_0, tpx_zz_xxxy_0, \
                                     tpx_zz_xxxz_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, tpx_zz_xxzz_0, tpx_zz_xyyy_0, \
                                     tpx_zz_xyyz_0, tpx_zz_xyzz_0, tpx_zz_xzzz_0, tpx_zz_yyyy_0, tpx_zz_yyyz_0, \
                                     tpx_zz_yyzz_0, tpy_yyyz_yyzz_0, tpy_yyyz_yzzz_0, tpy_yyyz_zzzz_0, tpy_yyz_yyzz_0, \
                                     tpy_yyz_yzz_0, tpy_yyz_yzzz_0, tpy_yyz_zzz_0, tpy_yyz_zzzz_0, tpy_yyzz_xxxx_0, \
                                     tpy_yyzz_xxxy_0, tpy_yyzz_xxxz_0, tpy_yyzz_xxyy_0, tpy_yyzz_xxyz_0, tpy_yyzz_xxzz_0, \
                                     tpy_yyzz_xyyy_0, tpy_yyzz_xyyz_0, tpy_yyzz_xyzz_0, tpy_yyzz_xzzz_0, tpy_yyzz_yyyy_0, \
                                     tpy_yyzz_yyyz_0, tpy_yyzz_yyzz_0, tpy_yz_yyzz_0, tpy_yz_yzzz_0, tpy_yz_zzzz_0, \
                                     tpy_yzz_xxx_0, tpy_yzz_xxxx_0, tpy_yzz_xxxy_0, tpy_yzz_xxxz_0, tpy_yzz_xxy_0, \
                                     tpy_yzz_xxyy_0, tpy_yzz_xxyz_0, tpy_yzz_xxz_0, tpy_yzz_xxzz_0, tpy_yzz_xyy_0, \
                                     tpy_yzz_xyyy_0, tpy_yzz_xyyz_0, tpy_yzz_xyz_0, tpy_yzz_xyzz_0, tpy_yzz_xzz_0, \
                                     tpy_yzz_xzzz_0, tpy_yzz_yyy_0, tpy_yzz_yyyy_0, tpy_yzz_yyyz_0, tpy_yzz_yyz_0, \
                                     tpy_yzz_yyzz_0, tpy_yzz_yzz_0, tpy_zz_xxxx_0, tpy_zz_xxxy_0, tpy_zz_xxxz_0, \
                                     tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxzz_0, tpy_zz_xyyy_0, tpy_zz_xyyz_0, \
                                     tpy_zz_xyzz_0, tpy_zz_xzzz_0, tpy_zz_yyyy_0, tpy_zz_yyyz_0, tpy_zz_yyzz_0, \
                                     tpz_yyyz_yyzz_0, tpz_yyyz_yzzz_0, tpz_yyyz_zzzz_0, tpz_yyz_yyzz_0, tpz_yyz_yzz_0, \
                                     tpz_yyz_yzzz_0, tpz_yyz_zzz_0, tpz_yyz_zzzz_0, tpz_yyzz_xxxx_0, tpz_yyzz_xxxy_0, \
                                     tpz_yyzz_xxxz_0, tpz_yyzz_xxyy_0, tpz_yyzz_xxyz_0, tpz_yyzz_xxzz_0, tpz_yyzz_xyyy_0, \
                                     tpz_yyzz_xyyz_0, tpz_yyzz_xyzz_0, tpz_yyzz_xzzz_0, tpz_yyzz_yyyy_0, tpz_yyzz_yyyz_0, \
                                     tpz_yyzz_yyzz_0, tpz_yz_yyzz_0, tpz_yz_yzzz_0, tpz_yz_zzzz_0, tpz_yzz_xxx_0, \
                                     tpz_yzz_xxxx_0, tpz_yzz_xxxy_0, tpz_yzz_xxxz_0, tpz_yzz_xxy_0, tpz_yzz_xxyy_0, \
                                     tpz_yzz_xxyz_0, tpz_yzz_xxz_0, tpz_yzz_xxzz_0, tpz_yzz_xyy_0, tpz_yzz_xyyy_0, \
                                     tpz_yzz_xyyz_0, tpz_yzz_xyz_0, tpz_yzz_xyzz_0, tpz_yzz_xzz_0, tpz_yzz_xzzz_0, \
                                     tpz_yzz_yyy_0, tpz_yzz_yyyy_0, tpz_yzz_yyyz_0, tpz_yzz_yyz_0, tpz_yzz_yyzz_0, \
                                     tpz_yzz_yzz_0, tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, tpz_zz_xxyy_0, \
                                     tpz_zz_xxyz_0, tpz_zz_xxzz_0, tpz_zz_xyyy_0, tpz_zz_xyyz_0, tpz_zz_xyzz_0, \
                                     tpz_zz_xzzz_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, tpz_zz_yyzz_0, ts_yyz_yyzz_0, \
                                     ts_yyz_yzzz_0, ts_yyz_zzzz_0, ts_yzz_xxxx_0, ts_yzz_xxxy_0, ts_yzz_xxxz_0, \
                                     ts_yzz_xxyy_0, ts_yzz_xxyz_0, ts_yzz_xxzz_0, ts_yzz_xyyy_0, ts_yzz_xyyz_0, \
                                     ts_yzz_xyzz_0, ts_yzz_xzzz_0, ts_yzz_yyyy_0, ts_yzz_yyyz_0, ts_yzz_yyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyyz_yyzz_0[j] = pa_y[j] * tpx_yyz_yyzz_0[j] + fl1_fx * tpx_yz_yyzz_0[j] + fl1_fx * tpx_yyz_yzz_0[j];

            tpy_yyyz_yyzz_0[j] =
                pa_y[j] * tpy_yyz_yyzz_0[j] + fl1_fx * tpy_yz_yyzz_0[j] + fl1_fx * tpy_yyz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyzz_0[j];

            tpz_yyyz_yyzz_0[j] = pa_y[j] * tpz_yyz_yyzz_0[j] + fl1_fx * tpz_yz_yyzz_0[j] + fl1_fx * tpz_yyz_yzz_0[j];

            tpx_yyyz_yzzz_0[j] = pa_y[j] * tpx_yyz_yzzz_0[j] + fl1_fx * tpx_yz_yzzz_0[j] + 0.5 * fl1_fx * tpx_yyz_zzz_0[j];

            tpy_yyyz_yzzz_0[j] =
                pa_y[j] * tpy_yyz_yzzz_0[j] + fl1_fx * tpy_yz_yzzz_0[j] + 0.5 * fl1_fx * tpy_yyz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yzzz_0[j];

            tpz_yyyz_yzzz_0[j] = pa_y[j] * tpz_yyz_yzzz_0[j] + fl1_fx * tpz_yz_yzzz_0[j] + 0.5 * fl1_fx * tpz_yyz_zzz_0[j];

            tpx_yyyz_zzzz_0[j] = pa_y[j] * tpx_yyz_zzzz_0[j] + fl1_fx * tpx_yz_zzzz_0[j];

            tpy_yyyz_zzzz_0[j] = pa_y[j] * tpy_yyz_zzzz_0[j] + fl1_fx * tpy_yz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_zzzz_0[j];

            tpz_yyyz_zzzz_0[j] = pa_y[j] * tpz_yyz_zzzz_0[j] + fl1_fx * tpz_yz_zzzz_0[j];

            tpx_yyzz_xxxx_0[j] = pa_y[j] * tpx_yzz_xxxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxxx_0[j];

            tpy_yyzz_xxxx_0[j] = pa_y[j] * tpy_yzz_xxxx_0[j] + 0.5 * fl1_fx * tpy_zz_xxxx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxxx_0[j];

            tpz_yyzz_xxxx_0[j] = pa_y[j] * tpz_yzz_xxxx_0[j] + 0.5 * fl1_fx * tpz_zz_xxxx_0[j];

            tpx_yyzz_xxxy_0[j] = pa_y[j] * tpx_yzz_xxxy_0[j] + 0.5 * fl1_fx * tpx_zz_xxxy_0[j] + 0.5 * fl1_fx * tpx_yzz_xxx_0[j];

            tpy_yyzz_xxxy_0[j] =
                pa_y[j] * tpy_yzz_xxxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxxy_0[j] + 0.5 * fl1_fx * tpy_yzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxxy_0[j];

            tpz_yyzz_xxxy_0[j] = pa_y[j] * tpz_yzz_xxxy_0[j] + 0.5 * fl1_fx * tpz_zz_xxxy_0[j] + 0.5 * fl1_fx * tpz_yzz_xxx_0[j];

            tpx_yyzz_xxxz_0[j] = pa_y[j] * tpx_yzz_xxxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxxz_0[j];

            tpy_yyzz_xxxz_0[j] = pa_y[j] * tpy_yzz_xxxz_0[j] + 0.5 * fl1_fx * tpy_zz_xxxz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxxz_0[j];

            tpz_yyzz_xxxz_0[j] = pa_y[j] * tpz_yzz_xxxz_0[j] + 0.5 * fl1_fx * tpz_zz_xxxz_0[j];

            tpx_yyzz_xxyy_0[j] = pa_y[j] * tpx_yzz_xxyy_0[j] + 0.5 * fl1_fx * tpx_zz_xxyy_0[j] + fl1_fx * tpx_yzz_xxy_0[j];

            tpy_yyzz_xxyy_0[j] =
                pa_y[j] * tpy_yzz_xxyy_0[j] + 0.5 * fl1_fx * tpy_zz_xxyy_0[j] + fl1_fx * tpy_yzz_xxy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxyy_0[j];

            tpz_yyzz_xxyy_0[j] = pa_y[j] * tpz_yzz_xxyy_0[j] + 0.5 * fl1_fx * tpz_zz_xxyy_0[j] + fl1_fx * tpz_yzz_xxy_0[j];

            tpx_yyzz_xxyz_0[j] = pa_y[j] * tpx_yzz_xxyz_0[j] + 0.5 * fl1_fx * tpx_zz_xxyz_0[j] + 0.5 * fl1_fx * tpx_yzz_xxz_0[j];

            tpy_yyzz_xxyz_0[j] =
                pa_y[j] * tpy_yzz_xxyz_0[j] + 0.5 * fl1_fx * tpy_zz_xxyz_0[j] + 0.5 * fl1_fx * tpy_yzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxyz_0[j];

            tpz_yyzz_xxyz_0[j] = pa_y[j] * tpz_yzz_xxyz_0[j] + 0.5 * fl1_fx * tpz_zz_xxyz_0[j] + 0.5 * fl1_fx * tpz_yzz_xxz_0[j];

            tpx_yyzz_xxzz_0[j] = pa_y[j] * tpx_yzz_xxzz_0[j] + 0.5 * fl1_fx * tpx_zz_xxzz_0[j];

            tpy_yyzz_xxzz_0[j] = pa_y[j] * tpy_yzz_xxzz_0[j] + 0.5 * fl1_fx * tpy_zz_xxzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxzz_0[j];

            tpz_yyzz_xxzz_0[j] = pa_y[j] * tpz_yzz_xxzz_0[j] + 0.5 * fl1_fx * tpz_zz_xxzz_0[j];

            tpx_yyzz_xyyy_0[j] = pa_y[j] * tpx_yzz_xyyy_0[j] + 0.5 * fl1_fx * tpx_zz_xyyy_0[j] + 1.5 * fl1_fx * tpx_yzz_xyy_0[j];

            tpy_yyzz_xyyy_0[j] =
                pa_y[j] * tpy_yzz_xyyy_0[j] + 0.5 * fl1_fx * tpy_zz_xyyy_0[j] + 1.5 * fl1_fx * tpy_yzz_xyy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyyy_0[j];

            tpz_yyzz_xyyy_0[j] = pa_y[j] * tpz_yzz_xyyy_0[j] + 0.5 * fl1_fx * tpz_zz_xyyy_0[j] + 1.5 * fl1_fx * tpz_yzz_xyy_0[j];

            tpx_yyzz_xyyz_0[j] = pa_y[j] * tpx_yzz_xyyz_0[j] + 0.5 * fl1_fx * tpx_zz_xyyz_0[j] + fl1_fx * tpx_yzz_xyz_0[j];

            tpy_yyzz_xyyz_0[j] =
                pa_y[j] * tpy_yzz_xyyz_0[j] + 0.5 * fl1_fx * tpy_zz_xyyz_0[j] + fl1_fx * tpy_yzz_xyz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyyz_0[j];

            tpz_yyzz_xyyz_0[j] = pa_y[j] * tpz_yzz_xyyz_0[j] + 0.5 * fl1_fx * tpz_zz_xyyz_0[j] + fl1_fx * tpz_yzz_xyz_0[j];

            tpx_yyzz_xyzz_0[j] = pa_y[j] * tpx_yzz_xyzz_0[j] + 0.5 * fl1_fx * tpx_zz_xyzz_0[j] + 0.5 * fl1_fx * tpx_yzz_xzz_0[j];

            tpy_yyzz_xyzz_0[j] =
                pa_y[j] * tpy_yzz_xyzz_0[j] + 0.5 * fl1_fx * tpy_zz_xyzz_0[j] + 0.5 * fl1_fx * tpy_yzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyzz_0[j];

            tpz_yyzz_xyzz_0[j] = pa_y[j] * tpz_yzz_xyzz_0[j] + 0.5 * fl1_fx * tpz_zz_xyzz_0[j] + 0.5 * fl1_fx * tpz_yzz_xzz_0[j];

            tpx_yyzz_xzzz_0[j] = pa_y[j] * tpx_yzz_xzzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzzz_0[j];

            tpy_yyzz_xzzz_0[j] = pa_y[j] * tpy_yzz_xzzz_0[j] + 0.5 * fl1_fx * tpy_zz_xzzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xzzz_0[j];

            tpz_yyzz_xzzz_0[j] = pa_y[j] * tpz_yzz_xzzz_0[j] + 0.5 * fl1_fx * tpz_zz_xzzz_0[j];

            tpx_yyzz_yyyy_0[j] = pa_y[j] * tpx_yzz_yyyy_0[j] + 0.5 * fl1_fx * tpx_zz_yyyy_0[j] + 2.0 * fl1_fx * tpx_yzz_yyy_0[j];

            tpy_yyzz_yyyy_0[j] =
                pa_y[j] * tpy_yzz_yyyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyyy_0[j] + 2.0 * fl1_fx * tpy_yzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyyy_0[j];

            tpz_yyzz_yyyy_0[j] = pa_y[j] * tpz_yzz_yyyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyyy_0[j] + 2.0 * fl1_fx * tpz_yzz_yyy_0[j];

            tpx_yyzz_yyyz_0[j] = pa_y[j] * tpx_yzz_yyyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyyz_0[j] + 1.5 * fl1_fx * tpx_yzz_yyz_0[j];

            tpy_yyzz_yyyz_0[j] =
                pa_y[j] * tpy_yzz_yyyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyyz_0[j] + 1.5 * fl1_fx * tpy_yzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyyz_0[j];

            tpz_yyzz_yyyz_0[j] = pa_y[j] * tpz_yzz_yyyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyyz_0[j] + 1.5 * fl1_fx * tpz_yzz_yyz_0[j];

            tpx_yyzz_yyzz_0[j] = pa_y[j] * tpx_yzz_yyzz_0[j] + 0.5 * fl1_fx * tpx_zz_yyzz_0[j] + fl1_fx * tpx_yzz_yzz_0[j];

            tpy_yyzz_yyzz_0[j] =
                pa_y[j] * tpy_yzz_yyzz_0[j] + 0.5 * fl1_fx * tpy_zz_yyzz_0[j] + fl1_fx * tpy_yzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyzz_0[j];

            tpz_yyzz_yyzz_0[j] = pa_y[j] * tpz_yzz_yyzz_0[j] + 0.5 * fl1_fx * tpz_zz_yyzz_0[j] + fl1_fx * tpz_yzz_yzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_579_627(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (579,627)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tpx_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 133);

        auto tpy_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tpz_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tpx_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 134);

        auto tpy_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tpz_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tpx_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 135);

        auto tpy_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tpz_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tpx_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 136);

        auto tpy_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tpz_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tpx_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 137);

        auto tpy_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tpz_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tpx_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 138);

        auto tpy_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tpz_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tpx_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 139);

        auto tpy_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tpz_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tpx_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 140);

        auto tpy_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tpz_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tpx_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 141);

        auto tpy_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tpz_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tpx_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 142);

        auto tpy_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tpz_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tpx_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 143);

        auto tpy_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tpz_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tpx_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 144);

        auto tpy_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tpz_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tpx_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 145);

        auto tpy_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tpz_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tpx_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 146);

        auto tpy_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tpz_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tpx_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 147);

        auto tpy_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tpz_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tpx_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 148);

        auto tpy_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tpz_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tpx_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 88);

        auto tpy_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tpz_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tpx_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 89);

        auto tpy_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tpx_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 89);

        auto tpy_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tpz_yzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tpx_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 90);

        auto tpy_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tpz_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tpx_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 91);

        auto tpy_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tpz_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tpx_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 92);

        auto tpy_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tpz_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tpx_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 93);

        auto tpy_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tpz_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tpx_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 94);

        auto tpy_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tpz_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tpx_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 95);

        auto tpy_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tpz_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tpx_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 96);

        auto tpy_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tpz_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tpx_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 97);

        auto tpy_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tpz_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tpx_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 98);

        auto tpy_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tpz_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tpx_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 99);

        auto tpy_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tpz_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto ts_yzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 133);

        auto ts_yzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 134);

        auto ts_zzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 135);

        auto ts_zzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 136);

        auto ts_zzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 137);

        auto ts_zzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 138);

        auto ts_zzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 139);

        auto ts_zzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 140);

        auto ts_zzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 141);

        auto ts_zzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 142);

        auto ts_zzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 143);

        auto ts_zzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 144);

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        // set up pointers to integrals

        auto tpx_yyzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 193);

        auto tpy_yyzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 193);

        auto tpz_yyzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 193);

        auto tpx_yyzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 194);

        auto tpy_yyzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 194);

        auto tpz_yyzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 194);

        auto tpx_yzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 195);

        auto tpy_yzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 195);

        auto tpz_yzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 195);

        auto tpx_yzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 196);

        auto tpy_yzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 196);

        auto tpz_yzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 196);

        auto tpx_yzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 197);

        auto tpy_yzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 197);

        auto tpz_yzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 197);

        auto tpx_yzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 198);

        auto tpy_yzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 198);

        auto tpz_yzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 198);

        auto tpx_yzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 199);

        auto tpy_yzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 199);

        auto tpz_yzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 199);

        auto tpx_yzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 200);

        auto tpy_yzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 200);

        auto tpz_yzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 200);

        auto tpx_yzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 201);

        auto tpy_yzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 201);

        auto tpz_yzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 201);

        auto tpx_yzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 202);

        auto tpy_yzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 202);

        auto tpz_yzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 202);

        auto tpx_yzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 203);

        auto tpy_yzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 203);

        auto tpz_yzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 203);

        auto tpx_yzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 204);

        auto tpy_yzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 204);

        auto tpz_yzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 204);

        auto tpx_yzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 205);

        auto tpy_yzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 205);

        auto tpz_yzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 205);

        auto tpx_yzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 206);

        auto tpy_yzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 206);

        auto tpz_yzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 206);

        auto tpx_yzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 207);

        auto tpy_yzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 207);

        auto tpz_yzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 207);

        auto tpx_yzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 208);

        auto tpy_yzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 208);

        auto tpz_yzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 208);

        // Batch of Integrals (579,627)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yyzz_yzzz_0, tpx_yyzz_zzzz_0, tpx_yzz_yzzz_0, \
                                     tpx_yzz_zzz_0, tpx_yzz_zzzz_0, tpx_yzzz_xxxx_0, tpx_yzzz_xxxy_0, tpx_yzzz_xxxz_0, \
                                     tpx_yzzz_xxyy_0, tpx_yzzz_xxyz_0, tpx_yzzz_xxzz_0, tpx_yzzz_xyyy_0, tpx_yzzz_xyyz_0, \
                                     tpx_yzzz_xyzz_0, tpx_yzzz_xzzz_0, tpx_yzzz_yyyy_0, tpx_yzzz_yyyz_0, tpx_yzzz_yyzz_0, \
                                     tpx_yzzz_yzzz_0, tpx_zz_yzzz_0, tpx_zz_zzzz_0, tpx_zzz_xxx_0, tpx_zzz_xxxx_0, \
                                     tpx_zzz_xxxy_0, tpx_zzz_xxxz_0, tpx_zzz_xxy_0, tpx_zzz_xxyy_0, tpx_zzz_xxyz_0, \
                                     tpx_zzz_xxz_0, tpx_zzz_xxzz_0, tpx_zzz_xyy_0, tpx_zzz_xyyy_0, tpx_zzz_xyyz_0, \
                                     tpx_zzz_xyz_0, tpx_zzz_xyzz_0, tpx_zzz_xzz_0, tpx_zzz_xzzz_0, tpx_zzz_yyy_0, \
                                     tpx_zzz_yyyy_0, tpx_zzz_yyyz_0, tpx_zzz_yyz_0, tpx_zzz_yyzz_0, tpx_zzz_yzz_0, \
                                     tpx_zzz_yzzz_0, tpx_zzz_zzz_0, tpy_yyzz_yzzz_0, tpy_yyzz_zzzz_0, tpy_yzz_yzzz_0, \
                                     tpy_yzz_zzz_0, tpy_yzz_zzzz_0, tpy_yzzz_xxxx_0, tpy_yzzz_xxxy_0, tpy_yzzz_xxxz_0, \
                                     tpy_yzzz_xxyy_0, tpy_yzzz_xxyz_0, tpy_yzzz_xxzz_0, tpy_yzzz_xyyy_0, tpy_yzzz_xyyz_0, \
                                     tpy_yzzz_xyzz_0, tpy_yzzz_xzzz_0, tpy_yzzz_yyyy_0, tpy_yzzz_yyyz_0, tpy_yzzz_yyzz_0, \
                                     tpy_yzzz_yzzz_0, tpy_zz_yzzz_0, tpy_zz_zzzz_0, tpy_zzz_xxx_0, tpy_zzz_xxxx_0, \
                                     tpy_zzz_xxxy_0, tpy_zzz_xxxz_0, tpy_zzz_xxy_0, tpy_zzz_xxyy_0, tpy_zzz_xxyz_0, \
                                     tpy_zzz_xxz_0, tpy_zzz_xxzz_0, tpy_zzz_xyy_0, tpy_zzz_xyyy_0, tpy_zzz_xyyz_0, \
                                     tpy_zzz_xyz_0, tpy_zzz_xyzz_0, tpy_zzz_xzz_0, tpy_zzz_xzzz_0, tpy_zzz_yyy_0, \
                                     tpy_zzz_yyyy_0, tpy_zzz_yyyz_0, tpy_zzz_yyz_0, tpy_zzz_yyzz_0, tpy_zzz_yzz_0, \
                                     tpy_zzz_yzzz_0, tpy_zzz_zzz_0, tpz_yyzz_yzzz_0, tpz_yyzz_zzzz_0, tpz_yzz_yzzz_0, \
                                     tpz_yzz_zzz_0, tpz_yzz_zzzz_0, tpz_yzzz_xxxx_0, tpz_yzzz_xxxy_0, tpz_yzzz_xxxz_0, \
                                     tpz_yzzz_xxyy_0, tpz_yzzz_xxyz_0, tpz_yzzz_xxzz_0, tpz_yzzz_xyyy_0, tpz_yzzz_xyyz_0, \
                                     tpz_yzzz_xyzz_0, tpz_yzzz_xzzz_0, tpz_yzzz_yyyy_0, tpz_yzzz_yyyz_0, tpz_yzzz_yyzz_0, \
                                     tpz_yzzz_yzzz_0, tpz_zz_yzzz_0, tpz_zz_zzzz_0, tpz_zzz_xxx_0, tpz_zzz_xxxx_0, \
                                     tpz_zzz_xxxy_0, tpz_zzz_xxxz_0, tpz_zzz_xxy_0, tpz_zzz_xxyy_0, tpz_zzz_xxyz_0, \
                                     tpz_zzz_xxz_0, tpz_zzz_xxzz_0, tpz_zzz_xyy_0, tpz_zzz_xyyy_0, tpz_zzz_xyyz_0, \
                                     tpz_zzz_xyz_0, tpz_zzz_xyzz_0, tpz_zzz_xzz_0, tpz_zzz_xzzz_0, tpz_zzz_yyy_0, \
                                     tpz_zzz_yyyy_0, tpz_zzz_yyyz_0, tpz_zzz_yyz_0, tpz_zzz_yyzz_0, tpz_zzz_yzz_0, \
                                     tpz_zzz_yzzz_0, tpz_zzz_zzz_0, ts_yzz_yzzz_0, ts_yzz_zzzz_0, ts_zzz_xxxx_0, \
                                     ts_zzz_xxxy_0, ts_zzz_xxxz_0, ts_zzz_xxyy_0, ts_zzz_xxyz_0, ts_zzz_xxzz_0, \
                                     ts_zzz_xyyy_0, ts_zzz_xyyz_0, ts_zzz_xyzz_0, ts_zzz_xzzz_0, ts_zzz_yyyy_0, \
                                     ts_zzz_yyyz_0, ts_zzz_yyzz_0, ts_zzz_yzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyzz_yzzz_0[j] = pa_y[j] * tpx_yzz_yzzz_0[j] + 0.5 * fl1_fx * tpx_zz_yzzz_0[j] + 0.5 * fl1_fx * tpx_yzz_zzz_0[j];

            tpy_yyzz_yzzz_0[j] =
                pa_y[j] * tpy_yzz_yzzz_0[j] + 0.5 * fl1_fx * tpy_zz_yzzz_0[j] + 0.5 * fl1_fx * tpy_yzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yzzz_0[j];

            tpz_yyzz_yzzz_0[j] = pa_y[j] * tpz_yzz_yzzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzzz_0[j] + 0.5 * fl1_fx * tpz_yzz_zzz_0[j];

            tpx_yyzz_zzzz_0[j] = pa_y[j] * tpx_yzz_zzzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzzz_0[j];

            tpy_yyzz_zzzz_0[j] = pa_y[j] * tpy_yzz_zzzz_0[j] + 0.5 * fl1_fx * tpy_zz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_zzzz_0[j];

            tpz_yyzz_zzzz_0[j] = pa_y[j] * tpz_yzz_zzzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzzz_0[j];

            tpx_yzzz_xxxx_0[j] = pa_y[j] * tpx_zzz_xxxx_0[j];

            tpy_yzzz_xxxx_0[j] = pa_y[j] * tpy_zzz_xxxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxx_0[j];

            tpz_yzzz_xxxx_0[j] = pa_y[j] * tpz_zzz_xxxx_0[j];

            tpx_yzzz_xxxy_0[j] = pa_y[j] * tpx_zzz_xxxy_0[j] + 0.5 * fl1_fx * tpx_zzz_xxx_0[j];

            tpy_yzzz_xxxy_0[j] = pa_y[j] * tpy_zzz_xxxy_0[j] + 0.5 * fl1_fx * tpy_zzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxy_0[j];

            tpz_yzzz_xxxy_0[j] = pa_y[j] * tpz_zzz_xxxy_0[j] + 0.5 * fl1_fx * tpz_zzz_xxx_0[j];

            tpx_yzzz_xxxz_0[j] = pa_y[j] * tpx_zzz_xxxz_0[j];

            tpy_yzzz_xxxz_0[j] = pa_y[j] * tpy_zzz_xxxz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxz_0[j];

            tpz_yzzz_xxxz_0[j] = pa_y[j] * tpz_zzz_xxxz_0[j];

            tpx_yzzz_xxyy_0[j] = pa_y[j] * tpx_zzz_xxyy_0[j] + fl1_fx * tpx_zzz_xxy_0[j];

            tpy_yzzz_xxyy_0[j] = pa_y[j] * tpy_zzz_xxyy_0[j] + fl1_fx * tpy_zzz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxyy_0[j];

            tpz_yzzz_xxyy_0[j] = pa_y[j] * tpz_zzz_xxyy_0[j] + fl1_fx * tpz_zzz_xxy_0[j];

            tpx_yzzz_xxyz_0[j] = pa_y[j] * tpx_zzz_xxyz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxz_0[j];

            tpy_yzzz_xxyz_0[j] = pa_y[j] * tpy_zzz_xxyz_0[j] + 0.5 * fl1_fx * tpy_zzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxyz_0[j];

            tpz_yzzz_xxyz_0[j] = pa_y[j] * tpz_zzz_xxyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxz_0[j];

            tpx_yzzz_xxzz_0[j] = pa_y[j] * tpx_zzz_xxzz_0[j];

            tpy_yzzz_xxzz_0[j] = pa_y[j] * tpy_zzz_xxzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxzz_0[j];

            tpz_yzzz_xxzz_0[j] = pa_y[j] * tpz_zzz_xxzz_0[j];

            tpx_yzzz_xyyy_0[j] = pa_y[j] * tpx_zzz_xyyy_0[j] + 1.5 * fl1_fx * tpx_zzz_xyy_0[j];

            tpy_yzzz_xyyy_0[j] = pa_y[j] * tpy_zzz_xyyy_0[j] + 1.5 * fl1_fx * tpy_zzz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyyy_0[j];

            tpz_yzzz_xyyy_0[j] = pa_y[j] * tpz_zzz_xyyy_0[j] + 1.5 * fl1_fx * tpz_zzz_xyy_0[j];

            tpx_yzzz_xyyz_0[j] = pa_y[j] * tpx_zzz_xyyz_0[j] + fl1_fx * tpx_zzz_xyz_0[j];

            tpy_yzzz_xyyz_0[j] = pa_y[j] * tpy_zzz_xyyz_0[j] + fl1_fx * tpy_zzz_xyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyyz_0[j];

            tpz_yzzz_xyyz_0[j] = pa_y[j] * tpz_zzz_xyyz_0[j] + fl1_fx * tpz_zzz_xyz_0[j];

            tpx_yzzz_xyzz_0[j] = pa_y[j] * tpx_zzz_xyzz_0[j] + 0.5 * fl1_fx * tpx_zzz_xzz_0[j];

            tpy_yzzz_xyzz_0[j] = pa_y[j] * tpy_zzz_xyzz_0[j] + 0.5 * fl1_fx * tpy_zzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyzz_0[j];

            tpz_yzzz_xyzz_0[j] = pa_y[j] * tpz_zzz_xyzz_0[j] + 0.5 * fl1_fx * tpz_zzz_xzz_0[j];

            tpx_yzzz_xzzz_0[j] = pa_y[j] * tpx_zzz_xzzz_0[j];

            tpy_yzzz_xzzz_0[j] = pa_y[j] * tpy_zzz_xzzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xzzz_0[j];

            tpz_yzzz_xzzz_0[j] = pa_y[j] * tpz_zzz_xzzz_0[j];

            tpx_yzzz_yyyy_0[j] = pa_y[j] * tpx_zzz_yyyy_0[j] + 2.0 * fl1_fx * tpx_zzz_yyy_0[j];

            tpy_yzzz_yyyy_0[j] = pa_y[j] * tpy_zzz_yyyy_0[j] + 2.0 * fl1_fx * tpy_zzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyyy_0[j];

            tpz_yzzz_yyyy_0[j] = pa_y[j] * tpz_zzz_yyyy_0[j] + 2.0 * fl1_fx * tpz_zzz_yyy_0[j];

            tpx_yzzz_yyyz_0[j] = pa_y[j] * tpx_zzz_yyyz_0[j] + 1.5 * fl1_fx * tpx_zzz_yyz_0[j];

            tpy_yzzz_yyyz_0[j] = pa_y[j] * tpy_zzz_yyyz_0[j] + 1.5 * fl1_fx * tpy_zzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyyz_0[j];

            tpz_yzzz_yyyz_0[j] = pa_y[j] * tpz_zzz_yyyz_0[j] + 1.5 * fl1_fx * tpz_zzz_yyz_0[j];

            tpx_yzzz_yyzz_0[j] = pa_y[j] * tpx_zzz_yyzz_0[j] + fl1_fx * tpx_zzz_yzz_0[j];

            tpy_yzzz_yyzz_0[j] = pa_y[j] * tpy_zzz_yyzz_0[j] + fl1_fx * tpy_zzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyzz_0[j];

            tpz_yzzz_yyzz_0[j] = pa_y[j] * tpz_zzz_yyzz_0[j] + fl1_fx * tpz_zzz_yzz_0[j];

            tpx_yzzz_yzzz_0[j] = pa_y[j] * tpx_zzz_yzzz_0[j] + 0.5 * fl1_fx * tpx_zzz_zzz_0[j];

            tpy_yzzz_yzzz_0[j] = pa_y[j] * tpy_zzz_yzzz_0[j] + 0.5 * fl1_fx * tpy_zzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yzzz_0[j];

            tpz_yzzz_yzzz_0[j] = pa_y[j] * tpz_zzz_yzzz_0[j] + 0.5 * fl1_fx * tpz_zzz_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGG_627_675(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (627,675)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tpx_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 135);

        auto tpy_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tpz_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tpx_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 136);

        auto tpy_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tpz_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tpx_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 137);

        auto tpy_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tpz_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tpx_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 138);

        auto tpy_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tpz_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tpx_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 139);

        auto tpy_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tpz_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tpx_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 140);

        auto tpy_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tpz_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tpx_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 141);

        auto tpy_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tpz_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tpx_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 142);

        auto tpy_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tpz_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tpx_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 143);

        auto tpy_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tpz_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tpx_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 144);

        auto tpy_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tpz_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tpx_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 145);

        auto tpy_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tpz_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tpx_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 146);

        auto tpy_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tpz_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tpx_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 147);

        auto tpy_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tpz_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tpx_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 148);

        auto tpy_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tpz_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tpx_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 149);

        auto tpy_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tpz_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 149);

        auto tpx_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 75);

        auto tpy_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tpz_zz_xxxx_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tpx_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 76);

        auto tpy_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tpz_zz_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tpx_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 77);

        auto tpy_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tpz_zz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tpx_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 78);

        auto tpy_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tpz_zz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tpx_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 79);

        auto tpy_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tpz_zz_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tpx_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 80);

        auto tpy_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tpz_zz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpy_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tpz_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tpx_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 84);

        auto tpy_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tpz_zz_xzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tpx_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 85);

        auto tpy_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tpz_zz_yyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tpx_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 86);

        auto tpy_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tpz_zz_yyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tpx_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 87);

        auto tpy_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tpz_zz_yyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tpx_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 88);

        auto tpy_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tpz_zz_yzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tpx_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 89);

        auto tpy_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tpz_zz_zzzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tpx_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 90);

        auto tpy_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tpz_zzz_xxx_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tpx_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 91);

        auto tpy_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tpz_zzz_xxy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tpx_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 92);

        auto tpy_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tpz_zzz_xxz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tpx_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 93);

        auto tpy_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tpz_zzz_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tpx_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 94);

        auto tpy_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tpz_zzz_xyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tpx_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 95);

        auto tpy_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tpz_zzz_xzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tpx_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 96);

        auto tpy_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tpz_zzz_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tpx_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 97);

        auto tpy_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tpz_zzz_yyz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tpx_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 98);

        auto tpy_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tpz_zzz_yzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tpx_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * idx + 99);

        auto tpy_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tpz_zzz_zzz_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto ts_zzz_xxxx_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 135);

        auto ts_zzz_xxxy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 136);

        auto ts_zzz_xxxz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 137);

        auto ts_zzz_xxyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 138);

        auto ts_zzz_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 139);

        auto ts_zzz_xxzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 140);

        auto ts_zzz_xyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 141);

        auto ts_zzz_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 142);

        auto ts_zzz_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 143);

        auto ts_zzz_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 144);

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        auto ts_zzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 149);

        // set up pointers to integrals

        auto tpx_yzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 209);

        auto tpy_yzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 209);

        auto tpz_yzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 209);

        auto tpx_zzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 210);

        auto tpy_zzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 210);

        auto tpz_zzzz_xxxx_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 210);

        auto tpx_zzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 211);

        auto tpy_zzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 211);

        auto tpz_zzzz_xxxy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 211);

        auto tpx_zzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 212);

        auto tpy_zzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 212);

        auto tpz_zzzz_xxxz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 212);

        auto tpx_zzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 213);

        auto tpy_zzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 213);

        auto tpz_zzzz_xxyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 213);

        auto tpx_zzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 214);

        auto tpy_zzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 214);

        auto tpz_zzzz_xxyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 214);

        auto tpx_zzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 215);

        auto tpy_zzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 215);

        auto tpz_zzzz_xxzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 215);

        auto tpx_zzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 216);

        auto tpy_zzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 216);

        auto tpz_zzzz_xyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 216);

        auto tpx_zzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 217);

        auto tpy_zzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 217);

        auto tpz_zzzz_xyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 217);

        auto tpx_zzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 218);

        auto tpy_zzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 218);

        auto tpz_zzzz_xyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 218);

        auto tpx_zzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 219);

        auto tpy_zzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 219);

        auto tpz_zzzz_xzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 219);

        auto tpx_zzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 220);

        auto tpy_zzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 220);

        auto tpz_zzzz_yyyy_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 220);

        auto tpx_zzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 221);

        auto tpy_zzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 221);

        auto tpz_zzzz_yyyz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 221);

        auto tpx_zzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 222);

        auto tpy_zzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 222);

        auto tpz_zzzz_yyzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 222);

        auto tpx_zzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 223);

        auto tpy_zzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 223);

        auto tpz_zzzz_yzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 223);

        auto tpx_zzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * idx + 224);

        auto tpy_zzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 225 * bdim + 225 * idx + 224);

        auto tpz_zzzz_zzzz_0 = primBuffer.data(pidx_p_4_4_m0 + 450 * bdim + 225 * idx + 224);

        // Batch of Integrals (627,675)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yzzz_zzzz_0, tpx_zz_xxxx_0, tpx_zz_xxxy_0, \
                                     tpx_zz_xxxz_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, tpx_zz_xxzz_0, tpx_zz_xyyy_0, \
                                     tpx_zz_xyyz_0, tpx_zz_xyzz_0, tpx_zz_xzzz_0, tpx_zz_yyyy_0, tpx_zz_yyyz_0, \
                                     tpx_zz_yyzz_0, tpx_zz_yzzz_0, tpx_zz_zzzz_0, tpx_zzz_xxx_0, tpx_zzz_xxxx_0, \
                                     tpx_zzz_xxxy_0, tpx_zzz_xxxz_0, tpx_zzz_xxy_0, tpx_zzz_xxyy_0, tpx_zzz_xxyz_0, \
                                     tpx_zzz_xxz_0, tpx_zzz_xxzz_0, tpx_zzz_xyy_0, tpx_zzz_xyyy_0, tpx_zzz_xyyz_0, \
                                     tpx_zzz_xyz_0, tpx_zzz_xyzz_0, tpx_zzz_xzz_0, tpx_zzz_xzzz_0, tpx_zzz_yyy_0, \
                                     tpx_zzz_yyyy_0, tpx_zzz_yyyz_0, tpx_zzz_yyz_0, tpx_zzz_yyzz_0, tpx_zzz_yzz_0, \
                                     tpx_zzz_yzzz_0, tpx_zzz_zzz_0, tpx_zzz_zzzz_0, tpx_zzzz_xxxx_0, tpx_zzzz_xxxy_0, \
                                     tpx_zzzz_xxxz_0, tpx_zzzz_xxyy_0, tpx_zzzz_xxyz_0, tpx_zzzz_xxzz_0, tpx_zzzz_xyyy_0, \
                                     tpx_zzzz_xyyz_0, tpx_zzzz_xyzz_0, tpx_zzzz_xzzz_0, tpx_zzzz_yyyy_0, tpx_zzzz_yyyz_0, \
                                     tpx_zzzz_yyzz_0, tpx_zzzz_yzzz_0, tpx_zzzz_zzzz_0, tpy_yzzz_zzzz_0, tpy_zz_xxxx_0, \
                                     tpy_zz_xxxy_0, tpy_zz_xxxz_0, tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxzz_0, \
                                     tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpy_zz_xyzz_0, tpy_zz_xzzz_0, tpy_zz_yyyy_0, \
                                     tpy_zz_yyyz_0, tpy_zz_yyzz_0, tpy_zz_yzzz_0, tpy_zz_zzzz_0, tpy_zzz_xxx_0, \
                                     tpy_zzz_xxxx_0, tpy_zzz_xxxy_0, tpy_zzz_xxxz_0, tpy_zzz_xxy_0, tpy_zzz_xxyy_0, \
                                     tpy_zzz_xxyz_0, tpy_zzz_xxz_0, tpy_zzz_xxzz_0, tpy_zzz_xyy_0, tpy_zzz_xyyy_0, \
                                     tpy_zzz_xyyz_0, tpy_zzz_xyz_0, tpy_zzz_xyzz_0, tpy_zzz_xzz_0, tpy_zzz_xzzz_0, \
                                     tpy_zzz_yyy_0, tpy_zzz_yyyy_0, tpy_zzz_yyyz_0, tpy_zzz_yyz_0, tpy_zzz_yyzz_0, \
                                     tpy_zzz_yzz_0, tpy_zzz_yzzz_0, tpy_zzz_zzz_0, tpy_zzz_zzzz_0, tpy_zzzz_xxxx_0, \
                                     tpy_zzzz_xxxy_0, tpy_zzzz_xxxz_0, tpy_zzzz_xxyy_0, tpy_zzzz_xxyz_0, tpy_zzzz_xxzz_0, \
                                     tpy_zzzz_xyyy_0, tpy_zzzz_xyyz_0, tpy_zzzz_xyzz_0, tpy_zzzz_xzzz_0, tpy_zzzz_yyyy_0, \
                                     tpy_zzzz_yyyz_0, tpy_zzzz_yyzz_0, tpy_zzzz_yzzz_0, tpy_zzzz_zzzz_0, tpz_yzzz_zzzz_0, \
                                     tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, tpz_zz_xxyy_0, tpz_zz_xxyz_0, \
                                     tpz_zz_xxzz_0, tpz_zz_xyyy_0, tpz_zz_xyyz_0, tpz_zz_xyzz_0, tpz_zz_xzzz_0, \
                                     tpz_zz_yyyy_0, tpz_zz_yyyz_0, tpz_zz_yyzz_0, tpz_zz_yzzz_0, tpz_zz_zzzz_0, \
                                     tpz_zzz_xxx_0, tpz_zzz_xxxx_0, tpz_zzz_xxxy_0, tpz_zzz_xxxz_0, tpz_zzz_xxy_0, \
                                     tpz_zzz_xxyy_0, tpz_zzz_xxyz_0, tpz_zzz_xxz_0, tpz_zzz_xxzz_0, tpz_zzz_xyy_0, \
                                     tpz_zzz_xyyy_0, tpz_zzz_xyyz_0, tpz_zzz_xyz_0, tpz_zzz_xyzz_0, tpz_zzz_xzz_0, \
                                     tpz_zzz_xzzz_0, tpz_zzz_yyy_0, tpz_zzz_yyyy_0, tpz_zzz_yyyz_0, tpz_zzz_yyz_0, \
                                     tpz_zzz_yyzz_0, tpz_zzz_yzz_0, tpz_zzz_yzzz_0, tpz_zzz_zzz_0, tpz_zzz_zzzz_0, \
                                     tpz_zzzz_xxxx_0, tpz_zzzz_xxxy_0, tpz_zzzz_xxxz_0, tpz_zzzz_xxyy_0, tpz_zzzz_xxyz_0, \
                                     tpz_zzzz_xxzz_0, tpz_zzzz_xyyy_0, tpz_zzzz_xyyz_0, tpz_zzzz_xyzz_0, tpz_zzzz_xzzz_0, \
                                     tpz_zzzz_yyyy_0, tpz_zzzz_yyyz_0, tpz_zzzz_yyzz_0, tpz_zzzz_yzzz_0, tpz_zzzz_zzzz_0, \
                                     ts_zzz_xxxx_0, ts_zzz_xxxy_0, ts_zzz_xxxz_0, ts_zzz_xxyy_0, ts_zzz_xxyz_0, \
                                     ts_zzz_xxzz_0, ts_zzz_xyyy_0, ts_zzz_xyyz_0, ts_zzz_xyzz_0, ts_zzz_xzzz_0, \
                                     ts_zzz_yyyy_0, ts_zzz_yyyz_0, ts_zzz_yyzz_0, ts_zzz_yzzz_0, ts_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yzzz_zzzz_0[j] = pa_y[j] * tpx_zzz_zzzz_0[j];

            tpy_yzzz_zzzz_0[j] = pa_y[j] * tpy_zzz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zzzz_0[j];

            tpz_yzzz_zzzz_0[j] = pa_y[j] * tpz_zzz_zzzz_0[j];

            tpx_zzzz_xxxx_0[j] = pa_z[j] * tpx_zzz_xxxx_0[j] + 1.5 * fl1_fx * tpx_zz_xxxx_0[j];

            tpy_zzzz_xxxx_0[j] = pa_z[j] * tpy_zzz_xxxx_0[j] + 1.5 * fl1_fx * tpy_zz_xxxx_0[j];

            tpz_zzzz_xxxx_0[j] = pa_z[j] * tpz_zzz_xxxx_0[j] + 1.5 * fl1_fx * tpz_zz_xxxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxx_0[j];

            tpx_zzzz_xxxy_0[j] = pa_z[j] * tpx_zzz_xxxy_0[j] + 1.5 * fl1_fx * tpx_zz_xxxy_0[j];

            tpy_zzzz_xxxy_0[j] = pa_z[j] * tpy_zzz_xxxy_0[j] + 1.5 * fl1_fx * tpy_zz_xxxy_0[j];

            tpz_zzzz_xxxy_0[j] = pa_z[j] * tpz_zzz_xxxy_0[j] + 1.5 * fl1_fx * tpz_zz_xxxy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxy_0[j];

            tpx_zzzz_xxxz_0[j] = pa_z[j] * tpx_zzz_xxxz_0[j] + 1.5 * fl1_fx * tpx_zz_xxxz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxx_0[j];

            tpy_zzzz_xxxz_0[j] = pa_z[j] * tpy_zzz_xxxz_0[j] + 1.5 * fl1_fx * tpy_zz_xxxz_0[j] + 0.5 * fl1_fx * tpy_zzz_xxx_0[j];

            tpz_zzzz_xxxz_0[j] =
                pa_z[j] * tpz_zzz_xxxz_0[j] + 1.5 * fl1_fx * tpz_zz_xxxz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxxz_0[j];

            tpx_zzzz_xxyy_0[j] = pa_z[j] * tpx_zzz_xxyy_0[j] + 1.5 * fl1_fx * tpx_zz_xxyy_0[j];

            tpy_zzzz_xxyy_0[j] = pa_z[j] * tpy_zzz_xxyy_0[j] + 1.5 * fl1_fx * tpy_zz_xxyy_0[j];

            tpz_zzzz_xxyy_0[j] = pa_z[j] * tpz_zzz_xxyy_0[j] + 1.5 * fl1_fx * tpz_zz_xxyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxyy_0[j];

            tpx_zzzz_xxyz_0[j] = pa_z[j] * tpx_zzz_xxyz_0[j] + 1.5 * fl1_fx * tpx_zz_xxyz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxy_0[j];

            tpy_zzzz_xxyz_0[j] = pa_z[j] * tpy_zzz_xxyz_0[j] + 1.5 * fl1_fx * tpy_zz_xxyz_0[j] + 0.5 * fl1_fx * tpy_zzz_xxy_0[j];

            tpz_zzzz_xxyz_0[j] =
                pa_z[j] * tpz_zzz_xxyz_0[j] + 1.5 * fl1_fx * tpz_zz_xxyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxyz_0[j];

            tpx_zzzz_xxzz_0[j] = pa_z[j] * tpx_zzz_xxzz_0[j] + 1.5 * fl1_fx * tpx_zz_xxzz_0[j] + fl1_fx * tpx_zzz_xxz_0[j];

            tpy_zzzz_xxzz_0[j] = pa_z[j] * tpy_zzz_xxzz_0[j] + 1.5 * fl1_fx * tpy_zz_xxzz_0[j] + fl1_fx * tpy_zzz_xxz_0[j];

            tpz_zzzz_xxzz_0[j] =
                pa_z[j] * tpz_zzz_xxzz_0[j] + 1.5 * fl1_fx * tpz_zz_xxzz_0[j] + fl1_fx * tpz_zzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxzz_0[j];

            tpx_zzzz_xyyy_0[j] = pa_z[j] * tpx_zzz_xyyy_0[j] + 1.5 * fl1_fx * tpx_zz_xyyy_0[j];

            tpy_zzzz_xyyy_0[j] = pa_z[j] * tpy_zzz_xyyy_0[j] + 1.5 * fl1_fx * tpy_zz_xyyy_0[j];

            tpz_zzzz_xyyy_0[j] = pa_z[j] * tpz_zzz_xyyy_0[j] + 1.5 * fl1_fx * tpz_zz_xyyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyyy_0[j];

            tpx_zzzz_xyyz_0[j] = pa_z[j] * tpx_zzz_xyyz_0[j] + 1.5 * fl1_fx * tpx_zz_xyyz_0[j] + 0.5 * fl1_fx * tpx_zzz_xyy_0[j];

            tpy_zzzz_xyyz_0[j] = pa_z[j] * tpy_zzz_xyyz_0[j] + 1.5 * fl1_fx * tpy_zz_xyyz_0[j] + 0.5 * fl1_fx * tpy_zzz_xyy_0[j];

            tpz_zzzz_xyyz_0[j] =
                pa_z[j] * tpz_zzz_xyyz_0[j] + 1.5 * fl1_fx * tpz_zz_xyyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyyz_0[j];

            tpx_zzzz_xyzz_0[j] = pa_z[j] * tpx_zzz_xyzz_0[j] + 1.5 * fl1_fx * tpx_zz_xyzz_0[j] + fl1_fx * tpx_zzz_xyz_0[j];

            tpy_zzzz_xyzz_0[j] = pa_z[j] * tpy_zzz_xyzz_0[j] + 1.5 * fl1_fx * tpy_zz_xyzz_0[j] + fl1_fx * tpy_zzz_xyz_0[j];

            tpz_zzzz_xyzz_0[j] =
                pa_z[j] * tpz_zzz_xyzz_0[j] + 1.5 * fl1_fx * tpz_zz_xyzz_0[j] + fl1_fx * tpz_zzz_xyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyzz_0[j];

            tpx_zzzz_xzzz_0[j] = pa_z[j] * tpx_zzz_xzzz_0[j] + 1.5 * fl1_fx * tpx_zz_xzzz_0[j] + 1.5 * fl1_fx * tpx_zzz_xzz_0[j];

            tpy_zzzz_xzzz_0[j] = pa_z[j] * tpy_zzz_xzzz_0[j] + 1.5 * fl1_fx * tpy_zz_xzzz_0[j] + 1.5 * fl1_fx * tpy_zzz_xzz_0[j];

            tpz_zzzz_xzzz_0[j] =
                pa_z[j] * tpz_zzz_xzzz_0[j] + 1.5 * fl1_fx * tpz_zz_xzzz_0[j] + 1.5 * fl1_fx * tpz_zzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xzzz_0[j];

            tpx_zzzz_yyyy_0[j] = pa_z[j] * tpx_zzz_yyyy_0[j] + 1.5 * fl1_fx * tpx_zz_yyyy_0[j];

            tpy_zzzz_yyyy_0[j] = pa_z[j] * tpy_zzz_yyyy_0[j] + 1.5 * fl1_fx * tpy_zz_yyyy_0[j];

            tpz_zzzz_yyyy_0[j] = pa_z[j] * tpz_zzz_yyyy_0[j] + 1.5 * fl1_fx * tpz_zz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyyy_0[j];

            tpx_zzzz_yyyz_0[j] = pa_z[j] * tpx_zzz_yyyz_0[j] + 1.5 * fl1_fx * tpx_zz_yyyz_0[j] + 0.5 * fl1_fx * tpx_zzz_yyy_0[j];

            tpy_zzzz_yyyz_0[j] = pa_z[j] * tpy_zzz_yyyz_0[j] + 1.5 * fl1_fx * tpy_zz_yyyz_0[j] + 0.5 * fl1_fx * tpy_zzz_yyy_0[j];

            tpz_zzzz_yyyz_0[j] =
                pa_z[j] * tpz_zzz_yyyz_0[j] + 1.5 * fl1_fx * tpz_zz_yyyz_0[j] + 0.5 * fl1_fx * tpz_zzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyyz_0[j];

            tpx_zzzz_yyzz_0[j] = pa_z[j] * tpx_zzz_yyzz_0[j] + 1.5 * fl1_fx * tpx_zz_yyzz_0[j] + fl1_fx * tpx_zzz_yyz_0[j];

            tpy_zzzz_yyzz_0[j] = pa_z[j] * tpy_zzz_yyzz_0[j] + 1.5 * fl1_fx * tpy_zz_yyzz_0[j] + fl1_fx * tpy_zzz_yyz_0[j];

            tpz_zzzz_yyzz_0[j] =
                pa_z[j] * tpz_zzz_yyzz_0[j] + 1.5 * fl1_fx * tpz_zz_yyzz_0[j] + fl1_fx * tpz_zzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyzz_0[j];

            tpx_zzzz_yzzz_0[j] = pa_z[j] * tpx_zzz_yzzz_0[j] + 1.5 * fl1_fx * tpx_zz_yzzz_0[j] + 1.5 * fl1_fx * tpx_zzz_yzz_0[j];

            tpy_zzzz_yzzz_0[j] = pa_z[j] * tpy_zzz_yzzz_0[j] + 1.5 * fl1_fx * tpy_zz_yzzz_0[j] + 1.5 * fl1_fx * tpy_zzz_yzz_0[j];

            tpz_zzzz_yzzz_0[j] =
                pa_z[j] * tpz_zzz_yzzz_0[j] + 1.5 * fl1_fx * tpz_zz_yzzz_0[j] + 1.5 * fl1_fx * tpz_zzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yzzz_0[j];

            tpx_zzzz_zzzz_0[j] = pa_z[j] * tpx_zzz_zzzz_0[j] + 1.5 * fl1_fx * tpx_zz_zzzz_0[j] + 2.0 * fl1_fx * tpx_zzz_zzz_0[j];

            tpy_zzzz_zzzz_0[j] = pa_z[j] * tpy_zzz_zzzz_0[j] + 1.5 * fl1_fx * tpy_zz_zzzz_0[j] + 2.0 * fl1_fx * tpy_zzz_zzz_0[j];

            tpz_zzzz_zzzz_0[j] =
                pa_z[j] * tpz_zzz_zzzz_0[j] + 1.5 * fl1_fx * tpz_zz_zzzz_0[j] + 2.0 * fl1_fx * tpz_zzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zzzz_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
