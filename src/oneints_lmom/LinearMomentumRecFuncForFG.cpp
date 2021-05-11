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

#include "LinearMomentumRecFuncForFG.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForFG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForFG_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_150_200(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_200_250(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_250_300(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_300_350(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_350_400(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFG_400_450(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForFG_0_50(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
{
    // Batch of Integrals (0,50)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpy_xy_xxxy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tpx_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx);

        auto tpy_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx);

        auto tpz_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx);

        auto tpx_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 1);

        auto tpy_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 1);

        auto tpz_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 1);

        auto tpx_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 2);

        auto tpy_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 2);

        auto tpz_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 2);

        auto tpx_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 3);

        auto tpy_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 3);

        auto tpz_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 3);

        auto tpx_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 4);

        auto tpy_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 4);

        auto tpz_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 4);

        auto tpx_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 5);

        auto tpy_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 5);

        auto tpz_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 5);

        auto tpx_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 6);

        auto tpy_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 6);

        auto tpz_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 6);

        auto tpx_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 7);

        auto tpy_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 7);

        auto tpz_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 7);

        auto tpx_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 8);

        auto tpy_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 8);

        auto tpz_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 8);

        auto tpx_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 9);

        auto tpy_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 9);

        auto tpz_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 9);

        auto tpx_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 10);

        auto tpy_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 10);

        auto tpz_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 10);

        auto tpx_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 11);

        auto tpy_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 11);

        auto tpz_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 11);

        auto tpx_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 12);

        auto tpy_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 12);

        auto tpz_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 12);

        auto tpx_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 13);

        auto tpy_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 13);

        auto tpz_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 13);

        auto tpx_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 14);

        auto tpy_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 14);

        auto tpz_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 14);

        auto tpx_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 15);

        auto tpy_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tpx_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 16);

        auto tpy_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tpx_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx);

        auto tpy_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx);

        auto tpz_xx_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx);

        auto tpx_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 1);

        auto tpy_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 1);

        auto tpz_xx_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 1);

        auto tpx_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 2);

        auto tpy_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 2);

        auto tpz_xx_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 2);

        auto tpx_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 3);

        auto tpy_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 3);

        auto tpz_xx_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 3);

        auto tpx_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 4);

        auto tpy_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 4);

        auto tpz_xx_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 4);

        auto tpx_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 5);

        auto tpy_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 5);

        auto tpz_xx_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 5);

        auto tpx_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 6);

        auto tpy_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 6);

        auto tpz_xx_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 6);

        auto tpx_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 7);

        auto tpy_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 7);

        auto tpz_xx_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 7);

        auto tpx_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 8);

        auto tpy_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 8);

        auto tpz_xx_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 8);

        auto tpx_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 9);

        auto tpy_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 9);

        auto tpz_xx_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 9);

        auto tpx_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 10);

        auto tpy_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 10);

        auto tpz_xy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 10);

        auto tpx_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 11);

        auto tpy_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 11);

        auto ts_xx_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx);

        auto ts_xx_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 1);

        auto ts_xx_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 2);

        auto ts_xx_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 3);

        auto ts_xx_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 4);

        auto ts_xx_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 5);

        auto ts_xx_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 6);

        auto ts_xx_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 7);

        auto ts_xx_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 8);

        auto ts_xx_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 9);

        auto ts_xx_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 10);

        auto ts_xx_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 11);

        auto ts_xx_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 12);

        auto ts_xx_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 13);

        auto ts_xx_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 14);

        auto ts_xy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 15);

        auto ts_xy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 16);

        // set up pointers to integrals

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

        auto tpy_xxy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 16);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_x_xxxx_0, tpx_x_xxxy_0, tpx_x_xxxz_0, tpx_x_xxyy_0, \
                                     tpx_x_xxyz_0, tpx_x_xxzz_0, tpx_x_xyyy_0, tpx_x_xyyz_0, tpx_x_xyzz_0, tpx_x_xzzz_0, \
                                     tpx_x_yyyy_0, tpx_x_yyyz_0, tpx_x_yyzz_0, tpx_x_yzzz_0, tpx_x_zzzz_0, tpx_xx_xxx_0, \
                                     tpx_xx_xxxx_0, tpx_xx_xxxy_0, tpx_xx_xxxz_0, tpx_xx_xxy_0, tpx_xx_xxyy_0, \
                                     tpx_xx_xxyz_0, tpx_xx_xxz_0, tpx_xx_xxzz_0, tpx_xx_xyy_0, tpx_xx_xyyy_0, \
                                     tpx_xx_xyyz_0, tpx_xx_xyz_0, tpx_xx_xyzz_0, tpx_xx_xzz_0, tpx_xx_xzzz_0, \
                                     tpx_xx_yyy_0, tpx_xx_yyyy_0, tpx_xx_yyyz_0, tpx_xx_yyz_0, tpx_xx_yyzz_0, \
                                     tpx_xx_yzz_0, tpx_xx_yzzz_0, tpx_xx_zzz_0, tpx_xx_zzzz_0, tpx_xxx_xxxx_0, \
                                     tpx_xxx_xxxy_0, tpx_xxx_xxxz_0, tpx_xxx_xxyy_0, tpx_xxx_xxyz_0, tpx_xxx_xxzz_0, \
                                     tpx_xxx_xyyy_0, tpx_xxx_xyyz_0, tpx_xxx_xyzz_0, tpx_xxx_xzzz_0, tpx_xxx_yyyy_0, \
                                     tpx_xxx_yyyz_0, tpx_xxx_yyzz_0, tpx_xxx_yzzz_0, tpx_xxx_zzzz_0, tpx_xxy_xxxx_0, \
                                     tpx_xxy_xxxy_0, tpx_xy_xxx_0, tpx_xy_xxxx_0, tpx_xy_xxxy_0, tpx_xy_xxy_0, \
                                     tpx_y_xxxx_0, tpx_y_xxxy_0, tpy_x_xxxx_0, tpy_x_xxxy_0, tpy_x_xxxz_0, tpy_x_xxyy_0, \
                                     tpy_x_xxyz_0, tpy_x_xxzz_0, tpy_x_xyyy_0, tpy_x_xyyz_0, tpy_x_xyzz_0, tpy_x_xzzz_0, \
                                     tpy_x_yyyy_0, tpy_x_yyyz_0, tpy_x_yyzz_0, tpy_x_yzzz_0, tpy_x_zzzz_0, tpy_xx_xxx_0, \
                                     tpy_xx_xxxx_0, tpy_xx_xxxy_0, tpy_xx_xxxz_0, tpy_xx_xxy_0, tpy_xx_xxyy_0, \
                                     tpy_xx_xxyz_0, tpy_xx_xxz_0, tpy_xx_xxzz_0, tpy_xx_xyy_0, tpy_xx_xyyy_0, \
                                     tpy_xx_xyyz_0, tpy_xx_xyz_0, tpy_xx_xyzz_0, tpy_xx_xzz_0, tpy_xx_xzzz_0, \
                                     tpy_xx_yyy_0, tpy_xx_yyyy_0, tpy_xx_yyyz_0, tpy_xx_yyz_0, tpy_xx_yyzz_0, \
                                     tpy_xx_yzz_0, tpy_xx_yzzz_0, tpy_xx_zzz_0, tpy_xx_zzzz_0, tpy_xxx_xxxx_0, \
                                     tpy_xxx_xxxy_0, tpy_xxx_xxxz_0, tpy_xxx_xxyy_0, tpy_xxx_xxyz_0, tpy_xxx_xxzz_0, \
                                     tpy_xxx_xyyy_0, tpy_xxx_xyyz_0, tpy_xxx_xyzz_0, tpy_xxx_xzzz_0, tpy_xxx_yyyy_0, \
                                     tpy_xxx_yyyz_0, tpy_xxx_yyzz_0, tpy_xxx_yzzz_0, tpy_xxx_zzzz_0, tpy_xxy_xxxx_0, \
                                     tpy_xxy_xxxy_0, tpy_xy_xxx_0, tpy_xy_xxxx_0, tpy_xy_xxxy_0, tpy_xy_xxy_0, \
                                     tpy_y_xxxx_0, tpy_y_xxxy_0, tpz_x_xxxx_0, tpz_x_xxxy_0, tpz_x_xxxz_0, tpz_x_xxyy_0, \
                                     tpz_x_xxyz_0, tpz_x_xxzz_0, tpz_x_xyyy_0, tpz_x_xyyz_0, tpz_x_xyzz_0, tpz_x_xzzz_0, \
                                     tpz_x_yyyy_0, tpz_x_yyyz_0, tpz_x_yyzz_0, tpz_x_yzzz_0, tpz_x_zzzz_0, tpz_xx_xxx_0, \
                                     tpz_xx_xxxx_0, tpz_xx_xxxy_0, tpz_xx_xxxz_0, tpz_xx_xxy_0, tpz_xx_xxyy_0, \
                                     tpz_xx_xxyz_0, tpz_xx_xxz_0, tpz_xx_xxzz_0, tpz_xx_xyy_0, tpz_xx_xyyy_0, \
                                     tpz_xx_xyyz_0, tpz_xx_xyz_0, tpz_xx_xyzz_0, tpz_xx_xzz_0, tpz_xx_xzzz_0, \
                                     tpz_xx_yyy_0, tpz_xx_yyyy_0, tpz_xx_yyyz_0, tpz_xx_yyz_0, tpz_xx_yyzz_0, \
                                     tpz_xx_yzz_0, tpz_xx_yzzz_0, tpz_xx_zzz_0, tpz_xx_zzzz_0, tpz_xxx_xxxx_0, \
                                     tpz_xxx_xxxy_0, tpz_xxx_xxxz_0, tpz_xxx_xxyy_0, tpz_xxx_xxyz_0, tpz_xxx_xxzz_0, \
                                     tpz_xxx_xyyy_0, tpz_xxx_xyyz_0, tpz_xxx_xyzz_0, tpz_xxx_xzzz_0, tpz_xxx_yyyy_0, \
                                     tpz_xxx_yyyz_0, tpz_xxx_yyzz_0, tpz_xxx_yzzz_0, tpz_xxx_zzzz_0, tpz_xxy_xxxx_0, \
                                     tpz_xy_xxx_0, tpz_xy_xxxx_0, tpz_y_xxxx_0, ts_xx_xxxx_0, ts_xx_xxxy_0, \
                                     ts_xx_xxxz_0, ts_xx_xxyy_0, ts_xx_xxyz_0, ts_xx_xxzz_0, ts_xx_xyyy_0, ts_xx_xyyz_0, \
                                     ts_xx_xyzz_0, ts_xx_xzzz_0, ts_xx_yyyy_0, ts_xx_yyyz_0, ts_xx_yyzz_0, ts_xx_yzzz_0, \
                                     ts_xx_zzzz_0, ts_xy_xxxx_0, ts_xy_xxxy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxx_xxxx_0[j] =
                pa_x[j] * tpx_xx_xxxx_0[j] + fl1_fx * tpx_x_xxxx_0[j] + 2.0 * fl1_fx * tpx_xx_xxx_0[j] - fl1_fgb * fl1_fx * ts_xx_xxxx_0[j];

            tpy_xxx_xxxx_0[j] = pa_x[j] * tpy_xx_xxxx_0[j] + fl1_fx * tpy_x_xxxx_0[j] + 2.0 * fl1_fx * tpy_xx_xxx_0[j];

            tpz_xxx_xxxx_0[j] = pa_x[j] * tpz_xx_xxxx_0[j] + fl1_fx * tpz_x_xxxx_0[j] + 2.0 * fl1_fx * tpz_xx_xxx_0[j];

            tpx_xxx_xxxy_0[j] =
                pa_x[j] * tpx_xx_xxxy_0[j] + fl1_fx * tpx_x_xxxy_0[j] + 1.5 * fl1_fx * tpx_xx_xxy_0[j] - fl1_fgb * fl1_fx * ts_xx_xxxy_0[j];

            tpy_xxx_xxxy_0[j] = pa_x[j] * tpy_xx_xxxy_0[j] + fl1_fx * tpy_x_xxxy_0[j] + 1.5 * fl1_fx * tpy_xx_xxy_0[j];

            tpz_xxx_xxxy_0[j] = pa_x[j] * tpz_xx_xxxy_0[j] + fl1_fx * tpz_x_xxxy_0[j] + 1.5 * fl1_fx * tpz_xx_xxy_0[j];

            tpx_xxx_xxxz_0[j] =
                pa_x[j] * tpx_xx_xxxz_0[j] + fl1_fx * tpx_x_xxxz_0[j] + 1.5 * fl1_fx * tpx_xx_xxz_0[j] - fl1_fgb * fl1_fx * ts_xx_xxxz_0[j];

            tpy_xxx_xxxz_0[j] = pa_x[j] * tpy_xx_xxxz_0[j] + fl1_fx * tpy_x_xxxz_0[j] + 1.5 * fl1_fx * tpy_xx_xxz_0[j];

            tpz_xxx_xxxz_0[j] = pa_x[j] * tpz_xx_xxxz_0[j] + fl1_fx * tpz_x_xxxz_0[j] + 1.5 * fl1_fx * tpz_xx_xxz_0[j];

            tpx_xxx_xxyy_0[j] = pa_x[j] * tpx_xx_xxyy_0[j] + fl1_fx * tpx_x_xxyy_0[j] + fl1_fx * tpx_xx_xyy_0[j] - fl1_fgb * fl1_fx * ts_xx_xxyy_0[j];

            tpy_xxx_xxyy_0[j] = pa_x[j] * tpy_xx_xxyy_0[j] + fl1_fx * tpy_x_xxyy_0[j] + fl1_fx * tpy_xx_xyy_0[j];

            tpz_xxx_xxyy_0[j] = pa_x[j] * tpz_xx_xxyy_0[j] + fl1_fx * tpz_x_xxyy_0[j] + fl1_fx * tpz_xx_xyy_0[j];

            tpx_xxx_xxyz_0[j] = pa_x[j] * tpx_xx_xxyz_0[j] + fl1_fx * tpx_x_xxyz_0[j] + fl1_fx * tpx_xx_xyz_0[j] - fl1_fgb * fl1_fx * ts_xx_xxyz_0[j];

            tpy_xxx_xxyz_0[j] = pa_x[j] * tpy_xx_xxyz_0[j] + fl1_fx * tpy_x_xxyz_0[j] + fl1_fx * tpy_xx_xyz_0[j];

            tpz_xxx_xxyz_0[j] = pa_x[j] * tpz_xx_xxyz_0[j] + fl1_fx * tpz_x_xxyz_0[j] + fl1_fx * tpz_xx_xyz_0[j];

            tpx_xxx_xxzz_0[j] = pa_x[j] * tpx_xx_xxzz_0[j] + fl1_fx * tpx_x_xxzz_0[j] + fl1_fx * tpx_xx_xzz_0[j] - fl1_fgb * fl1_fx * ts_xx_xxzz_0[j];

            tpy_xxx_xxzz_0[j] = pa_x[j] * tpy_xx_xxzz_0[j] + fl1_fx * tpy_x_xxzz_0[j] + fl1_fx * tpy_xx_xzz_0[j];

            tpz_xxx_xxzz_0[j] = pa_x[j] * tpz_xx_xxzz_0[j] + fl1_fx * tpz_x_xxzz_0[j] + fl1_fx * tpz_xx_xzz_0[j];

            tpx_xxx_xyyy_0[j] =
                pa_x[j] * tpx_xx_xyyy_0[j] + fl1_fx * tpx_x_xyyy_0[j] + 0.5 * fl1_fx * tpx_xx_yyy_0[j] - fl1_fgb * fl1_fx * ts_xx_xyyy_0[j];

            tpy_xxx_xyyy_0[j] = pa_x[j] * tpy_xx_xyyy_0[j] + fl1_fx * tpy_x_xyyy_0[j] + 0.5 * fl1_fx * tpy_xx_yyy_0[j];

            tpz_xxx_xyyy_0[j] = pa_x[j] * tpz_xx_xyyy_0[j] + fl1_fx * tpz_x_xyyy_0[j] + 0.5 * fl1_fx * tpz_xx_yyy_0[j];

            tpx_xxx_xyyz_0[j] =
                pa_x[j] * tpx_xx_xyyz_0[j] + fl1_fx * tpx_x_xyyz_0[j] + 0.5 * fl1_fx * tpx_xx_yyz_0[j] - fl1_fgb * fl1_fx * ts_xx_xyyz_0[j];

            tpy_xxx_xyyz_0[j] = pa_x[j] * tpy_xx_xyyz_0[j] + fl1_fx * tpy_x_xyyz_0[j] + 0.5 * fl1_fx * tpy_xx_yyz_0[j];

            tpz_xxx_xyyz_0[j] = pa_x[j] * tpz_xx_xyyz_0[j] + fl1_fx * tpz_x_xyyz_0[j] + 0.5 * fl1_fx * tpz_xx_yyz_0[j];

            tpx_xxx_xyzz_0[j] =
                pa_x[j] * tpx_xx_xyzz_0[j] + fl1_fx * tpx_x_xyzz_0[j] + 0.5 * fl1_fx * tpx_xx_yzz_0[j] - fl1_fgb * fl1_fx * ts_xx_xyzz_0[j];

            tpy_xxx_xyzz_0[j] = pa_x[j] * tpy_xx_xyzz_0[j] + fl1_fx * tpy_x_xyzz_0[j] + 0.5 * fl1_fx * tpy_xx_yzz_0[j];

            tpz_xxx_xyzz_0[j] = pa_x[j] * tpz_xx_xyzz_0[j] + fl1_fx * tpz_x_xyzz_0[j] + 0.5 * fl1_fx * tpz_xx_yzz_0[j];

            tpx_xxx_xzzz_0[j] =
                pa_x[j] * tpx_xx_xzzz_0[j] + fl1_fx * tpx_x_xzzz_0[j] + 0.5 * fl1_fx * tpx_xx_zzz_0[j] - fl1_fgb * fl1_fx * ts_xx_xzzz_0[j];

            tpy_xxx_xzzz_0[j] = pa_x[j] * tpy_xx_xzzz_0[j] + fl1_fx * tpy_x_xzzz_0[j] + 0.5 * fl1_fx * tpy_xx_zzz_0[j];

            tpz_xxx_xzzz_0[j] = pa_x[j] * tpz_xx_xzzz_0[j] + fl1_fx * tpz_x_xzzz_0[j] + 0.5 * fl1_fx * tpz_xx_zzz_0[j];

            tpx_xxx_yyyy_0[j] = pa_x[j] * tpx_xx_yyyy_0[j] + fl1_fx * tpx_x_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xx_yyyy_0[j];

            tpy_xxx_yyyy_0[j] = pa_x[j] * tpy_xx_yyyy_0[j] + fl1_fx * tpy_x_yyyy_0[j];

            tpz_xxx_yyyy_0[j] = pa_x[j] * tpz_xx_yyyy_0[j] + fl1_fx * tpz_x_yyyy_0[j];

            tpx_xxx_yyyz_0[j] = pa_x[j] * tpx_xx_yyyz_0[j] + fl1_fx * tpx_x_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xx_yyyz_0[j];

            tpy_xxx_yyyz_0[j] = pa_x[j] * tpy_xx_yyyz_0[j] + fl1_fx * tpy_x_yyyz_0[j];

            tpz_xxx_yyyz_0[j] = pa_x[j] * tpz_xx_yyyz_0[j] + fl1_fx * tpz_x_yyyz_0[j];

            tpx_xxx_yyzz_0[j] = pa_x[j] * tpx_xx_yyzz_0[j] + fl1_fx * tpx_x_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xx_yyzz_0[j];

            tpy_xxx_yyzz_0[j] = pa_x[j] * tpy_xx_yyzz_0[j] + fl1_fx * tpy_x_yyzz_0[j];

            tpz_xxx_yyzz_0[j] = pa_x[j] * tpz_xx_yyzz_0[j] + fl1_fx * tpz_x_yyzz_0[j];

            tpx_xxx_yzzz_0[j] = pa_x[j] * tpx_xx_yzzz_0[j] + fl1_fx * tpx_x_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xx_yzzz_0[j];

            tpy_xxx_yzzz_0[j] = pa_x[j] * tpy_xx_yzzz_0[j] + fl1_fx * tpy_x_yzzz_0[j];

            tpz_xxx_yzzz_0[j] = pa_x[j] * tpz_xx_yzzz_0[j] + fl1_fx * tpz_x_yzzz_0[j];

            tpx_xxx_zzzz_0[j] = pa_x[j] * tpx_xx_zzzz_0[j] + fl1_fx * tpx_x_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xx_zzzz_0[j];

            tpy_xxx_zzzz_0[j] = pa_x[j] * tpy_xx_zzzz_0[j] + fl1_fx * tpy_x_zzzz_0[j];

            tpz_xxx_zzzz_0[j] = pa_x[j] * tpz_xx_zzzz_0[j] + fl1_fx * tpz_x_zzzz_0[j];

            tpx_xxy_xxxx_0[j] =
                pa_x[j] * tpx_xy_xxxx_0[j] + 0.5 * fl1_fx * tpx_y_xxxx_0[j] + 2.0 * fl1_fx * tpx_xy_xxx_0[j] - fl1_fgb * fl1_fx * ts_xy_xxxx_0[j];

            tpy_xxy_xxxx_0[j] = pa_x[j] * tpy_xy_xxxx_0[j] + 0.5 * fl1_fx * tpy_y_xxxx_0[j] + 2.0 * fl1_fx * tpy_xy_xxx_0[j];

            tpz_xxy_xxxx_0[j] = pa_x[j] * tpz_xy_xxxx_0[j] + 0.5 * fl1_fx * tpz_y_xxxx_0[j] + 2.0 * fl1_fx * tpz_xy_xxx_0[j];

            tpx_xxy_xxxy_0[j] =
                pa_x[j] * tpx_xy_xxxy_0[j] + 0.5 * fl1_fx * tpx_y_xxxy_0[j] + 1.5 * fl1_fx * tpx_xy_xxy_0[j] - fl1_fgb * fl1_fx * ts_xy_xxxy_0[j];

            tpy_xxy_xxxy_0[j] = pa_x[j] * tpy_xy_xxxy_0[j] + 0.5 * fl1_fx * tpy_y_xxxy_0[j] + 1.5 * fl1_fx * tpy_xy_xxy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_50_100(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (50,100)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpz_xz_xxxz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tpx_xz_xxyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 33);

        auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tpx_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 17);

        auto tpy_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tpx_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 18);

        auto tpy_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tpx_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 19);

        auto tpy_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tpx_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 20);

        auto tpy_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tpx_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 21);

        auto tpy_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tpx_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 22);

        auto tpy_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tpx_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 23);

        auto tpy_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tpx_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 24);

        auto tpy_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tpx_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 25);

        auto tpy_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tpx_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 26);

        auto tpy_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tpx_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 27);

        auto tpy_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tpx_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 28);

        auto tpy_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tpx_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 29);

        auto tpy_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29);

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpz_xy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 11);

        auto tpx_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 12);

        auto tpy_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 12);

        auto tpz_xy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 12);

        auto tpx_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 13);

        auto tpy_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 13);

        auto tpz_xy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 13);

        auto tpx_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 14);

        auto tpy_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 14);

        auto tpz_xy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 14);

        auto tpx_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 15);

        auto tpy_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 15);

        auto tpz_xy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 15);

        auto tpx_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 16);

        auto tpy_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 16);

        auto tpz_xy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 16);

        auto tpx_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 17);

        auto tpy_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 17);

        auto tpz_xy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 17);

        auto tpx_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 18);

        auto tpy_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 18);

        auto tpz_xy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 18);

        auto tpx_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 19);

        auto tpy_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 19);

        auto tpz_xy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 19);

        auto tpx_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 20);

        auto tpy_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 20);

        auto tpz_xz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 20);

        auto tpx_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 21);

        auto tpy_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 21);

        auto tpz_xz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 21);

        auto tpx_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 22);

        auto tpy_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 22);

        auto tpz_xz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 22);

        auto tpx_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 23);

        auto ts_xy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 17);

        auto ts_xy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 18);

        auto ts_xy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 19);

        auto ts_xy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 20);

        auto ts_xy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 21);

        auto ts_xy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 22);

        auto ts_xy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 23);

        auto ts_xy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 24);

        auto ts_xy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 25);

        auto ts_xy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 26);

        auto ts_xy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 27);

        auto ts_xy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 28);

        auto ts_xy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 29);

        auto ts_xz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 30);

        auto ts_xz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 31);

        auto ts_xz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 32);

        auto ts_xz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 33);

        // set up pointers to integrals

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

        auto tpz_xxz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tpx_xxz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 33);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxy_xxxz_0, tpx_xxy_xxyy_0, tpx_xxy_xxyz_0, \
                                     tpx_xxy_xxzz_0, tpx_xxy_xyyy_0, tpx_xxy_xyyz_0, tpx_xxy_xyzz_0, tpx_xxy_xzzz_0, \
                                     tpx_xxy_yyyy_0, tpx_xxy_yyyz_0, tpx_xxy_yyzz_0, tpx_xxy_yzzz_0, tpx_xxy_zzzz_0, \
                                     tpx_xxz_xxxx_0, tpx_xxz_xxxy_0, tpx_xxz_xxxz_0, tpx_xxz_xxyy_0, tpx_xy_xxxz_0, \
                                     tpx_xy_xxyy_0, tpx_xy_xxyz_0, tpx_xy_xxz_0, tpx_xy_xxzz_0, tpx_xy_xyy_0, \
                                     tpx_xy_xyyy_0, tpx_xy_xyyz_0, tpx_xy_xyz_0, tpx_xy_xyzz_0, tpx_xy_xzz_0, \
                                     tpx_xy_xzzz_0, tpx_xy_yyy_0, tpx_xy_yyyy_0, tpx_xy_yyyz_0, tpx_xy_yyz_0, \
                                     tpx_xy_yyzz_0, tpx_xy_yzz_0, tpx_xy_yzzz_0, tpx_xy_zzz_0, tpx_xy_zzzz_0, \
                                     tpx_xz_xxx_0, tpx_xz_xxxx_0, tpx_xz_xxxy_0, tpx_xz_xxxz_0, tpx_xz_xxy_0, \
                                     tpx_xz_xxyy_0, tpx_xz_xxz_0, tpx_xz_xyy_0, tpx_y_xxxz_0, tpx_y_xxyy_0, tpx_y_xxyz_0, \
                                     tpx_y_xxzz_0, tpx_y_xyyy_0, tpx_y_xyyz_0, tpx_y_xyzz_0, tpx_y_xzzz_0, tpx_y_yyyy_0, \
                                     tpx_y_yyyz_0, tpx_y_yyzz_0, tpx_y_yzzz_0, tpx_y_zzzz_0, tpx_z_xxxx_0, tpx_z_xxxy_0, \
                                     tpx_z_xxxz_0, tpx_z_xxyy_0, tpy_xxy_xxxz_0, tpy_xxy_xxyy_0, tpy_xxy_xxyz_0, \
                                     tpy_xxy_xxzz_0, tpy_xxy_xyyy_0, tpy_xxy_xyyz_0, tpy_xxy_xyzz_0, tpy_xxy_xzzz_0, \
                                     tpy_xxy_yyyy_0, tpy_xxy_yyyz_0, tpy_xxy_yyzz_0, tpy_xxy_yzzz_0, tpy_xxy_zzzz_0, \
                                     tpy_xxz_xxxx_0, tpy_xxz_xxxy_0, tpy_xxz_xxxz_0, tpy_xy_xxxz_0, tpy_xy_xxyy_0, \
                                     tpy_xy_xxyz_0, tpy_xy_xxz_0, tpy_xy_xxzz_0, tpy_xy_xyy_0, tpy_xy_xyyy_0, \
                                     tpy_xy_xyyz_0, tpy_xy_xyz_0, tpy_xy_xyzz_0, tpy_xy_xzz_0, tpy_xy_xzzz_0, \
                                     tpy_xy_yyy_0, tpy_xy_yyyy_0, tpy_xy_yyyz_0, tpy_xy_yyz_0, tpy_xy_yyzz_0, \
                                     tpy_xy_yzz_0, tpy_xy_yzzz_0, tpy_xy_zzz_0, tpy_xy_zzzz_0, tpy_xz_xxx_0, \
                                     tpy_xz_xxxx_0, tpy_xz_xxxy_0, tpy_xz_xxxz_0, tpy_xz_xxy_0, tpy_xz_xxz_0, \
                                     tpy_y_xxxz_0, tpy_y_xxyy_0, tpy_y_xxyz_0, tpy_y_xxzz_0, tpy_y_xyyy_0, tpy_y_xyyz_0, \
                                     tpy_y_xyzz_0, tpy_y_xzzz_0, tpy_y_yyyy_0, tpy_y_yyyz_0, tpy_y_yyzz_0, tpy_y_yzzz_0, \
                                     tpy_y_zzzz_0, tpy_z_xxxx_0, tpy_z_xxxy_0, tpy_z_xxxz_0, tpz_xxy_xxxy_0, \
                                     tpz_xxy_xxxz_0, tpz_xxy_xxyy_0, tpz_xxy_xxyz_0, tpz_xxy_xxzz_0, tpz_xxy_xyyy_0, \
                                     tpz_xxy_xyyz_0, tpz_xxy_xyzz_0, tpz_xxy_xzzz_0, tpz_xxy_yyyy_0, tpz_xxy_yyyz_0, \
                                     tpz_xxy_yyzz_0, tpz_xxy_yzzz_0, tpz_xxy_zzzz_0, tpz_xxz_xxxx_0, tpz_xxz_xxxy_0, \
                                     tpz_xxz_xxxz_0, tpz_xy_xxxy_0, tpz_xy_xxxz_0, tpz_xy_xxy_0, tpz_xy_xxyy_0, \
                                     tpz_xy_xxyz_0, tpz_xy_xxz_0, tpz_xy_xxzz_0, tpz_xy_xyy_0, tpz_xy_xyyy_0, \
                                     tpz_xy_xyyz_0, tpz_xy_xyz_0, tpz_xy_xyzz_0, tpz_xy_xzz_0, tpz_xy_xzzz_0, \
                                     tpz_xy_yyy_0, tpz_xy_yyyy_0, tpz_xy_yyyz_0, tpz_xy_yyz_0, tpz_xy_yyzz_0, \
                                     tpz_xy_yzz_0, tpz_xy_yzzz_0, tpz_xy_zzz_0, tpz_xy_zzzz_0, tpz_xz_xxx_0, \
                                     tpz_xz_xxxx_0, tpz_xz_xxxy_0, tpz_xz_xxxz_0, tpz_xz_xxy_0, tpz_xz_xxz_0, \
                                     tpz_y_xxxy_0, tpz_y_xxxz_0, tpz_y_xxyy_0, tpz_y_xxyz_0, tpz_y_xxzz_0, tpz_y_xyyy_0, \
                                     tpz_y_xyyz_0, tpz_y_xyzz_0, tpz_y_xzzz_0, tpz_y_yyyy_0, tpz_y_yyyz_0, tpz_y_yyzz_0, \
                                     tpz_y_yzzz_0, tpz_y_zzzz_0, tpz_z_xxxx_0, tpz_z_xxxy_0, tpz_z_xxxz_0, ts_xy_xxxz_0, \
                                     ts_xy_xxyy_0, ts_xy_xxyz_0, ts_xy_xxzz_0, ts_xy_xyyy_0, ts_xy_xyyz_0, ts_xy_xyzz_0, \
                                     ts_xy_xzzz_0, ts_xy_yyyy_0, ts_xy_yyyz_0, ts_xy_yyzz_0, ts_xy_yzzz_0, ts_xy_zzzz_0, \
                                     ts_xz_xxxx_0, ts_xz_xxxy_0, ts_xz_xxxz_0, ts_xz_xxyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_xxy_xxxy_0[j] = pa_x[j] * tpz_xy_xxxy_0[j] + 0.5 * fl1_fx * tpz_y_xxxy_0[j] + 1.5 * fl1_fx * tpz_xy_xxy_0[j];

            tpx_xxy_xxxz_0[j] =
                pa_x[j] * tpx_xy_xxxz_0[j] + 0.5 * fl1_fx * tpx_y_xxxz_0[j] + 1.5 * fl1_fx * tpx_xy_xxz_0[j] - fl1_fgb * fl1_fx * ts_xy_xxxz_0[j];

            tpy_xxy_xxxz_0[j] = pa_x[j] * tpy_xy_xxxz_0[j] + 0.5 * fl1_fx * tpy_y_xxxz_0[j] + 1.5 * fl1_fx * tpy_xy_xxz_0[j];

            tpz_xxy_xxxz_0[j] = pa_x[j] * tpz_xy_xxxz_0[j] + 0.5 * fl1_fx * tpz_y_xxxz_0[j] + 1.5 * fl1_fx * tpz_xy_xxz_0[j];

            tpx_xxy_xxyy_0[j] =
                pa_x[j] * tpx_xy_xxyy_0[j] + 0.5 * fl1_fx * tpx_y_xxyy_0[j] + fl1_fx * tpx_xy_xyy_0[j] - fl1_fgb * fl1_fx * ts_xy_xxyy_0[j];

            tpy_xxy_xxyy_0[j] = pa_x[j] * tpy_xy_xxyy_0[j] + 0.5 * fl1_fx * tpy_y_xxyy_0[j] + fl1_fx * tpy_xy_xyy_0[j];

            tpz_xxy_xxyy_0[j] = pa_x[j] * tpz_xy_xxyy_0[j] + 0.5 * fl1_fx * tpz_y_xxyy_0[j] + fl1_fx * tpz_xy_xyy_0[j];

            tpx_xxy_xxyz_0[j] =
                pa_x[j] * tpx_xy_xxyz_0[j] + 0.5 * fl1_fx * tpx_y_xxyz_0[j] + fl1_fx * tpx_xy_xyz_0[j] - fl1_fgb * fl1_fx * ts_xy_xxyz_0[j];

            tpy_xxy_xxyz_0[j] = pa_x[j] * tpy_xy_xxyz_0[j] + 0.5 * fl1_fx * tpy_y_xxyz_0[j] + fl1_fx * tpy_xy_xyz_0[j];

            tpz_xxy_xxyz_0[j] = pa_x[j] * tpz_xy_xxyz_0[j] + 0.5 * fl1_fx * tpz_y_xxyz_0[j] + fl1_fx * tpz_xy_xyz_0[j];

            tpx_xxy_xxzz_0[j] =
                pa_x[j] * tpx_xy_xxzz_0[j] + 0.5 * fl1_fx * tpx_y_xxzz_0[j] + fl1_fx * tpx_xy_xzz_0[j] - fl1_fgb * fl1_fx * ts_xy_xxzz_0[j];

            tpy_xxy_xxzz_0[j] = pa_x[j] * tpy_xy_xxzz_0[j] + 0.5 * fl1_fx * tpy_y_xxzz_0[j] + fl1_fx * tpy_xy_xzz_0[j];

            tpz_xxy_xxzz_0[j] = pa_x[j] * tpz_xy_xxzz_0[j] + 0.5 * fl1_fx * tpz_y_xxzz_0[j] + fl1_fx * tpz_xy_xzz_0[j];

            tpx_xxy_xyyy_0[j] =
                pa_x[j] * tpx_xy_xyyy_0[j] + 0.5 * fl1_fx * tpx_y_xyyy_0[j] + 0.5 * fl1_fx * tpx_xy_yyy_0[j] - fl1_fgb * fl1_fx * ts_xy_xyyy_0[j];

            tpy_xxy_xyyy_0[j] = pa_x[j] * tpy_xy_xyyy_0[j] + 0.5 * fl1_fx * tpy_y_xyyy_0[j] + 0.5 * fl1_fx * tpy_xy_yyy_0[j];

            tpz_xxy_xyyy_0[j] = pa_x[j] * tpz_xy_xyyy_0[j] + 0.5 * fl1_fx * tpz_y_xyyy_0[j] + 0.5 * fl1_fx * tpz_xy_yyy_0[j];

            tpx_xxy_xyyz_0[j] =
                pa_x[j] * tpx_xy_xyyz_0[j] + 0.5 * fl1_fx * tpx_y_xyyz_0[j] + 0.5 * fl1_fx * tpx_xy_yyz_0[j] - fl1_fgb * fl1_fx * ts_xy_xyyz_0[j];

            tpy_xxy_xyyz_0[j] = pa_x[j] * tpy_xy_xyyz_0[j] + 0.5 * fl1_fx * tpy_y_xyyz_0[j] + 0.5 * fl1_fx * tpy_xy_yyz_0[j];

            tpz_xxy_xyyz_0[j] = pa_x[j] * tpz_xy_xyyz_0[j] + 0.5 * fl1_fx * tpz_y_xyyz_0[j] + 0.5 * fl1_fx * tpz_xy_yyz_0[j];

            tpx_xxy_xyzz_0[j] =
                pa_x[j] * tpx_xy_xyzz_0[j] + 0.5 * fl1_fx * tpx_y_xyzz_0[j] + 0.5 * fl1_fx * tpx_xy_yzz_0[j] - fl1_fgb * fl1_fx * ts_xy_xyzz_0[j];

            tpy_xxy_xyzz_0[j] = pa_x[j] * tpy_xy_xyzz_0[j] + 0.5 * fl1_fx * tpy_y_xyzz_0[j] + 0.5 * fl1_fx * tpy_xy_yzz_0[j];

            tpz_xxy_xyzz_0[j] = pa_x[j] * tpz_xy_xyzz_0[j] + 0.5 * fl1_fx * tpz_y_xyzz_0[j] + 0.5 * fl1_fx * tpz_xy_yzz_0[j];

            tpx_xxy_xzzz_0[j] =
                pa_x[j] * tpx_xy_xzzz_0[j] + 0.5 * fl1_fx * tpx_y_xzzz_0[j] + 0.5 * fl1_fx * tpx_xy_zzz_0[j] - fl1_fgb * fl1_fx * ts_xy_xzzz_0[j];

            tpy_xxy_xzzz_0[j] = pa_x[j] * tpy_xy_xzzz_0[j] + 0.5 * fl1_fx * tpy_y_xzzz_0[j] + 0.5 * fl1_fx * tpy_xy_zzz_0[j];

            tpz_xxy_xzzz_0[j] = pa_x[j] * tpz_xy_xzzz_0[j] + 0.5 * fl1_fx * tpz_y_xzzz_0[j] + 0.5 * fl1_fx * tpz_xy_zzz_0[j];

            tpx_xxy_yyyy_0[j] = pa_x[j] * tpx_xy_yyyy_0[j] + 0.5 * fl1_fx * tpx_y_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xy_yyyy_0[j];

            tpy_xxy_yyyy_0[j] = pa_x[j] * tpy_xy_yyyy_0[j] + 0.5 * fl1_fx * tpy_y_yyyy_0[j];

            tpz_xxy_yyyy_0[j] = pa_x[j] * tpz_xy_yyyy_0[j] + 0.5 * fl1_fx * tpz_y_yyyy_0[j];

            tpx_xxy_yyyz_0[j] = pa_x[j] * tpx_xy_yyyz_0[j] + 0.5 * fl1_fx * tpx_y_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xy_yyyz_0[j];

            tpy_xxy_yyyz_0[j] = pa_x[j] * tpy_xy_yyyz_0[j] + 0.5 * fl1_fx * tpy_y_yyyz_0[j];

            tpz_xxy_yyyz_0[j] = pa_x[j] * tpz_xy_yyyz_0[j] + 0.5 * fl1_fx * tpz_y_yyyz_0[j];

            tpx_xxy_yyzz_0[j] = pa_x[j] * tpx_xy_yyzz_0[j] + 0.5 * fl1_fx * tpx_y_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xy_yyzz_0[j];

            tpy_xxy_yyzz_0[j] = pa_x[j] * tpy_xy_yyzz_0[j] + 0.5 * fl1_fx * tpy_y_yyzz_0[j];

            tpz_xxy_yyzz_0[j] = pa_x[j] * tpz_xy_yyzz_0[j] + 0.5 * fl1_fx * tpz_y_yyzz_0[j];

            tpx_xxy_yzzz_0[j] = pa_x[j] * tpx_xy_yzzz_0[j] + 0.5 * fl1_fx * tpx_y_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xy_yzzz_0[j];

            tpy_xxy_yzzz_0[j] = pa_x[j] * tpy_xy_yzzz_0[j] + 0.5 * fl1_fx * tpy_y_yzzz_0[j];

            tpz_xxy_yzzz_0[j] = pa_x[j] * tpz_xy_yzzz_0[j] + 0.5 * fl1_fx * tpz_y_yzzz_0[j];

            tpx_xxy_zzzz_0[j] = pa_x[j] * tpx_xy_zzzz_0[j] + 0.5 * fl1_fx * tpx_y_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xy_zzzz_0[j];

            tpy_xxy_zzzz_0[j] = pa_x[j] * tpy_xy_zzzz_0[j] + 0.5 * fl1_fx * tpy_y_zzzz_0[j];

            tpz_xxy_zzzz_0[j] = pa_x[j] * tpz_xy_zzzz_0[j] + 0.5 * fl1_fx * tpz_y_zzzz_0[j];

            tpx_xxz_xxxx_0[j] =
                pa_x[j] * tpx_xz_xxxx_0[j] + 0.5 * fl1_fx * tpx_z_xxxx_0[j] + 2.0 * fl1_fx * tpx_xz_xxx_0[j] - fl1_fgb * fl1_fx * ts_xz_xxxx_0[j];

            tpy_xxz_xxxx_0[j] = pa_x[j] * tpy_xz_xxxx_0[j] + 0.5 * fl1_fx * tpy_z_xxxx_0[j] + 2.0 * fl1_fx * tpy_xz_xxx_0[j];

            tpz_xxz_xxxx_0[j] = pa_x[j] * tpz_xz_xxxx_0[j] + 0.5 * fl1_fx * tpz_z_xxxx_0[j] + 2.0 * fl1_fx * tpz_xz_xxx_0[j];

            tpx_xxz_xxxy_0[j] =
                pa_x[j] * tpx_xz_xxxy_0[j] + 0.5 * fl1_fx * tpx_z_xxxy_0[j] + 1.5 * fl1_fx * tpx_xz_xxy_0[j] - fl1_fgb * fl1_fx * ts_xz_xxxy_0[j];

            tpy_xxz_xxxy_0[j] = pa_x[j] * tpy_xz_xxxy_0[j] + 0.5 * fl1_fx * tpy_z_xxxy_0[j] + 1.5 * fl1_fx * tpy_xz_xxy_0[j];

            tpz_xxz_xxxy_0[j] = pa_x[j] * tpz_xz_xxxy_0[j] + 0.5 * fl1_fx * tpz_z_xxxy_0[j] + 1.5 * fl1_fx * tpz_xz_xxy_0[j];

            tpx_xxz_xxxz_0[j] =
                pa_x[j] * tpx_xz_xxxz_0[j] + 0.5 * fl1_fx * tpx_z_xxxz_0[j] + 1.5 * fl1_fx * tpx_xz_xxz_0[j] - fl1_fgb * fl1_fx * ts_xz_xxxz_0[j];

            tpy_xxz_xxxz_0[j] = pa_x[j] * tpy_xz_xxxz_0[j] + 0.5 * fl1_fx * tpy_z_xxxz_0[j] + 1.5 * fl1_fx * tpy_xz_xxz_0[j];

            tpz_xxz_xxxz_0[j] = pa_x[j] * tpz_xz_xxxz_0[j] + 0.5 * fl1_fx * tpz_z_xxxz_0[j] + 1.5 * fl1_fx * tpz_xz_xxz_0[j];

            tpx_xxz_xxyy_0[j] =
                pa_x[j] * tpx_xz_xxyy_0[j] + 0.5 * fl1_fx * tpx_z_xxyy_0[j] + fl1_fx * tpx_xz_xyy_0[j] - fl1_fgb * fl1_fx * ts_xz_xxyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_100_150(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (100,150)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 49);

        auto tpy_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tpz_yy_xxyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tpy_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 23);

        auto tpz_xz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 23);

        auto tpx_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 24);

        auto tpy_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 24);

        auto tpz_xz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 24);

        auto tpx_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 25);

        auto tpy_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 25);

        auto tpz_xz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 25);

        auto tpx_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 26);

        auto tpy_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 26);

        auto tpz_xz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 26);

        auto tpx_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 27);

        auto tpy_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 27);

        auto tpz_xz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 27);

        auto tpx_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 28);

        auto tpy_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 28);

        auto tpz_xz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 28);

        auto tpx_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 29);

        auto tpy_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 29);

        auto tpz_xz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 29);

        auto tpx_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 30);

        auto tpy_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 31);

        auto tpy_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 32);

        auto tpy_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 33);

        auto tpy_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 34);

        auto tpy_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto ts_xz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 34);

        auto ts_xz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 35);

        auto ts_xz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 36);

        auto ts_xz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 37);

        auto ts_xz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 38);

        auto ts_xz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 39);

        auto ts_xz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 40);

        auto ts_xz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 41);

        auto ts_xz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 42);

        auto ts_xz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 43);

        auto ts_xz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 44);

        auto ts_yy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 45);

        auto ts_yy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 46);

        auto ts_yy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 47);

        auto ts_yy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 48);

        auto ts_yy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 49);

        // set up pointers to integrals

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

        auto tpx_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 49);

        auto tpy_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tpz_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 49);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxz_xxyz_0, tpx_xxz_xxzz_0, tpx_xxz_xyyy_0, \
                                     tpx_xxz_xyyz_0, tpx_xxz_xyzz_0, tpx_xxz_xzzz_0, tpx_xxz_yyyy_0, tpx_xxz_yyyz_0, \
                                     tpx_xxz_yyzz_0, tpx_xxz_yzzz_0, tpx_xxz_zzzz_0, tpx_xyy_xxxx_0, tpx_xyy_xxxy_0, \
                                     tpx_xyy_xxxz_0, tpx_xyy_xxyy_0, tpx_xyy_xxyz_0, tpx_xz_xxyz_0, tpx_xz_xxzz_0, \
                                     tpx_xz_xyyy_0, tpx_xz_xyyz_0, tpx_xz_xyz_0, tpx_xz_xyzz_0, tpx_xz_xzz_0, \
                                     tpx_xz_xzzz_0, tpx_xz_yyy_0, tpx_xz_yyyy_0, tpx_xz_yyyz_0, tpx_xz_yyz_0, \
                                     tpx_xz_yyzz_0, tpx_xz_yzz_0, tpx_xz_yzzz_0, tpx_xz_zzz_0, tpx_xz_zzzz_0, \
                                     tpx_yy_xxx_0, tpx_yy_xxxx_0, tpx_yy_xxxy_0, tpx_yy_xxxz_0, tpx_yy_xxy_0, \
                                     tpx_yy_xxyy_0, tpx_yy_xxyz_0, tpx_yy_xxz_0, tpx_yy_xyy_0, tpx_yy_xyz_0, \
                                     tpx_z_xxyz_0, tpx_z_xxzz_0, tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyzz_0, tpx_z_xzzz_0, \
                                     tpx_z_yyyy_0, tpx_z_yyyz_0, tpx_z_yyzz_0, tpx_z_yzzz_0, tpx_z_zzzz_0, \
                                     tpy_xxz_xxyy_0, tpy_xxz_xxyz_0, tpy_xxz_xxzz_0, tpy_xxz_xyyy_0, tpy_xxz_xyyz_0, \
                                     tpy_xxz_xyzz_0, tpy_xxz_xzzz_0, tpy_xxz_yyyy_0, tpy_xxz_yyyz_0, tpy_xxz_yyzz_0, \
                                     tpy_xxz_yzzz_0, tpy_xxz_zzzz_0, tpy_xyy_xxxx_0, tpy_xyy_xxxy_0, tpy_xyy_xxxz_0, \
                                     tpy_xyy_xxyy_0, tpy_xyy_xxyz_0, tpy_xz_xxyy_0, tpy_xz_xxyz_0, tpy_xz_xxzz_0, \
                                     tpy_xz_xyy_0, tpy_xz_xyyy_0, tpy_xz_xyyz_0, tpy_xz_xyz_0, tpy_xz_xyzz_0, \
                                     tpy_xz_xzz_0, tpy_xz_xzzz_0, tpy_xz_yyy_0, tpy_xz_yyyy_0, tpy_xz_yyyz_0, \
                                     tpy_xz_yyz_0, tpy_xz_yyzz_0, tpy_xz_yzz_0, tpy_xz_yzzz_0, tpy_xz_zzz_0, \
                                     tpy_xz_zzzz_0, tpy_yy_xxx_0, tpy_yy_xxxx_0, tpy_yy_xxxy_0, tpy_yy_xxxz_0, \
                                     tpy_yy_xxy_0, tpy_yy_xxyy_0, tpy_yy_xxyz_0, tpy_yy_xxz_0, tpy_yy_xyy_0, \
                                     tpy_yy_xyz_0, tpy_z_xxyy_0, tpy_z_xxyz_0, tpy_z_xxzz_0, tpy_z_xyyy_0, tpy_z_xyyz_0, \
                                     tpy_z_xyzz_0, tpy_z_xzzz_0, tpy_z_yyyy_0, tpy_z_yyyz_0, tpy_z_yyzz_0, tpy_z_yzzz_0, \
                                     tpy_z_zzzz_0, tpz_xxz_xxyy_0, tpz_xxz_xxyz_0, tpz_xxz_xxzz_0, tpz_xxz_xyyy_0, \
                                     tpz_xxz_xyyz_0, tpz_xxz_xyzz_0, tpz_xxz_xzzz_0, tpz_xxz_yyyy_0, tpz_xxz_yyyz_0, \
                                     tpz_xxz_yyzz_0, tpz_xxz_yzzz_0, tpz_xxz_zzzz_0, tpz_xyy_xxxx_0, tpz_xyy_xxxy_0, \
                                     tpz_xyy_xxxz_0, tpz_xyy_xxyy_0, tpz_xyy_xxyz_0, tpz_xz_xxyy_0, tpz_xz_xxyz_0, \
                                     tpz_xz_xxzz_0, tpz_xz_xyy_0, tpz_xz_xyyy_0, tpz_xz_xyyz_0, tpz_xz_xyz_0, \
                                     tpz_xz_xyzz_0, tpz_xz_xzz_0, tpz_xz_xzzz_0, tpz_xz_yyy_0, tpz_xz_yyyy_0, \
                                     tpz_xz_yyyz_0, tpz_xz_yyz_0, tpz_xz_yyzz_0, tpz_xz_yzz_0, tpz_xz_yzzz_0, \
                                     tpz_xz_zzz_0, tpz_xz_zzzz_0, tpz_yy_xxx_0, tpz_yy_xxxx_0, tpz_yy_xxxy_0, \
                                     tpz_yy_xxxz_0, tpz_yy_xxy_0, tpz_yy_xxyy_0, tpz_yy_xxyz_0, tpz_yy_xxz_0, \
                                     tpz_yy_xyy_0, tpz_yy_xyz_0, tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxzz_0, tpz_z_xyyy_0, \
                                     tpz_z_xyyz_0, tpz_z_xyzz_0, tpz_z_xzzz_0, tpz_z_yyyy_0, tpz_z_yyyz_0, tpz_z_yyzz_0, \
                                     tpz_z_yzzz_0, tpz_z_zzzz_0, ts_xz_xxyz_0, ts_xz_xxzz_0, ts_xz_xyyy_0, ts_xz_xyyz_0, \
                                     ts_xz_xyzz_0, ts_xz_xzzz_0, ts_xz_yyyy_0, ts_xz_yyyz_0, ts_xz_yyzz_0, ts_xz_yzzz_0, \
                                     ts_xz_zzzz_0, ts_yy_xxxx_0, ts_yy_xxxy_0, ts_yy_xxxz_0, ts_yy_xxyy_0, ts_yy_xxyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_xxz_xxyy_0[j] = pa_x[j] * tpy_xz_xxyy_0[j] + 0.5 * fl1_fx * tpy_z_xxyy_0[j] + fl1_fx * tpy_xz_xyy_0[j];

            tpz_xxz_xxyy_0[j] = pa_x[j] * tpz_xz_xxyy_0[j] + 0.5 * fl1_fx * tpz_z_xxyy_0[j] + fl1_fx * tpz_xz_xyy_0[j];

            tpx_xxz_xxyz_0[j] =
                pa_x[j] * tpx_xz_xxyz_0[j] + 0.5 * fl1_fx * tpx_z_xxyz_0[j] + fl1_fx * tpx_xz_xyz_0[j] - fl1_fgb * fl1_fx * ts_xz_xxyz_0[j];

            tpy_xxz_xxyz_0[j] = pa_x[j] * tpy_xz_xxyz_0[j] + 0.5 * fl1_fx * tpy_z_xxyz_0[j] + fl1_fx * tpy_xz_xyz_0[j];

            tpz_xxz_xxyz_0[j] = pa_x[j] * tpz_xz_xxyz_0[j] + 0.5 * fl1_fx * tpz_z_xxyz_0[j] + fl1_fx * tpz_xz_xyz_0[j];

            tpx_xxz_xxzz_0[j] =
                pa_x[j] * tpx_xz_xxzz_0[j] + 0.5 * fl1_fx * tpx_z_xxzz_0[j] + fl1_fx * tpx_xz_xzz_0[j] - fl1_fgb * fl1_fx * ts_xz_xxzz_0[j];

            tpy_xxz_xxzz_0[j] = pa_x[j] * tpy_xz_xxzz_0[j] + 0.5 * fl1_fx * tpy_z_xxzz_0[j] + fl1_fx * tpy_xz_xzz_0[j];

            tpz_xxz_xxzz_0[j] = pa_x[j] * tpz_xz_xxzz_0[j] + 0.5 * fl1_fx * tpz_z_xxzz_0[j] + fl1_fx * tpz_xz_xzz_0[j];

            tpx_xxz_xyyy_0[j] =
                pa_x[j] * tpx_xz_xyyy_0[j] + 0.5 * fl1_fx * tpx_z_xyyy_0[j] + 0.5 * fl1_fx * tpx_xz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xz_xyyy_0[j];

            tpy_xxz_xyyy_0[j] = pa_x[j] * tpy_xz_xyyy_0[j] + 0.5 * fl1_fx * tpy_z_xyyy_0[j] + 0.5 * fl1_fx * tpy_xz_yyy_0[j];

            tpz_xxz_xyyy_0[j] = pa_x[j] * tpz_xz_xyyy_0[j] + 0.5 * fl1_fx * tpz_z_xyyy_0[j] + 0.5 * fl1_fx * tpz_xz_yyy_0[j];

            tpx_xxz_xyyz_0[j] =
                pa_x[j] * tpx_xz_xyyz_0[j] + 0.5 * fl1_fx * tpx_z_xyyz_0[j] + 0.5 * fl1_fx * tpx_xz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xz_xyyz_0[j];

            tpy_xxz_xyyz_0[j] = pa_x[j] * tpy_xz_xyyz_0[j] + 0.5 * fl1_fx * tpy_z_xyyz_0[j] + 0.5 * fl1_fx * tpy_xz_yyz_0[j];

            tpz_xxz_xyyz_0[j] = pa_x[j] * tpz_xz_xyyz_0[j] + 0.5 * fl1_fx * tpz_z_xyyz_0[j] + 0.5 * fl1_fx * tpz_xz_yyz_0[j];

            tpx_xxz_xyzz_0[j] =
                pa_x[j] * tpx_xz_xyzz_0[j] + 0.5 * fl1_fx * tpx_z_xyzz_0[j] + 0.5 * fl1_fx * tpx_xz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xz_xyzz_0[j];

            tpy_xxz_xyzz_0[j] = pa_x[j] * tpy_xz_xyzz_0[j] + 0.5 * fl1_fx * tpy_z_xyzz_0[j] + 0.5 * fl1_fx * tpy_xz_yzz_0[j];

            tpz_xxz_xyzz_0[j] = pa_x[j] * tpz_xz_xyzz_0[j] + 0.5 * fl1_fx * tpz_z_xyzz_0[j] + 0.5 * fl1_fx * tpz_xz_yzz_0[j];

            tpx_xxz_xzzz_0[j] =
                pa_x[j] * tpx_xz_xzzz_0[j] + 0.5 * fl1_fx * tpx_z_xzzz_0[j] + 0.5 * fl1_fx * tpx_xz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xz_xzzz_0[j];

            tpy_xxz_xzzz_0[j] = pa_x[j] * tpy_xz_xzzz_0[j] + 0.5 * fl1_fx * tpy_z_xzzz_0[j] + 0.5 * fl1_fx * tpy_xz_zzz_0[j];

            tpz_xxz_xzzz_0[j] = pa_x[j] * tpz_xz_xzzz_0[j] + 0.5 * fl1_fx * tpz_z_xzzz_0[j] + 0.5 * fl1_fx * tpz_xz_zzz_0[j];

            tpx_xxz_yyyy_0[j] = pa_x[j] * tpx_xz_yyyy_0[j] + 0.5 * fl1_fx * tpx_z_yyyy_0[j] - fl1_fgb * fl1_fx * ts_xz_yyyy_0[j];

            tpy_xxz_yyyy_0[j] = pa_x[j] * tpy_xz_yyyy_0[j] + 0.5 * fl1_fx * tpy_z_yyyy_0[j];

            tpz_xxz_yyyy_0[j] = pa_x[j] * tpz_xz_yyyy_0[j] + 0.5 * fl1_fx * tpz_z_yyyy_0[j];

            tpx_xxz_yyyz_0[j] = pa_x[j] * tpx_xz_yyyz_0[j] + 0.5 * fl1_fx * tpx_z_yyyz_0[j] - fl1_fgb * fl1_fx * ts_xz_yyyz_0[j];

            tpy_xxz_yyyz_0[j] = pa_x[j] * tpy_xz_yyyz_0[j] + 0.5 * fl1_fx * tpy_z_yyyz_0[j];

            tpz_xxz_yyyz_0[j] = pa_x[j] * tpz_xz_yyyz_0[j] + 0.5 * fl1_fx * tpz_z_yyyz_0[j];

            tpx_xxz_yyzz_0[j] = pa_x[j] * tpx_xz_yyzz_0[j] + 0.5 * fl1_fx * tpx_z_yyzz_0[j] - fl1_fgb * fl1_fx * ts_xz_yyzz_0[j];

            tpy_xxz_yyzz_0[j] = pa_x[j] * tpy_xz_yyzz_0[j] + 0.5 * fl1_fx * tpy_z_yyzz_0[j];

            tpz_xxz_yyzz_0[j] = pa_x[j] * tpz_xz_yyzz_0[j] + 0.5 * fl1_fx * tpz_z_yyzz_0[j];

            tpx_xxz_yzzz_0[j] = pa_x[j] * tpx_xz_yzzz_0[j] + 0.5 * fl1_fx * tpx_z_yzzz_0[j] - fl1_fgb * fl1_fx * ts_xz_yzzz_0[j];

            tpy_xxz_yzzz_0[j] = pa_x[j] * tpy_xz_yzzz_0[j] + 0.5 * fl1_fx * tpy_z_yzzz_0[j];

            tpz_xxz_yzzz_0[j] = pa_x[j] * tpz_xz_yzzz_0[j] + 0.5 * fl1_fx * tpz_z_yzzz_0[j];

            tpx_xxz_zzzz_0[j] = pa_x[j] * tpx_xz_zzzz_0[j] + 0.5 * fl1_fx * tpx_z_zzzz_0[j] - fl1_fgb * fl1_fx * ts_xz_zzzz_0[j];

            tpy_xxz_zzzz_0[j] = pa_x[j] * tpy_xz_zzzz_0[j] + 0.5 * fl1_fx * tpy_z_zzzz_0[j];

            tpz_xxz_zzzz_0[j] = pa_x[j] * tpz_xz_zzzz_0[j] + 0.5 * fl1_fx * tpz_z_zzzz_0[j];

            tpx_xyy_xxxx_0[j] = pa_x[j] * tpx_yy_xxxx_0[j] + 2.0 * fl1_fx * tpx_yy_xxx_0[j] - fl1_fgb * fl1_fx * ts_yy_xxxx_0[j];

            tpy_xyy_xxxx_0[j] = pa_x[j] * tpy_yy_xxxx_0[j] + 2.0 * fl1_fx * tpy_yy_xxx_0[j];

            tpz_xyy_xxxx_0[j] = pa_x[j] * tpz_yy_xxxx_0[j] + 2.0 * fl1_fx * tpz_yy_xxx_0[j];

            tpx_xyy_xxxy_0[j] = pa_x[j] * tpx_yy_xxxy_0[j] + 1.5 * fl1_fx * tpx_yy_xxy_0[j] - fl1_fgb * fl1_fx * ts_yy_xxxy_0[j];

            tpy_xyy_xxxy_0[j] = pa_x[j] * tpy_yy_xxxy_0[j] + 1.5 * fl1_fx * tpy_yy_xxy_0[j];

            tpz_xyy_xxxy_0[j] = pa_x[j] * tpz_yy_xxxy_0[j] + 1.5 * fl1_fx * tpz_yy_xxy_0[j];

            tpx_xyy_xxxz_0[j] = pa_x[j] * tpx_yy_xxxz_0[j] + 1.5 * fl1_fx * tpx_yy_xxz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxxz_0[j];

            tpy_xyy_xxxz_0[j] = pa_x[j] * tpy_yy_xxxz_0[j] + 1.5 * fl1_fx * tpy_yy_xxz_0[j];

            tpz_xyy_xxxz_0[j] = pa_x[j] * tpz_yy_xxxz_0[j] + 1.5 * fl1_fx * tpz_yy_xxz_0[j];

            tpx_xyy_xxyy_0[j] = pa_x[j] * tpx_yy_xxyy_0[j] + fl1_fx * tpx_yy_xyy_0[j] - fl1_fgb * fl1_fx * ts_yy_xxyy_0[j];

            tpy_xyy_xxyy_0[j] = pa_x[j] * tpy_yy_xxyy_0[j] + fl1_fx * tpy_yy_xyy_0[j];

            tpz_xyy_xxyy_0[j] = pa_x[j] * tpz_yy_xxyy_0[j] + fl1_fx * tpz_yy_xyy_0[j];

            tpx_xyy_xxyz_0[j] = pa_x[j] * tpx_yy_xxyz_0[j] + fl1_fx * tpx_yy_xyz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxyz_0[j];

            tpy_xyy_xxyz_0[j] = pa_x[j] * tpy_yy_xxyz_0[j] + fl1_fx * tpy_yy_xyz_0[j];

            tpz_xyy_xxyz_0[j] = pa_x[j] * tpz_yy_xxyz_0[j] + fl1_fx * tpz_yy_xyz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_150_200(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (150,200)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 65);

        auto tpy_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tpz_yz_xxzz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tpx_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 66);

        auto tpy_yz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tpx_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 35);

        auto tpy_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tpx_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 36);

        auto tpy_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 37);

        auto tpy_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 38);

        auto tpy_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 39);

        auto tpy_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 40);

        auto tpy_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 41);

        auto tpy_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 42);

        auto tpy_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 43);

        auto tpy_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 44);

        auto tpy_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 45);

        auto tpy_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 46);

        auto tpy_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto ts_yy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 50);

        auto ts_yy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 51);

        auto ts_yy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 52);

        auto ts_yy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 53);

        auto ts_yy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 54);

        auto ts_yy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 55);

        auto ts_yy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 56);

        auto ts_yy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 57);

        auto ts_yy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 58);

        auto ts_yy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 59);

        auto ts_yz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 60);

        auto ts_yz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 61);

        auto ts_yz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 62);

        auto ts_yz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 63);

        auto ts_yz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 64);

        auto ts_yz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 65);

        auto ts_yz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 66);

        // set up pointers to integrals

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

        auto tpx_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 65);

        auto tpy_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tpz_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tpx_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 66);

        auto tpy_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 66);

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyy_xxzz_0, tpx_xyy_xyyy_0, tpx_xyy_xyyz_0, \
                                     tpx_xyy_xyzz_0, tpx_xyy_xzzz_0, tpx_xyy_yyyy_0, tpx_xyy_yyyz_0, tpx_xyy_yyzz_0, \
                                     tpx_xyy_yzzz_0, tpx_xyy_zzzz_0, tpx_xyz_xxxx_0, tpx_xyz_xxxy_0, tpx_xyz_xxxz_0, \
                                     tpx_xyz_xxyy_0, tpx_xyz_xxyz_0, tpx_xyz_xxzz_0, tpx_xyz_xyyy_0, tpx_yy_xxzz_0, \
                                     tpx_yy_xyyy_0, tpx_yy_xyyz_0, tpx_yy_xyzz_0, tpx_yy_xzz_0, tpx_yy_xzzz_0, \
                                     tpx_yy_yyy_0, tpx_yy_yyyy_0, tpx_yy_yyyz_0, tpx_yy_yyz_0, tpx_yy_yyzz_0, \
                                     tpx_yy_yzz_0, tpx_yy_yzzz_0, tpx_yy_zzz_0, tpx_yy_zzzz_0, tpx_yz_xxx_0, \
                                     tpx_yz_xxxx_0, tpx_yz_xxxy_0, tpx_yz_xxxz_0, tpx_yz_xxy_0, tpx_yz_xxyy_0, \
                                     tpx_yz_xxyz_0, tpx_yz_xxz_0, tpx_yz_xxzz_0, tpx_yz_xyy_0, tpx_yz_xyyy_0, \
                                     tpx_yz_xyz_0, tpx_yz_xzz_0, tpx_yz_yyy_0, tpy_xyy_xxzz_0, tpy_xyy_xyyy_0, \
                                     tpy_xyy_xyyz_0, tpy_xyy_xyzz_0, tpy_xyy_xzzz_0, tpy_xyy_yyyy_0, tpy_xyy_yyyz_0, \
                                     tpy_xyy_yyzz_0, tpy_xyy_yzzz_0, tpy_xyy_zzzz_0, tpy_xyz_xxxx_0, tpy_xyz_xxxy_0, \
                                     tpy_xyz_xxxz_0, tpy_xyz_xxyy_0, tpy_xyz_xxyz_0, tpy_xyz_xxzz_0, tpy_xyz_xyyy_0, \
                                     tpy_yy_xxzz_0, tpy_yy_xyyy_0, tpy_yy_xyyz_0, tpy_yy_xyzz_0, tpy_yy_xzz_0, \
                                     tpy_yy_xzzz_0, tpy_yy_yyy_0, tpy_yy_yyyy_0, tpy_yy_yyyz_0, tpy_yy_yyz_0, \
                                     tpy_yy_yyzz_0, tpy_yy_yzz_0, tpy_yy_yzzz_0, tpy_yy_zzz_0, tpy_yy_zzzz_0, \
                                     tpy_yz_xxx_0, tpy_yz_xxxx_0, tpy_yz_xxxy_0, tpy_yz_xxxz_0, tpy_yz_xxy_0, \
                                     tpy_yz_xxyy_0, tpy_yz_xxyz_0, tpy_yz_xxz_0, tpy_yz_xxzz_0, tpy_yz_xyy_0, \
                                     tpy_yz_xyyy_0, tpy_yz_xyz_0, tpy_yz_xzz_0, tpy_yz_yyy_0, tpz_xyy_xxzz_0, \
                                     tpz_xyy_xyyy_0, tpz_xyy_xyyz_0, tpz_xyy_xyzz_0, tpz_xyy_xzzz_0, tpz_xyy_yyyy_0, \
                                     tpz_xyy_yyyz_0, tpz_xyy_yyzz_0, tpz_xyy_yzzz_0, tpz_xyy_zzzz_0, tpz_xyz_xxxx_0, \
                                     tpz_xyz_xxxy_0, tpz_xyz_xxxz_0, tpz_xyz_xxyy_0, tpz_xyz_xxyz_0, tpz_xyz_xxzz_0, \
                                     tpz_yy_xxzz_0, tpz_yy_xyyy_0, tpz_yy_xyyz_0, tpz_yy_xyzz_0, tpz_yy_xzz_0, \
                                     tpz_yy_xzzz_0, tpz_yy_yyy_0, tpz_yy_yyyy_0, tpz_yy_yyyz_0, tpz_yy_yyz_0, \
                                     tpz_yy_yyzz_0, tpz_yy_yzz_0, tpz_yy_yzzz_0, tpz_yy_zzz_0, tpz_yy_zzzz_0, \
                                     tpz_yz_xxx_0, tpz_yz_xxxx_0, tpz_yz_xxxy_0, tpz_yz_xxxz_0, tpz_yz_xxy_0, \
                                     tpz_yz_xxyy_0, tpz_yz_xxyz_0, tpz_yz_xxz_0, tpz_yz_xxzz_0, tpz_yz_xyy_0, \
                                     tpz_yz_xyz_0, tpz_yz_xzz_0, ts_yy_xxzz_0, ts_yy_xyyy_0, ts_yy_xyyz_0, ts_yy_xyzz_0, \
                                     ts_yy_xzzz_0, ts_yy_yyyy_0, ts_yy_yyyz_0, ts_yy_yyzz_0, ts_yy_yzzz_0, ts_yy_zzzz_0, \
                                     ts_yz_xxxx_0, ts_yz_xxxy_0, ts_yz_xxxz_0, ts_yz_xxyy_0, ts_yz_xxyz_0, ts_yz_xxzz_0, \
                                     ts_yz_xyyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xyy_xxzz_0[j] = pa_x[j] * tpx_yy_xxzz_0[j] + fl1_fx * tpx_yy_xzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxzz_0[j];

            tpy_xyy_xxzz_0[j] = pa_x[j] * tpy_yy_xxzz_0[j] + fl1_fx * tpy_yy_xzz_0[j];

            tpz_xyy_xxzz_0[j] = pa_x[j] * tpz_yy_xxzz_0[j] + fl1_fx * tpz_yy_xzz_0[j];

            tpx_xyy_xyyy_0[j] = pa_x[j] * tpx_yy_xyyy_0[j] + 0.5 * fl1_fx * tpx_yy_yyy_0[j] - fl1_fgb * fl1_fx * ts_yy_xyyy_0[j];

            tpy_xyy_xyyy_0[j] = pa_x[j] * tpy_yy_xyyy_0[j] + 0.5 * fl1_fx * tpy_yy_yyy_0[j];

            tpz_xyy_xyyy_0[j] = pa_x[j] * tpz_yy_xyyy_0[j] + 0.5 * fl1_fx * tpz_yy_yyy_0[j];

            tpx_xyy_xyyz_0[j] = pa_x[j] * tpx_yy_xyyz_0[j] + 0.5 * fl1_fx * tpx_yy_yyz_0[j] - fl1_fgb * fl1_fx * ts_yy_xyyz_0[j];

            tpy_xyy_xyyz_0[j] = pa_x[j] * tpy_yy_xyyz_0[j] + 0.5 * fl1_fx * tpy_yy_yyz_0[j];

            tpz_xyy_xyyz_0[j] = pa_x[j] * tpz_yy_xyyz_0[j] + 0.5 * fl1_fx * tpz_yy_yyz_0[j];

            tpx_xyy_xyzz_0[j] = pa_x[j] * tpx_yy_xyzz_0[j] + 0.5 * fl1_fx * tpx_yy_yzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xyzz_0[j];

            tpy_xyy_xyzz_0[j] = pa_x[j] * tpy_yy_xyzz_0[j] + 0.5 * fl1_fx * tpy_yy_yzz_0[j];

            tpz_xyy_xyzz_0[j] = pa_x[j] * tpz_yy_xyzz_0[j] + 0.5 * fl1_fx * tpz_yy_yzz_0[j];

            tpx_xyy_xzzz_0[j] = pa_x[j] * tpx_yy_xzzz_0[j] + 0.5 * fl1_fx * tpx_yy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xzzz_0[j];

            tpy_xyy_xzzz_0[j] = pa_x[j] * tpy_yy_xzzz_0[j] + 0.5 * fl1_fx * tpy_yy_zzz_0[j];

            tpz_xyy_xzzz_0[j] = pa_x[j] * tpz_yy_xzzz_0[j] + 0.5 * fl1_fx * tpz_yy_zzz_0[j];

            tpx_xyy_yyyy_0[j] = pa_x[j] * tpx_yy_yyyy_0[j] - fl1_fgb * fl1_fx * ts_yy_yyyy_0[j];

            tpy_xyy_yyyy_0[j] = pa_x[j] * tpy_yy_yyyy_0[j];

            tpz_xyy_yyyy_0[j] = pa_x[j] * tpz_yy_yyyy_0[j];

            tpx_xyy_yyyz_0[j] = pa_x[j] * tpx_yy_yyyz_0[j] - fl1_fgb * fl1_fx * ts_yy_yyyz_0[j];

            tpy_xyy_yyyz_0[j] = pa_x[j] * tpy_yy_yyyz_0[j];

            tpz_xyy_yyyz_0[j] = pa_x[j] * tpz_yy_yyyz_0[j];

            tpx_xyy_yyzz_0[j] = pa_x[j] * tpx_yy_yyzz_0[j] - fl1_fgb * fl1_fx * ts_yy_yyzz_0[j];

            tpy_xyy_yyzz_0[j] = pa_x[j] * tpy_yy_yyzz_0[j];

            tpz_xyy_yyzz_0[j] = pa_x[j] * tpz_yy_yyzz_0[j];

            tpx_xyy_yzzz_0[j] = pa_x[j] * tpx_yy_yzzz_0[j] - fl1_fgb * fl1_fx * ts_yy_yzzz_0[j];

            tpy_xyy_yzzz_0[j] = pa_x[j] * tpy_yy_yzzz_0[j];

            tpz_xyy_yzzz_0[j] = pa_x[j] * tpz_yy_yzzz_0[j];

            tpx_xyy_zzzz_0[j] = pa_x[j] * tpx_yy_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yy_zzzz_0[j];

            tpy_xyy_zzzz_0[j] = pa_x[j] * tpy_yy_zzzz_0[j];

            tpz_xyy_zzzz_0[j] = pa_x[j] * tpz_yy_zzzz_0[j];

            tpx_xyz_xxxx_0[j] = pa_x[j] * tpx_yz_xxxx_0[j] + 2.0 * fl1_fx * tpx_yz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yz_xxxx_0[j];

            tpy_xyz_xxxx_0[j] = pa_x[j] * tpy_yz_xxxx_0[j] + 2.0 * fl1_fx * tpy_yz_xxx_0[j];

            tpz_xyz_xxxx_0[j] = pa_x[j] * tpz_yz_xxxx_0[j] + 2.0 * fl1_fx * tpz_yz_xxx_0[j];

            tpx_xyz_xxxy_0[j] = pa_x[j] * tpx_yz_xxxy_0[j] + 1.5 * fl1_fx * tpx_yz_xxy_0[j] - fl1_fgb * fl1_fx * ts_yz_xxxy_0[j];

            tpy_xyz_xxxy_0[j] = pa_x[j] * tpy_yz_xxxy_0[j] + 1.5 * fl1_fx * tpy_yz_xxy_0[j];

            tpz_xyz_xxxy_0[j] = pa_x[j] * tpz_yz_xxxy_0[j] + 1.5 * fl1_fx * tpz_yz_xxy_0[j];

            tpx_xyz_xxxz_0[j] = pa_x[j] * tpx_yz_xxxz_0[j] + 1.5 * fl1_fx * tpx_yz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxxz_0[j];

            tpy_xyz_xxxz_0[j] = pa_x[j] * tpy_yz_xxxz_0[j] + 1.5 * fl1_fx * tpy_yz_xxz_0[j];

            tpz_xyz_xxxz_0[j] = pa_x[j] * tpz_yz_xxxz_0[j] + 1.5 * fl1_fx * tpz_yz_xxz_0[j];

            tpx_xyz_xxyy_0[j] = pa_x[j] * tpx_yz_xxyy_0[j] + fl1_fx * tpx_yz_xyy_0[j] - fl1_fgb * fl1_fx * ts_yz_xxyy_0[j];

            tpy_xyz_xxyy_0[j] = pa_x[j] * tpy_yz_xxyy_0[j] + fl1_fx * tpy_yz_xyy_0[j];

            tpz_xyz_xxyy_0[j] = pa_x[j] * tpz_yz_xxyy_0[j] + fl1_fx * tpz_yz_xyy_0[j];

            tpx_xyz_xxyz_0[j] = pa_x[j] * tpx_yz_xxyz_0[j] + fl1_fx * tpx_yz_xyz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxyz_0[j];

            tpy_xyz_xxyz_0[j] = pa_x[j] * tpy_yz_xxyz_0[j] + fl1_fx * tpy_yz_xyz_0[j];

            tpz_xyz_xxyz_0[j] = pa_x[j] * tpz_yz_xxyz_0[j] + fl1_fx * tpz_yz_xyz_0[j];

            tpx_xyz_xxzz_0[j] = pa_x[j] * tpx_yz_xxzz_0[j] + fl1_fx * tpx_yz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxzz_0[j];

            tpy_xyz_xxzz_0[j] = pa_x[j] * tpy_yz_xxzz_0[j] + fl1_fx * tpy_yz_xzz_0[j];

            tpz_xyz_xxzz_0[j] = pa_x[j] * tpz_yz_xxzz_0[j] + fl1_fx * tpz_yz_xzz_0[j];

            tpx_xyz_xyyy_0[j] = pa_x[j] * tpx_yz_xyyy_0[j] + 0.5 * fl1_fx * tpx_yz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yz_xyyy_0[j];

            tpy_xyz_xyyy_0[j] = pa_x[j] * tpy_yz_xyyy_0[j] + 0.5 * fl1_fx * tpy_yz_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_200_250(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (200,250)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 81);

        auto tpy_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tpz_zz_xyyy_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tpx_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 82);

        auto tpy_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tpz_zz_xyyz_0 = primBuffer.data(pidx_p_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tpx_zz_xyzz_0 = primBuffer.data(pidx_p_2_4_m0 + 90 * idx + 83);

        auto tpz_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 47);

        auto tpy_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 48);

        auto tpy_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 49);

        auto tpy_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 50);

        auto tpy_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 51);

        auto tpy_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 52);

        auto tpy_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 53);

        auto tpy_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 54);

        auto tpy_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 55);

        auto tpy_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 56);

        auto tpy_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 57);

        auto tpy_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 58);

        auto ts_yz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 67);

        auto ts_yz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 68);

        auto ts_yz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 69);

        auto ts_yz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 70);

        auto ts_yz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 71);

        auto ts_yz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 72);

        auto ts_yz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 73);

        auto ts_yz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 74);

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        // set up pointers to integrals

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

        auto tpx_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 81);

        auto tpy_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tpz_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tpx_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 82);

        auto tpy_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tpz_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tpx_xzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 83);

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyz_xyyz_0, tpx_xyz_xyzz_0, tpx_xyz_xzzz_0, \
                                     tpx_xyz_yyyy_0, tpx_xyz_yyyz_0, tpx_xyz_yyzz_0, tpx_xyz_yzzz_0, tpx_xyz_zzzz_0, \
                                     tpx_xzz_xxxx_0, tpx_xzz_xxxy_0, tpx_xzz_xxxz_0, tpx_xzz_xxyy_0, tpx_xzz_xxyz_0, \
                                     tpx_xzz_xxzz_0, tpx_xzz_xyyy_0, tpx_xzz_xyyz_0, tpx_xzz_xyzz_0, tpx_yz_xyyz_0, \
                                     tpx_yz_xyzz_0, tpx_yz_xzzz_0, tpx_yz_yyyy_0, tpx_yz_yyyz_0, tpx_yz_yyz_0, \
                                     tpx_yz_yyzz_0, tpx_yz_yzz_0, tpx_yz_yzzz_0, tpx_yz_zzz_0, tpx_yz_zzzz_0, \
                                     tpx_zz_xxx_0, tpx_zz_xxxx_0, tpx_zz_xxxy_0, tpx_zz_xxxz_0, tpx_zz_xxy_0, \
                                     tpx_zz_xxyy_0, tpx_zz_xxyz_0, tpx_zz_xxz_0, tpx_zz_xxzz_0, tpx_zz_xyy_0, \
                                     tpx_zz_xyyy_0, tpx_zz_xyyz_0, tpx_zz_xyz_0, tpx_zz_xyzz_0, tpx_zz_xzz_0, \
                                     tpx_zz_yyy_0, tpx_zz_yyz_0, tpx_zz_yzz_0, tpy_xyz_xyyz_0, tpy_xyz_xyzz_0, \
                                     tpy_xyz_xzzz_0, tpy_xyz_yyyy_0, tpy_xyz_yyyz_0, tpy_xyz_yyzz_0, tpy_xyz_yzzz_0, \
                                     tpy_xyz_zzzz_0, tpy_xzz_xxxx_0, tpy_xzz_xxxy_0, tpy_xzz_xxxz_0, tpy_xzz_xxyy_0, \
                                     tpy_xzz_xxyz_0, tpy_xzz_xxzz_0, tpy_xzz_xyyy_0, tpy_xzz_xyyz_0, tpy_yz_xyyz_0, \
                                     tpy_yz_xyzz_0, tpy_yz_xzzz_0, tpy_yz_yyyy_0, tpy_yz_yyyz_0, tpy_yz_yyz_0, \
                                     tpy_yz_yyzz_0, tpy_yz_yzz_0, tpy_yz_yzzz_0, tpy_yz_zzz_0, tpy_yz_zzzz_0, \
                                     tpy_zz_xxx_0, tpy_zz_xxxx_0, tpy_zz_xxxy_0, tpy_zz_xxxz_0, tpy_zz_xxy_0, \
                                     tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxz_0, tpy_zz_xxzz_0, tpy_zz_xyy_0, \
                                     tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpy_zz_xyz_0, tpy_zz_xzz_0, tpy_zz_yyy_0, \
                                     tpy_zz_yyz_0, tpz_xyz_xyyy_0, tpz_xyz_xyyz_0, tpz_xyz_xyzz_0, tpz_xyz_xzzz_0, \
                                     tpz_xyz_yyyy_0, tpz_xyz_yyyz_0, tpz_xyz_yyzz_0, tpz_xyz_yzzz_0, tpz_xyz_zzzz_0, \
                                     tpz_xzz_xxxx_0, tpz_xzz_xxxy_0, tpz_xzz_xxxz_0, tpz_xzz_xxyy_0, tpz_xzz_xxyz_0, \
                                     tpz_xzz_xxzz_0, tpz_xzz_xyyy_0, tpz_xzz_xyyz_0, tpz_yz_xyyy_0, tpz_yz_xyyz_0, \
                                     tpz_yz_xyzz_0, tpz_yz_xzzz_0, tpz_yz_yyy_0, tpz_yz_yyyy_0, tpz_yz_yyyz_0, \
                                     tpz_yz_yyz_0, tpz_yz_yyzz_0, tpz_yz_yzz_0, tpz_yz_yzzz_0, tpz_yz_zzz_0, \
                                     tpz_yz_zzzz_0, tpz_zz_xxx_0, tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, \
                                     tpz_zz_xxy_0, tpz_zz_xxyy_0, tpz_zz_xxyz_0, tpz_zz_xxz_0, tpz_zz_xxzz_0, \
                                     tpz_zz_xyy_0, tpz_zz_xyyy_0, tpz_zz_xyyz_0, tpz_zz_xyz_0, tpz_zz_xzz_0, \
                                     tpz_zz_yyy_0, tpz_zz_yyz_0, ts_yz_xyyz_0, ts_yz_xyzz_0, ts_yz_xzzz_0, ts_yz_yyyy_0, \
                                     ts_yz_yyyz_0, ts_yz_yyzz_0, ts_yz_yzzz_0, ts_yz_zzzz_0, ts_zz_xxxx_0, ts_zz_xxxy_0, \
                                     ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, ts_zz_xxzz_0, ts_zz_xyyy_0, ts_zz_xyyz_0, \
                                     ts_zz_xyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_xyz_xyyy_0[j] = pa_x[j] * tpz_yz_xyyy_0[j] + 0.5 * fl1_fx * tpz_yz_yyy_0[j];

            tpx_xyz_xyyz_0[j] = pa_x[j] * tpx_yz_xyyz_0[j] + 0.5 * fl1_fx * tpx_yz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yz_xyyz_0[j];

            tpy_xyz_xyyz_0[j] = pa_x[j] * tpy_yz_xyyz_0[j] + 0.5 * fl1_fx * tpy_yz_yyz_0[j];

            tpz_xyz_xyyz_0[j] = pa_x[j] * tpz_yz_xyyz_0[j] + 0.5 * fl1_fx * tpz_yz_yyz_0[j];

            tpx_xyz_xyzz_0[j] = pa_x[j] * tpx_yz_xyzz_0[j] + 0.5 * fl1_fx * tpx_yz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xyzz_0[j];

            tpy_xyz_xyzz_0[j] = pa_x[j] * tpy_yz_xyzz_0[j] + 0.5 * fl1_fx * tpy_yz_yzz_0[j];

            tpz_xyz_xyzz_0[j] = pa_x[j] * tpz_yz_xyzz_0[j] + 0.5 * fl1_fx * tpz_yz_yzz_0[j];

            tpx_xyz_xzzz_0[j] = pa_x[j] * tpx_yz_xzzz_0[j] + 0.5 * fl1_fx * tpx_yz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xzzz_0[j];

            tpy_xyz_xzzz_0[j] = pa_x[j] * tpy_yz_xzzz_0[j] + 0.5 * fl1_fx * tpy_yz_zzz_0[j];

            tpz_xyz_xzzz_0[j] = pa_x[j] * tpz_yz_xzzz_0[j] + 0.5 * fl1_fx * tpz_yz_zzz_0[j];

            tpx_xyz_yyyy_0[j] = pa_x[j] * tpx_yz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_yz_yyyy_0[j];

            tpy_xyz_yyyy_0[j] = pa_x[j] * tpy_yz_yyyy_0[j];

            tpz_xyz_yyyy_0[j] = pa_x[j] * tpz_yz_yyyy_0[j];

            tpx_xyz_yyyz_0[j] = pa_x[j] * tpx_yz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_yz_yyyz_0[j];

            tpy_xyz_yyyz_0[j] = pa_x[j] * tpy_yz_yyyz_0[j];

            tpz_xyz_yyyz_0[j] = pa_x[j] * tpz_yz_yyyz_0[j];

            tpx_xyz_yyzz_0[j] = pa_x[j] * tpx_yz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_yz_yyzz_0[j];

            tpy_xyz_yyzz_0[j] = pa_x[j] * tpy_yz_yyzz_0[j];

            tpz_xyz_yyzz_0[j] = pa_x[j] * tpz_yz_yyzz_0[j];

            tpx_xyz_yzzz_0[j] = pa_x[j] * tpx_yz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_yz_yzzz_0[j];

            tpy_xyz_yzzz_0[j] = pa_x[j] * tpy_yz_yzzz_0[j];

            tpz_xyz_yzzz_0[j] = pa_x[j] * tpz_yz_yzzz_0[j];

            tpx_xyz_zzzz_0[j] = pa_x[j] * tpx_yz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yz_zzzz_0[j];

            tpy_xyz_zzzz_0[j] = pa_x[j] * tpy_yz_zzzz_0[j];

            tpz_xyz_zzzz_0[j] = pa_x[j] * tpz_yz_zzzz_0[j];

            tpx_xzz_xxxx_0[j] = pa_x[j] * tpx_zz_xxxx_0[j] + 2.0 * fl1_fx * tpx_zz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxx_0[j];

            tpy_xzz_xxxx_0[j] = pa_x[j] * tpy_zz_xxxx_0[j] + 2.0 * fl1_fx * tpy_zz_xxx_0[j];

            tpz_xzz_xxxx_0[j] = pa_x[j] * tpz_zz_xxxx_0[j] + 2.0 * fl1_fx * tpz_zz_xxx_0[j];

            tpx_xzz_xxxy_0[j] = pa_x[j] * tpx_zz_xxxy_0[j] + 1.5 * fl1_fx * tpx_zz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxy_0[j];

            tpy_xzz_xxxy_0[j] = pa_x[j] * tpy_zz_xxxy_0[j] + 1.5 * fl1_fx * tpy_zz_xxy_0[j];

            tpz_xzz_xxxy_0[j] = pa_x[j] * tpz_zz_xxxy_0[j] + 1.5 * fl1_fx * tpz_zz_xxy_0[j];

            tpx_xzz_xxxz_0[j] = pa_x[j] * tpx_zz_xxxz_0[j] + 1.5 * fl1_fx * tpx_zz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxz_0[j];

            tpy_xzz_xxxz_0[j] = pa_x[j] * tpy_zz_xxxz_0[j] + 1.5 * fl1_fx * tpy_zz_xxz_0[j];

            tpz_xzz_xxxz_0[j] = pa_x[j] * tpz_zz_xxxz_0[j] + 1.5 * fl1_fx * tpz_zz_xxz_0[j];

            tpx_xzz_xxyy_0[j] = pa_x[j] * tpx_zz_xxyy_0[j] + fl1_fx * tpx_zz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxyy_0[j];

            tpy_xzz_xxyy_0[j] = pa_x[j] * tpy_zz_xxyy_0[j] + fl1_fx * tpy_zz_xyy_0[j];

            tpz_xzz_xxyy_0[j] = pa_x[j] * tpz_zz_xxyy_0[j] + fl1_fx * tpz_zz_xyy_0[j];

            tpx_xzz_xxyz_0[j] = pa_x[j] * tpx_zz_xxyz_0[j] + fl1_fx * tpx_zz_xyz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxyz_0[j];

            tpy_xzz_xxyz_0[j] = pa_x[j] * tpy_zz_xxyz_0[j] + fl1_fx * tpy_zz_xyz_0[j];

            tpz_xzz_xxyz_0[j] = pa_x[j] * tpz_zz_xxyz_0[j] + fl1_fx * tpz_zz_xyz_0[j];

            tpx_xzz_xxzz_0[j] = pa_x[j] * tpx_zz_xxzz_0[j] + fl1_fx * tpx_zz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxzz_0[j];

            tpy_xzz_xxzz_0[j] = pa_x[j] * tpy_zz_xxzz_0[j] + fl1_fx * tpy_zz_xzz_0[j];

            tpz_xzz_xxzz_0[j] = pa_x[j] * tpz_zz_xxzz_0[j] + fl1_fx * tpz_zz_xzz_0[j];

            tpx_xzz_xyyy_0[j] = pa_x[j] * tpx_zz_xyyy_0[j] + 0.5 * fl1_fx * tpx_zz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyyy_0[j];

            tpy_xzz_xyyy_0[j] = pa_x[j] * tpy_zz_xyyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyy_0[j];

            tpz_xzz_xyyy_0[j] = pa_x[j] * tpz_zz_xyyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyy_0[j];

            tpx_xzz_xyyz_0[j] = pa_x[j] * tpx_zz_xyyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyyz_0[j];

            tpy_xzz_xyyz_0[j] = pa_x[j] * tpy_zz_xyyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyz_0[j];

            tpz_xzz_xyyz_0[j] = pa_x[j] * tpz_zz_xyyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyz_0[j];

            tpx_xzz_xyzz_0[j] = pa_x[j] * tpx_zz_xyzz_0[j] + 0.5 * fl1_fx * tpx_zz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_250_300(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (250,300)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 15);

        auto tpy_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tpx_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 16);

        auto tpy_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tpx_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 17);

        auto tpy_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tpx_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 18);

        auto tpy_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tpx_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 19);

        auto tpy_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tpx_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 20);

        auto tpy_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tpx_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 21);

        auto tpy_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tpx_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 22);

        auto tpy_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tpx_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 23);

        auto tpy_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tpx_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 24);

        auto tpy_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tpx_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 30);

        auto tpy_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_yy_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 31);

        auto tpy_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_yy_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 32);

        auto tpy_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_yy_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 33);

        auto tpy_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 34);

        auto tpy_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 34);

        auto tpx_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 35);

        auto tpy_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_yy_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 35);

        auto tpy_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 59);

        auto tpy_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto ts_yy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 45);

        auto ts_yy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 46);

        auto ts_yy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 47);

        auto ts_yy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 48);

        auto ts_yy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 49);

        auto ts_yy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 50);

        auto ts_yy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 51);

        auto ts_yy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 52);

        auto ts_yy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 53);

        auto ts_yy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 54);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        auto ts_zz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 88);

        auto ts_zz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 89);

        // set up pointers to integrals

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

        auto tpx_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 97);

        auto tpy_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tpz_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tpx_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 98);

        auto tpy_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tpz_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tpx_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 99);

        auto tpy_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tpz_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 99);

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_xzz_xzzz_0, tpx_xzz_yyyy_0, tpx_xzz_yyyz_0, \
                                     tpx_xzz_yyzz_0, tpx_xzz_yzzz_0, tpx_xzz_zzzz_0, tpx_y_xxxx_0, tpx_y_xxxy_0, \
                                     tpx_y_xxxz_0, tpx_y_xxyy_0, tpx_y_xxyz_0, tpx_y_xxzz_0, tpx_y_xyyy_0, tpx_y_xyyz_0, \
                                     tpx_y_xyzz_0, tpx_y_xzzz_0, tpx_yy_xxx_0, tpx_yy_xxxx_0, tpx_yy_xxxy_0, \
                                     tpx_yy_xxxz_0, tpx_yy_xxy_0, tpx_yy_xxyy_0, tpx_yy_xxyz_0, tpx_yy_xxz_0, \
                                     tpx_yy_xxzz_0, tpx_yy_xyy_0, tpx_yy_xyyy_0, tpx_yy_xyyz_0, tpx_yy_xyz_0, \
                                     tpx_yy_xyzz_0, tpx_yy_xzz_0, tpx_yy_xzzz_0, tpx_yyy_xxxx_0, tpx_yyy_xxxy_0, \
                                     tpx_yyy_xxxz_0, tpx_yyy_xxyy_0, tpx_yyy_xxyz_0, tpx_yyy_xxzz_0, tpx_yyy_xyyy_0, \
                                     tpx_yyy_xyyz_0, tpx_yyy_xyzz_0, tpx_yyy_xzzz_0, tpx_zz_xzzz_0, tpx_zz_yyyy_0, \
                                     tpx_zz_yyyz_0, tpx_zz_yyzz_0, tpx_zz_yzzz_0, tpx_zz_zzz_0, tpx_zz_zzzz_0, \
                                     tpy_xzz_xyzz_0, tpy_xzz_xzzz_0, tpy_xzz_yyyy_0, tpy_xzz_yyyz_0, tpy_xzz_yyzz_0, \
                                     tpy_xzz_yzzz_0, tpy_xzz_zzzz_0, tpy_y_xxxx_0, tpy_y_xxxy_0, tpy_y_xxxz_0, \
                                     tpy_y_xxyy_0, tpy_y_xxyz_0, tpy_y_xxzz_0, tpy_y_xyyy_0, tpy_y_xyyz_0, tpy_y_xyzz_0, \
                                     tpy_y_xzzz_0, tpy_yy_xxx_0, tpy_yy_xxxx_0, tpy_yy_xxxy_0, tpy_yy_xxxz_0, \
                                     tpy_yy_xxy_0, tpy_yy_xxyy_0, tpy_yy_xxyz_0, tpy_yy_xxz_0, tpy_yy_xxzz_0, \
                                     tpy_yy_xyy_0, tpy_yy_xyyy_0, tpy_yy_xyyz_0, tpy_yy_xyz_0, tpy_yy_xyzz_0, \
                                     tpy_yy_xzz_0, tpy_yy_xzzz_0, tpy_yyy_xxxx_0, tpy_yyy_xxxy_0, tpy_yyy_xxxz_0, \
                                     tpy_yyy_xxyy_0, tpy_yyy_xxyz_0, tpy_yyy_xxzz_0, tpy_yyy_xyyy_0, tpy_yyy_xyyz_0, \
                                     tpy_yyy_xyzz_0, tpy_yyy_xzzz_0, tpy_zz_xyzz_0, tpy_zz_xzzz_0, tpy_zz_yyyy_0, \
                                     tpy_zz_yyyz_0, tpy_zz_yyzz_0, tpy_zz_yzz_0, tpy_zz_yzzz_0, tpy_zz_zzz_0, \
                                     tpy_zz_zzzz_0, tpz_xzz_xyzz_0, tpz_xzz_xzzz_0, tpz_xzz_yyyy_0, tpz_xzz_yyyz_0, \
                                     tpz_xzz_yyzz_0, tpz_xzz_yzzz_0, tpz_xzz_zzzz_0, tpz_y_xxxx_0, tpz_y_xxxy_0, \
                                     tpz_y_xxxz_0, tpz_y_xxyy_0, tpz_y_xxyz_0, tpz_y_xxzz_0, tpz_y_xyyy_0, tpz_y_xyyz_0, \
                                     tpz_y_xyzz_0, tpz_y_xzzz_0, tpz_yy_xxx_0, tpz_yy_xxxx_0, tpz_yy_xxxy_0, \
                                     tpz_yy_xxxz_0, tpz_yy_xxy_0, tpz_yy_xxyy_0, tpz_yy_xxyz_0, tpz_yy_xxz_0, \
                                     tpz_yy_xxzz_0, tpz_yy_xyy_0, tpz_yy_xyyy_0, tpz_yy_xyyz_0, tpz_yy_xyz_0, \
                                     tpz_yy_xyzz_0, tpz_yy_xzz_0, tpz_yy_xzzz_0, tpz_yyy_xxxx_0, tpz_yyy_xxxy_0, \
                                     tpz_yyy_xxxz_0, tpz_yyy_xxyy_0, tpz_yyy_xxyz_0, tpz_yyy_xxzz_0, tpz_yyy_xyyy_0, \
                                     tpz_yyy_xyyz_0, tpz_yyy_xyzz_0, tpz_yyy_xzzz_0, tpz_zz_xyzz_0, tpz_zz_xzzz_0, \
                                     tpz_zz_yyyy_0, tpz_zz_yyyz_0, tpz_zz_yyzz_0, tpz_zz_yzz_0, tpz_zz_yzzz_0, \
                                     tpz_zz_zzz_0, tpz_zz_zzzz_0, ts_yy_xxxx_0, ts_yy_xxxy_0, ts_yy_xxxz_0, \
                                     ts_yy_xxyy_0, ts_yy_xxyz_0, ts_yy_xxzz_0, ts_yy_xyyy_0, ts_yy_xyyz_0, ts_yy_xyzz_0, \
                                     ts_yy_xzzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, ts_zz_yyyz_0, ts_zz_yyzz_0, ts_zz_yzzz_0, \
                                     ts_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_xzz_xyzz_0[j] = pa_x[j] * tpy_zz_xyzz_0[j] + 0.5 * fl1_fx * tpy_zz_yzz_0[j];

            tpz_xzz_xyzz_0[j] = pa_x[j] * tpz_zz_xyzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzz_0[j];

            tpx_xzz_xzzz_0[j] = pa_x[j] * tpx_zz_xzzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xzzz_0[j];

            tpy_xzz_xzzz_0[j] = pa_x[j] * tpy_zz_xzzz_0[j] + 0.5 * fl1_fx * tpy_zz_zzz_0[j];

            tpz_xzz_xzzz_0[j] = pa_x[j] * tpz_zz_xzzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzz_0[j];

            tpx_xzz_yyyy_0[j] = pa_x[j] * tpx_zz_yyyy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyyy_0[j];

            tpy_xzz_yyyy_0[j] = pa_x[j] * tpy_zz_yyyy_0[j];

            tpz_xzz_yyyy_0[j] = pa_x[j] * tpz_zz_yyyy_0[j];

            tpx_xzz_yyyz_0[j] = pa_x[j] * tpx_zz_yyyz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyyz_0[j];

            tpy_xzz_yyyz_0[j] = pa_x[j] * tpy_zz_yyyz_0[j];

            tpz_xzz_yyyz_0[j] = pa_x[j] * tpz_zz_yyyz_0[j];

            tpx_xzz_yyzz_0[j] = pa_x[j] * tpx_zz_yyzz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyzz_0[j];

            tpy_xzz_yyzz_0[j] = pa_x[j] * tpy_zz_yyzz_0[j];

            tpz_xzz_yyzz_0[j] = pa_x[j] * tpz_zz_yyzz_0[j];

            tpx_xzz_yzzz_0[j] = pa_x[j] * tpx_zz_yzzz_0[j] - fl1_fgb * fl1_fx * ts_zz_yzzz_0[j];

            tpy_xzz_yzzz_0[j] = pa_x[j] * tpy_zz_yzzz_0[j];

            tpz_xzz_yzzz_0[j] = pa_x[j] * tpz_zz_yzzz_0[j];

            tpx_xzz_zzzz_0[j] = pa_x[j] * tpx_zz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_zz_zzzz_0[j];

            tpy_xzz_zzzz_0[j] = pa_x[j] * tpy_zz_zzzz_0[j];

            tpz_xzz_zzzz_0[j] = pa_x[j] * tpz_zz_zzzz_0[j];

            tpx_yyy_xxxx_0[j] = pa_y[j] * tpx_yy_xxxx_0[j] + fl1_fx * tpx_y_xxxx_0[j];

            tpy_yyy_xxxx_0[j] = pa_y[j] * tpy_yy_xxxx_0[j] + fl1_fx * tpy_y_xxxx_0[j] - fl1_fgb * fl1_fx * ts_yy_xxxx_0[j];

            tpz_yyy_xxxx_0[j] = pa_y[j] * tpz_yy_xxxx_0[j] + fl1_fx * tpz_y_xxxx_0[j];

            tpx_yyy_xxxy_0[j] = pa_y[j] * tpx_yy_xxxy_0[j] + fl1_fx * tpx_y_xxxy_0[j] + 0.5 * fl1_fx * tpx_yy_xxx_0[j];

            tpy_yyy_xxxy_0[j] =
                pa_y[j] * tpy_yy_xxxy_0[j] + fl1_fx * tpy_y_xxxy_0[j] + 0.5 * fl1_fx * tpy_yy_xxx_0[j] - fl1_fgb * fl1_fx * ts_yy_xxxy_0[j];

            tpz_yyy_xxxy_0[j] = pa_y[j] * tpz_yy_xxxy_0[j] + fl1_fx * tpz_y_xxxy_0[j] + 0.5 * fl1_fx * tpz_yy_xxx_0[j];

            tpx_yyy_xxxz_0[j] = pa_y[j] * tpx_yy_xxxz_0[j] + fl1_fx * tpx_y_xxxz_0[j];

            tpy_yyy_xxxz_0[j] = pa_y[j] * tpy_yy_xxxz_0[j] + fl1_fx * tpy_y_xxxz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxxz_0[j];

            tpz_yyy_xxxz_0[j] = pa_y[j] * tpz_yy_xxxz_0[j] + fl1_fx * tpz_y_xxxz_0[j];

            tpx_yyy_xxyy_0[j] = pa_y[j] * tpx_yy_xxyy_0[j] + fl1_fx * tpx_y_xxyy_0[j] + fl1_fx * tpx_yy_xxy_0[j];

            tpy_yyy_xxyy_0[j] = pa_y[j] * tpy_yy_xxyy_0[j] + fl1_fx * tpy_y_xxyy_0[j] + fl1_fx * tpy_yy_xxy_0[j] - fl1_fgb * fl1_fx * ts_yy_xxyy_0[j];

            tpz_yyy_xxyy_0[j] = pa_y[j] * tpz_yy_xxyy_0[j] + fl1_fx * tpz_y_xxyy_0[j] + fl1_fx * tpz_yy_xxy_0[j];

            tpx_yyy_xxyz_0[j] = pa_y[j] * tpx_yy_xxyz_0[j] + fl1_fx * tpx_y_xxyz_0[j] + 0.5 * fl1_fx * tpx_yy_xxz_0[j];

            tpy_yyy_xxyz_0[j] =
                pa_y[j] * tpy_yy_xxyz_0[j] + fl1_fx * tpy_y_xxyz_0[j] + 0.5 * fl1_fx * tpy_yy_xxz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxyz_0[j];

            tpz_yyy_xxyz_0[j] = pa_y[j] * tpz_yy_xxyz_0[j] + fl1_fx * tpz_y_xxyz_0[j] + 0.5 * fl1_fx * tpz_yy_xxz_0[j];

            tpx_yyy_xxzz_0[j] = pa_y[j] * tpx_yy_xxzz_0[j] + fl1_fx * tpx_y_xxzz_0[j];

            tpy_yyy_xxzz_0[j] = pa_y[j] * tpy_yy_xxzz_0[j] + fl1_fx * tpy_y_xxzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxzz_0[j];

            tpz_yyy_xxzz_0[j] = pa_y[j] * tpz_yy_xxzz_0[j] + fl1_fx * tpz_y_xxzz_0[j];

            tpx_yyy_xyyy_0[j] = pa_y[j] * tpx_yy_xyyy_0[j] + fl1_fx * tpx_y_xyyy_0[j] + 1.5 * fl1_fx * tpx_yy_xyy_0[j];

            tpy_yyy_xyyy_0[j] =
                pa_y[j] * tpy_yy_xyyy_0[j] + fl1_fx * tpy_y_xyyy_0[j] + 1.5 * fl1_fx * tpy_yy_xyy_0[j] - fl1_fgb * fl1_fx * ts_yy_xyyy_0[j];

            tpz_yyy_xyyy_0[j] = pa_y[j] * tpz_yy_xyyy_0[j] + fl1_fx * tpz_y_xyyy_0[j] + 1.5 * fl1_fx * tpz_yy_xyy_0[j];

            tpx_yyy_xyyz_0[j] = pa_y[j] * tpx_yy_xyyz_0[j] + fl1_fx * tpx_y_xyyz_0[j] + fl1_fx * tpx_yy_xyz_0[j];

            tpy_yyy_xyyz_0[j] = pa_y[j] * tpy_yy_xyyz_0[j] + fl1_fx * tpy_y_xyyz_0[j] + fl1_fx * tpy_yy_xyz_0[j] - fl1_fgb * fl1_fx * ts_yy_xyyz_0[j];

            tpz_yyy_xyyz_0[j] = pa_y[j] * tpz_yy_xyyz_0[j] + fl1_fx * tpz_y_xyyz_0[j] + fl1_fx * tpz_yy_xyz_0[j];

            tpx_yyy_xyzz_0[j] = pa_y[j] * tpx_yy_xyzz_0[j] + fl1_fx * tpx_y_xyzz_0[j] + 0.5 * fl1_fx * tpx_yy_xzz_0[j];

            tpy_yyy_xyzz_0[j] =
                pa_y[j] * tpy_yy_xyzz_0[j] + fl1_fx * tpy_y_xyzz_0[j] + 0.5 * fl1_fx * tpy_yy_xzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xyzz_0[j];

            tpz_yyy_xyzz_0[j] = pa_y[j] * tpz_yy_xyzz_0[j] + fl1_fx * tpz_y_xyzz_0[j] + 0.5 * fl1_fx * tpz_yy_xzz_0[j];

            tpx_yyy_xzzz_0[j] = pa_y[j] * tpx_yy_xzzz_0[j] + fl1_fx * tpx_y_xzzz_0[j];

            tpy_yyy_xzzz_0[j] = pa_y[j] * tpy_yy_xzzz_0[j] + fl1_fx * tpy_y_xzzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xzzz_0[j];

            tpz_yyy_xzzz_0[j] = pa_y[j] * tpz_yy_xzzz_0[j] + fl1_fx * tpz_y_xzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_300_350(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (300,350)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 25);

        auto tpy_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tpx_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 26);

        auto tpy_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tpx_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 27);

        auto tpy_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tpx_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 28);

        auto tpy_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tpx_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 29);

        auto tpy_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29);

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpx_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 36);

        auto tpy_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 37);

        auto tpy_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yy_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 38);

        auto tpy_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yy_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 39);

        auto tpy_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yy_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 40);

        auto tpy_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 41);

        auto tpy_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 42);

        auto tpy_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 43);

        auto tpy_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 44);

        auto tpy_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 45);

        auto tpy_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 46);

        auto tpy_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 47);

        auto tpy_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 47);

        auto ts_yy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 55);

        auto ts_yy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 56);

        auto ts_yy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 57);

        auto ts_yy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 58);

        auto ts_yy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 59);

        auto ts_yz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 60);

        auto ts_yz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 61);

        auto ts_yz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 62);

        auto ts_yz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 63);

        auto ts_yz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 64);

        auto ts_yz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 65);

        auto ts_yz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 66);

        auto ts_yz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 67);

        auto ts_yz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 68);

        auto ts_yz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 69);

        auto ts_yz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 70);

        auto ts_yz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 71);

        // set up pointers to integrals

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

        // Batch of Integrals (300,350)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_y_yyyy_0, tpx_y_yyyz_0, tpx_y_yyzz_0, tpx_y_yzzz_0, \
                                     tpx_y_zzzz_0, tpx_yy_yyy_0, tpx_yy_yyyy_0, tpx_yy_yyyz_0, tpx_yy_yyz_0, \
                                     tpx_yy_yyzz_0, tpx_yy_yzz_0, tpx_yy_yzzz_0, tpx_yy_zzz_0, tpx_yy_zzzz_0, \
                                     tpx_yyy_yyyy_0, tpx_yyy_yyyz_0, tpx_yyy_yyzz_0, tpx_yyy_yzzz_0, tpx_yyy_zzzz_0, \
                                     tpx_yyz_xxxx_0, tpx_yyz_xxxy_0, tpx_yyz_xxxz_0, tpx_yyz_xxyy_0, tpx_yyz_xxyz_0, \
                                     tpx_yyz_xxzz_0, tpx_yyz_xyyy_0, tpx_yyz_xyyz_0, tpx_yyz_xyzz_0, tpx_yyz_xzzz_0, \
                                     tpx_yyz_yyyy_0, tpx_yyz_yyyz_0, tpx_yz_xxx_0, tpx_yz_xxxx_0, tpx_yz_xxxy_0, \
                                     tpx_yz_xxxz_0, tpx_yz_xxy_0, tpx_yz_xxyy_0, tpx_yz_xxyz_0, tpx_yz_xxz_0, \
                                     tpx_yz_xxzz_0, tpx_yz_xyy_0, tpx_yz_xyyy_0, tpx_yz_xyyz_0, tpx_yz_xyz_0, \
                                     tpx_yz_xyzz_0, tpx_yz_xzz_0, tpx_yz_xzzz_0, tpx_yz_yyy_0, tpx_yz_yyyy_0, \
                                     tpx_yz_yyyz_0, tpx_yz_yyz_0, tpx_z_xxxx_0, tpx_z_xxxy_0, tpx_z_xxxz_0, tpx_z_xxyy_0, \
                                     tpx_z_xxyz_0, tpx_z_xxzz_0, tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyzz_0, tpx_z_xzzz_0, \
                                     tpx_z_yyyy_0, tpx_z_yyyz_0, tpy_y_yyyy_0, tpy_y_yyyz_0, tpy_y_yyzz_0, tpy_y_yzzz_0, \
                                     tpy_y_zzzz_0, tpy_yy_yyy_0, tpy_yy_yyyy_0, tpy_yy_yyyz_0, tpy_yy_yyz_0, \
                                     tpy_yy_yyzz_0, tpy_yy_yzz_0, tpy_yy_yzzz_0, tpy_yy_zzz_0, tpy_yy_zzzz_0, \
                                     tpy_yyy_yyyy_0, tpy_yyy_yyyz_0, tpy_yyy_yyzz_0, tpy_yyy_yzzz_0, tpy_yyy_zzzz_0, \
                                     tpy_yyz_xxxx_0, tpy_yyz_xxxy_0, tpy_yyz_xxxz_0, tpy_yyz_xxyy_0, tpy_yyz_xxyz_0, \
                                     tpy_yyz_xxzz_0, tpy_yyz_xyyy_0, tpy_yyz_xyyz_0, tpy_yyz_xyzz_0, tpy_yyz_xzzz_0, \
                                     tpy_yyz_yyyy_0, tpy_yyz_yyyz_0, tpy_yz_xxx_0, tpy_yz_xxxx_0, tpy_yz_xxxy_0, \
                                     tpy_yz_xxxz_0, tpy_yz_xxy_0, tpy_yz_xxyy_0, tpy_yz_xxyz_0, tpy_yz_xxz_0, \
                                     tpy_yz_xxzz_0, tpy_yz_xyy_0, tpy_yz_xyyy_0, tpy_yz_xyyz_0, tpy_yz_xyz_0, \
                                     tpy_yz_xyzz_0, tpy_yz_xzz_0, tpy_yz_xzzz_0, tpy_yz_yyy_0, tpy_yz_yyyy_0, \
                                     tpy_yz_yyyz_0, tpy_yz_yyz_0, tpy_z_xxxx_0, tpy_z_xxxy_0, tpy_z_xxxz_0, tpy_z_xxyy_0, \
                                     tpy_z_xxyz_0, tpy_z_xxzz_0, tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyzz_0, tpy_z_xzzz_0, \
                                     tpy_z_yyyy_0, tpy_z_yyyz_0, tpz_y_yyyy_0, tpz_y_yyyz_0, tpz_y_yyzz_0, tpz_y_yzzz_0, \
                                     tpz_y_zzzz_0, tpz_yy_yyy_0, tpz_yy_yyyy_0, tpz_yy_yyyz_0, tpz_yy_yyz_0, \
                                     tpz_yy_yyzz_0, tpz_yy_yzz_0, tpz_yy_yzzz_0, tpz_yy_zzz_0, tpz_yy_zzzz_0, \
                                     tpz_yyy_yyyy_0, tpz_yyy_yyyz_0, tpz_yyy_yyzz_0, tpz_yyy_yzzz_0, tpz_yyy_zzzz_0, \
                                     tpz_yyz_xxxx_0, tpz_yyz_xxxy_0, tpz_yyz_xxxz_0, tpz_yyz_xxyy_0, tpz_yyz_xxyz_0, \
                                     tpz_yyz_xxzz_0, tpz_yyz_xyyy_0, tpz_yyz_xyyz_0, tpz_yyz_xyzz_0, tpz_yyz_xzzz_0, \
                                     tpz_yyz_yyyy_0, tpz_yz_xxx_0, tpz_yz_xxxx_0, tpz_yz_xxxy_0, tpz_yz_xxxz_0, \
                                     tpz_yz_xxy_0, tpz_yz_xxyy_0, tpz_yz_xxyz_0, tpz_yz_xxz_0, tpz_yz_xxzz_0, \
                                     tpz_yz_xyy_0, tpz_yz_xyyy_0, tpz_yz_xyyz_0, tpz_yz_xyz_0, tpz_yz_xyzz_0, \
                                     tpz_yz_xzz_0, tpz_yz_xzzz_0, tpz_yz_yyy_0, tpz_yz_yyyy_0, tpz_z_xxxx_0, \
                                     tpz_z_xxxy_0, tpz_z_xxxz_0, tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxzz_0, tpz_z_xyyy_0, \
                                     tpz_z_xyyz_0, tpz_z_xyzz_0, tpz_z_xzzz_0, tpz_z_yyyy_0, ts_yy_yyyy_0, ts_yy_yyyz_0, \
                                     ts_yy_yyzz_0, ts_yy_yzzz_0, ts_yy_zzzz_0, ts_yz_xxxx_0, ts_yz_xxxy_0, ts_yz_xxxz_0, \
                                     ts_yz_xxyy_0, ts_yz_xxyz_0, ts_yz_xxzz_0, ts_yz_xyyy_0, ts_yz_xyyz_0, ts_yz_xyzz_0, \
                                     ts_yz_xzzz_0, ts_yz_yyyy_0, ts_yz_yyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyy_yyyy_0[j] = pa_y[j] * tpx_yy_yyyy_0[j] + fl1_fx * tpx_y_yyyy_0[j] + 2.0 * fl1_fx * tpx_yy_yyy_0[j];

            tpy_yyy_yyyy_0[j] =
                pa_y[j] * tpy_yy_yyyy_0[j] + fl1_fx * tpy_y_yyyy_0[j] + 2.0 * fl1_fx * tpy_yy_yyy_0[j] - fl1_fgb * fl1_fx * ts_yy_yyyy_0[j];

            tpz_yyy_yyyy_0[j] = pa_y[j] * tpz_yy_yyyy_0[j] + fl1_fx * tpz_y_yyyy_0[j] + 2.0 * fl1_fx * tpz_yy_yyy_0[j];

            tpx_yyy_yyyz_0[j] = pa_y[j] * tpx_yy_yyyz_0[j] + fl1_fx * tpx_y_yyyz_0[j] + 1.5 * fl1_fx * tpx_yy_yyz_0[j];

            tpy_yyy_yyyz_0[j] =
                pa_y[j] * tpy_yy_yyyz_0[j] + fl1_fx * tpy_y_yyyz_0[j] + 1.5 * fl1_fx * tpy_yy_yyz_0[j] - fl1_fgb * fl1_fx * ts_yy_yyyz_0[j];

            tpz_yyy_yyyz_0[j] = pa_y[j] * tpz_yy_yyyz_0[j] + fl1_fx * tpz_y_yyyz_0[j] + 1.5 * fl1_fx * tpz_yy_yyz_0[j];

            tpx_yyy_yyzz_0[j] = pa_y[j] * tpx_yy_yyzz_0[j] + fl1_fx * tpx_y_yyzz_0[j] + fl1_fx * tpx_yy_yzz_0[j];

            tpy_yyy_yyzz_0[j] = pa_y[j] * tpy_yy_yyzz_0[j] + fl1_fx * tpy_y_yyzz_0[j] + fl1_fx * tpy_yy_yzz_0[j] - fl1_fgb * fl1_fx * ts_yy_yyzz_0[j];

            tpz_yyy_yyzz_0[j] = pa_y[j] * tpz_yy_yyzz_0[j] + fl1_fx * tpz_y_yyzz_0[j] + fl1_fx * tpz_yy_yzz_0[j];

            tpx_yyy_yzzz_0[j] = pa_y[j] * tpx_yy_yzzz_0[j] + fl1_fx * tpx_y_yzzz_0[j] + 0.5 * fl1_fx * tpx_yy_zzz_0[j];

            tpy_yyy_yzzz_0[j] =
                pa_y[j] * tpy_yy_yzzz_0[j] + fl1_fx * tpy_y_yzzz_0[j] + 0.5 * fl1_fx * tpy_yy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yy_yzzz_0[j];

            tpz_yyy_yzzz_0[j] = pa_y[j] * tpz_yy_yzzz_0[j] + fl1_fx * tpz_y_yzzz_0[j] + 0.5 * fl1_fx * tpz_yy_zzz_0[j];

            tpx_yyy_zzzz_0[j] = pa_y[j] * tpx_yy_zzzz_0[j] + fl1_fx * tpx_y_zzzz_0[j];

            tpy_yyy_zzzz_0[j] = pa_y[j] * tpy_yy_zzzz_0[j] + fl1_fx * tpy_y_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yy_zzzz_0[j];

            tpz_yyy_zzzz_0[j] = pa_y[j] * tpz_yy_zzzz_0[j] + fl1_fx * tpz_y_zzzz_0[j];

            tpx_yyz_xxxx_0[j] = pa_y[j] * tpx_yz_xxxx_0[j] + 0.5 * fl1_fx * tpx_z_xxxx_0[j];

            tpy_yyz_xxxx_0[j] = pa_y[j] * tpy_yz_xxxx_0[j] + 0.5 * fl1_fx * tpy_z_xxxx_0[j] - fl1_fgb * fl1_fx * ts_yz_xxxx_0[j];

            tpz_yyz_xxxx_0[j] = pa_y[j] * tpz_yz_xxxx_0[j] + 0.5 * fl1_fx * tpz_z_xxxx_0[j];

            tpx_yyz_xxxy_0[j] = pa_y[j] * tpx_yz_xxxy_0[j] + 0.5 * fl1_fx * tpx_z_xxxy_0[j] + 0.5 * fl1_fx * tpx_yz_xxx_0[j];

            tpy_yyz_xxxy_0[j] =
                pa_y[j] * tpy_yz_xxxy_0[j] + 0.5 * fl1_fx * tpy_z_xxxy_0[j] + 0.5 * fl1_fx * tpy_yz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yz_xxxy_0[j];

            tpz_yyz_xxxy_0[j] = pa_y[j] * tpz_yz_xxxy_0[j] + 0.5 * fl1_fx * tpz_z_xxxy_0[j] + 0.5 * fl1_fx * tpz_yz_xxx_0[j];

            tpx_yyz_xxxz_0[j] = pa_y[j] * tpx_yz_xxxz_0[j] + 0.5 * fl1_fx * tpx_z_xxxz_0[j];

            tpy_yyz_xxxz_0[j] = pa_y[j] * tpy_yz_xxxz_0[j] + 0.5 * fl1_fx * tpy_z_xxxz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxxz_0[j];

            tpz_yyz_xxxz_0[j] = pa_y[j] * tpz_yz_xxxz_0[j] + 0.5 * fl1_fx * tpz_z_xxxz_0[j];

            tpx_yyz_xxyy_0[j] = pa_y[j] * tpx_yz_xxyy_0[j] + 0.5 * fl1_fx * tpx_z_xxyy_0[j] + fl1_fx * tpx_yz_xxy_0[j];

            tpy_yyz_xxyy_0[j] =
                pa_y[j] * tpy_yz_xxyy_0[j] + 0.5 * fl1_fx * tpy_z_xxyy_0[j] + fl1_fx * tpy_yz_xxy_0[j] - fl1_fgb * fl1_fx * ts_yz_xxyy_0[j];

            tpz_yyz_xxyy_0[j] = pa_y[j] * tpz_yz_xxyy_0[j] + 0.5 * fl1_fx * tpz_z_xxyy_0[j] + fl1_fx * tpz_yz_xxy_0[j];

            tpx_yyz_xxyz_0[j] = pa_y[j] * tpx_yz_xxyz_0[j] + 0.5 * fl1_fx * tpx_z_xxyz_0[j] + 0.5 * fl1_fx * tpx_yz_xxz_0[j];

            tpy_yyz_xxyz_0[j] =
                pa_y[j] * tpy_yz_xxyz_0[j] + 0.5 * fl1_fx * tpy_z_xxyz_0[j] + 0.5 * fl1_fx * tpy_yz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxyz_0[j];

            tpz_yyz_xxyz_0[j] = pa_y[j] * tpz_yz_xxyz_0[j] + 0.5 * fl1_fx * tpz_z_xxyz_0[j] + 0.5 * fl1_fx * tpz_yz_xxz_0[j];

            tpx_yyz_xxzz_0[j] = pa_y[j] * tpx_yz_xxzz_0[j] + 0.5 * fl1_fx * tpx_z_xxzz_0[j];

            tpy_yyz_xxzz_0[j] = pa_y[j] * tpy_yz_xxzz_0[j] + 0.5 * fl1_fx * tpy_z_xxzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxzz_0[j];

            tpz_yyz_xxzz_0[j] = pa_y[j] * tpz_yz_xxzz_0[j] + 0.5 * fl1_fx * tpz_z_xxzz_0[j];

            tpx_yyz_xyyy_0[j] = pa_y[j] * tpx_yz_xyyy_0[j] + 0.5 * fl1_fx * tpx_z_xyyy_0[j] + 1.5 * fl1_fx * tpx_yz_xyy_0[j];

            tpy_yyz_xyyy_0[j] =
                pa_y[j] * tpy_yz_xyyy_0[j] + 0.5 * fl1_fx * tpy_z_xyyy_0[j] + 1.5 * fl1_fx * tpy_yz_xyy_0[j] - fl1_fgb * fl1_fx * ts_yz_xyyy_0[j];

            tpz_yyz_xyyy_0[j] = pa_y[j] * tpz_yz_xyyy_0[j] + 0.5 * fl1_fx * tpz_z_xyyy_0[j] + 1.5 * fl1_fx * tpz_yz_xyy_0[j];

            tpx_yyz_xyyz_0[j] = pa_y[j] * tpx_yz_xyyz_0[j] + 0.5 * fl1_fx * tpx_z_xyyz_0[j] + fl1_fx * tpx_yz_xyz_0[j];

            tpy_yyz_xyyz_0[j] =
                pa_y[j] * tpy_yz_xyyz_0[j] + 0.5 * fl1_fx * tpy_z_xyyz_0[j] + fl1_fx * tpy_yz_xyz_0[j] - fl1_fgb * fl1_fx * ts_yz_xyyz_0[j];

            tpz_yyz_xyyz_0[j] = pa_y[j] * tpz_yz_xyyz_0[j] + 0.5 * fl1_fx * tpz_z_xyyz_0[j] + fl1_fx * tpz_yz_xyz_0[j];

            tpx_yyz_xyzz_0[j] = pa_y[j] * tpx_yz_xyzz_0[j] + 0.5 * fl1_fx * tpx_z_xyzz_0[j] + 0.5 * fl1_fx * tpx_yz_xzz_0[j];

            tpy_yyz_xyzz_0[j] =
                pa_y[j] * tpy_yz_xyzz_0[j] + 0.5 * fl1_fx * tpy_z_xyzz_0[j] + 0.5 * fl1_fx * tpy_yz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xyzz_0[j];

            tpz_yyz_xyzz_0[j] = pa_y[j] * tpz_yz_xyzz_0[j] + 0.5 * fl1_fx * tpz_z_xyzz_0[j] + 0.5 * fl1_fx * tpz_yz_xzz_0[j];

            tpx_yyz_xzzz_0[j] = pa_y[j] * tpx_yz_xzzz_0[j] + 0.5 * fl1_fx * tpx_z_xzzz_0[j];

            tpy_yyz_xzzz_0[j] = pa_y[j] * tpy_yz_xzzz_0[j] + 0.5 * fl1_fx * tpy_z_xzzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xzzz_0[j];

            tpz_yyz_xzzz_0[j] = pa_y[j] * tpz_yz_xzzz_0[j] + 0.5 * fl1_fx * tpz_z_xzzz_0[j];

            tpx_yyz_yyyy_0[j] = pa_y[j] * tpx_yz_yyyy_0[j] + 0.5 * fl1_fx * tpx_z_yyyy_0[j] + 2.0 * fl1_fx * tpx_yz_yyy_0[j];

            tpy_yyz_yyyy_0[j] =
                pa_y[j] * tpy_yz_yyyy_0[j] + 0.5 * fl1_fx * tpy_z_yyyy_0[j] + 2.0 * fl1_fx * tpy_yz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yz_yyyy_0[j];

            tpz_yyz_yyyy_0[j] = pa_y[j] * tpz_yz_yyyy_0[j] + 0.5 * fl1_fx * tpz_z_yyyy_0[j] + 2.0 * fl1_fx * tpz_yz_yyy_0[j];

            tpx_yyz_yyyz_0[j] = pa_y[j] * tpx_yz_yyyz_0[j] + 0.5 * fl1_fx * tpx_z_yyyz_0[j] + 1.5 * fl1_fx * tpx_yz_yyz_0[j];

            tpy_yyz_yyyz_0[j] =
                pa_y[j] * tpy_yz_yyyz_0[j] + 0.5 * fl1_fx * tpy_z_yyyz_0[j] + 1.5 * fl1_fx * tpy_yz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yz_yyyz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_350_400(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (350,400)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tpz_yz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 48);

        auto tpy_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 49);

        auto tpy_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 50);

        auto tpy_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 51);

        auto tpy_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 52);

        auto tpy_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 53);

        auto tpy_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 54);

        auto tpy_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 55);

        auto tpy_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 56);

        auto tpy_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 57);

        auto tpy_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 58);

        auto tpy_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 59);

        auto ts_yz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 72);

        auto ts_yz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 73);

        auto ts_yz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 74);

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        // set up pointers to integrals

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

        // Batch of Integrals (350,400)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yyz_yyzz_0, tpx_yyz_yzzz_0, tpx_yyz_zzzz_0, \
                                     tpx_yz_yyzz_0, tpx_yz_yzz_0, tpx_yz_yzzz_0, tpx_yz_zzz_0, tpx_yz_zzzz_0, \
                                     tpx_yzz_xxxx_0, tpx_yzz_xxxy_0, tpx_yzz_xxxz_0, tpx_yzz_xxyy_0, tpx_yzz_xxyz_0, \
                                     tpx_yzz_xxzz_0, tpx_yzz_xyyy_0, tpx_yzz_xyyz_0, tpx_yzz_xyzz_0, tpx_yzz_xzzz_0, \
                                     tpx_yzz_yyyy_0, tpx_yzz_yyyz_0, tpx_yzz_yyzz_0, tpx_yzz_yzzz_0, tpx_z_yyzz_0, \
                                     tpx_z_yzzz_0, tpx_z_zzzz_0, tpx_zz_xxx_0, tpx_zz_xxxx_0, tpx_zz_xxxy_0, \
                                     tpx_zz_xxxz_0, tpx_zz_xxy_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, tpx_zz_xxz_0, \
                                     tpx_zz_xxzz_0, tpx_zz_xyy_0, tpx_zz_xyyy_0, tpx_zz_xyyz_0, tpx_zz_xyz_0, \
                                     tpx_zz_xyzz_0, tpx_zz_xzz_0, tpx_zz_xzzz_0, tpx_zz_yyy_0, tpx_zz_yyyy_0, \
                                     tpx_zz_yyyz_0, tpx_zz_yyz_0, tpx_zz_yyzz_0, tpx_zz_yzz_0, tpx_zz_yzzz_0, \
                                     tpx_zz_zzz_0, tpy_yyz_yyzz_0, tpy_yyz_yzzz_0, tpy_yyz_zzzz_0, tpy_yz_yyzz_0, \
                                     tpy_yz_yzz_0, tpy_yz_yzzz_0, tpy_yz_zzz_0, tpy_yz_zzzz_0, tpy_yzz_xxxx_0, \
                                     tpy_yzz_xxxy_0, tpy_yzz_xxxz_0, tpy_yzz_xxyy_0, tpy_yzz_xxyz_0, tpy_yzz_xxzz_0, \
                                     tpy_yzz_xyyy_0, tpy_yzz_xyyz_0, tpy_yzz_xyzz_0, tpy_yzz_xzzz_0, tpy_yzz_yyyy_0, \
                                     tpy_yzz_yyyz_0, tpy_yzz_yyzz_0, tpy_z_yyzz_0, tpy_z_yzzz_0, tpy_z_zzzz_0, \
                                     tpy_zz_xxx_0, tpy_zz_xxxx_0, tpy_zz_xxxy_0, tpy_zz_xxxz_0, tpy_zz_xxy_0, \
                                     tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxz_0, tpy_zz_xxzz_0, tpy_zz_xyy_0, \
                                     tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpy_zz_xyz_0, tpy_zz_xyzz_0, tpy_zz_xzz_0, \
                                     tpy_zz_xzzz_0, tpy_zz_yyy_0, tpy_zz_yyyy_0, tpy_zz_yyyz_0, tpy_zz_yyz_0, \
                                     tpy_zz_yyzz_0, tpy_zz_yzz_0, tpz_yyz_yyyz_0, tpz_yyz_yyzz_0, tpz_yyz_yzzz_0, \
                                     tpz_yyz_zzzz_0, tpz_yz_yyyz_0, tpz_yz_yyz_0, tpz_yz_yyzz_0, tpz_yz_yzz_0, \
                                     tpz_yz_yzzz_0, tpz_yz_zzz_0, tpz_yz_zzzz_0, tpz_yzz_xxxx_0, tpz_yzz_xxxy_0, \
                                     tpz_yzz_xxxz_0, tpz_yzz_xxyy_0, tpz_yzz_xxyz_0, tpz_yzz_xxzz_0, tpz_yzz_xyyy_0, \
                                     tpz_yzz_xyyz_0, tpz_yzz_xyzz_0, tpz_yzz_xzzz_0, tpz_yzz_yyyy_0, tpz_yzz_yyyz_0, \
                                     tpz_yzz_yyzz_0, tpz_z_yyyz_0, tpz_z_yyzz_0, tpz_z_yzzz_0, tpz_z_zzzz_0, tpz_zz_xxx_0, \
                                     tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, tpz_zz_xxy_0, tpz_zz_xxyy_0, \
                                     tpz_zz_xxyz_0, tpz_zz_xxz_0, tpz_zz_xxzz_0, tpz_zz_xyy_0, tpz_zz_xyyy_0, \
                                     tpz_zz_xyyz_0, tpz_zz_xyz_0, tpz_zz_xyzz_0, tpz_zz_xzz_0, tpz_zz_xzzz_0, \
                                     tpz_zz_yyy_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, tpz_zz_yyz_0, tpz_zz_yyzz_0, \
                                     tpz_zz_yzz_0, ts_yz_yyzz_0, ts_yz_yzzz_0, ts_yz_zzzz_0, ts_zz_xxxx_0, ts_zz_xxxy_0, \
                                     ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, ts_zz_xxzz_0, ts_zz_xyyy_0, ts_zz_xyyz_0, \
                                     ts_zz_xyzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, ts_zz_yyyz_0, ts_zz_yyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_yyz_yyyz_0[j] = pa_y[j] * tpz_yz_yyyz_0[j] + 0.5 * fl1_fx * tpz_z_yyyz_0[j] + 1.5 * fl1_fx * tpz_yz_yyz_0[j];

            tpx_yyz_yyzz_0[j] = pa_y[j] * tpx_yz_yyzz_0[j] + 0.5 * fl1_fx * tpx_z_yyzz_0[j] + fl1_fx * tpx_yz_yzz_0[j];

            tpy_yyz_yyzz_0[j] =
                pa_y[j] * tpy_yz_yyzz_0[j] + 0.5 * fl1_fx * tpy_z_yyzz_0[j] + fl1_fx * tpy_yz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yz_yyzz_0[j];

            tpz_yyz_yyzz_0[j] = pa_y[j] * tpz_yz_yyzz_0[j] + 0.5 * fl1_fx * tpz_z_yyzz_0[j] + fl1_fx * tpz_yz_yzz_0[j];

            tpx_yyz_yzzz_0[j] = pa_y[j] * tpx_yz_yzzz_0[j] + 0.5 * fl1_fx * tpx_z_yzzz_0[j] + 0.5 * fl1_fx * tpx_yz_zzz_0[j];

            tpy_yyz_yzzz_0[j] =
                pa_y[j] * tpy_yz_yzzz_0[j] + 0.5 * fl1_fx * tpy_z_yzzz_0[j] + 0.5 * fl1_fx * tpy_yz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yz_yzzz_0[j];

            tpz_yyz_yzzz_0[j] = pa_y[j] * tpz_yz_yzzz_0[j] + 0.5 * fl1_fx * tpz_z_yzzz_0[j] + 0.5 * fl1_fx * tpz_yz_zzz_0[j];

            tpx_yyz_zzzz_0[j] = pa_y[j] * tpx_yz_zzzz_0[j] + 0.5 * fl1_fx * tpx_z_zzzz_0[j];

            tpy_yyz_zzzz_0[j] = pa_y[j] * tpy_yz_zzzz_0[j] + 0.5 * fl1_fx * tpy_z_zzzz_0[j] - fl1_fgb * fl1_fx * ts_yz_zzzz_0[j];

            tpz_yyz_zzzz_0[j] = pa_y[j] * tpz_yz_zzzz_0[j] + 0.5 * fl1_fx * tpz_z_zzzz_0[j];

            tpx_yzz_xxxx_0[j] = pa_y[j] * tpx_zz_xxxx_0[j];

            tpy_yzz_xxxx_0[j] = pa_y[j] * tpy_zz_xxxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxx_0[j];

            tpz_yzz_xxxx_0[j] = pa_y[j] * tpz_zz_xxxx_0[j];

            tpx_yzz_xxxy_0[j] = pa_y[j] * tpx_zz_xxxy_0[j] + 0.5 * fl1_fx * tpx_zz_xxx_0[j];

            tpy_yzz_xxxy_0[j] = pa_y[j] * tpy_zz_xxxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxy_0[j];

            tpz_yzz_xxxy_0[j] = pa_y[j] * tpz_zz_xxxy_0[j] + 0.5 * fl1_fx * tpz_zz_xxx_0[j];

            tpx_yzz_xxxz_0[j] = pa_y[j] * tpx_zz_xxxz_0[j];

            tpy_yzz_xxxz_0[j] = pa_y[j] * tpy_zz_xxxz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxz_0[j];

            tpz_yzz_xxxz_0[j] = pa_y[j] * tpz_zz_xxxz_0[j];

            tpx_yzz_xxyy_0[j] = pa_y[j] * tpx_zz_xxyy_0[j] + fl1_fx * tpx_zz_xxy_0[j];

            tpy_yzz_xxyy_0[j] = pa_y[j] * tpy_zz_xxyy_0[j] + fl1_fx * tpy_zz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxyy_0[j];

            tpz_yzz_xxyy_0[j] = pa_y[j] * tpz_zz_xxyy_0[j] + fl1_fx * tpz_zz_xxy_0[j];

            tpx_yzz_xxyz_0[j] = pa_y[j] * tpx_zz_xxyz_0[j] + 0.5 * fl1_fx * tpx_zz_xxz_0[j];

            tpy_yzz_xxyz_0[j] = pa_y[j] * tpy_zz_xxyz_0[j] + 0.5 * fl1_fx * tpy_zz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxyz_0[j];

            tpz_yzz_xxyz_0[j] = pa_y[j] * tpz_zz_xxyz_0[j] + 0.5 * fl1_fx * tpz_zz_xxz_0[j];

            tpx_yzz_xxzz_0[j] = pa_y[j] * tpx_zz_xxzz_0[j];

            tpy_yzz_xxzz_0[j] = pa_y[j] * tpy_zz_xxzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxzz_0[j];

            tpz_yzz_xxzz_0[j] = pa_y[j] * tpz_zz_xxzz_0[j];

            tpx_yzz_xyyy_0[j] = pa_y[j] * tpx_zz_xyyy_0[j] + 1.5 * fl1_fx * tpx_zz_xyy_0[j];

            tpy_yzz_xyyy_0[j] = pa_y[j] * tpy_zz_xyyy_0[j] + 1.5 * fl1_fx * tpy_zz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyyy_0[j];

            tpz_yzz_xyyy_0[j] = pa_y[j] * tpz_zz_xyyy_0[j] + 1.5 * fl1_fx * tpz_zz_xyy_0[j];

            tpx_yzz_xyyz_0[j] = pa_y[j] * tpx_zz_xyyz_0[j] + fl1_fx * tpx_zz_xyz_0[j];

            tpy_yzz_xyyz_0[j] = pa_y[j] * tpy_zz_xyyz_0[j] + fl1_fx * tpy_zz_xyz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyyz_0[j];

            tpz_yzz_xyyz_0[j] = pa_y[j] * tpz_zz_xyyz_0[j] + fl1_fx * tpz_zz_xyz_0[j];

            tpx_yzz_xyzz_0[j] = pa_y[j] * tpx_zz_xyzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzz_0[j];

            tpy_yzz_xyzz_0[j] = pa_y[j] * tpy_zz_xyzz_0[j] + 0.5 * fl1_fx * tpy_zz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyzz_0[j];

            tpz_yzz_xyzz_0[j] = pa_y[j] * tpz_zz_xyzz_0[j] + 0.5 * fl1_fx * tpz_zz_xzz_0[j];

            tpx_yzz_xzzz_0[j] = pa_y[j] * tpx_zz_xzzz_0[j];

            tpy_yzz_xzzz_0[j] = pa_y[j] * tpy_zz_xzzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xzzz_0[j];

            tpz_yzz_xzzz_0[j] = pa_y[j] * tpz_zz_xzzz_0[j];

            tpx_yzz_yyyy_0[j] = pa_y[j] * tpx_zz_yyyy_0[j] + 2.0 * fl1_fx * tpx_zz_yyy_0[j];

            tpy_yzz_yyyy_0[j] = pa_y[j] * tpy_zz_yyyy_0[j] + 2.0 * fl1_fx * tpy_zz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyyy_0[j];

            tpz_yzz_yyyy_0[j] = pa_y[j] * tpz_zz_yyyy_0[j] + 2.0 * fl1_fx * tpz_zz_yyy_0[j];

            tpx_yzz_yyyz_0[j] = pa_y[j] * tpx_zz_yyyz_0[j] + 1.5 * fl1_fx * tpx_zz_yyz_0[j];

            tpy_yzz_yyyz_0[j] = pa_y[j] * tpy_zz_yyyz_0[j] + 1.5 * fl1_fx * tpy_zz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyyz_0[j];

            tpz_yzz_yyyz_0[j] = pa_y[j] * tpz_zz_yyyz_0[j] + 1.5 * fl1_fx * tpz_zz_yyz_0[j];

            tpx_yzz_yyzz_0[j] = pa_y[j] * tpx_zz_yyzz_0[j] + fl1_fx * tpx_zz_yzz_0[j];

            tpy_yzz_yyzz_0[j] = pa_y[j] * tpy_zz_yyzz_0[j] + fl1_fx * tpy_zz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyzz_0[j];

            tpz_yzz_yyzz_0[j] = pa_y[j] * tpz_zz_yyzz_0[j] + fl1_fx * tpz_zz_yzz_0[j];

            tpx_yzz_yzzz_0[j] = pa_y[j] * tpx_zz_yzzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFG_400_450(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (400,450)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        auto tpx_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 50);

        auto tpy_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_zz_xxx_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 51);

        auto tpy_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_zz_xxy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 52);

        auto tpy_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_zz_xxz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 53);

        auto tpy_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_zz_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 54);

        auto tpy_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zz_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 55);

        auto tpy_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zz_xzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 56);

        auto tpy_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zz_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 57);

        auto tpy_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zz_yyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 58);

        auto tpy_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 59);

        auto tpy_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        auto ts_zz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 88);

        auto ts_zz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 89);

        // set up pointers to integrals

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

        auto tpx_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 149);

        auto tpy_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tpz_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 149);

        // Batch of Integrals (400,450)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yzz_zzzz_0, tpx_z_xxxx_0, tpx_z_xxxy_0, \
                                     tpx_z_xxxz_0, tpx_z_xxyy_0, tpx_z_xxyz_0, tpx_z_xxzz_0, tpx_z_xyyy_0, tpx_z_xyyz_0, \
                                     tpx_z_xyzz_0, tpx_z_xzzz_0, tpx_z_yyyy_0, tpx_z_yyyz_0, tpx_z_yyzz_0, tpx_z_yzzz_0, \
                                     tpx_z_zzzz_0, tpx_zz_xxx_0, tpx_zz_xxxx_0, tpx_zz_xxxy_0, tpx_zz_xxxz_0, \
                                     tpx_zz_xxy_0, tpx_zz_xxyy_0, tpx_zz_xxyz_0, tpx_zz_xxz_0, tpx_zz_xxzz_0, \
                                     tpx_zz_xyy_0, tpx_zz_xyyy_0, tpx_zz_xyyz_0, tpx_zz_xyz_0, tpx_zz_xyzz_0, \
                                     tpx_zz_xzz_0, tpx_zz_xzzz_0, tpx_zz_yyy_0, tpx_zz_yyyy_0, tpx_zz_yyyz_0, \
                                     tpx_zz_yyz_0, tpx_zz_yyzz_0, tpx_zz_yzz_0, tpx_zz_yzzz_0, tpx_zz_zzz_0, \
                                     tpx_zz_zzzz_0, tpx_zzz_xxxx_0, tpx_zzz_xxxy_0, tpx_zzz_xxxz_0, tpx_zzz_xxyy_0, \
                                     tpx_zzz_xxyz_0, tpx_zzz_xxzz_0, tpx_zzz_xyyy_0, tpx_zzz_xyyz_0, tpx_zzz_xyzz_0, \
                                     tpx_zzz_xzzz_0, tpx_zzz_yyyy_0, tpx_zzz_yyyz_0, tpx_zzz_yyzz_0, tpx_zzz_yzzz_0, \
                                     tpx_zzz_zzzz_0, tpy_yzz_yzzz_0, tpy_yzz_zzzz_0, tpy_z_xxxx_0, tpy_z_xxxy_0, \
                                     tpy_z_xxxz_0, tpy_z_xxyy_0, tpy_z_xxyz_0, tpy_z_xxzz_0, tpy_z_xyyy_0, tpy_z_xyyz_0, \
                                     tpy_z_xyzz_0, tpy_z_xzzz_0, tpy_z_yyyy_0, tpy_z_yyyz_0, tpy_z_yyzz_0, tpy_z_yzzz_0, \
                                     tpy_z_zzzz_0, tpy_zz_xxx_0, tpy_zz_xxxx_0, tpy_zz_xxxy_0, tpy_zz_xxxz_0, \
                                     tpy_zz_xxy_0, tpy_zz_xxyy_0, tpy_zz_xxyz_0, tpy_zz_xxz_0, tpy_zz_xxzz_0, \
                                     tpy_zz_xyy_0, tpy_zz_xyyy_0, tpy_zz_xyyz_0, tpy_zz_xyz_0, tpy_zz_xyzz_0, \
                                     tpy_zz_xzz_0, tpy_zz_xzzz_0, tpy_zz_yyy_0, tpy_zz_yyyy_0, tpy_zz_yyyz_0, \
                                     tpy_zz_yyz_0, tpy_zz_yyzz_0, tpy_zz_yzz_0, tpy_zz_yzzz_0, tpy_zz_zzz_0, \
                                     tpy_zz_zzzz_0, tpy_zzz_xxxx_0, tpy_zzz_xxxy_0, tpy_zzz_xxxz_0, tpy_zzz_xxyy_0, \
                                     tpy_zzz_xxyz_0, tpy_zzz_xxzz_0, tpy_zzz_xyyy_0, tpy_zzz_xyyz_0, tpy_zzz_xyzz_0, \
                                     tpy_zzz_xzzz_0, tpy_zzz_yyyy_0, tpy_zzz_yyyz_0, tpy_zzz_yyzz_0, tpy_zzz_yzzz_0, \
                                     tpy_zzz_zzzz_0, tpz_yzz_yzzz_0, tpz_yzz_zzzz_0, tpz_z_xxxx_0, tpz_z_xxxy_0, \
                                     tpz_z_xxxz_0, tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxzz_0, tpz_z_xyyy_0, tpz_z_xyyz_0, \
                                     tpz_z_xyzz_0, tpz_z_xzzz_0, tpz_z_yyyy_0, tpz_z_yyyz_0, tpz_z_yyzz_0, tpz_z_yzzz_0, \
                                     tpz_z_zzzz_0, tpz_zz_xxx_0, tpz_zz_xxxx_0, tpz_zz_xxxy_0, tpz_zz_xxxz_0, \
                                     tpz_zz_xxy_0, tpz_zz_xxyy_0, tpz_zz_xxyz_0, tpz_zz_xxz_0, tpz_zz_xxzz_0, \
                                     tpz_zz_xyy_0, tpz_zz_xyyy_0, tpz_zz_xyyz_0, tpz_zz_xyz_0, tpz_zz_xyzz_0, \
                                     tpz_zz_xzz_0, tpz_zz_xzzz_0, tpz_zz_yyy_0, tpz_zz_yyyy_0, tpz_zz_yyyz_0, \
                                     tpz_zz_yyz_0, tpz_zz_yyzz_0, tpz_zz_yzz_0, tpz_zz_yzzz_0, tpz_zz_zzz_0, \
                                     tpz_zz_zzzz_0, tpz_zzz_xxxx_0, tpz_zzz_xxxy_0, tpz_zzz_xxxz_0, tpz_zzz_xxyy_0, \
                                     tpz_zzz_xxyz_0, tpz_zzz_xxzz_0, tpz_zzz_xyyy_0, tpz_zzz_xyyz_0, tpz_zzz_xyzz_0, \
                                     tpz_zzz_xzzz_0, tpz_zzz_yyyy_0, tpz_zzz_yyyz_0, tpz_zzz_yyzz_0, tpz_zzz_yzzz_0, \
                                     tpz_zzz_zzzz_0, ts_zz_xxxx_0, ts_zz_xxxy_0, ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, \
                                     ts_zz_xxzz_0, ts_zz_xyyy_0, ts_zz_xyyz_0, ts_zz_xyzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, \
                                     ts_zz_yyyz_0, ts_zz_yyzz_0, ts_zz_yzzz_0, ts_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_yzz_yzzz_0[j] = pa_y[j] * tpy_zz_yzzz_0[j] + 0.5 * fl1_fx * tpy_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zz_yzzz_0[j];

            tpz_yzz_yzzz_0[j] = pa_y[j] * tpz_zz_yzzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzz_0[j];

            tpx_yzz_zzzz_0[j] = pa_y[j] * tpx_zz_zzzz_0[j];

            tpy_yzz_zzzz_0[j] = pa_y[j] * tpy_zz_zzzz_0[j] - fl1_fgb * fl1_fx * ts_zz_zzzz_0[j];

            tpz_yzz_zzzz_0[j] = pa_y[j] * tpz_zz_zzzz_0[j];

            tpx_zzz_xxxx_0[j] = pa_z[j] * tpx_zz_xxxx_0[j] + fl1_fx * tpx_z_xxxx_0[j];

            tpy_zzz_xxxx_0[j] = pa_z[j] * tpy_zz_xxxx_0[j] + fl1_fx * tpy_z_xxxx_0[j];

            tpz_zzz_xxxx_0[j] = pa_z[j] * tpz_zz_xxxx_0[j] + fl1_fx * tpz_z_xxxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxx_0[j];

            tpx_zzz_xxxy_0[j] = pa_z[j] * tpx_zz_xxxy_0[j] + fl1_fx * tpx_z_xxxy_0[j];

            tpy_zzz_xxxy_0[j] = pa_z[j] * tpy_zz_xxxy_0[j] + fl1_fx * tpy_z_xxxy_0[j];

            tpz_zzz_xxxy_0[j] = pa_z[j] * tpz_zz_xxxy_0[j] + fl1_fx * tpz_z_xxxy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxy_0[j];

            tpx_zzz_xxxz_0[j] = pa_z[j] * tpx_zz_xxxz_0[j] + fl1_fx * tpx_z_xxxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxx_0[j];

            tpy_zzz_xxxz_0[j] = pa_z[j] * tpy_zz_xxxz_0[j] + fl1_fx * tpy_z_xxxz_0[j] + 0.5 * fl1_fx * tpy_zz_xxx_0[j];

            tpz_zzz_xxxz_0[j] =
                pa_z[j] * tpz_zz_xxxz_0[j] + fl1_fx * tpz_z_xxxz_0[j] + 0.5 * fl1_fx * tpz_zz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxxz_0[j];

            tpx_zzz_xxyy_0[j] = pa_z[j] * tpx_zz_xxyy_0[j] + fl1_fx * tpx_z_xxyy_0[j];

            tpy_zzz_xxyy_0[j] = pa_z[j] * tpy_zz_xxyy_0[j] + fl1_fx * tpy_z_xxyy_0[j];

            tpz_zzz_xxyy_0[j] = pa_z[j] * tpz_zz_xxyy_0[j] + fl1_fx * tpz_z_xxyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxyy_0[j];

            tpx_zzz_xxyz_0[j] = pa_z[j] * tpx_zz_xxyz_0[j] + fl1_fx * tpx_z_xxyz_0[j] + 0.5 * fl1_fx * tpx_zz_xxy_0[j];

            tpy_zzz_xxyz_0[j] = pa_z[j] * tpy_zz_xxyz_0[j] + fl1_fx * tpy_z_xxyz_0[j] + 0.5 * fl1_fx * tpy_zz_xxy_0[j];

            tpz_zzz_xxyz_0[j] =
                pa_z[j] * tpz_zz_xxyz_0[j] + fl1_fx * tpz_z_xxyz_0[j] + 0.5 * fl1_fx * tpz_zz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxyz_0[j];

            tpx_zzz_xxzz_0[j] = pa_z[j] * tpx_zz_xxzz_0[j] + fl1_fx * tpx_z_xxzz_0[j] + fl1_fx * tpx_zz_xxz_0[j];

            tpy_zzz_xxzz_0[j] = pa_z[j] * tpy_zz_xxzz_0[j] + fl1_fx * tpy_z_xxzz_0[j] + fl1_fx * tpy_zz_xxz_0[j];

            tpz_zzz_xxzz_0[j] = pa_z[j] * tpz_zz_xxzz_0[j] + fl1_fx * tpz_z_xxzz_0[j] + fl1_fx * tpz_zz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxzz_0[j];

            tpx_zzz_xyyy_0[j] = pa_z[j] * tpx_zz_xyyy_0[j] + fl1_fx * tpx_z_xyyy_0[j];

            tpy_zzz_xyyy_0[j] = pa_z[j] * tpy_zz_xyyy_0[j] + fl1_fx * tpy_z_xyyy_0[j];

            tpz_zzz_xyyy_0[j] = pa_z[j] * tpz_zz_xyyy_0[j] + fl1_fx * tpz_z_xyyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyyy_0[j];

            tpx_zzz_xyyz_0[j] = pa_z[j] * tpx_zz_xyyz_0[j] + fl1_fx * tpx_z_xyyz_0[j] + 0.5 * fl1_fx * tpx_zz_xyy_0[j];

            tpy_zzz_xyyz_0[j] = pa_z[j] * tpy_zz_xyyz_0[j] + fl1_fx * tpy_z_xyyz_0[j] + 0.5 * fl1_fx * tpy_zz_xyy_0[j];

            tpz_zzz_xyyz_0[j] =
                pa_z[j] * tpz_zz_xyyz_0[j] + fl1_fx * tpz_z_xyyz_0[j] + 0.5 * fl1_fx * tpz_zz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyyz_0[j];

            tpx_zzz_xyzz_0[j] = pa_z[j] * tpx_zz_xyzz_0[j] + fl1_fx * tpx_z_xyzz_0[j] + fl1_fx * tpx_zz_xyz_0[j];

            tpy_zzz_xyzz_0[j] = pa_z[j] * tpy_zz_xyzz_0[j] + fl1_fx * tpy_z_xyzz_0[j] + fl1_fx * tpy_zz_xyz_0[j];

            tpz_zzz_xyzz_0[j] = pa_z[j] * tpz_zz_xyzz_0[j] + fl1_fx * tpz_z_xyzz_0[j] + fl1_fx * tpz_zz_xyz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyzz_0[j];

            tpx_zzz_xzzz_0[j] = pa_z[j] * tpx_zz_xzzz_0[j] + fl1_fx * tpx_z_xzzz_0[j] + 1.5 * fl1_fx * tpx_zz_xzz_0[j];

            tpy_zzz_xzzz_0[j] = pa_z[j] * tpy_zz_xzzz_0[j] + fl1_fx * tpy_z_xzzz_0[j] + 1.5 * fl1_fx * tpy_zz_xzz_0[j];

            tpz_zzz_xzzz_0[j] =
                pa_z[j] * tpz_zz_xzzz_0[j] + fl1_fx * tpz_z_xzzz_0[j] + 1.5 * fl1_fx * tpz_zz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xzzz_0[j];

            tpx_zzz_yyyy_0[j] = pa_z[j] * tpx_zz_yyyy_0[j] + fl1_fx * tpx_z_yyyy_0[j];

            tpy_zzz_yyyy_0[j] = pa_z[j] * tpy_zz_yyyy_0[j] + fl1_fx * tpy_z_yyyy_0[j];

            tpz_zzz_yyyy_0[j] = pa_z[j] * tpz_zz_yyyy_0[j] + fl1_fx * tpz_z_yyyy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyyy_0[j];

            tpx_zzz_yyyz_0[j] = pa_z[j] * tpx_zz_yyyz_0[j] + fl1_fx * tpx_z_yyyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyy_0[j];

            tpy_zzz_yyyz_0[j] = pa_z[j] * tpy_zz_yyyz_0[j] + fl1_fx * tpy_z_yyyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyy_0[j];

            tpz_zzz_yyyz_0[j] =
                pa_z[j] * tpz_zz_yyyz_0[j] + fl1_fx * tpz_z_yyyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyyz_0[j];

            tpx_zzz_yyzz_0[j] = pa_z[j] * tpx_zz_yyzz_0[j] + fl1_fx * tpx_z_yyzz_0[j] + fl1_fx * tpx_zz_yyz_0[j];

            tpy_zzz_yyzz_0[j] = pa_z[j] * tpy_zz_yyzz_0[j] + fl1_fx * tpy_z_yyzz_0[j] + fl1_fx * tpy_zz_yyz_0[j];

            tpz_zzz_yyzz_0[j] = pa_z[j] * tpz_zz_yyzz_0[j] + fl1_fx * tpz_z_yyzz_0[j] + fl1_fx * tpz_zz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyzz_0[j];

            tpx_zzz_yzzz_0[j] = pa_z[j] * tpx_zz_yzzz_0[j] + fl1_fx * tpx_z_yzzz_0[j] + 1.5 * fl1_fx * tpx_zz_yzz_0[j];

            tpy_zzz_yzzz_0[j] = pa_z[j] * tpy_zz_yzzz_0[j] + fl1_fx * tpy_z_yzzz_0[j] + 1.5 * fl1_fx * tpy_zz_yzz_0[j];

            tpz_zzz_yzzz_0[j] =
                pa_z[j] * tpz_zz_yzzz_0[j] + fl1_fx * tpz_z_yzzz_0[j] + 1.5 * fl1_fx * tpz_zz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zz_yzzz_0[j];

            tpx_zzz_zzzz_0[j] = pa_z[j] * tpx_zz_zzzz_0[j] + fl1_fx * tpx_z_zzzz_0[j] + 2.0 * fl1_fx * tpx_zz_zzz_0[j];

            tpy_zzz_zzzz_0[j] = pa_z[j] * tpy_zz_zzzz_0[j] + fl1_fx * tpy_z_zzzz_0[j] + 2.0 * fl1_fx * tpy_zz_zzz_0[j];

            tpz_zzz_zzzz_0[j] =
                pa_z[j] * tpz_zz_zzzz_0[j] + fl1_fx * tpz_z_zzzz_0[j] + 2.0 * fl1_fx * tpz_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zz_zzzz_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
