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

#include "LinearMomentumRecFuncForFF.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForFF(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForFF_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFF_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFF_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFF_150_200(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFF_200_250(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFF_250_300(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForFF_0_50(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx);

        auto tpy_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx);

        auto tpz_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx);

        auto tpx_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 1);

        auto tpy_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 2);

        auto tpy_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 3);

        auto tpy_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 4);

        auto tpy_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 5);

        auto tpy_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 6);

        auto tpy_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 7);

        auto tpy_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 8);

        auto tpy_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 8);

        auto tpx_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 9);

        auto tpy_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 9);

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpx_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx);

        auto tpy_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx);

        auto tpz_xx_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx);

        auto tpx_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 1);

        auto tpy_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 1);

        auto tpz_xx_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 1);

        auto tpx_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 2);

        auto tpy_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 2);

        auto tpz_xx_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 2);

        auto tpx_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 3);

        auto tpy_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 3);

        auto tpz_xx_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 3);

        auto tpx_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 4);

        auto tpy_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 4);

        auto tpz_xx_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 4);

        auto tpx_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 5);

        auto tpy_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 5);

        auto tpz_xx_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 5);

        auto tpx_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 6);

        auto tpy_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 6);

        auto tpz_xy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 6);

        auto tpx_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 7);

        auto tpy_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 7);

        auto tpz_xy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 7);

        auto tpx_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 8);

        auto tpy_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 8);

        auto tpz_xy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 8);

        auto tpx_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 9);

        auto tpy_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 9);

        auto tpz_xy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 9);

        auto tpx_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 10);

        auto tpy_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 10);

        auto tpz_xy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 10);

        auto tpx_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 11);

        auto tpy_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 11);

        auto tpz_xy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 11);

        auto ts_xx_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx);

        auto ts_xx_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 1);

        auto ts_xx_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 2);

        auto ts_xx_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 3);

        auto ts_xx_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 4);

        auto ts_xx_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 5);

        auto ts_xx_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 6);

        auto ts_xx_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 7);

        auto ts_xx_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 8);

        auto ts_xx_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 9);

        auto ts_xy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 10);

        auto ts_xy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 11);

        auto ts_xy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 12);

        auto ts_xy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 13);

        auto ts_xy_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 14);

        auto ts_xy_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 15);

        auto ts_xy_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 16);

        // set up pointers to integrals

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

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_x_xxx_0, tpx_x_xxy_0, tpx_x_xxz_0, tpx_x_xyy_0, \
                                     tpx_x_xyz_0, tpx_x_xzz_0, tpx_x_yyy_0, tpx_x_yyz_0, tpx_x_yzz_0, tpx_x_zzz_0, \
                                     tpx_xx_xx_0, tpx_xx_xxx_0, tpx_xx_xxy_0, tpx_xx_xxz_0, tpx_xx_xy_0, tpx_xx_xyy_0, \
                                     tpx_xx_xyz_0, tpx_xx_xz_0, tpx_xx_xzz_0, tpx_xx_yy_0, tpx_xx_yyy_0, tpx_xx_yyz_0, \
                                     tpx_xx_yz_0, tpx_xx_yzz_0, tpx_xx_zz_0, tpx_xx_zzz_0, tpx_xxx_xxx_0, \
                                     tpx_xxx_xxy_0, tpx_xxx_xxz_0, tpx_xxx_xyy_0, tpx_xxx_xyz_0, tpx_xxx_xzz_0, \
                                     tpx_xxx_yyy_0, tpx_xxx_yyz_0, tpx_xxx_yzz_0, tpx_xxx_zzz_0, tpx_xxy_xxx_0, \
                                     tpx_xxy_xxy_0, tpx_xxy_xxz_0, tpx_xxy_xyy_0, tpx_xxy_xyz_0, tpx_xxy_xzz_0, \
                                     tpx_xxy_yyy_0, tpx_xy_xx_0, tpx_xy_xxx_0, tpx_xy_xxy_0, tpx_xy_xxz_0, tpx_xy_xy_0, \
                                     tpx_xy_xyy_0, tpx_xy_xyz_0, tpx_xy_xz_0, tpx_xy_xzz_0, tpx_xy_yy_0, tpx_xy_yyy_0, \
                                     tpx_xy_yz_0, tpx_xy_zz_0, tpx_y_xxx_0, tpx_y_xxy_0, tpx_y_xxz_0, tpx_y_xyy_0, \
                                     tpx_y_xyz_0, tpx_y_xzz_0, tpx_y_yyy_0, tpy_x_xxx_0, tpy_x_xxy_0, tpy_x_xxz_0, \
                                     tpy_x_xyy_0, tpy_x_xyz_0, tpy_x_xzz_0, tpy_x_yyy_0, tpy_x_yyz_0, tpy_x_yzz_0, \
                                     tpy_x_zzz_0, tpy_xx_xx_0, tpy_xx_xxx_0, tpy_xx_xxy_0, tpy_xx_xxz_0, tpy_xx_xy_0, \
                                     tpy_xx_xyy_0, tpy_xx_xyz_0, tpy_xx_xz_0, tpy_xx_xzz_0, tpy_xx_yy_0, tpy_xx_yyy_0, \
                                     tpy_xx_yyz_0, tpy_xx_yz_0, tpy_xx_yzz_0, tpy_xx_zz_0, tpy_xx_zzz_0, tpy_xxx_xxx_0, \
                                     tpy_xxx_xxy_0, tpy_xxx_xxz_0, tpy_xxx_xyy_0, tpy_xxx_xyz_0, tpy_xxx_xzz_0, \
                                     tpy_xxx_yyy_0, tpy_xxx_yyz_0, tpy_xxx_yzz_0, tpy_xxx_zzz_0, tpy_xxy_xxx_0, \
                                     tpy_xxy_xxy_0, tpy_xxy_xxz_0, tpy_xxy_xyy_0, tpy_xxy_xyz_0, tpy_xxy_xzz_0, \
                                     tpy_xxy_yyy_0, tpy_xy_xx_0, tpy_xy_xxx_0, tpy_xy_xxy_0, tpy_xy_xxz_0, tpy_xy_xy_0, \
                                     tpy_xy_xyy_0, tpy_xy_xyz_0, tpy_xy_xz_0, tpy_xy_xzz_0, tpy_xy_yy_0, tpy_xy_yyy_0, \
                                     tpy_xy_yz_0, tpy_xy_zz_0, tpy_y_xxx_0, tpy_y_xxy_0, tpy_y_xxz_0, tpy_y_xyy_0, \
                                     tpy_y_xyz_0, tpy_y_xzz_0, tpy_y_yyy_0, tpz_x_xxx_0, tpz_x_xxy_0, tpz_x_xxz_0, \
                                     tpz_x_xyy_0, tpz_x_xyz_0, tpz_x_xzz_0, tpz_x_yyy_0, tpz_x_yyz_0, tpz_x_yzz_0, \
                                     tpz_x_zzz_0, tpz_xx_xx_0, tpz_xx_xxx_0, tpz_xx_xxy_0, tpz_xx_xxz_0, tpz_xx_xy_0, \
                                     tpz_xx_xyy_0, tpz_xx_xyz_0, tpz_xx_xz_0, tpz_xx_xzz_0, tpz_xx_yy_0, tpz_xx_yyy_0, \
                                     tpz_xx_yyz_0, tpz_xx_yz_0, tpz_xx_yzz_0, tpz_xx_zz_0, tpz_xx_zzz_0, tpz_xxx_xxx_0, \
                                     tpz_xxx_xxy_0, tpz_xxx_xxz_0, tpz_xxx_xyy_0, tpz_xxx_xyz_0, tpz_xxx_xzz_0, \
                                     tpz_xxx_yyy_0, tpz_xxx_yyz_0, tpz_xxx_yzz_0, tpz_xxx_zzz_0, tpz_xxy_xxx_0, \
                                     tpz_xxy_xxy_0, tpz_xxy_xxz_0, tpz_xxy_xyy_0, tpz_xxy_xyz_0, tpz_xxy_xzz_0, \
                                     tpz_xy_xx_0, tpz_xy_xxx_0, tpz_xy_xxy_0, tpz_xy_xxz_0, tpz_xy_xy_0, tpz_xy_xyy_0, \
                                     tpz_xy_xyz_0, tpz_xy_xz_0, tpz_xy_xzz_0, tpz_xy_yy_0, tpz_xy_yz_0, tpz_xy_zz_0, \
                                     tpz_y_xxx_0, tpz_y_xxy_0, tpz_y_xxz_0, tpz_y_xyy_0, tpz_y_xyz_0, tpz_y_xzz_0, \
                                     ts_xx_xxx_0, ts_xx_xxy_0, ts_xx_xxz_0, ts_xx_xyy_0, ts_xx_xyz_0, ts_xx_xzz_0, \
                                     ts_xx_yyy_0, ts_xx_yyz_0, ts_xx_yzz_0, ts_xx_zzz_0, ts_xy_xxx_0, ts_xy_xxy_0, \
                                     ts_xy_xxz_0, ts_xy_xyy_0, ts_xy_xyz_0, ts_xy_xzz_0, ts_xy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxx_xxx_0[j] =
                pa_x[j] * tpx_xx_xxx_0[j] + fl1_fx * tpx_x_xxx_0[j] + 1.5 * fl1_fx * tpx_xx_xx_0[j] - fl1_fgb * fl1_fx * ts_xx_xxx_0[j];

            tpy_xxx_xxx_0[j] = pa_x[j] * tpy_xx_xxx_0[j] + fl1_fx * tpy_x_xxx_0[j] + 1.5 * fl1_fx * tpy_xx_xx_0[j];

            tpz_xxx_xxx_0[j] = pa_x[j] * tpz_xx_xxx_0[j] + fl1_fx * tpz_x_xxx_0[j] + 1.5 * fl1_fx * tpz_xx_xx_0[j];

            tpx_xxx_xxy_0[j] = pa_x[j] * tpx_xx_xxy_0[j] + fl1_fx * tpx_x_xxy_0[j] + fl1_fx * tpx_xx_xy_0[j] - fl1_fgb * fl1_fx * ts_xx_xxy_0[j];

            tpy_xxx_xxy_0[j] = pa_x[j] * tpy_xx_xxy_0[j] + fl1_fx * tpy_x_xxy_0[j] + fl1_fx * tpy_xx_xy_0[j];

            tpz_xxx_xxy_0[j] = pa_x[j] * tpz_xx_xxy_0[j] + fl1_fx * tpz_x_xxy_0[j] + fl1_fx * tpz_xx_xy_0[j];

            tpx_xxx_xxz_0[j] = pa_x[j] * tpx_xx_xxz_0[j] + fl1_fx * tpx_x_xxz_0[j] + fl1_fx * tpx_xx_xz_0[j] - fl1_fgb * fl1_fx * ts_xx_xxz_0[j];

            tpy_xxx_xxz_0[j] = pa_x[j] * tpy_xx_xxz_0[j] + fl1_fx * tpy_x_xxz_0[j] + fl1_fx * tpy_xx_xz_0[j];

            tpz_xxx_xxz_0[j] = pa_x[j] * tpz_xx_xxz_0[j] + fl1_fx * tpz_x_xxz_0[j] + fl1_fx * tpz_xx_xz_0[j];

            tpx_xxx_xyy_0[j] =
                pa_x[j] * tpx_xx_xyy_0[j] + fl1_fx * tpx_x_xyy_0[j] + 0.5 * fl1_fx * tpx_xx_yy_0[j] - fl1_fgb * fl1_fx * ts_xx_xyy_0[j];

            tpy_xxx_xyy_0[j] = pa_x[j] * tpy_xx_xyy_0[j] + fl1_fx * tpy_x_xyy_0[j] + 0.5 * fl1_fx * tpy_xx_yy_0[j];

            tpz_xxx_xyy_0[j] = pa_x[j] * tpz_xx_xyy_0[j] + fl1_fx * tpz_x_xyy_0[j] + 0.5 * fl1_fx * tpz_xx_yy_0[j];

            tpx_xxx_xyz_0[j] =
                pa_x[j] * tpx_xx_xyz_0[j] + fl1_fx * tpx_x_xyz_0[j] + 0.5 * fl1_fx * tpx_xx_yz_0[j] - fl1_fgb * fl1_fx * ts_xx_xyz_0[j];

            tpy_xxx_xyz_0[j] = pa_x[j] * tpy_xx_xyz_0[j] + fl1_fx * tpy_x_xyz_0[j] + 0.5 * fl1_fx * tpy_xx_yz_0[j];

            tpz_xxx_xyz_0[j] = pa_x[j] * tpz_xx_xyz_0[j] + fl1_fx * tpz_x_xyz_0[j] + 0.5 * fl1_fx * tpz_xx_yz_0[j];

            tpx_xxx_xzz_0[j] =
                pa_x[j] * tpx_xx_xzz_0[j] + fl1_fx * tpx_x_xzz_0[j] + 0.5 * fl1_fx * tpx_xx_zz_0[j] - fl1_fgb * fl1_fx * ts_xx_xzz_0[j];

            tpy_xxx_xzz_0[j] = pa_x[j] * tpy_xx_xzz_0[j] + fl1_fx * tpy_x_xzz_0[j] + 0.5 * fl1_fx * tpy_xx_zz_0[j];

            tpz_xxx_xzz_0[j] = pa_x[j] * tpz_xx_xzz_0[j] + fl1_fx * tpz_x_xzz_0[j] + 0.5 * fl1_fx * tpz_xx_zz_0[j];

            tpx_xxx_yyy_0[j] = pa_x[j] * tpx_xx_yyy_0[j] + fl1_fx * tpx_x_yyy_0[j] - fl1_fgb * fl1_fx * ts_xx_yyy_0[j];

            tpy_xxx_yyy_0[j] = pa_x[j] * tpy_xx_yyy_0[j] + fl1_fx * tpy_x_yyy_0[j];

            tpz_xxx_yyy_0[j] = pa_x[j] * tpz_xx_yyy_0[j] + fl1_fx * tpz_x_yyy_0[j];

            tpx_xxx_yyz_0[j] = pa_x[j] * tpx_xx_yyz_0[j] + fl1_fx * tpx_x_yyz_0[j] - fl1_fgb * fl1_fx * ts_xx_yyz_0[j];

            tpy_xxx_yyz_0[j] = pa_x[j] * tpy_xx_yyz_0[j] + fl1_fx * tpy_x_yyz_0[j];

            tpz_xxx_yyz_0[j] = pa_x[j] * tpz_xx_yyz_0[j] + fl1_fx * tpz_x_yyz_0[j];

            tpx_xxx_yzz_0[j] = pa_x[j] * tpx_xx_yzz_0[j] + fl1_fx * tpx_x_yzz_0[j] - fl1_fgb * fl1_fx * ts_xx_yzz_0[j];

            tpy_xxx_yzz_0[j] = pa_x[j] * tpy_xx_yzz_0[j] + fl1_fx * tpy_x_yzz_0[j];

            tpz_xxx_yzz_0[j] = pa_x[j] * tpz_xx_yzz_0[j] + fl1_fx * tpz_x_yzz_0[j];

            tpx_xxx_zzz_0[j] = pa_x[j] * tpx_xx_zzz_0[j] + fl1_fx * tpx_x_zzz_0[j] - fl1_fgb * fl1_fx * ts_xx_zzz_0[j];

            tpy_xxx_zzz_0[j] = pa_x[j] * tpy_xx_zzz_0[j] + fl1_fx * tpy_x_zzz_0[j];

            tpz_xxx_zzz_0[j] = pa_x[j] * tpz_xx_zzz_0[j] + fl1_fx * tpz_x_zzz_0[j];

            tpx_xxy_xxx_0[j] =
                pa_x[j] * tpx_xy_xxx_0[j] + 0.5 * fl1_fx * tpx_y_xxx_0[j] + 1.5 * fl1_fx * tpx_xy_xx_0[j] - fl1_fgb * fl1_fx * ts_xy_xxx_0[j];

            tpy_xxy_xxx_0[j] = pa_x[j] * tpy_xy_xxx_0[j] + 0.5 * fl1_fx * tpy_y_xxx_0[j] + 1.5 * fl1_fx * tpy_xy_xx_0[j];

            tpz_xxy_xxx_0[j] = pa_x[j] * tpz_xy_xxx_0[j] + 0.5 * fl1_fx * tpz_y_xxx_0[j] + 1.5 * fl1_fx * tpz_xy_xx_0[j];

            tpx_xxy_xxy_0[j] =
                pa_x[j] * tpx_xy_xxy_0[j] + 0.5 * fl1_fx * tpx_y_xxy_0[j] + fl1_fx * tpx_xy_xy_0[j] - fl1_fgb * fl1_fx * ts_xy_xxy_0[j];

            tpy_xxy_xxy_0[j] = pa_x[j] * tpy_xy_xxy_0[j] + 0.5 * fl1_fx * tpy_y_xxy_0[j] + fl1_fx * tpy_xy_xy_0[j];

            tpz_xxy_xxy_0[j] = pa_x[j] * tpz_xy_xxy_0[j] + 0.5 * fl1_fx * tpz_y_xxy_0[j] + fl1_fx * tpz_xy_xy_0[j];

            tpx_xxy_xxz_0[j] =
                pa_x[j] * tpx_xy_xxz_0[j] + 0.5 * fl1_fx * tpx_y_xxz_0[j] + fl1_fx * tpx_xy_xz_0[j] - fl1_fgb * fl1_fx * ts_xy_xxz_0[j];

            tpy_xxy_xxz_0[j] = pa_x[j] * tpy_xy_xxz_0[j] + 0.5 * fl1_fx * tpy_y_xxz_0[j] + fl1_fx * tpy_xy_xz_0[j];

            tpz_xxy_xxz_0[j] = pa_x[j] * tpz_xy_xxz_0[j] + 0.5 * fl1_fx * tpz_y_xxz_0[j] + fl1_fx * tpz_xy_xz_0[j];

            tpx_xxy_xyy_0[j] =
                pa_x[j] * tpx_xy_xyy_0[j] + 0.5 * fl1_fx * tpx_y_xyy_0[j] + 0.5 * fl1_fx * tpx_xy_yy_0[j] - fl1_fgb * fl1_fx * ts_xy_xyy_0[j];

            tpy_xxy_xyy_0[j] = pa_x[j] * tpy_xy_xyy_0[j] + 0.5 * fl1_fx * tpy_y_xyy_0[j] + 0.5 * fl1_fx * tpy_xy_yy_0[j];

            tpz_xxy_xyy_0[j] = pa_x[j] * tpz_xy_xyy_0[j] + 0.5 * fl1_fx * tpz_y_xyy_0[j] + 0.5 * fl1_fx * tpz_xy_yy_0[j];

            tpx_xxy_xyz_0[j] =
                pa_x[j] * tpx_xy_xyz_0[j] + 0.5 * fl1_fx * tpx_y_xyz_0[j] + 0.5 * fl1_fx * tpx_xy_yz_0[j] - fl1_fgb * fl1_fx * ts_xy_xyz_0[j];

            tpy_xxy_xyz_0[j] = pa_x[j] * tpy_xy_xyz_0[j] + 0.5 * fl1_fx * tpy_y_xyz_0[j] + 0.5 * fl1_fx * tpy_xy_yz_0[j];

            tpz_xxy_xyz_0[j] = pa_x[j] * tpz_xy_xyz_0[j] + 0.5 * fl1_fx * tpz_y_xyz_0[j] + 0.5 * fl1_fx * tpz_xy_yz_0[j];

            tpx_xxy_xzz_0[j] =
                pa_x[j] * tpx_xy_xzz_0[j] + 0.5 * fl1_fx * tpx_y_xzz_0[j] + 0.5 * fl1_fx * tpx_xy_zz_0[j] - fl1_fgb * fl1_fx * ts_xy_xzz_0[j];

            tpy_xxy_xzz_0[j] = pa_x[j] * tpy_xy_xzz_0[j] + 0.5 * fl1_fx * tpy_y_xzz_0[j] + 0.5 * fl1_fx * tpy_xy_zz_0[j];

            tpz_xxy_xzz_0[j] = pa_x[j] * tpz_xy_xzz_0[j] + 0.5 * fl1_fx * tpz_y_xzz_0[j] + 0.5 * fl1_fx * tpz_xy_zz_0[j];

            tpx_xxy_yyy_0[j] = pa_x[j] * tpx_xy_yyy_0[j] + 0.5 * fl1_fx * tpx_y_yyy_0[j] - fl1_fgb * fl1_fx * ts_xy_yyy_0[j];

            tpy_xxy_yyy_0[j] = pa_x[j] * tpy_xy_yyy_0[j] + 0.5 * fl1_fx * tpy_y_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFF_50_100(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tpx_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 12);

        auto tpy_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 12);

        auto tpz_xz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 12);

        auto tpx_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 13);

        auto tpy_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 13);

        auto tpz_xz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 13);

        auto tpx_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 14);

        auto tpy_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 14);

        auto tpz_xz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 14);

        auto tpx_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 15);

        auto tpy_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 15);

        auto tpz_xz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 15);

        auto tpx_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 16);

        auto tpy_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 16);

        auto tpz_xz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 16);

        auto tpx_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 17);

        auto tpy_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 17);

        auto tpz_xz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 17);

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto ts_xy_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 17);

        auto ts_xy_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 18);

        auto ts_xy_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 19);

        auto ts_xz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 20);

        auto ts_xz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 21);

        auto ts_xz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 22);

        auto ts_xz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 23);

        auto ts_xz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 24);

        auto ts_xz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 25);

        auto ts_xz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 26);

        auto ts_xz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 27);

        auto ts_xz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 28);

        auto ts_xz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 29);

        auto ts_yy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 30);

        auto ts_yy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 31);

        auto ts_yy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 32);

        auto ts_yy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 33);

        // set up pointers to integrals

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

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxy_yyz_0, tpx_xxy_yzz_0, tpx_xxy_zzz_0, \
                                     tpx_xxz_xxx_0, tpx_xxz_xxy_0, tpx_xxz_xxz_0, tpx_xxz_xyy_0, tpx_xxz_xyz_0, \
                                     tpx_xxz_xzz_0, tpx_xxz_yyy_0, tpx_xxz_yyz_0, tpx_xxz_yzz_0, tpx_xxz_zzz_0, \
                                     tpx_xy_yyz_0, tpx_xy_yzz_0, tpx_xy_zzz_0, tpx_xyy_xxx_0, tpx_xyy_xxy_0, \
                                     tpx_xyy_xxz_0, tpx_xyy_xyy_0, tpx_xz_xx_0, tpx_xz_xxx_0, tpx_xz_xxy_0, tpx_xz_xxz_0, \
                                     tpx_xz_xy_0, tpx_xz_xyy_0, tpx_xz_xyz_0, tpx_xz_xz_0, tpx_xz_xzz_0, tpx_xz_yy_0, \
                                     tpx_xz_yyy_0, tpx_xz_yyz_0, tpx_xz_yz_0, tpx_xz_yzz_0, tpx_xz_zz_0, tpx_xz_zzz_0, \
                                     tpx_y_yyz_0, tpx_y_yzz_0, tpx_y_zzz_0, tpx_yy_xx_0, tpx_yy_xxx_0, tpx_yy_xxy_0, \
                                     tpx_yy_xxz_0, tpx_yy_xy_0, tpx_yy_xyy_0, tpx_yy_xz_0, tpx_yy_yy_0, tpx_z_xxx_0, \
                                     tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xyy_0, tpx_z_xyz_0, tpx_z_xzz_0, tpx_z_yyy_0, \
                                     tpx_z_yyz_0, tpx_z_yzz_0, tpx_z_zzz_0, tpy_xxy_yyz_0, tpy_xxy_yzz_0, \
                                     tpy_xxy_zzz_0, tpy_xxz_xxx_0, tpy_xxz_xxy_0, tpy_xxz_xxz_0, tpy_xxz_xyy_0, \
                                     tpy_xxz_xyz_0, tpy_xxz_xzz_0, tpy_xxz_yyy_0, tpy_xxz_yyz_0, tpy_xxz_yzz_0, \
                                     tpy_xxz_zzz_0, tpy_xy_yyz_0, tpy_xy_yzz_0, tpy_xy_zzz_0, tpy_xyy_xxx_0, \
                                     tpy_xyy_xxy_0, tpy_xyy_xxz_0, tpy_xz_xx_0, tpy_xz_xxx_0, tpy_xz_xxy_0, tpy_xz_xxz_0, \
                                     tpy_xz_xy_0, tpy_xz_xyy_0, tpy_xz_xyz_0, tpy_xz_xz_0, tpy_xz_xzz_0, tpy_xz_yy_0, \
                                     tpy_xz_yyy_0, tpy_xz_yyz_0, tpy_xz_yz_0, tpy_xz_yzz_0, tpy_xz_zz_0, tpy_xz_zzz_0, \
                                     tpy_y_yyz_0, tpy_y_yzz_0, tpy_y_zzz_0, tpy_yy_xx_0, tpy_yy_xxx_0, tpy_yy_xxy_0, \
                                     tpy_yy_xxz_0, tpy_yy_xy_0, tpy_yy_xz_0, tpy_z_xxx_0, tpy_z_xxy_0, tpy_z_xxz_0, \
                                     tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xzz_0, tpy_z_yyy_0, tpy_z_yyz_0, tpy_z_yzz_0, \
                                     tpy_z_zzz_0, tpz_xxy_yyy_0, tpz_xxy_yyz_0, tpz_xxy_yzz_0, tpz_xxy_zzz_0, \
                                     tpz_xxz_xxx_0, tpz_xxz_xxy_0, tpz_xxz_xxz_0, tpz_xxz_xyy_0, tpz_xxz_xyz_0, \
                                     tpz_xxz_xzz_0, tpz_xxz_yyy_0, tpz_xxz_yyz_0, tpz_xxz_yzz_0, tpz_xxz_zzz_0, \
                                     tpz_xy_yyy_0, tpz_xy_yyz_0, tpz_xy_yzz_0, tpz_xy_zzz_0, tpz_xyy_xxx_0, \
                                     tpz_xyy_xxy_0, tpz_xyy_xxz_0, tpz_xz_xx_0, tpz_xz_xxx_0, tpz_xz_xxy_0, tpz_xz_xxz_0, \
                                     tpz_xz_xy_0, tpz_xz_xyy_0, tpz_xz_xyz_0, tpz_xz_xz_0, tpz_xz_xzz_0, tpz_xz_yy_0, \
                                     tpz_xz_yyy_0, tpz_xz_yyz_0, tpz_xz_yz_0, tpz_xz_yzz_0, tpz_xz_zz_0, tpz_xz_zzz_0, \
                                     tpz_y_yyy_0, tpz_y_yyz_0, tpz_y_yzz_0, tpz_y_zzz_0, tpz_yy_xx_0, tpz_yy_xxx_0, \
                                     tpz_yy_xxy_0, tpz_yy_xxz_0, tpz_yy_xy_0, tpz_yy_xz_0, tpz_z_xxx_0, tpz_z_xxy_0, \
                                     tpz_z_xxz_0, tpz_z_xyy_0, tpz_z_xyz_0, tpz_z_xzz_0, tpz_z_yyy_0, tpz_z_yyz_0, \
                                     tpz_z_yzz_0, tpz_z_zzz_0, ts_xy_yyz_0, ts_xy_yzz_0, ts_xy_zzz_0, ts_xz_xxx_0, \
                                     ts_xz_xxy_0, ts_xz_xxz_0, ts_xz_xyy_0, ts_xz_xyz_0, ts_xz_xzz_0, ts_xz_yyy_0, \
                                     ts_xz_yyz_0, ts_xz_yzz_0, ts_xz_zzz_0, ts_yy_xxx_0, ts_yy_xxy_0, ts_yy_xxz_0, \
                                     ts_yy_xyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_xxy_yyy_0[j] = pa_x[j] * tpz_xy_yyy_0[j] + 0.5 * fl1_fx * tpz_y_yyy_0[j];

            tpx_xxy_yyz_0[j] = pa_x[j] * tpx_xy_yyz_0[j] + 0.5 * fl1_fx * tpx_y_yyz_0[j] - fl1_fgb * fl1_fx * ts_xy_yyz_0[j];

            tpy_xxy_yyz_0[j] = pa_x[j] * tpy_xy_yyz_0[j] + 0.5 * fl1_fx * tpy_y_yyz_0[j];

            tpz_xxy_yyz_0[j] = pa_x[j] * tpz_xy_yyz_0[j] + 0.5 * fl1_fx * tpz_y_yyz_0[j];

            tpx_xxy_yzz_0[j] = pa_x[j] * tpx_xy_yzz_0[j] + 0.5 * fl1_fx * tpx_y_yzz_0[j] - fl1_fgb * fl1_fx * ts_xy_yzz_0[j];

            tpy_xxy_yzz_0[j] = pa_x[j] * tpy_xy_yzz_0[j] + 0.5 * fl1_fx * tpy_y_yzz_0[j];

            tpz_xxy_yzz_0[j] = pa_x[j] * tpz_xy_yzz_0[j] + 0.5 * fl1_fx * tpz_y_yzz_0[j];

            tpx_xxy_zzz_0[j] = pa_x[j] * tpx_xy_zzz_0[j] + 0.5 * fl1_fx * tpx_y_zzz_0[j] - fl1_fgb * fl1_fx * ts_xy_zzz_0[j];

            tpy_xxy_zzz_0[j] = pa_x[j] * tpy_xy_zzz_0[j] + 0.5 * fl1_fx * tpy_y_zzz_0[j];

            tpz_xxy_zzz_0[j] = pa_x[j] * tpz_xy_zzz_0[j] + 0.5 * fl1_fx * tpz_y_zzz_0[j];

            tpx_xxz_xxx_0[j] =
                pa_x[j] * tpx_xz_xxx_0[j] + 0.5 * fl1_fx * tpx_z_xxx_0[j] + 1.5 * fl1_fx * tpx_xz_xx_0[j] - fl1_fgb * fl1_fx * ts_xz_xxx_0[j];

            tpy_xxz_xxx_0[j] = pa_x[j] * tpy_xz_xxx_0[j] + 0.5 * fl1_fx * tpy_z_xxx_0[j] + 1.5 * fl1_fx * tpy_xz_xx_0[j];

            tpz_xxz_xxx_0[j] = pa_x[j] * tpz_xz_xxx_0[j] + 0.5 * fl1_fx * tpz_z_xxx_0[j] + 1.5 * fl1_fx * tpz_xz_xx_0[j];

            tpx_xxz_xxy_0[j] =
                pa_x[j] * tpx_xz_xxy_0[j] + 0.5 * fl1_fx * tpx_z_xxy_0[j] + fl1_fx * tpx_xz_xy_0[j] - fl1_fgb * fl1_fx * ts_xz_xxy_0[j];

            tpy_xxz_xxy_0[j] = pa_x[j] * tpy_xz_xxy_0[j] + 0.5 * fl1_fx * tpy_z_xxy_0[j] + fl1_fx * tpy_xz_xy_0[j];

            tpz_xxz_xxy_0[j] = pa_x[j] * tpz_xz_xxy_0[j] + 0.5 * fl1_fx * tpz_z_xxy_0[j] + fl1_fx * tpz_xz_xy_0[j];

            tpx_xxz_xxz_0[j] =
                pa_x[j] * tpx_xz_xxz_0[j] + 0.5 * fl1_fx * tpx_z_xxz_0[j] + fl1_fx * tpx_xz_xz_0[j] - fl1_fgb * fl1_fx * ts_xz_xxz_0[j];

            tpy_xxz_xxz_0[j] = pa_x[j] * tpy_xz_xxz_0[j] + 0.5 * fl1_fx * tpy_z_xxz_0[j] + fl1_fx * tpy_xz_xz_0[j];

            tpz_xxz_xxz_0[j] = pa_x[j] * tpz_xz_xxz_0[j] + 0.5 * fl1_fx * tpz_z_xxz_0[j] + fl1_fx * tpz_xz_xz_0[j];

            tpx_xxz_xyy_0[j] =
                pa_x[j] * tpx_xz_xyy_0[j] + 0.5 * fl1_fx * tpx_z_xyy_0[j] + 0.5 * fl1_fx * tpx_xz_yy_0[j] - fl1_fgb * fl1_fx * ts_xz_xyy_0[j];

            tpy_xxz_xyy_0[j] = pa_x[j] * tpy_xz_xyy_0[j] + 0.5 * fl1_fx * tpy_z_xyy_0[j] + 0.5 * fl1_fx * tpy_xz_yy_0[j];

            tpz_xxz_xyy_0[j] = pa_x[j] * tpz_xz_xyy_0[j] + 0.5 * fl1_fx * tpz_z_xyy_0[j] + 0.5 * fl1_fx * tpz_xz_yy_0[j];

            tpx_xxz_xyz_0[j] =
                pa_x[j] * tpx_xz_xyz_0[j] + 0.5 * fl1_fx * tpx_z_xyz_0[j] + 0.5 * fl1_fx * tpx_xz_yz_0[j] - fl1_fgb * fl1_fx * ts_xz_xyz_0[j];

            tpy_xxz_xyz_0[j] = pa_x[j] * tpy_xz_xyz_0[j] + 0.5 * fl1_fx * tpy_z_xyz_0[j] + 0.5 * fl1_fx * tpy_xz_yz_0[j];

            tpz_xxz_xyz_0[j] = pa_x[j] * tpz_xz_xyz_0[j] + 0.5 * fl1_fx * tpz_z_xyz_0[j] + 0.5 * fl1_fx * tpz_xz_yz_0[j];

            tpx_xxz_xzz_0[j] =
                pa_x[j] * tpx_xz_xzz_0[j] + 0.5 * fl1_fx * tpx_z_xzz_0[j] + 0.5 * fl1_fx * tpx_xz_zz_0[j] - fl1_fgb * fl1_fx * ts_xz_xzz_0[j];

            tpy_xxz_xzz_0[j] = pa_x[j] * tpy_xz_xzz_0[j] + 0.5 * fl1_fx * tpy_z_xzz_0[j] + 0.5 * fl1_fx * tpy_xz_zz_0[j];

            tpz_xxz_xzz_0[j] = pa_x[j] * tpz_xz_xzz_0[j] + 0.5 * fl1_fx * tpz_z_xzz_0[j] + 0.5 * fl1_fx * tpz_xz_zz_0[j];

            tpx_xxz_yyy_0[j] = pa_x[j] * tpx_xz_yyy_0[j] + 0.5 * fl1_fx * tpx_z_yyy_0[j] - fl1_fgb * fl1_fx * ts_xz_yyy_0[j];

            tpy_xxz_yyy_0[j] = pa_x[j] * tpy_xz_yyy_0[j] + 0.5 * fl1_fx * tpy_z_yyy_0[j];

            tpz_xxz_yyy_0[j] = pa_x[j] * tpz_xz_yyy_0[j] + 0.5 * fl1_fx * tpz_z_yyy_0[j];

            tpx_xxz_yyz_0[j] = pa_x[j] * tpx_xz_yyz_0[j] + 0.5 * fl1_fx * tpx_z_yyz_0[j] - fl1_fgb * fl1_fx * ts_xz_yyz_0[j];

            tpy_xxz_yyz_0[j] = pa_x[j] * tpy_xz_yyz_0[j] + 0.5 * fl1_fx * tpy_z_yyz_0[j];

            tpz_xxz_yyz_0[j] = pa_x[j] * tpz_xz_yyz_0[j] + 0.5 * fl1_fx * tpz_z_yyz_0[j];

            tpx_xxz_yzz_0[j] = pa_x[j] * tpx_xz_yzz_0[j] + 0.5 * fl1_fx * tpx_z_yzz_0[j] - fl1_fgb * fl1_fx * ts_xz_yzz_0[j];

            tpy_xxz_yzz_0[j] = pa_x[j] * tpy_xz_yzz_0[j] + 0.5 * fl1_fx * tpy_z_yzz_0[j];

            tpz_xxz_yzz_0[j] = pa_x[j] * tpz_xz_yzz_0[j] + 0.5 * fl1_fx * tpz_z_yzz_0[j];

            tpx_xxz_zzz_0[j] = pa_x[j] * tpx_xz_zzz_0[j] + 0.5 * fl1_fx * tpx_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_xz_zzz_0[j];

            tpy_xxz_zzz_0[j] = pa_x[j] * tpy_xz_zzz_0[j] + 0.5 * fl1_fx * tpy_z_zzz_0[j];

            tpz_xxz_zzz_0[j] = pa_x[j] * tpz_xz_zzz_0[j] + 0.5 * fl1_fx * tpz_z_zzz_0[j];

            tpx_xyy_xxx_0[j] = pa_x[j] * tpx_yy_xxx_0[j] + 1.5 * fl1_fx * tpx_yy_xx_0[j] - fl1_fgb * fl1_fx * ts_yy_xxx_0[j];

            tpy_xyy_xxx_0[j] = pa_x[j] * tpy_yy_xxx_0[j] + 1.5 * fl1_fx * tpy_yy_xx_0[j];

            tpz_xyy_xxx_0[j] = pa_x[j] * tpz_yy_xxx_0[j] + 1.5 * fl1_fx * tpz_yy_xx_0[j];

            tpx_xyy_xxy_0[j] = pa_x[j] * tpx_yy_xxy_0[j] + fl1_fx * tpx_yy_xy_0[j] - fl1_fgb * fl1_fx * ts_yy_xxy_0[j];

            tpy_xyy_xxy_0[j] = pa_x[j] * tpy_yy_xxy_0[j] + fl1_fx * tpy_yy_xy_0[j];

            tpz_xyy_xxy_0[j] = pa_x[j] * tpz_yy_xxy_0[j] + fl1_fx * tpz_yy_xy_0[j];

            tpx_xyy_xxz_0[j] = pa_x[j] * tpx_yy_xxz_0[j] + fl1_fx * tpx_yy_xz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxz_0[j];

            tpy_xyy_xxz_0[j] = pa_x[j] * tpy_yy_xxz_0[j] + fl1_fx * tpy_yy_xz_0[j];

            tpz_xyy_xxz_0[j] = pa_x[j] * tpz_yy_xxz_0[j] + fl1_fx * tpz_yy_xz_0[j];

            tpx_xyy_xyy_0[j] = pa_x[j] * tpx_yy_xyy_0[j] + 0.5 * fl1_fx * tpx_yy_yy_0[j] - fl1_fgb * fl1_fx * ts_yy_xyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFF_100_150(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpy_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_yy_xyy_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 34);

        auto tpy_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_yy_xyz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 34);

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

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto ts_yy_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 34);

        auto ts_yy_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 35);

        auto ts_yy_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 36);

        auto ts_yy_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 37);

        auto ts_yy_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 38);

        auto ts_yy_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 39);

        auto ts_yz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 40);

        auto ts_yz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 41);

        auto ts_yz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 42);

        auto ts_yz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 43);

        auto ts_yz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 44);

        auto ts_yz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 45);

        auto ts_yz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 46);

        auto ts_yz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 47);

        auto ts_yz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 48);

        auto ts_yz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 49);

        // set up pointers to integrals

        auto tpy_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tpz_xyy_xyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 33);

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

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyy_xyz_0, tpx_xyy_xzz_0, tpx_xyy_yyy_0, \
                                     tpx_xyy_yyz_0, tpx_xyy_yzz_0, tpx_xyy_zzz_0, tpx_xyz_xxx_0, tpx_xyz_xxy_0, \
                                     tpx_xyz_xxz_0, tpx_xyz_xyy_0, tpx_xyz_xyz_0, tpx_xyz_xzz_0, tpx_xyz_yyy_0, \
                                     tpx_xyz_yyz_0, tpx_xyz_yzz_0, tpx_xyz_zzz_0, tpx_yy_xyz_0, tpx_yy_xzz_0, \
                                     tpx_yy_yyy_0, tpx_yy_yyz_0, tpx_yy_yz_0, tpx_yy_yzz_0, tpx_yy_zz_0, tpx_yy_zzz_0, \
                                     tpx_yz_xx_0, tpx_yz_xxx_0, tpx_yz_xxy_0, tpx_yz_xxz_0, tpx_yz_xy_0, tpx_yz_xyy_0, \
                                     tpx_yz_xyz_0, tpx_yz_xz_0, tpx_yz_xzz_0, tpx_yz_yy_0, tpx_yz_yyy_0, tpx_yz_yyz_0, \
                                     tpx_yz_yz_0, tpx_yz_yzz_0, tpx_yz_zz_0, tpx_yz_zzz_0, tpy_xyy_xyy_0, \
                                     tpy_xyy_xyz_0, tpy_xyy_xzz_0, tpy_xyy_yyy_0, tpy_xyy_yyz_0, tpy_xyy_yzz_0, \
                                     tpy_xyy_zzz_0, tpy_xyz_xxx_0, tpy_xyz_xxy_0, tpy_xyz_xxz_0, tpy_xyz_xyy_0, \
                                     tpy_xyz_xyz_0, tpy_xyz_xzz_0, tpy_xyz_yyy_0, tpy_xyz_yyz_0, tpy_xyz_yzz_0, \
                                     tpy_xyz_zzz_0, tpy_yy_xyy_0, tpy_yy_xyz_0, tpy_yy_xzz_0, tpy_yy_yy_0, tpy_yy_yyy_0, \
                                     tpy_yy_yyz_0, tpy_yy_yz_0, tpy_yy_yzz_0, tpy_yy_zz_0, tpy_yy_zzz_0, tpy_yz_xx_0, \
                                     tpy_yz_xxx_0, tpy_yz_xxy_0, tpy_yz_xxz_0, tpy_yz_xy_0, tpy_yz_xyy_0, tpy_yz_xyz_0, \
                                     tpy_yz_xz_0, tpy_yz_xzz_0, tpy_yz_yy_0, tpy_yz_yyy_0, tpy_yz_yyz_0, tpy_yz_yz_0, \
                                     tpy_yz_yzz_0, tpy_yz_zz_0, tpy_yz_zzz_0, tpz_xyy_xyy_0, tpz_xyy_xyz_0, \
                                     tpz_xyy_xzz_0, tpz_xyy_yyy_0, tpz_xyy_yyz_0, tpz_xyy_yzz_0, tpz_xyy_zzz_0, \
                                     tpz_xyz_xxx_0, tpz_xyz_xxy_0, tpz_xyz_xxz_0, tpz_xyz_xyy_0, tpz_xyz_xyz_0, \
                                     tpz_xyz_xzz_0, tpz_xyz_yyy_0, tpz_xyz_yyz_0, tpz_xyz_yzz_0, tpz_xyz_zzz_0, \
                                     tpz_yy_xyy_0, tpz_yy_xyz_0, tpz_yy_xzz_0, tpz_yy_yy_0, tpz_yy_yyy_0, tpz_yy_yyz_0, \
                                     tpz_yy_yz_0, tpz_yy_yzz_0, tpz_yy_zz_0, tpz_yy_zzz_0, tpz_yz_xx_0, tpz_yz_xxx_0, \
                                     tpz_yz_xxy_0, tpz_yz_xxz_0, tpz_yz_xy_0, tpz_yz_xyy_0, tpz_yz_xyz_0, tpz_yz_xz_0, \
                                     tpz_yz_xzz_0, tpz_yz_yy_0, tpz_yz_yyy_0, tpz_yz_yyz_0, tpz_yz_yz_0, tpz_yz_yzz_0, \
                                     tpz_yz_zz_0, tpz_yz_zzz_0, ts_yy_xyz_0, ts_yy_xzz_0, ts_yy_yyy_0, ts_yy_yyz_0, \
                                     ts_yy_yzz_0, ts_yy_zzz_0, ts_yz_xxx_0, ts_yz_xxy_0, ts_yz_xxz_0, ts_yz_xyy_0, \
                                     ts_yz_xyz_0, ts_yz_xzz_0, ts_yz_yyy_0, ts_yz_yyz_0, ts_yz_yzz_0, ts_yz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_xyy_xyy_0[j] = pa_x[j] * tpy_yy_xyy_0[j] + 0.5 * fl1_fx * tpy_yy_yy_0[j];

            tpz_xyy_xyy_0[j] = pa_x[j] * tpz_yy_xyy_0[j] + 0.5 * fl1_fx * tpz_yy_yy_0[j];

            tpx_xyy_xyz_0[j] = pa_x[j] * tpx_yy_xyz_0[j] + 0.5 * fl1_fx * tpx_yy_yz_0[j] - fl1_fgb * fl1_fx * ts_yy_xyz_0[j];

            tpy_xyy_xyz_0[j] = pa_x[j] * tpy_yy_xyz_0[j] + 0.5 * fl1_fx * tpy_yy_yz_0[j];

            tpz_xyy_xyz_0[j] = pa_x[j] * tpz_yy_xyz_0[j] + 0.5 * fl1_fx * tpz_yy_yz_0[j];

            tpx_xyy_xzz_0[j] = pa_x[j] * tpx_yy_xzz_0[j] + 0.5 * fl1_fx * tpx_yy_zz_0[j] - fl1_fgb * fl1_fx * ts_yy_xzz_0[j];

            tpy_xyy_xzz_0[j] = pa_x[j] * tpy_yy_xzz_0[j] + 0.5 * fl1_fx * tpy_yy_zz_0[j];

            tpz_xyy_xzz_0[j] = pa_x[j] * tpz_yy_xzz_0[j] + 0.5 * fl1_fx * tpz_yy_zz_0[j];

            tpx_xyy_yyy_0[j] = pa_x[j] * tpx_yy_yyy_0[j] - fl1_fgb * fl1_fx * ts_yy_yyy_0[j];

            tpy_xyy_yyy_0[j] = pa_x[j] * tpy_yy_yyy_0[j];

            tpz_xyy_yyy_0[j] = pa_x[j] * tpz_yy_yyy_0[j];

            tpx_xyy_yyz_0[j] = pa_x[j] * tpx_yy_yyz_0[j] - fl1_fgb * fl1_fx * ts_yy_yyz_0[j];

            tpy_xyy_yyz_0[j] = pa_x[j] * tpy_yy_yyz_0[j];

            tpz_xyy_yyz_0[j] = pa_x[j] * tpz_yy_yyz_0[j];

            tpx_xyy_yzz_0[j] = pa_x[j] * tpx_yy_yzz_0[j] - fl1_fgb * fl1_fx * ts_yy_yzz_0[j];

            tpy_xyy_yzz_0[j] = pa_x[j] * tpy_yy_yzz_0[j];

            tpz_xyy_yzz_0[j] = pa_x[j] * tpz_yy_yzz_0[j];

            tpx_xyy_zzz_0[j] = pa_x[j] * tpx_yy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yy_zzz_0[j];

            tpy_xyy_zzz_0[j] = pa_x[j] * tpy_yy_zzz_0[j];

            tpz_xyy_zzz_0[j] = pa_x[j] * tpz_yy_zzz_0[j];

            tpx_xyz_xxx_0[j] = pa_x[j] * tpx_yz_xxx_0[j] + 1.5 * fl1_fx * tpx_yz_xx_0[j] - fl1_fgb * fl1_fx * ts_yz_xxx_0[j];

            tpy_xyz_xxx_0[j] = pa_x[j] * tpy_yz_xxx_0[j] + 1.5 * fl1_fx * tpy_yz_xx_0[j];

            tpz_xyz_xxx_0[j] = pa_x[j] * tpz_yz_xxx_0[j] + 1.5 * fl1_fx * tpz_yz_xx_0[j];

            tpx_xyz_xxy_0[j] = pa_x[j] * tpx_yz_xxy_0[j] + fl1_fx * tpx_yz_xy_0[j] - fl1_fgb * fl1_fx * ts_yz_xxy_0[j];

            tpy_xyz_xxy_0[j] = pa_x[j] * tpy_yz_xxy_0[j] + fl1_fx * tpy_yz_xy_0[j];

            tpz_xyz_xxy_0[j] = pa_x[j] * tpz_yz_xxy_0[j] + fl1_fx * tpz_yz_xy_0[j];

            tpx_xyz_xxz_0[j] = pa_x[j] * tpx_yz_xxz_0[j] + fl1_fx * tpx_yz_xz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxz_0[j];

            tpy_xyz_xxz_0[j] = pa_x[j] * tpy_yz_xxz_0[j] + fl1_fx * tpy_yz_xz_0[j];

            tpz_xyz_xxz_0[j] = pa_x[j] * tpz_yz_xxz_0[j] + fl1_fx * tpz_yz_xz_0[j];

            tpx_xyz_xyy_0[j] = pa_x[j] * tpx_yz_xyy_0[j] + 0.5 * fl1_fx * tpx_yz_yy_0[j] - fl1_fgb * fl1_fx * ts_yz_xyy_0[j];

            tpy_xyz_xyy_0[j] = pa_x[j] * tpy_yz_xyy_0[j] + 0.5 * fl1_fx * tpy_yz_yy_0[j];

            tpz_xyz_xyy_0[j] = pa_x[j] * tpz_yz_xyy_0[j] + 0.5 * fl1_fx * tpz_yz_yy_0[j];

            tpx_xyz_xyz_0[j] = pa_x[j] * tpx_yz_xyz_0[j] + 0.5 * fl1_fx * tpx_yz_yz_0[j] - fl1_fgb * fl1_fx * ts_yz_xyz_0[j];

            tpy_xyz_xyz_0[j] = pa_x[j] * tpy_yz_xyz_0[j] + 0.5 * fl1_fx * tpy_yz_yz_0[j];

            tpz_xyz_xyz_0[j] = pa_x[j] * tpz_yz_xyz_0[j] + 0.5 * fl1_fx * tpz_yz_yz_0[j];

            tpx_xyz_xzz_0[j] = pa_x[j] * tpx_yz_xzz_0[j] + 0.5 * fl1_fx * tpx_yz_zz_0[j] - fl1_fgb * fl1_fx * ts_yz_xzz_0[j];

            tpy_xyz_xzz_0[j] = pa_x[j] * tpy_yz_xzz_0[j] + 0.5 * fl1_fx * tpy_yz_zz_0[j];

            tpz_xyz_xzz_0[j] = pa_x[j] * tpz_yz_xzz_0[j] + 0.5 * fl1_fx * tpz_yz_zz_0[j];

            tpx_xyz_yyy_0[j] = pa_x[j] * tpx_yz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yz_yyy_0[j];

            tpy_xyz_yyy_0[j] = pa_x[j] * tpy_yz_yyy_0[j];

            tpz_xyz_yyy_0[j] = pa_x[j] * tpz_yz_yyy_0[j];

            tpx_xyz_yyz_0[j] = pa_x[j] * tpx_yz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yz_yyz_0[j];

            tpy_xyz_yyz_0[j] = pa_x[j] * tpy_yz_yyz_0[j];

            tpz_xyz_yyz_0[j] = pa_x[j] * tpz_yz_yyz_0[j];

            tpx_xyz_yzz_0[j] = pa_x[j] * tpx_yz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yz_yzz_0[j];

            tpy_xyz_yzz_0[j] = pa_x[j] * tpy_yz_yzz_0[j];

            tpz_xyz_yzz_0[j] = pa_x[j] * tpz_yz_yzz_0[j];

            tpx_xyz_zzz_0[j] = pa_x[j] * tpx_yz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yz_zzz_0[j];

            tpy_xyz_zzz_0[j] = pa_x[j] * tpy_yz_zzz_0[j];

            tpz_xyz_zzz_0[j] = pa_x[j] * tpz_yz_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFF_150_200(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 36);

        auto tpy_yy_yyy_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 36);

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

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpx_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 18);

        auto tpy_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 18);

        auto tpz_yy_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 18);

        auto tpx_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 19);

        auto tpy_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 19);

        auto tpz_yy_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 19);

        auto tpx_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 20);

        auto tpy_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 20);

        auto tpz_yy_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 20);

        auto tpx_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 21);

        auto tpy_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 21);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto ts_yy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 30);

        auto ts_yy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 31);

        auto ts_yy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 32);

        auto ts_yy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 33);

        auto ts_yy_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 34);

        auto ts_yy_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 35);

        auto ts_yy_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 36);

        auto ts_zz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 50);

        auto ts_zz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 51);

        auto ts_zz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 52);

        auto ts_zz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 53);

        auto ts_zz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 54);

        auto ts_zz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 55);

        auto ts_zz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 56);

        auto ts_zz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 57);

        auto ts_zz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 58);

        auto ts_zz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 59);

        // set up pointers to integrals

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

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_xzz_xxx_0, tpx_xzz_xxy_0, tpx_xzz_xxz_0, \
                                     tpx_xzz_xyy_0, tpx_xzz_xyz_0, tpx_xzz_xzz_0, tpx_xzz_yyy_0, tpx_xzz_yyz_0, \
                                     tpx_xzz_yzz_0, tpx_xzz_zzz_0, tpx_y_xxx_0, tpx_y_xxy_0, tpx_y_xxz_0, tpx_y_xyy_0, \
                                     tpx_y_xyz_0, tpx_y_xzz_0, tpx_y_yyy_0, tpx_yy_xx_0, tpx_yy_xxx_0, tpx_yy_xxy_0, \
                                     tpx_yy_xxz_0, tpx_yy_xy_0, tpx_yy_xyy_0, tpx_yy_xyz_0, tpx_yy_xz_0, tpx_yy_xzz_0, \
                                     tpx_yy_yy_0, tpx_yy_yyy_0, tpx_yyy_xxx_0, tpx_yyy_xxy_0, tpx_yyy_xxz_0, \
                                     tpx_yyy_xyy_0, tpx_yyy_xyz_0, tpx_yyy_xzz_0, tpx_yyy_yyy_0, tpx_zz_xx_0, \
                                     tpx_zz_xxx_0, tpx_zz_xxy_0, tpx_zz_xxz_0, tpx_zz_xy_0, tpx_zz_xyy_0, tpx_zz_xyz_0, \
                                     tpx_zz_xz_0, tpx_zz_xzz_0, tpx_zz_yy_0, tpx_zz_yyy_0, tpx_zz_yyz_0, tpx_zz_yz_0, \
                                     tpx_zz_yzz_0, tpx_zz_zz_0, tpx_zz_zzz_0, tpy_xzz_xxx_0, tpy_xzz_xxy_0, \
                                     tpy_xzz_xxz_0, tpy_xzz_xyy_0, tpy_xzz_xyz_0, tpy_xzz_xzz_0, tpy_xzz_yyy_0, \
                                     tpy_xzz_yyz_0, tpy_xzz_yzz_0, tpy_xzz_zzz_0, tpy_y_xxx_0, tpy_y_xxy_0, tpy_y_xxz_0, \
                                     tpy_y_xyy_0, tpy_y_xyz_0, tpy_y_xzz_0, tpy_y_yyy_0, tpy_yy_xx_0, tpy_yy_xxx_0, \
                                     tpy_yy_xxy_0, tpy_yy_xxz_0, tpy_yy_xy_0, tpy_yy_xyy_0, tpy_yy_xyz_0, tpy_yy_xz_0, \
                                     tpy_yy_xzz_0, tpy_yy_yy_0, tpy_yy_yyy_0, tpy_yyy_xxx_0, tpy_yyy_xxy_0, \
                                     tpy_yyy_xxz_0, tpy_yyy_xyy_0, tpy_yyy_xyz_0, tpy_yyy_xzz_0, tpy_yyy_yyy_0, \
                                     tpy_zz_xx_0, tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, tpy_zz_xy_0, tpy_zz_xyy_0, \
                                     tpy_zz_xyz_0, tpy_zz_xz_0, tpy_zz_xzz_0, tpy_zz_yy_0, tpy_zz_yyy_0, tpy_zz_yyz_0, \
                                     tpy_zz_yz_0, tpy_zz_yzz_0, tpy_zz_zz_0, tpy_zz_zzz_0, tpz_xzz_xxx_0, \
                                     tpz_xzz_xxy_0, tpz_xzz_xxz_0, tpz_xzz_xyy_0, tpz_xzz_xyz_0, tpz_xzz_xzz_0, \
                                     tpz_xzz_yyy_0, tpz_xzz_yyz_0, tpz_xzz_yzz_0, tpz_xzz_zzz_0, tpz_y_xxx_0, \
                                     tpz_y_xxy_0, tpz_y_xxz_0, tpz_y_xyy_0, tpz_y_xyz_0, tpz_y_xzz_0, tpz_yy_xx_0, \
                                     tpz_yy_xxx_0, tpz_yy_xxy_0, tpz_yy_xxz_0, tpz_yy_xy_0, tpz_yy_xyy_0, tpz_yy_xyz_0, \
                                     tpz_yy_xz_0, tpz_yy_xzz_0, tpz_yyy_xxx_0, tpz_yyy_xxy_0, tpz_yyy_xxz_0, \
                                     tpz_yyy_xyy_0, tpz_yyy_xyz_0, tpz_yyy_xzz_0, tpz_zz_xx_0, tpz_zz_xxx_0, \
                                     tpz_zz_xxy_0, tpz_zz_xxz_0, tpz_zz_xy_0, tpz_zz_xyy_0, tpz_zz_xyz_0, tpz_zz_xz_0, \
                                     tpz_zz_xzz_0, tpz_zz_yy_0, tpz_zz_yyy_0, tpz_zz_yyz_0, tpz_zz_yz_0, tpz_zz_yzz_0, \
                                     tpz_zz_zz_0, tpz_zz_zzz_0, ts_yy_xxx_0, ts_yy_xxy_0, ts_yy_xxz_0, ts_yy_xyy_0, \
                                     ts_yy_xyz_0, ts_yy_xzz_0, ts_yy_yyy_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0, \
                                     ts_zz_xyy_0, ts_zz_xyz_0, ts_zz_xzz_0, ts_zz_yyy_0, ts_zz_yyz_0, ts_zz_yzz_0, \
                                     ts_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xzz_xxx_0[j] = pa_x[j] * tpx_zz_xxx_0[j] + 1.5 * fl1_fx * tpx_zz_xx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxx_0[j];

            tpy_xzz_xxx_0[j] = pa_x[j] * tpy_zz_xxx_0[j] + 1.5 * fl1_fx * tpy_zz_xx_0[j];

            tpz_xzz_xxx_0[j] = pa_x[j] * tpz_zz_xxx_0[j] + 1.5 * fl1_fx * tpz_zz_xx_0[j];

            tpx_xzz_xxy_0[j] = pa_x[j] * tpx_zz_xxy_0[j] + fl1_fx * tpx_zz_xy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxy_0[j];

            tpy_xzz_xxy_0[j] = pa_x[j] * tpy_zz_xxy_0[j] + fl1_fx * tpy_zz_xy_0[j];

            tpz_xzz_xxy_0[j] = pa_x[j] * tpz_zz_xxy_0[j] + fl1_fx * tpz_zz_xy_0[j];

            tpx_xzz_xxz_0[j] = pa_x[j] * tpx_zz_xxz_0[j] + fl1_fx * tpx_zz_xz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxz_0[j];

            tpy_xzz_xxz_0[j] = pa_x[j] * tpy_zz_xxz_0[j] + fl1_fx * tpy_zz_xz_0[j];

            tpz_xzz_xxz_0[j] = pa_x[j] * tpz_zz_xxz_0[j] + fl1_fx * tpz_zz_xz_0[j];

            tpx_xzz_xyy_0[j] = pa_x[j] * tpx_zz_xyy_0[j] + 0.5 * fl1_fx * tpx_zz_yy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyy_0[j];

            tpy_xzz_xyy_0[j] = pa_x[j] * tpy_zz_xyy_0[j] + 0.5 * fl1_fx * tpy_zz_yy_0[j];

            tpz_xzz_xyy_0[j] = pa_x[j] * tpz_zz_xyy_0[j] + 0.5 * fl1_fx * tpz_zz_yy_0[j];

            tpx_xzz_xyz_0[j] = pa_x[j] * tpx_zz_xyz_0[j] + 0.5 * fl1_fx * tpx_zz_yz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyz_0[j];

            tpy_xzz_xyz_0[j] = pa_x[j] * tpy_zz_xyz_0[j] + 0.5 * fl1_fx * tpy_zz_yz_0[j];

            tpz_xzz_xyz_0[j] = pa_x[j] * tpz_zz_xyz_0[j] + 0.5 * fl1_fx * tpz_zz_yz_0[j];

            tpx_xzz_xzz_0[j] = pa_x[j] * tpx_zz_xzz_0[j] + 0.5 * fl1_fx * tpx_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_zz_xzz_0[j];

            tpy_xzz_xzz_0[j] = pa_x[j] * tpy_zz_xzz_0[j] + 0.5 * fl1_fx * tpy_zz_zz_0[j];

            tpz_xzz_xzz_0[j] = pa_x[j] * tpz_zz_xzz_0[j] + 0.5 * fl1_fx * tpz_zz_zz_0[j];

            tpx_xzz_yyy_0[j] = pa_x[j] * tpx_zz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyy_0[j];

            tpy_xzz_yyy_0[j] = pa_x[j] * tpy_zz_yyy_0[j];

            tpz_xzz_yyy_0[j] = pa_x[j] * tpz_zz_yyy_0[j];

            tpx_xzz_yyz_0[j] = pa_x[j] * tpx_zz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyz_0[j];

            tpy_xzz_yyz_0[j] = pa_x[j] * tpy_zz_yyz_0[j];

            tpz_xzz_yyz_0[j] = pa_x[j] * tpz_zz_yyz_0[j];

            tpx_xzz_yzz_0[j] = pa_x[j] * tpx_zz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zz_yzz_0[j];

            tpy_xzz_yzz_0[j] = pa_x[j] * tpy_zz_yzz_0[j];

            tpz_xzz_yzz_0[j] = pa_x[j] * tpz_zz_yzz_0[j];

            tpx_xzz_zzz_0[j] = pa_x[j] * tpx_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zz_zzz_0[j];

            tpy_xzz_zzz_0[j] = pa_x[j] * tpy_zz_zzz_0[j];

            tpz_xzz_zzz_0[j] = pa_x[j] * tpz_zz_zzz_0[j];

            tpx_yyy_xxx_0[j] = pa_y[j] * tpx_yy_xxx_0[j] + fl1_fx * tpx_y_xxx_0[j];

            tpy_yyy_xxx_0[j] = pa_y[j] * tpy_yy_xxx_0[j] + fl1_fx * tpy_y_xxx_0[j] - fl1_fgb * fl1_fx * ts_yy_xxx_0[j];

            tpz_yyy_xxx_0[j] = pa_y[j] * tpz_yy_xxx_0[j] + fl1_fx * tpz_y_xxx_0[j];

            tpx_yyy_xxy_0[j] = pa_y[j] * tpx_yy_xxy_0[j] + fl1_fx * tpx_y_xxy_0[j] + 0.5 * fl1_fx * tpx_yy_xx_0[j];

            tpy_yyy_xxy_0[j] =
                pa_y[j] * tpy_yy_xxy_0[j] + fl1_fx * tpy_y_xxy_0[j] + 0.5 * fl1_fx * tpy_yy_xx_0[j] - fl1_fgb * fl1_fx * ts_yy_xxy_0[j];

            tpz_yyy_xxy_0[j] = pa_y[j] * tpz_yy_xxy_0[j] + fl1_fx * tpz_y_xxy_0[j] + 0.5 * fl1_fx * tpz_yy_xx_0[j];

            tpx_yyy_xxz_0[j] = pa_y[j] * tpx_yy_xxz_0[j] + fl1_fx * tpx_y_xxz_0[j];

            tpy_yyy_xxz_0[j] = pa_y[j] * tpy_yy_xxz_0[j] + fl1_fx * tpy_y_xxz_0[j] - fl1_fgb * fl1_fx * ts_yy_xxz_0[j];

            tpz_yyy_xxz_0[j] = pa_y[j] * tpz_yy_xxz_0[j] + fl1_fx * tpz_y_xxz_0[j];

            tpx_yyy_xyy_0[j] = pa_y[j] * tpx_yy_xyy_0[j] + fl1_fx * tpx_y_xyy_0[j] + fl1_fx * tpx_yy_xy_0[j];

            tpy_yyy_xyy_0[j] = pa_y[j] * tpy_yy_xyy_0[j] + fl1_fx * tpy_y_xyy_0[j] + fl1_fx * tpy_yy_xy_0[j] - fl1_fgb * fl1_fx * ts_yy_xyy_0[j];

            tpz_yyy_xyy_0[j] = pa_y[j] * tpz_yy_xyy_0[j] + fl1_fx * tpz_y_xyy_0[j] + fl1_fx * tpz_yy_xy_0[j];

            tpx_yyy_xyz_0[j] = pa_y[j] * tpx_yy_xyz_0[j] + fl1_fx * tpx_y_xyz_0[j] + 0.5 * fl1_fx * tpx_yy_xz_0[j];

            tpy_yyy_xyz_0[j] =
                pa_y[j] * tpy_yy_xyz_0[j] + fl1_fx * tpy_y_xyz_0[j] + 0.5 * fl1_fx * tpy_yy_xz_0[j] - fl1_fgb * fl1_fx * ts_yy_xyz_0[j];

            tpz_yyy_xyz_0[j] = pa_y[j] * tpz_yy_xyz_0[j] + fl1_fx * tpz_y_xyz_0[j] + 0.5 * fl1_fx * tpz_yy_xz_0[j];

            tpx_yyy_xzz_0[j] = pa_y[j] * tpx_yy_xzz_0[j] + fl1_fx * tpx_y_xzz_0[j];

            tpy_yyy_xzz_0[j] = pa_y[j] * tpy_yy_xzz_0[j] + fl1_fx * tpy_y_xzz_0[j] - fl1_fgb * fl1_fx * ts_yy_xzz_0[j];

            tpz_yyy_xzz_0[j] = pa_y[j] * tpz_yy_xzz_0[j] + fl1_fx * tpz_y_xzz_0[j];

            tpx_yyy_yyy_0[j] = pa_y[j] * tpx_yy_yyy_0[j] + fl1_fx * tpx_y_yyy_0[j] + 1.5 * fl1_fx * tpx_yy_yy_0[j];

            tpy_yyy_yyy_0[j] =
                pa_y[j] * tpy_yy_yyy_0[j] + fl1_fx * tpy_y_yyy_0[j] + 1.5 * fl1_fx * tpy_yy_yy_0[j] - fl1_fgb * fl1_fx * ts_yy_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFF_200_250(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tpz_yy_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 21);

        auto tpx_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 22);

        auto tpy_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 22);

        auto tpz_yy_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 22);

        auto tpx_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 23);

        auto tpy_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 23);

        auto tpz_yy_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 23);

        auto tpx_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 24);

        auto tpy_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 24);

        auto tpz_yz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 24);

        auto tpx_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 25);

        auto tpy_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 25);

        auto tpz_yz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 25);

        auto tpx_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 26);

        auto tpy_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 26);

        auto tpz_yz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 26);

        auto tpx_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 27);

        auto tpy_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 27);

        auto tpz_yz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 27);

        auto tpx_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 28);

        auto tpy_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 28);

        auto tpz_yz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 28);

        auto tpx_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 29);

        auto tpy_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 29);

        auto tpz_yz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 29);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto ts_yy_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 37);

        auto ts_yy_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 38);

        auto ts_yy_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 39);

        auto ts_yz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 40);

        auto ts_yz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 41);

        auto ts_yz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 42);

        auto ts_yz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 43);

        auto ts_yz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 44);

        auto ts_yz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 45);

        auto ts_yz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 46);

        auto ts_yz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 47);

        auto ts_yz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 48);

        auto ts_yz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 49);

        auto ts_zz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 50);

        auto ts_zz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 51);

        auto ts_zz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 52);

        // set up pointers to integrals

        auto tpz_yyy_yyy_0 = primBuffer.data(pidx_p_3_3_m0 + 200 * bdim + 100 * idx + 66);

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

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_y_yyz_0, tpx_y_yzz_0, tpx_y_zzz_0, tpx_yy_yyz_0, \
                                     tpx_yy_yz_0, tpx_yy_yzz_0, tpx_yy_zz_0, tpx_yy_zzz_0, tpx_yyy_yyz_0, \
                                     tpx_yyy_yzz_0, tpx_yyy_zzz_0, tpx_yyz_xxx_0, tpx_yyz_xxy_0, tpx_yyz_xxz_0, \
                                     tpx_yyz_xyy_0, tpx_yyz_xyz_0, tpx_yyz_xzz_0, tpx_yyz_yyy_0, tpx_yyz_yyz_0, \
                                     tpx_yyz_yzz_0, tpx_yyz_zzz_0, tpx_yz_xx_0, tpx_yz_xxx_0, tpx_yz_xxy_0, tpx_yz_xxz_0, \
                                     tpx_yz_xy_0, tpx_yz_xyy_0, tpx_yz_xyz_0, tpx_yz_xz_0, tpx_yz_xzz_0, tpx_yz_yy_0, \
                                     tpx_yz_yyy_0, tpx_yz_yyz_0, tpx_yz_yz_0, tpx_yz_yzz_0, tpx_yz_zz_0, tpx_yz_zzz_0, \
                                     tpx_yzz_xxx_0, tpx_yzz_xxy_0, tpx_yzz_xxz_0, tpx_yzz_xyy_0, tpx_z_xxx_0, \
                                     tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xyy_0, tpx_z_xyz_0, tpx_z_xzz_0, tpx_z_yyy_0, \
                                     tpx_z_yyz_0, tpx_z_yzz_0, tpx_z_zzz_0, tpx_zz_xx_0, tpx_zz_xxx_0, tpx_zz_xxy_0, \
                                     tpx_zz_xxz_0, tpx_zz_xy_0, tpx_zz_xyy_0, tpy_y_yyz_0, tpy_y_yzz_0, tpy_y_zzz_0, \
                                     tpy_yy_yyz_0, tpy_yy_yz_0, tpy_yy_yzz_0, tpy_yy_zz_0, tpy_yy_zzz_0, tpy_yyy_yyz_0, \
                                     tpy_yyy_yzz_0, tpy_yyy_zzz_0, tpy_yyz_xxx_0, tpy_yyz_xxy_0, tpy_yyz_xxz_0, \
                                     tpy_yyz_xyy_0, tpy_yyz_xyz_0, tpy_yyz_xzz_0, tpy_yyz_yyy_0, tpy_yyz_yyz_0, \
                                     tpy_yyz_yzz_0, tpy_yyz_zzz_0, tpy_yz_xx_0, tpy_yz_xxx_0, tpy_yz_xxy_0, tpy_yz_xxz_0, \
                                     tpy_yz_xy_0, tpy_yz_xyy_0, tpy_yz_xyz_0, tpy_yz_xz_0, tpy_yz_xzz_0, tpy_yz_yy_0, \
                                     tpy_yz_yyy_0, tpy_yz_yyz_0, tpy_yz_yz_0, tpy_yz_yzz_0, tpy_yz_zz_0, tpy_yz_zzz_0, \
                                     tpy_yzz_xxx_0, tpy_yzz_xxy_0, tpy_yzz_xxz_0, tpy_z_xxx_0, tpy_z_xxy_0, tpy_z_xxz_0, \
                                     tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xzz_0, tpy_z_yyy_0, tpy_z_yyz_0, tpy_z_yzz_0, \
                                     tpy_z_zzz_0, tpy_zz_xx_0, tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, tpz_y_yyy_0, \
                                     tpz_y_yyz_0, tpz_y_yzz_0, tpz_y_zzz_0, tpz_yy_yy_0, tpz_yy_yyy_0, tpz_yy_yyz_0, \
                                     tpz_yy_yz_0, tpz_yy_yzz_0, tpz_yy_zz_0, tpz_yy_zzz_0, tpz_yyy_yyy_0, \
                                     tpz_yyy_yyz_0, tpz_yyy_yzz_0, tpz_yyy_zzz_0, tpz_yyz_xxx_0, tpz_yyz_xxy_0, \
                                     tpz_yyz_xxz_0, tpz_yyz_xyy_0, tpz_yyz_xyz_0, tpz_yyz_xzz_0, tpz_yyz_yyy_0, \
                                     tpz_yyz_yyz_0, tpz_yyz_yzz_0, tpz_yyz_zzz_0, tpz_yz_xx_0, tpz_yz_xxx_0, \
                                     tpz_yz_xxy_0, tpz_yz_xxz_0, tpz_yz_xy_0, tpz_yz_xyy_0, tpz_yz_xyz_0, tpz_yz_xz_0, \
                                     tpz_yz_xzz_0, tpz_yz_yy_0, tpz_yz_yyy_0, tpz_yz_yyz_0, tpz_yz_yz_0, tpz_yz_yzz_0, \
                                     tpz_yz_zz_0, tpz_yz_zzz_0, tpz_yzz_xxx_0, tpz_yzz_xxy_0, tpz_yzz_xxz_0, \
                                     tpz_z_xxx_0, tpz_z_xxy_0, tpz_z_xxz_0, tpz_z_xyy_0, tpz_z_xyz_0, tpz_z_xzz_0, \
                                     tpz_z_yyy_0, tpz_z_yyz_0, tpz_z_yzz_0, tpz_z_zzz_0, tpz_zz_xx_0, tpz_zz_xxx_0, \
                                     tpz_zz_xxy_0, tpz_zz_xxz_0, ts_yy_yyz_0, ts_yy_yzz_0, ts_yy_zzz_0, ts_yz_xxx_0, \
                                     ts_yz_xxy_0, ts_yz_xxz_0, ts_yz_xyy_0, ts_yz_xyz_0, ts_yz_xzz_0, ts_yz_yyy_0, \
                                     ts_yz_yyz_0, ts_yz_yzz_0, ts_yz_zzz_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_yyy_yyy_0[j] = pa_y[j] * tpz_yy_yyy_0[j] + fl1_fx * tpz_y_yyy_0[j] + 1.5 * fl1_fx * tpz_yy_yy_0[j];

            tpx_yyy_yyz_0[j] = pa_y[j] * tpx_yy_yyz_0[j] + fl1_fx * tpx_y_yyz_0[j] + fl1_fx * tpx_yy_yz_0[j];

            tpy_yyy_yyz_0[j] = pa_y[j] * tpy_yy_yyz_0[j] + fl1_fx * tpy_y_yyz_0[j] + fl1_fx * tpy_yy_yz_0[j] - fl1_fgb * fl1_fx * ts_yy_yyz_0[j];

            tpz_yyy_yyz_0[j] = pa_y[j] * tpz_yy_yyz_0[j] + fl1_fx * tpz_y_yyz_0[j] + fl1_fx * tpz_yy_yz_0[j];

            tpx_yyy_yzz_0[j] = pa_y[j] * tpx_yy_yzz_0[j] + fl1_fx * tpx_y_yzz_0[j] + 0.5 * fl1_fx * tpx_yy_zz_0[j];

            tpy_yyy_yzz_0[j] =
                pa_y[j] * tpy_yy_yzz_0[j] + fl1_fx * tpy_y_yzz_0[j] + 0.5 * fl1_fx * tpy_yy_zz_0[j] - fl1_fgb * fl1_fx * ts_yy_yzz_0[j];

            tpz_yyy_yzz_0[j] = pa_y[j] * tpz_yy_yzz_0[j] + fl1_fx * tpz_y_yzz_0[j] + 0.5 * fl1_fx * tpz_yy_zz_0[j];

            tpx_yyy_zzz_0[j] = pa_y[j] * tpx_yy_zzz_0[j] + fl1_fx * tpx_y_zzz_0[j];

            tpy_yyy_zzz_0[j] = pa_y[j] * tpy_yy_zzz_0[j] + fl1_fx * tpy_y_zzz_0[j] - fl1_fgb * fl1_fx * ts_yy_zzz_0[j];

            tpz_yyy_zzz_0[j] = pa_y[j] * tpz_yy_zzz_0[j] + fl1_fx * tpz_y_zzz_0[j];

            tpx_yyz_xxx_0[j] = pa_y[j] * tpx_yz_xxx_0[j] + 0.5 * fl1_fx * tpx_z_xxx_0[j];

            tpy_yyz_xxx_0[j] = pa_y[j] * tpy_yz_xxx_0[j] + 0.5 * fl1_fx * tpy_z_xxx_0[j] - fl1_fgb * fl1_fx * ts_yz_xxx_0[j];

            tpz_yyz_xxx_0[j] = pa_y[j] * tpz_yz_xxx_0[j] + 0.5 * fl1_fx * tpz_z_xxx_0[j];

            tpx_yyz_xxy_0[j] = pa_y[j] * tpx_yz_xxy_0[j] + 0.5 * fl1_fx * tpx_z_xxy_0[j] + 0.5 * fl1_fx * tpx_yz_xx_0[j];

            tpy_yyz_xxy_0[j] =
                pa_y[j] * tpy_yz_xxy_0[j] + 0.5 * fl1_fx * tpy_z_xxy_0[j] + 0.5 * fl1_fx * tpy_yz_xx_0[j] - fl1_fgb * fl1_fx * ts_yz_xxy_0[j];

            tpz_yyz_xxy_0[j] = pa_y[j] * tpz_yz_xxy_0[j] + 0.5 * fl1_fx * tpz_z_xxy_0[j] + 0.5 * fl1_fx * tpz_yz_xx_0[j];

            tpx_yyz_xxz_0[j] = pa_y[j] * tpx_yz_xxz_0[j] + 0.5 * fl1_fx * tpx_z_xxz_0[j];

            tpy_yyz_xxz_0[j] = pa_y[j] * tpy_yz_xxz_0[j] + 0.5 * fl1_fx * tpy_z_xxz_0[j] - fl1_fgb * fl1_fx * ts_yz_xxz_0[j];

            tpz_yyz_xxz_0[j] = pa_y[j] * tpz_yz_xxz_0[j] + 0.5 * fl1_fx * tpz_z_xxz_0[j];

            tpx_yyz_xyy_0[j] = pa_y[j] * tpx_yz_xyy_0[j] + 0.5 * fl1_fx * tpx_z_xyy_0[j] + fl1_fx * tpx_yz_xy_0[j];

            tpy_yyz_xyy_0[j] =
                pa_y[j] * tpy_yz_xyy_0[j] + 0.5 * fl1_fx * tpy_z_xyy_0[j] + fl1_fx * tpy_yz_xy_0[j] - fl1_fgb * fl1_fx * ts_yz_xyy_0[j];

            tpz_yyz_xyy_0[j] = pa_y[j] * tpz_yz_xyy_0[j] + 0.5 * fl1_fx * tpz_z_xyy_0[j] + fl1_fx * tpz_yz_xy_0[j];

            tpx_yyz_xyz_0[j] = pa_y[j] * tpx_yz_xyz_0[j] + 0.5 * fl1_fx * tpx_z_xyz_0[j] + 0.5 * fl1_fx * tpx_yz_xz_0[j];

            tpy_yyz_xyz_0[j] =
                pa_y[j] * tpy_yz_xyz_0[j] + 0.5 * fl1_fx * tpy_z_xyz_0[j] + 0.5 * fl1_fx * tpy_yz_xz_0[j] - fl1_fgb * fl1_fx * ts_yz_xyz_0[j];

            tpz_yyz_xyz_0[j] = pa_y[j] * tpz_yz_xyz_0[j] + 0.5 * fl1_fx * tpz_z_xyz_0[j] + 0.5 * fl1_fx * tpz_yz_xz_0[j];

            tpx_yyz_xzz_0[j] = pa_y[j] * tpx_yz_xzz_0[j] + 0.5 * fl1_fx * tpx_z_xzz_0[j];

            tpy_yyz_xzz_0[j] = pa_y[j] * tpy_yz_xzz_0[j] + 0.5 * fl1_fx * tpy_z_xzz_0[j] - fl1_fgb * fl1_fx * ts_yz_xzz_0[j];

            tpz_yyz_xzz_0[j] = pa_y[j] * tpz_yz_xzz_0[j] + 0.5 * fl1_fx * tpz_z_xzz_0[j];

            tpx_yyz_yyy_0[j] = pa_y[j] * tpx_yz_yyy_0[j] + 0.5 * fl1_fx * tpx_z_yyy_0[j] + 1.5 * fl1_fx * tpx_yz_yy_0[j];

            tpy_yyz_yyy_0[j] =
                pa_y[j] * tpy_yz_yyy_0[j] + 0.5 * fl1_fx * tpy_z_yyy_0[j] + 1.5 * fl1_fx * tpy_yz_yy_0[j] - fl1_fgb * fl1_fx * ts_yz_yyy_0[j];

            tpz_yyz_yyy_0[j] = pa_y[j] * tpz_yz_yyy_0[j] + 0.5 * fl1_fx * tpz_z_yyy_0[j] + 1.5 * fl1_fx * tpz_yz_yy_0[j];

            tpx_yyz_yyz_0[j] = pa_y[j] * tpx_yz_yyz_0[j] + 0.5 * fl1_fx * tpx_z_yyz_0[j] + fl1_fx * tpx_yz_yz_0[j];

            tpy_yyz_yyz_0[j] =
                pa_y[j] * tpy_yz_yyz_0[j] + 0.5 * fl1_fx * tpy_z_yyz_0[j] + fl1_fx * tpy_yz_yz_0[j] - fl1_fgb * fl1_fx * ts_yz_yyz_0[j];

            tpz_yyz_yyz_0[j] = pa_y[j] * tpz_yz_yyz_0[j] + 0.5 * fl1_fx * tpz_z_yyz_0[j] + fl1_fx * tpz_yz_yz_0[j];

            tpx_yyz_yzz_0[j] = pa_y[j] * tpx_yz_yzz_0[j] + 0.5 * fl1_fx * tpx_z_yzz_0[j] + 0.5 * fl1_fx * tpx_yz_zz_0[j];

            tpy_yyz_yzz_0[j] =
                pa_y[j] * tpy_yz_yzz_0[j] + 0.5 * fl1_fx * tpy_z_yzz_0[j] + 0.5 * fl1_fx * tpy_yz_zz_0[j] - fl1_fgb * fl1_fx * ts_yz_yzz_0[j];

            tpz_yyz_yzz_0[j] = pa_y[j] * tpz_yz_yzz_0[j] + 0.5 * fl1_fx * tpz_z_yzz_0[j] + 0.5 * fl1_fx * tpz_yz_zz_0[j];

            tpx_yyz_zzz_0[j] = pa_y[j] * tpx_yz_zzz_0[j] + 0.5 * fl1_fx * tpx_z_zzz_0[j];

            tpy_yyz_zzz_0[j] = pa_y[j] * tpy_yz_zzz_0[j] + 0.5 * fl1_fx * tpy_z_zzz_0[j] - fl1_fgb * fl1_fx * ts_yz_zzz_0[j];

            tpz_yyz_zzz_0[j] = pa_y[j] * tpz_yz_zzz_0[j] + 0.5 * fl1_fx * tpz_z_zzz_0[j];

            tpx_yzz_xxx_0[j] = pa_y[j] * tpx_zz_xxx_0[j];

            tpy_yzz_xxx_0[j] = pa_y[j] * tpy_zz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxx_0[j];

            tpz_yzz_xxx_0[j] = pa_y[j] * tpz_zz_xxx_0[j];

            tpx_yzz_xxy_0[j] = pa_y[j] * tpx_zz_xxy_0[j] + 0.5 * fl1_fx * tpx_zz_xx_0[j];

            tpy_yzz_xxy_0[j] = pa_y[j] * tpy_zz_xxy_0[j] + 0.5 * fl1_fx * tpy_zz_xx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxy_0[j];

            tpz_yzz_xxy_0[j] = pa_y[j] * tpz_zz_xxy_0[j] + 0.5 * fl1_fx * tpz_zz_xx_0[j];

            tpx_yzz_xxz_0[j] = pa_y[j] * tpx_zz_xxz_0[j];

            tpy_yzz_xxz_0[j] = pa_y[j] * tpy_zz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zz_xxz_0[j];

            tpz_yzz_xxz_0[j] = pa_y[j] * tpz_zz_xxz_0[j];

            tpx_yzz_xyy_0[j] = pa_y[j] * tpx_zz_xyy_0[j] + fl1_fx * tpx_zz_xy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFF_250_300(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

        auto tpx_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 30);

        auto tpy_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 30);

        auto tpz_zz_xx_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 30);

        auto tpx_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 31);

        auto tpy_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 31);

        auto tpz_zz_xy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 31);

        auto tpx_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 32);

        auto tpy_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 32);

        auto tpz_zz_xz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 32);

        auto tpx_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 33);

        auto tpy_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 33);

        auto tpz_zz_yy_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 33);

        auto tpx_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 34);

        auto tpy_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 34);

        auto tpz_zz_yz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 34);

        auto tpx_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * idx + 35);

        auto tpy_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 36 * bdim + 36 * idx + 35);

        auto tpz_zz_zz_0 = primBuffer.data(pidx_p_2_2_m0 + 72 * bdim + 36 * idx + 35);

        auto ts_zz_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 50);

        auto ts_zz_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 51);

        auto ts_zz_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 52);

        auto ts_zz_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 53);

        auto ts_zz_xyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 54);

        auto ts_zz_xzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 55);

        auto ts_zz_yyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 56);

        auto ts_zz_yyz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 57);

        auto ts_zz_yzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 58);

        auto ts_zz_zzz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 59);

        // set up pointers to integrals

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

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yzz_xyz_0, tpx_yzz_xzz_0, tpx_yzz_yyy_0, \
                                     tpx_yzz_yyz_0, tpx_yzz_yzz_0, tpx_yzz_zzz_0, tpx_z_xxx_0, tpx_z_xxy_0, tpx_z_xxz_0, \
                                     tpx_z_xyy_0, tpx_z_xyz_0, tpx_z_xzz_0, tpx_z_yyy_0, tpx_z_yyz_0, tpx_z_yzz_0, \
                                     tpx_z_zzz_0, tpx_zz_xx_0, tpx_zz_xxx_0, tpx_zz_xxy_0, tpx_zz_xxz_0, tpx_zz_xy_0, \
                                     tpx_zz_xyy_0, tpx_zz_xyz_0, tpx_zz_xz_0, tpx_zz_xzz_0, tpx_zz_yy_0, tpx_zz_yyy_0, \
                                     tpx_zz_yyz_0, tpx_zz_yz_0, tpx_zz_yzz_0, tpx_zz_zz_0, tpx_zz_zzz_0, tpx_zzz_xxx_0, \
                                     tpx_zzz_xxy_0, tpx_zzz_xxz_0, tpx_zzz_xyy_0, tpx_zzz_xyz_0, tpx_zzz_xzz_0, \
                                     tpx_zzz_yyy_0, tpx_zzz_yyz_0, tpx_zzz_yzz_0, tpx_zzz_zzz_0, tpy_yzz_xyy_0, \
                                     tpy_yzz_xyz_0, tpy_yzz_xzz_0, tpy_yzz_yyy_0, tpy_yzz_yyz_0, tpy_yzz_yzz_0, \
                                     tpy_yzz_zzz_0, tpy_z_xxx_0, tpy_z_xxy_0, tpy_z_xxz_0, tpy_z_xyy_0, tpy_z_xyz_0, \
                                     tpy_z_xzz_0, tpy_z_yyy_0, tpy_z_yyz_0, tpy_z_yzz_0, tpy_z_zzz_0, tpy_zz_xx_0, \
                                     tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, tpy_zz_xy_0, tpy_zz_xyy_0, tpy_zz_xyz_0, \
                                     tpy_zz_xz_0, tpy_zz_xzz_0, tpy_zz_yy_0, tpy_zz_yyy_0, tpy_zz_yyz_0, tpy_zz_yz_0, \
                                     tpy_zz_yzz_0, tpy_zz_zz_0, tpy_zz_zzz_0, tpy_zzz_xxx_0, tpy_zzz_xxy_0, \
                                     tpy_zzz_xxz_0, tpy_zzz_xyy_0, tpy_zzz_xyz_0, tpy_zzz_xzz_0, tpy_zzz_yyy_0, \
                                     tpy_zzz_yyz_0, tpy_zzz_yzz_0, tpy_zzz_zzz_0, tpz_yzz_xyy_0, tpz_yzz_xyz_0, \
                                     tpz_yzz_xzz_0, tpz_yzz_yyy_0, tpz_yzz_yyz_0, tpz_yzz_yzz_0, tpz_yzz_zzz_0, \
                                     tpz_z_xxx_0, tpz_z_xxy_0, tpz_z_xxz_0, tpz_z_xyy_0, tpz_z_xyz_0, tpz_z_xzz_0, \
                                     tpz_z_yyy_0, tpz_z_yyz_0, tpz_z_yzz_0, tpz_z_zzz_0, tpz_zz_xx_0, tpz_zz_xxx_0, \
                                     tpz_zz_xxy_0, tpz_zz_xxz_0, tpz_zz_xy_0, tpz_zz_xyy_0, tpz_zz_xyz_0, tpz_zz_xz_0, \
                                     tpz_zz_xzz_0, tpz_zz_yy_0, tpz_zz_yyy_0, tpz_zz_yyz_0, tpz_zz_yz_0, tpz_zz_yzz_0, \
                                     tpz_zz_zz_0, tpz_zz_zzz_0, tpz_zzz_xxx_0, tpz_zzz_xxy_0, tpz_zzz_xxz_0, \
                                     tpz_zzz_xyy_0, tpz_zzz_xyz_0, tpz_zzz_xzz_0, tpz_zzz_yyy_0, tpz_zzz_yyz_0, \
                                     tpz_zzz_yzz_0, tpz_zzz_zzz_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0, ts_zz_xyy_0, \
                                     ts_zz_xyz_0, ts_zz_xzz_0, ts_zz_yyy_0, ts_zz_yyz_0, ts_zz_yzz_0, ts_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_yzz_xyy_0[j] = pa_y[j] * tpy_zz_xyy_0[j] + fl1_fx * tpy_zz_xy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyy_0[j];

            tpz_yzz_xyy_0[j] = pa_y[j] * tpz_zz_xyy_0[j] + fl1_fx * tpz_zz_xy_0[j];

            tpx_yzz_xyz_0[j] = pa_y[j] * tpx_zz_xyz_0[j] + 0.5 * fl1_fx * tpx_zz_xz_0[j];

            tpy_yzz_xyz_0[j] = pa_y[j] * tpy_zz_xyz_0[j] + 0.5 * fl1_fx * tpy_zz_xz_0[j] - fl1_fgb * fl1_fx * ts_zz_xyz_0[j];

            tpz_yzz_xyz_0[j] = pa_y[j] * tpz_zz_xyz_0[j] + 0.5 * fl1_fx * tpz_zz_xz_0[j];

            tpx_yzz_xzz_0[j] = pa_y[j] * tpx_zz_xzz_0[j];

            tpy_yzz_xzz_0[j] = pa_y[j] * tpy_zz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zz_xzz_0[j];

            tpz_yzz_xzz_0[j] = pa_y[j] * tpz_zz_xzz_0[j];

            tpx_yzz_yyy_0[j] = pa_y[j] * tpx_zz_yyy_0[j] + 1.5 * fl1_fx * tpx_zz_yy_0[j];

            tpy_yzz_yyy_0[j] = pa_y[j] * tpy_zz_yyy_0[j] + 1.5 * fl1_fx * tpy_zz_yy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyy_0[j];

            tpz_yzz_yyy_0[j] = pa_y[j] * tpz_zz_yyy_0[j] + 1.5 * fl1_fx * tpz_zz_yy_0[j];

            tpx_yzz_yyz_0[j] = pa_y[j] * tpx_zz_yyz_0[j] + fl1_fx * tpx_zz_yz_0[j];

            tpy_yzz_yyz_0[j] = pa_y[j] * tpy_zz_yyz_0[j] + fl1_fx * tpy_zz_yz_0[j] - fl1_fgb * fl1_fx * ts_zz_yyz_0[j];

            tpz_yzz_yyz_0[j] = pa_y[j] * tpz_zz_yyz_0[j] + fl1_fx * tpz_zz_yz_0[j];

            tpx_yzz_yzz_0[j] = pa_y[j] * tpx_zz_yzz_0[j] + 0.5 * fl1_fx * tpx_zz_zz_0[j];

            tpy_yzz_yzz_0[j] = pa_y[j] * tpy_zz_yzz_0[j] + 0.5 * fl1_fx * tpy_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_zz_yzz_0[j];

            tpz_yzz_yzz_0[j] = pa_y[j] * tpz_zz_yzz_0[j] + 0.5 * fl1_fx * tpz_zz_zz_0[j];

            tpx_yzz_zzz_0[j] = pa_y[j] * tpx_zz_zzz_0[j];

            tpy_yzz_zzz_0[j] = pa_y[j] * tpy_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zz_zzz_0[j];

            tpz_yzz_zzz_0[j] = pa_y[j] * tpz_zz_zzz_0[j];

            tpx_zzz_xxx_0[j] = pa_z[j] * tpx_zz_xxx_0[j] + fl1_fx * tpx_z_xxx_0[j];

            tpy_zzz_xxx_0[j] = pa_z[j] * tpy_zz_xxx_0[j] + fl1_fx * tpy_z_xxx_0[j];

            tpz_zzz_xxx_0[j] = pa_z[j] * tpz_zz_xxx_0[j] + fl1_fx * tpz_z_xxx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxx_0[j];

            tpx_zzz_xxy_0[j] = pa_z[j] * tpx_zz_xxy_0[j] + fl1_fx * tpx_z_xxy_0[j];

            tpy_zzz_xxy_0[j] = pa_z[j] * tpy_zz_xxy_0[j] + fl1_fx * tpy_z_xxy_0[j];

            tpz_zzz_xxy_0[j] = pa_z[j] * tpz_zz_xxy_0[j] + fl1_fx * tpz_z_xxy_0[j] - fl1_fgb * fl1_fx * ts_zz_xxy_0[j];

            tpx_zzz_xxz_0[j] = pa_z[j] * tpx_zz_xxz_0[j] + fl1_fx * tpx_z_xxz_0[j] + 0.5 * fl1_fx * tpx_zz_xx_0[j];

            tpy_zzz_xxz_0[j] = pa_z[j] * tpy_zz_xxz_0[j] + fl1_fx * tpy_z_xxz_0[j] + 0.5 * fl1_fx * tpy_zz_xx_0[j];

            tpz_zzz_xxz_0[j] =
                pa_z[j] * tpz_zz_xxz_0[j] + fl1_fx * tpz_z_xxz_0[j] + 0.5 * fl1_fx * tpz_zz_xx_0[j] - fl1_fgb * fl1_fx * ts_zz_xxz_0[j];

            tpx_zzz_xyy_0[j] = pa_z[j] * tpx_zz_xyy_0[j] + fl1_fx * tpx_z_xyy_0[j];

            tpy_zzz_xyy_0[j] = pa_z[j] * tpy_zz_xyy_0[j] + fl1_fx * tpy_z_xyy_0[j];

            tpz_zzz_xyy_0[j] = pa_z[j] * tpz_zz_xyy_0[j] + fl1_fx * tpz_z_xyy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyy_0[j];

            tpx_zzz_xyz_0[j] = pa_z[j] * tpx_zz_xyz_0[j] + fl1_fx * tpx_z_xyz_0[j] + 0.5 * fl1_fx * tpx_zz_xy_0[j];

            tpy_zzz_xyz_0[j] = pa_z[j] * tpy_zz_xyz_0[j] + fl1_fx * tpy_z_xyz_0[j] + 0.5 * fl1_fx * tpy_zz_xy_0[j];

            tpz_zzz_xyz_0[j] =
                pa_z[j] * tpz_zz_xyz_0[j] + fl1_fx * tpz_z_xyz_0[j] + 0.5 * fl1_fx * tpz_zz_xy_0[j] - fl1_fgb * fl1_fx * ts_zz_xyz_0[j];

            tpx_zzz_xzz_0[j] = pa_z[j] * tpx_zz_xzz_0[j] + fl1_fx * tpx_z_xzz_0[j] + fl1_fx * tpx_zz_xz_0[j];

            tpy_zzz_xzz_0[j] = pa_z[j] * tpy_zz_xzz_0[j] + fl1_fx * tpy_z_xzz_0[j] + fl1_fx * tpy_zz_xz_0[j];

            tpz_zzz_xzz_0[j] = pa_z[j] * tpz_zz_xzz_0[j] + fl1_fx * tpz_z_xzz_0[j] + fl1_fx * tpz_zz_xz_0[j] - fl1_fgb * fl1_fx * ts_zz_xzz_0[j];

            tpx_zzz_yyy_0[j] = pa_z[j] * tpx_zz_yyy_0[j] + fl1_fx * tpx_z_yyy_0[j];

            tpy_zzz_yyy_0[j] = pa_z[j] * tpy_zz_yyy_0[j] + fl1_fx * tpy_z_yyy_0[j];

            tpz_zzz_yyy_0[j] = pa_z[j] * tpz_zz_yyy_0[j] + fl1_fx * tpz_z_yyy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyy_0[j];

            tpx_zzz_yyz_0[j] = pa_z[j] * tpx_zz_yyz_0[j] + fl1_fx * tpx_z_yyz_0[j] + 0.5 * fl1_fx * tpx_zz_yy_0[j];

            tpy_zzz_yyz_0[j] = pa_z[j] * tpy_zz_yyz_0[j] + fl1_fx * tpy_z_yyz_0[j] + 0.5 * fl1_fx * tpy_zz_yy_0[j];

            tpz_zzz_yyz_0[j] =
                pa_z[j] * tpz_zz_yyz_0[j] + fl1_fx * tpz_z_yyz_0[j] + 0.5 * fl1_fx * tpz_zz_yy_0[j] - fl1_fgb * fl1_fx * ts_zz_yyz_0[j];

            tpx_zzz_yzz_0[j] = pa_z[j] * tpx_zz_yzz_0[j] + fl1_fx * tpx_z_yzz_0[j] + fl1_fx * tpx_zz_yz_0[j];

            tpy_zzz_yzz_0[j] = pa_z[j] * tpy_zz_yzz_0[j] + fl1_fx * tpy_z_yzz_0[j] + fl1_fx * tpy_zz_yz_0[j];

            tpz_zzz_yzz_0[j] = pa_z[j] * tpz_zz_yzz_0[j] + fl1_fx * tpz_z_yzz_0[j] + fl1_fx * tpz_zz_yz_0[j] - fl1_fgb * fl1_fx * ts_zz_yzz_0[j];

            tpx_zzz_zzz_0[j] = pa_z[j] * tpx_zz_zzz_0[j] + fl1_fx * tpx_z_zzz_0[j] + 1.5 * fl1_fx * tpx_zz_zz_0[j];

            tpy_zzz_zzz_0[j] = pa_z[j] * tpy_zz_zzz_0[j] + fl1_fx * tpy_z_zzz_0[j] + 1.5 * fl1_fx * tpy_zz_zz_0[j];

            tpz_zzz_zzz_0[j] =
                pa_z[j] * tpz_zz_zzz_0[j] + fl1_fx * tpz_z_zzz_0[j] + 1.5 * fl1_fx * tpz_zz_zz_0[j] - fl1_fgb * fl1_fx * ts_zz_zzz_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
