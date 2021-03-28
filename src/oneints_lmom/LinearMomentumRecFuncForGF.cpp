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

#include "LinearMomentumRecFuncForGF.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForGF(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForGF_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_150_200(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_200_250(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_250_300(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_300_350(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_350_400(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGF_400_450(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForGF_0_50(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx);

        auto tpy_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx);

        auto tpz_xxx_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx);

        auto tpx_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 1);

        auto tpy_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 1);

        auto tpz_xxx_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 1);

        auto tpx_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 2);

        auto tpy_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 2);

        auto tpz_xxx_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 2);

        auto tpx_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 3);

        auto tpy_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 3);

        auto tpz_xxx_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 3);

        auto tpx_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 4);

        auto tpy_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 4);

        auto tpz_xxx_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 4);

        auto tpx_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 5);

        auto tpy_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 5);

        auto tpz_xxx_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 5);

        auto tpx_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 6);

        auto tpy_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 6);

        auto tpz_xxy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 6);

        auto tpx_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 7);

        auto tpy_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 7);

        auto tpz_xxy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 7);

        auto tpx_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 8);

        auto tpy_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 8);

        auto tpz_xxy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 8);

        auto tpx_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 9);

        auto tpy_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 9);

        auto tpz_xxy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 9);

        auto tpx_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 10);

        auto tpy_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 10);

        auto tpz_xxy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 10);

        auto tpx_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 11);

        auto tpy_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 11);

        auto tpz_xxy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 11);

        auto ts_xxx_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx);

        auto ts_xxx_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 1);

        auto ts_xxx_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 2);

        auto ts_xxx_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 3);

        auto ts_xxx_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 4);

        auto ts_xxx_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 5);

        auto ts_xxx_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 6);

        auto ts_xxx_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 7);

        auto ts_xxx_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 8);

        auto ts_xxx_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 9);

        auto ts_xxy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 10);

        auto ts_xxy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 11);

        auto ts_xxy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 12);

        auto ts_xxy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 13);

        auto ts_xxy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 14);

        auto ts_xxy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 15);

        auto ts_xxy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 16);

        // set up pointers to integrals

        auto tpx_xxxx_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx);

        auto tpy_xxxx_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx);

        auto tpz_xxxx_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx);

        auto tpx_xxxx_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 1);

        auto tpy_xxxx_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 1);

        auto tpz_xxxx_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 1);

        auto tpx_xxxx_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 2);

        auto tpy_xxxx_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 2);

        auto tpz_xxxx_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 2);

        auto tpx_xxxx_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 3);

        auto tpy_xxxx_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 3);

        auto tpz_xxxx_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 3);

        auto tpx_xxxx_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 4);

        auto tpy_xxxx_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 4);

        auto tpz_xxxx_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 4);

        auto tpx_xxxx_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 5);

        auto tpy_xxxx_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 5);

        auto tpz_xxxx_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 5);

        auto tpx_xxxx_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 6);

        auto tpy_xxxx_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 6);

        auto tpz_xxxx_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 6);

        auto tpx_xxxx_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 7);

        auto tpy_xxxx_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 7);

        auto tpz_xxxx_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 7);

        auto tpx_xxxx_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 8);

        auto tpy_xxxx_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 8);

        auto tpz_xxxx_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 8);

        auto tpx_xxxx_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 9);

        auto tpy_xxxx_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 9);

        auto tpz_xxxx_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 9);

        auto tpx_xxxy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 10);

        auto tpy_xxxy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 10);

        auto tpz_xxxy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 10);

        auto tpx_xxxy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 11);

        auto tpy_xxxy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 11);

        auto tpz_xxxy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 11);

        auto tpx_xxxy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 12);

        auto tpy_xxxy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 12);

        auto tpz_xxxy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 12);

        auto tpx_xxxy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 13);

        auto tpy_xxxy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 13);

        auto tpz_xxxy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 13);

        auto tpx_xxxy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 14);

        auto tpy_xxxy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 14);

        auto tpz_xxxy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 14);

        auto tpx_xxxy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 15);

        auto tpy_xxxy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 15);

        auto tpz_xxxy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 15);

        auto tpx_xxxy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 16);

        auto tpy_xxxy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 16);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xx_xxx_0, tpx_xx_xxy_0, tpx_xx_xxz_0, tpx_xx_xyy_0, \
                                     tpx_xx_xyz_0, tpx_xx_xzz_0, tpx_xx_yyy_0, tpx_xx_yyz_0, tpx_xx_yzz_0, tpx_xx_zzz_0, \
                                     tpx_xxx_xx_0, tpx_xxx_xxx_0, tpx_xxx_xxy_0, tpx_xxx_xxz_0, tpx_xxx_xy_0, \
                                     tpx_xxx_xyy_0, tpx_xxx_xyz_0, tpx_xxx_xz_0, tpx_xxx_xzz_0, tpx_xxx_yy_0, \
                                     tpx_xxx_yyy_0, tpx_xxx_yyz_0, tpx_xxx_yz_0, tpx_xxx_yzz_0, tpx_xxx_zz_0, \
                                     tpx_xxx_zzz_0, tpx_xxxx_xxx_0, tpx_xxxx_xxy_0, tpx_xxxx_xxz_0, tpx_xxxx_xyy_0, \
                                     tpx_xxxx_xyz_0, tpx_xxxx_xzz_0, tpx_xxxx_yyy_0, tpx_xxxx_yyz_0, tpx_xxxx_yzz_0, \
                                     tpx_xxxx_zzz_0, tpx_xxxy_xxx_0, tpx_xxxy_xxy_0, tpx_xxxy_xxz_0, tpx_xxxy_xyy_0, \
                                     tpx_xxxy_xyz_0, tpx_xxxy_xzz_0, tpx_xxxy_yyy_0, tpx_xxy_xx_0, tpx_xxy_xxx_0, \
                                     tpx_xxy_xxy_0, tpx_xxy_xxz_0, tpx_xxy_xy_0, tpx_xxy_xyy_0, tpx_xxy_xyz_0, \
                                     tpx_xxy_xz_0, tpx_xxy_xzz_0, tpx_xxy_yy_0, tpx_xxy_yyy_0, tpx_xxy_yz_0, \
                                     tpx_xxy_zz_0, tpx_xy_xxx_0, tpx_xy_xxy_0, tpx_xy_xxz_0, tpx_xy_xyy_0, tpx_xy_xyz_0, \
                                     tpx_xy_xzz_0, tpx_xy_yyy_0, tpy_xx_xxx_0, tpy_xx_xxy_0, tpy_xx_xxz_0, tpy_xx_xyy_0, \
                                     tpy_xx_xyz_0, tpy_xx_xzz_0, tpy_xx_yyy_0, tpy_xx_yyz_0, tpy_xx_yzz_0, tpy_xx_zzz_0, \
                                     tpy_xxx_xx_0, tpy_xxx_xxx_0, tpy_xxx_xxy_0, tpy_xxx_xxz_0, tpy_xxx_xy_0, \
                                     tpy_xxx_xyy_0, tpy_xxx_xyz_0, tpy_xxx_xz_0, tpy_xxx_xzz_0, tpy_xxx_yy_0, \
                                     tpy_xxx_yyy_0, tpy_xxx_yyz_0, tpy_xxx_yz_0, tpy_xxx_yzz_0, tpy_xxx_zz_0, \
                                     tpy_xxx_zzz_0, tpy_xxxx_xxx_0, tpy_xxxx_xxy_0, tpy_xxxx_xxz_0, tpy_xxxx_xyy_0, \
                                     tpy_xxxx_xyz_0, tpy_xxxx_xzz_0, tpy_xxxx_yyy_0, tpy_xxxx_yyz_0, tpy_xxxx_yzz_0, \
                                     tpy_xxxx_zzz_0, tpy_xxxy_xxx_0, tpy_xxxy_xxy_0, tpy_xxxy_xxz_0, tpy_xxxy_xyy_0, \
                                     tpy_xxxy_xyz_0, tpy_xxxy_xzz_0, tpy_xxxy_yyy_0, tpy_xxy_xx_0, tpy_xxy_xxx_0, \
                                     tpy_xxy_xxy_0, tpy_xxy_xxz_0, tpy_xxy_xy_0, tpy_xxy_xyy_0, tpy_xxy_xyz_0, \
                                     tpy_xxy_xz_0, tpy_xxy_xzz_0, tpy_xxy_yy_0, tpy_xxy_yyy_0, tpy_xxy_yz_0, \
                                     tpy_xxy_zz_0, tpy_xy_xxx_0, tpy_xy_xxy_0, tpy_xy_xxz_0, tpy_xy_xyy_0, tpy_xy_xyz_0, \
                                     tpy_xy_xzz_0, tpy_xy_yyy_0, tpz_xx_xxx_0, tpz_xx_xxy_0, tpz_xx_xxz_0, tpz_xx_xyy_0, \
                                     tpz_xx_xyz_0, tpz_xx_xzz_0, tpz_xx_yyy_0, tpz_xx_yyz_0, tpz_xx_yzz_0, tpz_xx_zzz_0, \
                                     tpz_xxx_xx_0, tpz_xxx_xxx_0, tpz_xxx_xxy_0, tpz_xxx_xxz_0, tpz_xxx_xy_0, \
                                     tpz_xxx_xyy_0, tpz_xxx_xyz_0, tpz_xxx_xz_0, tpz_xxx_xzz_0, tpz_xxx_yy_0, \
                                     tpz_xxx_yyy_0, tpz_xxx_yyz_0, tpz_xxx_yz_0, tpz_xxx_yzz_0, tpz_xxx_zz_0, \
                                     tpz_xxx_zzz_0, tpz_xxxx_xxx_0, tpz_xxxx_xxy_0, tpz_xxxx_xxz_0, tpz_xxxx_xyy_0, \
                                     tpz_xxxx_xyz_0, tpz_xxxx_xzz_0, tpz_xxxx_yyy_0, tpz_xxxx_yyz_0, tpz_xxxx_yzz_0, \
                                     tpz_xxxx_zzz_0, tpz_xxxy_xxx_0, tpz_xxxy_xxy_0, tpz_xxxy_xxz_0, tpz_xxxy_xyy_0, \
                                     tpz_xxxy_xyz_0, tpz_xxxy_xzz_0, tpz_xxy_xx_0, tpz_xxy_xxx_0, tpz_xxy_xxy_0, \
                                     tpz_xxy_xxz_0, tpz_xxy_xy_0, tpz_xxy_xyy_0, tpz_xxy_xyz_0, tpz_xxy_xz_0, \
                                     tpz_xxy_xzz_0, tpz_xxy_yy_0, tpz_xxy_yz_0, tpz_xxy_zz_0, tpz_xy_xxx_0, tpz_xy_xxy_0, \
                                     tpz_xy_xxz_0, tpz_xy_xyy_0, tpz_xy_xyz_0, tpz_xy_xzz_0, ts_xxx_xxx_0, ts_xxx_xxy_0, \
                                     ts_xxx_xxz_0, ts_xxx_xyy_0, ts_xxx_xyz_0, ts_xxx_xzz_0, ts_xxx_yyy_0, ts_xxx_yyz_0, \
                                     ts_xxx_yzz_0, ts_xxx_zzz_0, ts_xxy_xxx_0, ts_xxy_xxy_0, ts_xxy_xxz_0, ts_xxy_xyy_0, \
                                     ts_xxy_xyz_0, ts_xxy_xzz_0, ts_xxy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxxx_xxx_0[j] =
                pa_x[j] * tpx_xxx_xxx_0[j] + 1.5 * fl1_fx * tpx_xx_xxx_0[j] + 1.5 * fl1_fx * tpx_xxx_xx_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxx_0[j];

            tpy_xxxx_xxx_0[j] = pa_x[j] * tpy_xxx_xxx_0[j] + 1.5 * fl1_fx * tpy_xx_xxx_0[j] + 1.5 * fl1_fx * tpy_xxx_xx_0[j];

            tpz_xxxx_xxx_0[j] = pa_x[j] * tpz_xxx_xxx_0[j] + 1.5 * fl1_fx * tpz_xx_xxx_0[j] + 1.5 * fl1_fx * tpz_xxx_xx_0[j];

            tpx_xxxx_xxy_0[j] =
                pa_x[j] * tpx_xxx_xxy_0[j] + 1.5 * fl1_fx * tpx_xx_xxy_0[j] + fl1_fx * tpx_xxx_xy_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxy_0[j];

            tpy_xxxx_xxy_0[j] = pa_x[j] * tpy_xxx_xxy_0[j] + 1.5 * fl1_fx * tpy_xx_xxy_0[j] + fl1_fx * tpy_xxx_xy_0[j];

            tpz_xxxx_xxy_0[j] = pa_x[j] * tpz_xxx_xxy_0[j] + 1.5 * fl1_fx * tpz_xx_xxy_0[j] + fl1_fx * tpz_xxx_xy_0[j];

            tpx_xxxx_xxz_0[j] =
                pa_x[j] * tpx_xxx_xxz_0[j] + 1.5 * fl1_fx * tpx_xx_xxz_0[j] + fl1_fx * tpx_xxx_xz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xxz_0[j];

            tpy_xxxx_xxz_0[j] = pa_x[j] * tpy_xxx_xxz_0[j] + 1.5 * fl1_fx * tpy_xx_xxz_0[j] + fl1_fx * tpy_xxx_xz_0[j];

            tpz_xxxx_xxz_0[j] = pa_x[j] * tpz_xxx_xxz_0[j] + 1.5 * fl1_fx * tpz_xx_xxz_0[j] + fl1_fx * tpz_xxx_xz_0[j];

            tpx_xxxx_xyy_0[j] =
                pa_x[j] * tpx_xxx_xyy_0[j] + 1.5 * fl1_fx * tpx_xx_xyy_0[j] + 0.5 * fl1_fx * tpx_xxx_yy_0[j] - fl1_fgb * fl1_fx * ts_xxx_xyy_0[j];

            tpy_xxxx_xyy_0[j] = pa_x[j] * tpy_xxx_xyy_0[j] + 1.5 * fl1_fx * tpy_xx_xyy_0[j] + 0.5 * fl1_fx * tpy_xxx_yy_0[j];

            tpz_xxxx_xyy_0[j] = pa_x[j] * tpz_xxx_xyy_0[j] + 1.5 * fl1_fx * tpz_xx_xyy_0[j] + 0.5 * fl1_fx * tpz_xxx_yy_0[j];

            tpx_xxxx_xyz_0[j] =
                pa_x[j] * tpx_xxx_xyz_0[j] + 1.5 * fl1_fx * tpx_xx_xyz_0[j] + 0.5 * fl1_fx * tpx_xxx_yz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xyz_0[j];

            tpy_xxxx_xyz_0[j] = pa_x[j] * tpy_xxx_xyz_0[j] + 1.5 * fl1_fx * tpy_xx_xyz_0[j] + 0.5 * fl1_fx * tpy_xxx_yz_0[j];

            tpz_xxxx_xyz_0[j] = pa_x[j] * tpz_xxx_xyz_0[j] + 1.5 * fl1_fx * tpz_xx_xyz_0[j] + 0.5 * fl1_fx * tpz_xxx_yz_0[j];

            tpx_xxxx_xzz_0[j] =
                pa_x[j] * tpx_xxx_xzz_0[j] + 1.5 * fl1_fx * tpx_xx_xzz_0[j] + 0.5 * fl1_fx * tpx_xxx_zz_0[j] - fl1_fgb * fl1_fx * ts_xxx_xzz_0[j];

            tpy_xxxx_xzz_0[j] = pa_x[j] * tpy_xxx_xzz_0[j] + 1.5 * fl1_fx * tpy_xx_xzz_0[j] + 0.5 * fl1_fx * tpy_xxx_zz_0[j];

            tpz_xxxx_xzz_0[j] = pa_x[j] * tpz_xxx_xzz_0[j] + 1.5 * fl1_fx * tpz_xx_xzz_0[j] + 0.5 * fl1_fx * tpz_xxx_zz_0[j];

            tpx_xxxx_yyy_0[j] = pa_x[j] * tpx_xxx_yyy_0[j] + 1.5 * fl1_fx * tpx_xx_yyy_0[j] - fl1_fgb * fl1_fx * ts_xxx_yyy_0[j];

            tpy_xxxx_yyy_0[j] = pa_x[j] * tpy_xxx_yyy_0[j] + 1.5 * fl1_fx * tpy_xx_yyy_0[j];

            tpz_xxxx_yyy_0[j] = pa_x[j] * tpz_xxx_yyy_0[j] + 1.5 * fl1_fx * tpz_xx_yyy_0[j];

            tpx_xxxx_yyz_0[j] = pa_x[j] * tpx_xxx_yyz_0[j] + 1.5 * fl1_fx * tpx_xx_yyz_0[j] - fl1_fgb * fl1_fx * ts_xxx_yyz_0[j];

            tpy_xxxx_yyz_0[j] = pa_x[j] * tpy_xxx_yyz_0[j] + 1.5 * fl1_fx * tpy_xx_yyz_0[j];

            tpz_xxxx_yyz_0[j] = pa_x[j] * tpz_xxx_yyz_0[j] + 1.5 * fl1_fx * tpz_xx_yyz_0[j];

            tpx_xxxx_yzz_0[j] = pa_x[j] * tpx_xxx_yzz_0[j] + 1.5 * fl1_fx * tpx_xx_yzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_yzz_0[j];

            tpy_xxxx_yzz_0[j] = pa_x[j] * tpy_xxx_yzz_0[j] + 1.5 * fl1_fx * tpy_xx_yzz_0[j];

            tpz_xxxx_yzz_0[j] = pa_x[j] * tpz_xxx_yzz_0[j] + 1.5 * fl1_fx * tpz_xx_yzz_0[j];

            tpx_xxxx_zzz_0[j] = pa_x[j] * tpx_xxx_zzz_0[j] + 1.5 * fl1_fx * tpx_xx_zzz_0[j] - fl1_fgb * fl1_fx * ts_xxx_zzz_0[j];

            tpy_xxxx_zzz_0[j] = pa_x[j] * tpy_xxx_zzz_0[j] + 1.5 * fl1_fx * tpy_xx_zzz_0[j];

            tpz_xxxx_zzz_0[j] = pa_x[j] * tpz_xxx_zzz_0[j] + 1.5 * fl1_fx * tpz_xx_zzz_0[j];

            tpx_xxxy_xxx_0[j] =
                pa_x[j] * tpx_xxy_xxx_0[j] + fl1_fx * tpx_xy_xxx_0[j] + 1.5 * fl1_fx * tpx_xxy_xx_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxx_0[j];

            tpy_xxxy_xxx_0[j] = pa_x[j] * tpy_xxy_xxx_0[j] + fl1_fx * tpy_xy_xxx_0[j] + 1.5 * fl1_fx * tpy_xxy_xx_0[j];

            tpz_xxxy_xxx_0[j] = pa_x[j] * tpz_xxy_xxx_0[j] + fl1_fx * tpz_xy_xxx_0[j] + 1.5 * fl1_fx * tpz_xxy_xx_0[j];

            tpx_xxxy_xxy_0[j] = pa_x[j] * tpx_xxy_xxy_0[j] + fl1_fx * tpx_xy_xxy_0[j] + fl1_fx * tpx_xxy_xy_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxy_0[j];

            tpy_xxxy_xxy_0[j] = pa_x[j] * tpy_xxy_xxy_0[j] + fl1_fx * tpy_xy_xxy_0[j] + fl1_fx * tpy_xxy_xy_0[j];

            tpz_xxxy_xxy_0[j] = pa_x[j] * tpz_xxy_xxy_0[j] + fl1_fx * tpz_xy_xxy_0[j] + fl1_fx * tpz_xxy_xy_0[j];

            tpx_xxxy_xxz_0[j] = pa_x[j] * tpx_xxy_xxz_0[j] + fl1_fx * tpx_xy_xxz_0[j] + fl1_fx * tpx_xxy_xz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xxz_0[j];

            tpy_xxxy_xxz_0[j] = pa_x[j] * tpy_xxy_xxz_0[j] + fl1_fx * tpy_xy_xxz_0[j] + fl1_fx * tpy_xxy_xz_0[j];

            tpz_xxxy_xxz_0[j] = pa_x[j] * tpz_xxy_xxz_0[j] + fl1_fx * tpz_xy_xxz_0[j] + fl1_fx * tpz_xxy_xz_0[j];

            tpx_xxxy_xyy_0[j] =
                pa_x[j] * tpx_xxy_xyy_0[j] + fl1_fx * tpx_xy_xyy_0[j] + 0.5 * fl1_fx * tpx_xxy_yy_0[j] - fl1_fgb * fl1_fx * ts_xxy_xyy_0[j];

            tpy_xxxy_xyy_0[j] = pa_x[j] * tpy_xxy_xyy_0[j] + fl1_fx * tpy_xy_xyy_0[j] + 0.5 * fl1_fx * tpy_xxy_yy_0[j];

            tpz_xxxy_xyy_0[j] = pa_x[j] * tpz_xxy_xyy_0[j] + fl1_fx * tpz_xy_xyy_0[j] + 0.5 * fl1_fx * tpz_xxy_yy_0[j];

            tpx_xxxy_xyz_0[j] =
                pa_x[j] * tpx_xxy_xyz_0[j] + fl1_fx * tpx_xy_xyz_0[j] + 0.5 * fl1_fx * tpx_xxy_yz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xyz_0[j];

            tpy_xxxy_xyz_0[j] = pa_x[j] * tpy_xxy_xyz_0[j] + fl1_fx * tpy_xy_xyz_0[j] + 0.5 * fl1_fx * tpy_xxy_yz_0[j];

            tpz_xxxy_xyz_0[j] = pa_x[j] * tpz_xxy_xyz_0[j] + fl1_fx * tpz_xy_xyz_0[j] + 0.5 * fl1_fx * tpz_xxy_yz_0[j];

            tpx_xxxy_xzz_0[j] =
                pa_x[j] * tpx_xxy_xzz_0[j] + fl1_fx * tpx_xy_xzz_0[j] + 0.5 * fl1_fx * tpx_xxy_zz_0[j] - fl1_fgb * fl1_fx * ts_xxy_xzz_0[j];

            tpy_xxxy_xzz_0[j] = pa_x[j] * tpy_xxy_xzz_0[j] + fl1_fx * tpy_xy_xzz_0[j] + 0.5 * fl1_fx * tpy_xxy_zz_0[j];

            tpz_xxxy_xzz_0[j] = pa_x[j] * tpz_xxy_xzz_0[j] + fl1_fx * tpz_xy_xzz_0[j] + 0.5 * fl1_fx * tpz_xxy_zz_0[j];

            tpx_xxxy_yyy_0[j] = pa_x[j] * tpx_xxy_yyy_0[j] + fl1_fx * tpx_xy_yyy_0[j] - fl1_fgb * fl1_fx * ts_xxy_yyy_0[j];

            tpy_xxxy_yyy_0[j] = pa_x[j] * tpy_xxy_yyy_0[j] + fl1_fx * tpy_xy_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_50_100(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 12);

        auto tpy_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 12);

        auto tpz_xxz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 12);

        auto tpx_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 13);

        auto tpy_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 13);

        auto tpz_xxz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 13);

        auto tpx_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 14);

        auto tpy_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 14);

        auto tpz_xxz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 14);

        auto tpx_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 15);

        auto tpy_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 15);

        auto tpz_xxz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 15);

        auto tpx_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 16);

        auto tpy_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 16);

        auto tpz_xxz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 16);

        auto tpx_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 17);

        auto tpy_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 17);

        auto tpz_xxz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 17);

        auto tpx_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 18);

        auto tpy_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 18);

        auto tpz_xyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 18);

        auto tpx_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 19);

        auto tpy_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 19);

        auto tpz_xyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 19);

        auto tpx_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 20);

        auto tpy_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 20);

        auto tpz_xyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 20);

        auto tpx_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 21);

        auto ts_xxy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 17);

        auto ts_xxy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 18);

        auto ts_xxy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 19);

        auto ts_xxz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 20);

        auto ts_xxz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 21);

        auto ts_xxz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 22);

        auto ts_xxz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 23);

        auto ts_xxz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 24);

        auto ts_xxz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 25);

        auto ts_xxz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 26);

        auto ts_xxz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 27);

        auto ts_xxz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 28);

        auto ts_xxz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 29);

        auto ts_xyy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 30);

        auto ts_xyy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 31);

        auto ts_xyy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 32);

        auto ts_xyy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 33);

        // set up pointers to integrals

        auto tpz_xxxy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 16);

        auto tpx_xxxy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 17);

        auto tpy_xxxy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 17);

        auto tpz_xxxy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 17);

        auto tpx_xxxy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 18);

        auto tpy_xxxy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 18);

        auto tpz_xxxy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 18);

        auto tpx_xxxy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 19);

        auto tpy_xxxy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 19);

        auto tpz_xxxy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 19);

        auto tpx_xxxz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 20);

        auto tpy_xxxz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 20);

        auto tpz_xxxz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 20);

        auto tpx_xxxz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 21);

        auto tpy_xxxz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 21);

        auto tpz_xxxz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 21);

        auto tpx_xxxz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 22);

        auto tpy_xxxz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 22);

        auto tpz_xxxz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 22);

        auto tpx_xxxz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 23);

        auto tpy_xxxz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 23);

        auto tpz_xxxz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 23);

        auto tpx_xxxz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 24);

        auto tpy_xxxz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 24);

        auto tpz_xxxz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 24);

        auto tpx_xxxz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 25);

        auto tpy_xxxz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 25);

        auto tpz_xxxz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 25);

        auto tpx_xxxz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 26);

        auto tpy_xxxz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 26);

        auto tpz_xxxz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 26);

        auto tpx_xxxz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 27);

        auto tpy_xxxz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 27);

        auto tpz_xxxz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 27);

        auto tpx_xxxz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 28);

        auto tpy_xxxz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 28);

        auto tpz_xxxz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 28);

        auto tpx_xxxz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 29);

        auto tpy_xxxz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 29);

        auto tpz_xxxz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 29);

        auto tpx_xxyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 30);

        auto tpy_xxyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 30);

        auto tpz_xxyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 30);

        auto tpx_xxyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 31);

        auto tpy_xxyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 31);

        auto tpz_xxyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 31);

        auto tpx_xxyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 32);

        auto tpy_xxyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 32);

        auto tpz_xxyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 32);

        auto tpx_xxyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 33);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxxy_yyz_0, tpx_xxxy_yzz_0, tpx_xxxy_zzz_0, \
                                     tpx_xxxz_xxx_0, tpx_xxxz_xxy_0, tpx_xxxz_xxz_0, tpx_xxxz_xyy_0, tpx_xxxz_xyz_0, \
                                     tpx_xxxz_xzz_0, tpx_xxxz_yyy_0, tpx_xxxz_yyz_0, tpx_xxxz_yzz_0, tpx_xxxz_zzz_0, \
                                     tpx_xxy_yyz_0, tpx_xxy_yzz_0, tpx_xxy_zzz_0, tpx_xxyy_xxx_0, tpx_xxyy_xxy_0, \
                                     tpx_xxyy_xxz_0, tpx_xxyy_xyy_0, tpx_xxz_xx_0, tpx_xxz_xxx_0, tpx_xxz_xxy_0, \
                                     tpx_xxz_xxz_0, tpx_xxz_xy_0, tpx_xxz_xyy_0, tpx_xxz_xyz_0, tpx_xxz_xz_0, \
                                     tpx_xxz_xzz_0, tpx_xxz_yy_0, tpx_xxz_yyy_0, tpx_xxz_yyz_0, tpx_xxz_yz_0, \
                                     tpx_xxz_yzz_0, tpx_xxz_zz_0, tpx_xxz_zzz_0, tpx_xy_yyz_0, tpx_xy_yzz_0, \
                                     tpx_xy_zzz_0, tpx_xyy_xx_0, tpx_xyy_xxx_0, tpx_xyy_xxy_0, tpx_xyy_xxz_0, \
                                     tpx_xyy_xy_0, tpx_xyy_xyy_0, tpx_xyy_xz_0, tpx_xyy_yy_0, tpx_xz_xxx_0, \
                                     tpx_xz_xxy_0, tpx_xz_xxz_0, tpx_xz_xyy_0, tpx_xz_xyz_0, tpx_xz_xzz_0, tpx_xz_yyy_0, \
                                     tpx_xz_yyz_0, tpx_xz_yzz_0, tpx_xz_zzz_0, tpx_yy_xxx_0, tpx_yy_xxy_0, tpx_yy_xxz_0, \
                                     tpx_yy_xyy_0, tpy_xxxy_yyz_0, tpy_xxxy_yzz_0, tpy_xxxy_zzz_0, tpy_xxxz_xxx_0, \
                                     tpy_xxxz_xxy_0, tpy_xxxz_xxz_0, tpy_xxxz_xyy_0, tpy_xxxz_xyz_0, tpy_xxxz_xzz_0, \
                                     tpy_xxxz_yyy_0, tpy_xxxz_yyz_0, tpy_xxxz_yzz_0, tpy_xxxz_zzz_0, tpy_xxy_yyz_0, \
                                     tpy_xxy_yzz_0, tpy_xxy_zzz_0, tpy_xxyy_xxx_0, tpy_xxyy_xxy_0, tpy_xxyy_xxz_0, \
                                     tpy_xxz_xx_0, tpy_xxz_xxx_0, tpy_xxz_xxy_0, tpy_xxz_xxz_0, tpy_xxz_xy_0, \
                                     tpy_xxz_xyy_0, tpy_xxz_xyz_0, tpy_xxz_xz_0, tpy_xxz_xzz_0, tpy_xxz_yy_0, \
                                     tpy_xxz_yyy_0, tpy_xxz_yyz_0, tpy_xxz_yz_0, tpy_xxz_yzz_0, tpy_xxz_zz_0, \
                                     tpy_xxz_zzz_0, tpy_xy_yyz_0, tpy_xy_yzz_0, tpy_xy_zzz_0, tpy_xyy_xx_0, \
                                     tpy_xyy_xxx_0, tpy_xyy_xxy_0, tpy_xyy_xxz_0, tpy_xyy_xy_0, tpy_xyy_xz_0, \
                                     tpy_xz_xxx_0, tpy_xz_xxy_0, tpy_xz_xxz_0, tpy_xz_xyy_0, tpy_xz_xyz_0, tpy_xz_xzz_0, \
                                     tpy_xz_yyy_0, tpy_xz_yyz_0, tpy_xz_yzz_0, tpy_xz_zzz_0, tpy_yy_xxx_0, tpy_yy_xxy_0, \
                                     tpy_yy_xxz_0, tpz_xxxy_yyy_0, tpz_xxxy_yyz_0, tpz_xxxy_yzz_0, tpz_xxxy_zzz_0, \
                                     tpz_xxxz_xxx_0, tpz_xxxz_xxy_0, tpz_xxxz_xxz_0, tpz_xxxz_xyy_0, tpz_xxxz_xyz_0, \
                                     tpz_xxxz_xzz_0, tpz_xxxz_yyy_0, tpz_xxxz_yyz_0, tpz_xxxz_yzz_0, tpz_xxxz_zzz_0, \
                                     tpz_xxy_yyy_0, tpz_xxy_yyz_0, tpz_xxy_yzz_0, tpz_xxy_zzz_0, tpz_xxyy_xxx_0, \
                                     tpz_xxyy_xxy_0, tpz_xxyy_xxz_0, tpz_xxz_xx_0, tpz_xxz_xxx_0, tpz_xxz_xxy_0, \
                                     tpz_xxz_xxz_0, tpz_xxz_xy_0, tpz_xxz_xyy_0, tpz_xxz_xyz_0, tpz_xxz_xz_0, \
                                     tpz_xxz_xzz_0, tpz_xxz_yy_0, tpz_xxz_yyy_0, tpz_xxz_yyz_0, tpz_xxz_yz_0, \
                                     tpz_xxz_yzz_0, tpz_xxz_zz_0, tpz_xxz_zzz_0, tpz_xy_yyy_0, tpz_xy_yyz_0, \
                                     tpz_xy_yzz_0, tpz_xy_zzz_0, tpz_xyy_xx_0, tpz_xyy_xxx_0, tpz_xyy_xxy_0, \
                                     tpz_xyy_xxz_0, tpz_xyy_xy_0, tpz_xyy_xz_0, tpz_xz_xxx_0, tpz_xz_xxy_0, tpz_xz_xxz_0, \
                                     tpz_xz_xyy_0, tpz_xz_xyz_0, tpz_xz_xzz_0, tpz_xz_yyy_0, tpz_xz_yyz_0, tpz_xz_yzz_0, \
                                     tpz_xz_zzz_0, tpz_yy_xxx_0, tpz_yy_xxy_0, tpz_yy_xxz_0, ts_xxy_yyz_0, ts_xxy_yzz_0, \
                                     ts_xxy_zzz_0, ts_xxz_xxx_0, ts_xxz_xxy_0, ts_xxz_xxz_0, ts_xxz_xyy_0, ts_xxz_xyz_0, \
                                     ts_xxz_xzz_0, ts_xxz_yyy_0, ts_xxz_yyz_0, ts_xxz_yzz_0, ts_xxz_zzz_0, ts_xyy_xxx_0, \
                                     ts_xyy_xxy_0, ts_xyy_xxz_0, ts_xyy_xyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_xxxy_yyy_0[j] = pa_x[j] * tpz_xxy_yyy_0[j] + fl1_fx * tpz_xy_yyy_0[j];

            tpx_xxxy_yyz_0[j] = pa_x[j] * tpx_xxy_yyz_0[j] + fl1_fx * tpx_xy_yyz_0[j] - fl1_fgb * fl1_fx * ts_xxy_yyz_0[j];

            tpy_xxxy_yyz_0[j] = pa_x[j] * tpy_xxy_yyz_0[j] + fl1_fx * tpy_xy_yyz_0[j];

            tpz_xxxy_yyz_0[j] = pa_x[j] * tpz_xxy_yyz_0[j] + fl1_fx * tpz_xy_yyz_0[j];

            tpx_xxxy_yzz_0[j] = pa_x[j] * tpx_xxy_yzz_0[j] + fl1_fx * tpx_xy_yzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_yzz_0[j];

            tpy_xxxy_yzz_0[j] = pa_x[j] * tpy_xxy_yzz_0[j] + fl1_fx * tpy_xy_yzz_0[j];

            tpz_xxxy_yzz_0[j] = pa_x[j] * tpz_xxy_yzz_0[j] + fl1_fx * tpz_xy_yzz_0[j];

            tpx_xxxy_zzz_0[j] = pa_x[j] * tpx_xxy_zzz_0[j] + fl1_fx * tpx_xy_zzz_0[j] - fl1_fgb * fl1_fx * ts_xxy_zzz_0[j];

            tpy_xxxy_zzz_0[j] = pa_x[j] * tpy_xxy_zzz_0[j] + fl1_fx * tpy_xy_zzz_0[j];

            tpz_xxxy_zzz_0[j] = pa_x[j] * tpz_xxy_zzz_0[j] + fl1_fx * tpz_xy_zzz_0[j];

            tpx_xxxz_xxx_0[j] =
                pa_x[j] * tpx_xxz_xxx_0[j] + fl1_fx * tpx_xz_xxx_0[j] + 1.5 * fl1_fx * tpx_xxz_xx_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxx_0[j];

            tpy_xxxz_xxx_0[j] = pa_x[j] * tpy_xxz_xxx_0[j] + fl1_fx * tpy_xz_xxx_0[j] + 1.5 * fl1_fx * tpy_xxz_xx_0[j];

            tpz_xxxz_xxx_0[j] = pa_x[j] * tpz_xxz_xxx_0[j] + fl1_fx * tpz_xz_xxx_0[j] + 1.5 * fl1_fx * tpz_xxz_xx_0[j];

            tpx_xxxz_xxy_0[j] = pa_x[j] * tpx_xxz_xxy_0[j] + fl1_fx * tpx_xz_xxy_0[j] + fl1_fx * tpx_xxz_xy_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxy_0[j];

            tpy_xxxz_xxy_0[j] = pa_x[j] * tpy_xxz_xxy_0[j] + fl1_fx * tpy_xz_xxy_0[j] + fl1_fx * tpy_xxz_xy_0[j];

            tpz_xxxz_xxy_0[j] = pa_x[j] * tpz_xxz_xxy_0[j] + fl1_fx * tpz_xz_xxy_0[j] + fl1_fx * tpz_xxz_xy_0[j];

            tpx_xxxz_xxz_0[j] = pa_x[j] * tpx_xxz_xxz_0[j] + fl1_fx * tpx_xz_xxz_0[j] + fl1_fx * tpx_xxz_xz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xxz_0[j];

            tpy_xxxz_xxz_0[j] = pa_x[j] * tpy_xxz_xxz_0[j] + fl1_fx * tpy_xz_xxz_0[j] + fl1_fx * tpy_xxz_xz_0[j];

            tpz_xxxz_xxz_0[j] = pa_x[j] * tpz_xxz_xxz_0[j] + fl1_fx * tpz_xz_xxz_0[j] + fl1_fx * tpz_xxz_xz_0[j];

            tpx_xxxz_xyy_0[j] =
                pa_x[j] * tpx_xxz_xyy_0[j] + fl1_fx * tpx_xz_xyy_0[j] + 0.5 * fl1_fx * tpx_xxz_yy_0[j] - fl1_fgb * fl1_fx * ts_xxz_xyy_0[j];

            tpy_xxxz_xyy_0[j] = pa_x[j] * tpy_xxz_xyy_0[j] + fl1_fx * tpy_xz_xyy_0[j] + 0.5 * fl1_fx * tpy_xxz_yy_0[j];

            tpz_xxxz_xyy_0[j] = pa_x[j] * tpz_xxz_xyy_0[j] + fl1_fx * tpz_xz_xyy_0[j] + 0.5 * fl1_fx * tpz_xxz_yy_0[j];

            tpx_xxxz_xyz_0[j] =
                pa_x[j] * tpx_xxz_xyz_0[j] + fl1_fx * tpx_xz_xyz_0[j] + 0.5 * fl1_fx * tpx_xxz_yz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xyz_0[j];

            tpy_xxxz_xyz_0[j] = pa_x[j] * tpy_xxz_xyz_0[j] + fl1_fx * tpy_xz_xyz_0[j] + 0.5 * fl1_fx * tpy_xxz_yz_0[j];

            tpz_xxxz_xyz_0[j] = pa_x[j] * tpz_xxz_xyz_0[j] + fl1_fx * tpz_xz_xyz_0[j] + 0.5 * fl1_fx * tpz_xxz_yz_0[j];

            tpx_xxxz_xzz_0[j] =
                pa_x[j] * tpx_xxz_xzz_0[j] + fl1_fx * tpx_xz_xzz_0[j] + 0.5 * fl1_fx * tpx_xxz_zz_0[j] - fl1_fgb * fl1_fx * ts_xxz_xzz_0[j];

            tpy_xxxz_xzz_0[j] = pa_x[j] * tpy_xxz_xzz_0[j] + fl1_fx * tpy_xz_xzz_0[j] + 0.5 * fl1_fx * tpy_xxz_zz_0[j];

            tpz_xxxz_xzz_0[j] = pa_x[j] * tpz_xxz_xzz_0[j] + fl1_fx * tpz_xz_xzz_0[j] + 0.5 * fl1_fx * tpz_xxz_zz_0[j];

            tpx_xxxz_yyy_0[j] = pa_x[j] * tpx_xxz_yyy_0[j] + fl1_fx * tpx_xz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xxz_yyy_0[j];

            tpy_xxxz_yyy_0[j] = pa_x[j] * tpy_xxz_yyy_0[j] + fl1_fx * tpy_xz_yyy_0[j];

            tpz_xxxz_yyy_0[j] = pa_x[j] * tpz_xxz_yyy_0[j] + fl1_fx * tpz_xz_yyy_0[j];

            tpx_xxxz_yyz_0[j] = pa_x[j] * tpx_xxz_yyz_0[j] + fl1_fx * tpx_xz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xxz_yyz_0[j];

            tpy_xxxz_yyz_0[j] = pa_x[j] * tpy_xxz_yyz_0[j] + fl1_fx * tpy_xz_yyz_0[j];

            tpz_xxxz_yyz_0[j] = pa_x[j] * tpz_xxz_yyz_0[j] + fl1_fx * tpz_xz_yyz_0[j];

            tpx_xxxz_yzz_0[j] = pa_x[j] * tpx_xxz_yzz_0[j] + fl1_fx * tpx_xz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_yzz_0[j];

            tpy_xxxz_yzz_0[j] = pa_x[j] * tpy_xxz_yzz_0[j] + fl1_fx * tpy_xz_yzz_0[j];

            tpz_xxxz_yzz_0[j] = pa_x[j] * tpz_xxz_yzz_0[j] + fl1_fx * tpz_xz_yzz_0[j];

            tpx_xxxz_zzz_0[j] = pa_x[j] * tpx_xxz_zzz_0[j] + fl1_fx * tpx_xz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xxz_zzz_0[j];

            tpy_xxxz_zzz_0[j] = pa_x[j] * tpy_xxz_zzz_0[j] + fl1_fx * tpy_xz_zzz_0[j];

            tpz_xxxz_zzz_0[j] = pa_x[j] * tpz_xxz_zzz_0[j] + fl1_fx * tpz_xz_zzz_0[j];

            tpx_xxyy_xxx_0[j] =
                pa_x[j] * tpx_xyy_xxx_0[j] + 0.5 * fl1_fx * tpx_yy_xxx_0[j] + 1.5 * fl1_fx * tpx_xyy_xx_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxx_0[j];

            tpy_xxyy_xxx_0[j] = pa_x[j] * tpy_xyy_xxx_0[j] + 0.5 * fl1_fx * tpy_yy_xxx_0[j] + 1.5 * fl1_fx * tpy_xyy_xx_0[j];

            tpz_xxyy_xxx_0[j] = pa_x[j] * tpz_xyy_xxx_0[j] + 0.5 * fl1_fx * tpz_yy_xxx_0[j] + 1.5 * fl1_fx * tpz_xyy_xx_0[j];

            tpx_xxyy_xxy_0[j] =
                pa_x[j] * tpx_xyy_xxy_0[j] + 0.5 * fl1_fx * tpx_yy_xxy_0[j] + fl1_fx * tpx_xyy_xy_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxy_0[j];

            tpy_xxyy_xxy_0[j] = pa_x[j] * tpy_xyy_xxy_0[j] + 0.5 * fl1_fx * tpy_yy_xxy_0[j] + fl1_fx * tpy_xyy_xy_0[j];

            tpz_xxyy_xxy_0[j] = pa_x[j] * tpz_xyy_xxy_0[j] + 0.5 * fl1_fx * tpz_yy_xxy_0[j] + fl1_fx * tpz_xyy_xy_0[j];

            tpx_xxyy_xxz_0[j] =
                pa_x[j] * tpx_xyy_xxz_0[j] + 0.5 * fl1_fx * tpx_yy_xxz_0[j] + fl1_fx * tpx_xyy_xz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xxz_0[j];

            tpy_xxyy_xxz_0[j] = pa_x[j] * tpy_xyy_xxz_0[j] + 0.5 * fl1_fx * tpy_yy_xxz_0[j] + fl1_fx * tpy_xyy_xz_0[j];

            tpz_xxyy_xxz_0[j] = pa_x[j] * tpz_xyy_xxz_0[j] + 0.5 * fl1_fx * tpz_yy_xxz_0[j] + fl1_fx * tpz_xyy_xz_0[j];

            tpx_xxyy_xyy_0[j] =
                pa_x[j] * tpx_xyy_xyy_0[j] + 0.5 * fl1_fx * tpx_yy_xyy_0[j] + 0.5 * fl1_fx * tpx_xyy_yy_0[j] - fl1_fgb * fl1_fx * ts_xyy_xyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_100_150(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpy_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 21);

        auto tpz_xyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 21);

        auto tpx_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 22);

        auto tpy_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 22);

        auto tpz_xyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 22);

        auto tpx_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 23);

        auto tpy_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 23);

        auto tpz_xyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 23);

        auto tpx_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 24);

        auto tpy_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 24);

        auto tpz_xyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 24);

        auto tpx_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 25);

        auto tpy_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 25);

        auto tpz_xyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 25);

        auto tpx_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 26);

        auto tpy_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 26);

        auto tpz_xyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 26);

        auto tpx_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 27);

        auto tpy_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 27);

        auto tpz_xyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 27);

        auto tpx_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 28);

        auto tpy_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 28);

        auto tpz_xyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 28);

        auto tpx_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 29);

        auto tpy_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 29);

        auto tpz_xyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 29);

        auto ts_xyy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 34);

        auto ts_xyy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 35);

        auto ts_xyy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 36);

        auto ts_xyy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 37);

        auto ts_xyy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 38);

        auto ts_xyy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 39);

        auto ts_xyz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 40);

        auto ts_xyz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 41);

        auto ts_xyz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 42);

        auto ts_xyz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 43);

        auto ts_xyz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 44);

        auto ts_xyz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 45);

        auto ts_xyz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 46);

        auto ts_xyz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 47);

        auto ts_xyz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 48);

        auto ts_xyz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 49);

        // set up pointers to integrals

        auto tpy_xxyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 33);

        auto tpz_xxyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 33);

        auto tpx_xxyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 34);

        auto tpy_xxyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 34);

        auto tpz_xxyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 34);

        auto tpx_xxyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 35);

        auto tpy_xxyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 35);

        auto tpz_xxyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 35);

        auto tpx_xxyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 36);

        auto tpy_xxyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 36);

        auto tpz_xxyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 36);

        auto tpx_xxyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 37);

        auto tpy_xxyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 37);

        auto tpz_xxyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 37);

        auto tpx_xxyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 38);

        auto tpy_xxyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 38);

        auto tpz_xxyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 38);

        auto tpx_xxyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 39);

        auto tpy_xxyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 39);

        auto tpz_xxyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 39);

        auto tpx_xxyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 40);

        auto tpy_xxyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 40);

        auto tpz_xxyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 40);

        auto tpx_xxyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 41);

        auto tpy_xxyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 41);

        auto tpz_xxyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 41);

        auto tpx_xxyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 42);

        auto tpy_xxyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 42);

        auto tpz_xxyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 42);

        auto tpx_xxyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 43);

        auto tpy_xxyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 43);

        auto tpz_xxyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 43);

        auto tpx_xxyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 44);

        auto tpy_xxyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 44);

        auto tpz_xxyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 44);

        auto tpx_xxyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 45);

        auto tpy_xxyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 45);

        auto tpz_xxyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 45);

        auto tpx_xxyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 46);

        auto tpy_xxyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 46);

        auto tpz_xxyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 46);

        auto tpx_xxyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 47);

        auto tpy_xxyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 47);

        auto tpz_xxyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 47);

        auto tpx_xxyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 48);

        auto tpy_xxyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 48);

        auto tpz_xxyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 48);

        auto tpx_xxyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 49);

        auto tpy_xxyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 49);

        auto tpz_xxyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 49);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxyy_xyz_0, tpx_xxyy_xzz_0, tpx_xxyy_yyy_0, \
                                     tpx_xxyy_yyz_0, tpx_xxyy_yzz_0, tpx_xxyy_zzz_0, tpx_xxyz_xxx_0, tpx_xxyz_xxy_0, \
                                     tpx_xxyz_xxz_0, tpx_xxyz_xyy_0, tpx_xxyz_xyz_0, tpx_xxyz_xzz_0, tpx_xxyz_yyy_0, \
                                     tpx_xxyz_yyz_0, tpx_xxyz_yzz_0, tpx_xxyz_zzz_0, tpx_xyy_xyz_0, tpx_xyy_xzz_0, \
                                     tpx_xyy_yyy_0, tpx_xyy_yyz_0, tpx_xyy_yz_0, tpx_xyy_yzz_0, tpx_xyy_zz_0, \
                                     tpx_xyy_zzz_0, tpx_xyz_xx_0, tpx_xyz_xxx_0, tpx_xyz_xxy_0, tpx_xyz_xxz_0, \
                                     tpx_xyz_xy_0, tpx_xyz_xyy_0, tpx_xyz_xyz_0, tpx_xyz_xz_0, tpx_xyz_xzz_0, \
                                     tpx_xyz_yy_0, tpx_xyz_yyy_0, tpx_xyz_yyz_0, tpx_xyz_yz_0, tpx_xyz_yzz_0, \
                                     tpx_xyz_zz_0, tpx_xyz_zzz_0, tpx_yy_xyz_0, tpx_yy_xzz_0, tpx_yy_yyy_0, \
                                     tpx_yy_yyz_0, tpx_yy_yzz_0, tpx_yy_zzz_0, tpx_yz_xxx_0, tpx_yz_xxy_0, tpx_yz_xxz_0, \
                                     tpx_yz_xyy_0, tpx_yz_xyz_0, tpx_yz_xzz_0, tpx_yz_yyy_0, tpx_yz_yyz_0, tpx_yz_yzz_0, \
                                     tpx_yz_zzz_0, tpy_xxyy_xyy_0, tpy_xxyy_xyz_0, tpy_xxyy_xzz_0, tpy_xxyy_yyy_0, \
                                     tpy_xxyy_yyz_0, tpy_xxyy_yzz_0, tpy_xxyy_zzz_0, tpy_xxyz_xxx_0, tpy_xxyz_xxy_0, \
                                     tpy_xxyz_xxz_0, tpy_xxyz_xyy_0, tpy_xxyz_xyz_0, tpy_xxyz_xzz_0, tpy_xxyz_yyy_0, \
                                     tpy_xxyz_yyz_0, tpy_xxyz_yzz_0, tpy_xxyz_zzz_0, tpy_xyy_xyy_0, tpy_xyy_xyz_0, \
                                     tpy_xyy_xzz_0, tpy_xyy_yy_0, tpy_xyy_yyy_0, tpy_xyy_yyz_0, tpy_xyy_yz_0, \
                                     tpy_xyy_yzz_0, tpy_xyy_zz_0, tpy_xyy_zzz_0, tpy_xyz_xx_0, tpy_xyz_xxx_0, \
                                     tpy_xyz_xxy_0, tpy_xyz_xxz_0, tpy_xyz_xy_0, tpy_xyz_xyy_0, tpy_xyz_xyz_0, \
                                     tpy_xyz_xz_0, tpy_xyz_xzz_0, tpy_xyz_yy_0, tpy_xyz_yyy_0, tpy_xyz_yyz_0, \
                                     tpy_xyz_yz_0, tpy_xyz_yzz_0, tpy_xyz_zz_0, tpy_xyz_zzz_0, tpy_yy_xyy_0, \
                                     tpy_yy_xyz_0, tpy_yy_xzz_0, tpy_yy_yyy_0, tpy_yy_yyz_0, tpy_yy_yzz_0, tpy_yy_zzz_0, \
                                     tpy_yz_xxx_0, tpy_yz_xxy_0, tpy_yz_xxz_0, tpy_yz_xyy_0, tpy_yz_xyz_0, tpy_yz_xzz_0, \
                                     tpy_yz_yyy_0, tpy_yz_yyz_0, tpy_yz_yzz_0, tpy_yz_zzz_0, tpz_xxyy_xyy_0, \
                                     tpz_xxyy_xyz_0, tpz_xxyy_xzz_0, tpz_xxyy_yyy_0, tpz_xxyy_yyz_0, tpz_xxyy_yzz_0, \
                                     tpz_xxyy_zzz_0, tpz_xxyz_xxx_0, tpz_xxyz_xxy_0, tpz_xxyz_xxz_0, tpz_xxyz_xyy_0, \
                                     tpz_xxyz_xyz_0, tpz_xxyz_xzz_0, tpz_xxyz_yyy_0, tpz_xxyz_yyz_0, tpz_xxyz_yzz_0, \
                                     tpz_xxyz_zzz_0, tpz_xyy_xyy_0, tpz_xyy_xyz_0, tpz_xyy_xzz_0, tpz_xyy_yy_0, \
                                     tpz_xyy_yyy_0, tpz_xyy_yyz_0, tpz_xyy_yz_0, tpz_xyy_yzz_0, tpz_xyy_zz_0, \
                                     tpz_xyy_zzz_0, tpz_xyz_xx_0, tpz_xyz_xxx_0, tpz_xyz_xxy_0, tpz_xyz_xxz_0, \
                                     tpz_xyz_xy_0, tpz_xyz_xyy_0, tpz_xyz_xyz_0, tpz_xyz_xz_0, tpz_xyz_xzz_0, \
                                     tpz_xyz_yy_0, tpz_xyz_yyy_0, tpz_xyz_yyz_0, tpz_xyz_yz_0, tpz_xyz_yzz_0, \
                                     tpz_xyz_zz_0, tpz_xyz_zzz_0, tpz_yy_xyy_0, tpz_yy_xyz_0, tpz_yy_xzz_0, \
                                     tpz_yy_yyy_0, tpz_yy_yyz_0, tpz_yy_yzz_0, tpz_yy_zzz_0, tpz_yz_xxx_0, tpz_yz_xxy_0, \
                                     tpz_yz_xxz_0, tpz_yz_xyy_0, tpz_yz_xyz_0, tpz_yz_xzz_0, tpz_yz_yyy_0, tpz_yz_yyz_0, \
                                     tpz_yz_yzz_0, tpz_yz_zzz_0, ts_xyy_xyz_0, ts_xyy_xzz_0, ts_xyy_yyy_0, ts_xyy_yyz_0, \
                                     ts_xyy_yzz_0, ts_xyy_zzz_0, ts_xyz_xxx_0, ts_xyz_xxy_0, ts_xyz_xxz_0, ts_xyz_xyy_0, \
                                     ts_xyz_xyz_0, ts_xyz_xzz_0, ts_xyz_yyy_0, ts_xyz_yyz_0, ts_xyz_yzz_0, ts_xyz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_xxyy_xyy_0[j] = pa_x[j] * tpy_xyy_xyy_0[j] + 0.5 * fl1_fx * tpy_yy_xyy_0[j] + 0.5 * fl1_fx * tpy_xyy_yy_0[j];

            tpz_xxyy_xyy_0[j] = pa_x[j] * tpz_xyy_xyy_0[j] + 0.5 * fl1_fx * tpz_yy_xyy_0[j] + 0.5 * fl1_fx * tpz_xyy_yy_0[j];

            tpx_xxyy_xyz_0[j] =
                pa_x[j] * tpx_xyy_xyz_0[j] + 0.5 * fl1_fx * tpx_yy_xyz_0[j] + 0.5 * fl1_fx * tpx_xyy_yz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xyz_0[j];

            tpy_xxyy_xyz_0[j] = pa_x[j] * tpy_xyy_xyz_0[j] + 0.5 * fl1_fx * tpy_yy_xyz_0[j] + 0.5 * fl1_fx * tpy_xyy_yz_0[j];

            tpz_xxyy_xyz_0[j] = pa_x[j] * tpz_xyy_xyz_0[j] + 0.5 * fl1_fx * tpz_yy_xyz_0[j] + 0.5 * fl1_fx * tpz_xyy_yz_0[j];

            tpx_xxyy_xzz_0[j] =
                pa_x[j] * tpx_xyy_xzz_0[j] + 0.5 * fl1_fx * tpx_yy_xzz_0[j] + 0.5 * fl1_fx * tpx_xyy_zz_0[j] - fl1_fgb * fl1_fx * ts_xyy_xzz_0[j];

            tpy_xxyy_xzz_0[j] = pa_x[j] * tpy_xyy_xzz_0[j] + 0.5 * fl1_fx * tpy_yy_xzz_0[j] + 0.5 * fl1_fx * tpy_xyy_zz_0[j];

            tpz_xxyy_xzz_0[j] = pa_x[j] * tpz_xyy_xzz_0[j] + 0.5 * fl1_fx * tpz_yy_xzz_0[j] + 0.5 * fl1_fx * tpz_xyy_zz_0[j];

            tpx_xxyy_yyy_0[j] = pa_x[j] * tpx_xyy_yyy_0[j] + 0.5 * fl1_fx * tpx_yy_yyy_0[j] - fl1_fgb * fl1_fx * ts_xyy_yyy_0[j];

            tpy_xxyy_yyy_0[j] = pa_x[j] * tpy_xyy_yyy_0[j] + 0.5 * fl1_fx * tpy_yy_yyy_0[j];

            tpz_xxyy_yyy_0[j] = pa_x[j] * tpz_xyy_yyy_0[j] + 0.5 * fl1_fx * tpz_yy_yyy_0[j];

            tpx_xxyy_yyz_0[j] = pa_x[j] * tpx_xyy_yyz_0[j] + 0.5 * fl1_fx * tpx_yy_yyz_0[j] - fl1_fgb * fl1_fx * ts_xyy_yyz_0[j];

            tpy_xxyy_yyz_0[j] = pa_x[j] * tpy_xyy_yyz_0[j] + 0.5 * fl1_fx * tpy_yy_yyz_0[j];

            tpz_xxyy_yyz_0[j] = pa_x[j] * tpz_xyy_yyz_0[j] + 0.5 * fl1_fx * tpz_yy_yyz_0[j];

            tpx_xxyy_yzz_0[j] = pa_x[j] * tpx_xyy_yzz_0[j] + 0.5 * fl1_fx * tpx_yy_yzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_yzz_0[j];

            tpy_xxyy_yzz_0[j] = pa_x[j] * tpy_xyy_yzz_0[j] + 0.5 * fl1_fx * tpy_yy_yzz_0[j];

            tpz_xxyy_yzz_0[j] = pa_x[j] * tpz_xyy_yzz_0[j] + 0.5 * fl1_fx * tpz_yy_yzz_0[j];

            tpx_xxyy_zzz_0[j] = pa_x[j] * tpx_xyy_zzz_0[j] + 0.5 * fl1_fx * tpx_yy_zzz_0[j] - fl1_fgb * fl1_fx * ts_xyy_zzz_0[j];

            tpy_xxyy_zzz_0[j] = pa_x[j] * tpy_xyy_zzz_0[j] + 0.5 * fl1_fx * tpy_yy_zzz_0[j];

            tpz_xxyy_zzz_0[j] = pa_x[j] * tpz_xyy_zzz_0[j] + 0.5 * fl1_fx * tpz_yy_zzz_0[j];

            tpx_xxyz_xxx_0[j] =
                pa_x[j] * tpx_xyz_xxx_0[j] + 0.5 * fl1_fx * tpx_yz_xxx_0[j] + 1.5 * fl1_fx * tpx_xyz_xx_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxx_0[j];

            tpy_xxyz_xxx_0[j] = pa_x[j] * tpy_xyz_xxx_0[j] + 0.5 * fl1_fx * tpy_yz_xxx_0[j] + 1.5 * fl1_fx * tpy_xyz_xx_0[j];

            tpz_xxyz_xxx_0[j] = pa_x[j] * tpz_xyz_xxx_0[j] + 0.5 * fl1_fx * tpz_yz_xxx_0[j] + 1.5 * fl1_fx * tpz_xyz_xx_0[j];

            tpx_xxyz_xxy_0[j] =
                pa_x[j] * tpx_xyz_xxy_0[j] + 0.5 * fl1_fx * tpx_yz_xxy_0[j] + fl1_fx * tpx_xyz_xy_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxy_0[j];

            tpy_xxyz_xxy_0[j] = pa_x[j] * tpy_xyz_xxy_0[j] + 0.5 * fl1_fx * tpy_yz_xxy_0[j] + fl1_fx * tpy_xyz_xy_0[j];

            tpz_xxyz_xxy_0[j] = pa_x[j] * tpz_xyz_xxy_0[j] + 0.5 * fl1_fx * tpz_yz_xxy_0[j] + fl1_fx * tpz_xyz_xy_0[j];

            tpx_xxyz_xxz_0[j] =
                pa_x[j] * tpx_xyz_xxz_0[j] + 0.5 * fl1_fx * tpx_yz_xxz_0[j] + fl1_fx * tpx_xyz_xz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xxz_0[j];

            tpy_xxyz_xxz_0[j] = pa_x[j] * tpy_xyz_xxz_0[j] + 0.5 * fl1_fx * tpy_yz_xxz_0[j] + fl1_fx * tpy_xyz_xz_0[j];

            tpz_xxyz_xxz_0[j] = pa_x[j] * tpz_xyz_xxz_0[j] + 0.5 * fl1_fx * tpz_yz_xxz_0[j] + fl1_fx * tpz_xyz_xz_0[j];

            tpx_xxyz_xyy_0[j] =
                pa_x[j] * tpx_xyz_xyy_0[j] + 0.5 * fl1_fx * tpx_yz_xyy_0[j] + 0.5 * fl1_fx * tpx_xyz_yy_0[j] - fl1_fgb * fl1_fx * ts_xyz_xyy_0[j];

            tpy_xxyz_xyy_0[j] = pa_x[j] * tpy_xyz_xyy_0[j] + 0.5 * fl1_fx * tpy_yz_xyy_0[j] + 0.5 * fl1_fx * tpy_xyz_yy_0[j];

            tpz_xxyz_xyy_0[j] = pa_x[j] * tpz_xyz_xyy_0[j] + 0.5 * fl1_fx * tpz_yz_xyy_0[j] + 0.5 * fl1_fx * tpz_xyz_yy_0[j];

            tpx_xxyz_xyz_0[j] =
                pa_x[j] * tpx_xyz_xyz_0[j] + 0.5 * fl1_fx * tpx_yz_xyz_0[j] + 0.5 * fl1_fx * tpx_xyz_yz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xyz_0[j];

            tpy_xxyz_xyz_0[j] = pa_x[j] * tpy_xyz_xyz_0[j] + 0.5 * fl1_fx * tpy_yz_xyz_0[j] + 0.5 * fl1_fx * tpy_xyz_yz_0[j];

            tpz_xxyz_xyz_0[j] = pa_x[j] * tpz_xyz_xyz_0[j] + 0.5 * fl1_fx * tpz_yz_xyz_0[j] + 0.5 * fl1_fx * tpz_xyz_yz_0[j];

            tpx_xxyz_xzz_0[j] =
                pa_x[j] * tpx_xyz_xzz_0[j] + 0.5 * fl1_fx * tpx_yz_xzz_0[j] + 0.5 * fl1_fx * tpx_xyz_zz_0[j] - fl1_fgb * fl1_fx * ts_xyz_xzz_0[j];

            tpy_xxyz_xzz_0[j] = pa_x[j] * tpy_xyz_xzz_0[j] + 0.5 * fl1_fx * tpy_yz_xzz_0[j] + 0.5 * fl1_fx * tpy_xyz_zz_0[j];

            tpz_xxyz_xzz_0[j] = pa_x[j] * tpz_xyz_xzz_0[j] + 0.5 * fl1_fx * tpz_yz_xzz_0[j] + 0.5 * fl1_fx * tpz_xyz_zz_0[j];

            tpx_xxyz_yyy_0[j] = pa_x[j] * tpx_xyz_yyy_0[j] + 0.5 * fl1_fx * tpx_yz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xyz_yyy_0[j];

            tpy_xxyz_yyy_0[j] = pa_x[j] * tpy_xyz_yyy_0[j] + 0.5 * fl1_fx * tpy_yz_yyy_0[j];

            tpz_xxyz_yyy_0[j] = pa_x[j] * tpz_xyz_yyy_0[j] + 0.5 * fl1_fx * tpz_yz_yyy_0[j];

            tpx_xxyz_yyz_0[j] = pa_x[j] * tpx_xyz_yyz_0[j] + 0.5 * fl1_fx * tpx_yz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xyz_yyz_0[j];

            tpy_xxyz_yyz_0[j] = pa_x[j] * tpy_xyz_yyz_0[j] + 0.5 * fl1_fx * tpy_yz_yyz_0[j];

            tpz_xxyz_yyz_0[j] = pa_x[j] * tpz_xyz_yyz_0[j] + 0.5 * fl1_fx * tpz_yz_yyz_0[j];

            tpx_xxyz_yzz_0[j] = pa_x[j] * tpx_xyz_yzz_0[j] + 0.5 * fl1_fx * tpx_yz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_yzz_0[j];

            tpy_xxyz_yzz_0[j] = pa_x[j] * tpy_xyz_yzz_0[j] + 0.5 * fl1_fx * tpy_yz_yzz_0[j];

            tpz_xxyz_yzz_0[j] = pa_x[j] * tpz_xyz_yzz_0[j] + 0.5 * fl1_fx * tpz_yz_yzz_0[j];

            tpx_xxyz_zzz_0[j] = pa_x[j] * tpx_xyz_zzz_0[j] + 0.5 * fl1_fx * tpx_yz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xyz_zzz_0[j];

            tpy_xxyz_zzz_0[j] = pa_x[j] * tpy_xyz_zzz_0[j] + 0.5 * fl1_fx * tpy_yz_zzz_0[j];

            tpz_xxyz_zzz_0[j] = pa_x[j] * tpz_xyz_zzz_0[j] + 0.5 * fl1_fx * tpz_yz_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_150_200(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 30);

        auto tpy_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 30);

        auto tpz_xzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 30);

        auto tpx_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 31);

        auto tpy_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 31);

        auto tpz_xzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 31);

        auto tpx_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 32);

        auto tpy_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 32);

        auto tpz_xzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 32);

        auto tpx_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 33);

        auto tpy_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 33);

        auto tpz_xzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 33);

        auto tpx_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 34);

        auto tpy_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 34);

        auto tpz_xzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 34);

        auto tpx_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 35);

        auto tpy_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 35);

        auto tpz_xzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 35);

        auto tpx_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 36);

        auto tpy_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 37);

        auto tpy_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 38);

        auto tpy_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 39);

        auto tpy_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 40);

        auto tpy_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 41);

        auto tpy_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto ts_xzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 50);

        auto ts_xzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 51);

        auto ts_xzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 52);

        auto ts_xzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 53);

        auto ts_xzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 54);

        auto ts_xzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 55);

        auto ts_xzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 56);

        auto ts_xzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 57);

        auto ts_xzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 58);

        auto ts_xzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 59);

        auto ts_yyy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 60);

        auto ts_yyy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 61);

        auto ts_yyy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 62);

        auto ts_yyy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 63);

        auto ts_yyy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 64);

        auto ts_yyy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 65);

        auto ts_yyy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 66);

        // set up pointers to integrals

        auto tpx_xxzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 50);

        auto tpy_xxzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 50);

        auto tpz_xxzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 50);

        auto tpx_xxzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 51);

        auto tpy_xxzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 51);

        auto tpz_xxzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 51);

        auto tpx_xxzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 52);

        auto tpy_xxzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 52);

        auto tpz_xxzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 52);

        auto tpx_xxzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 53);

        auto tpy_xxzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 53);

        auto tpz_xxzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 53);

        auto tpx_xxzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 54);

        auto tpy_xxzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 54);

        auto tpz_xxzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 54);

        auto tpx_xxzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 55);

        auto tpy_xxzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 55);

        auto tpz_xxzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 55);

        auto tpx_xxzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 56);

        auto tpy_xxzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 56);

        auto tpz_xxzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 56);

        auto tpx_xxzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 57);

        auto tpy_xxzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 57);

        auto tpz_xxzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 57);

        auto tpx_xxzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 58);

        auto tpy_xxzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 58);

        auto tpz_xxzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 58);

        auto tpx_xxzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 59);

        auto tpy_xxzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 59);

        auto tpz_xxzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 59);

        auto tpx_xyyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 60);

        auto tpy_xyyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 60);

        auto tpz_xyyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 60);

        auto tpx_xyyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 61);

        auto tpy_xyyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 61);

        auto tpz_xyyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 61);

        auto tpx_xyyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 62);

        auto tpy_xyyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 62);

        auto tpz_xyyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 62);

        auto tpx_xyyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 63);

        auto tpy_xyyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 63);

        auto tpz_xyyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 63);

        auto tpx_xyyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 64);

        auto tpy_xyyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 64);

        auto tpz_xyyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 64);

        auto tpx_xyyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 65);

        auto tpy_xyyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 65);

        auto tpz_xyyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 65);

        auto tpx_xyyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 66);

        auto tpy_xyyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 66);

        // Batch of Integrals (150,200)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxzz_xxx_0, tpx_xxzz_xxy_0, tpx_xxzz_xxz_0, \
                                     tpx_xxzz_xyy_0, tpx_xxzz_xyz_0, tpx_xxzz_xzz_0, tpx_xxzz_yyy_0, tpx_xxzz_yyz_0, \
                                     tpx_xxzz_yzz_0, tpx_xxzz_zzz_0, tpx_xyyy_xxx_0, tpx_xyyy_xxy_0, tpx_xyyy_xxz_0, \
                                     tpx_xyyy_xyy_0, tpx_xyyy_xyz_0, tpx_xyyy_xzz_0, tpx_xyyy_yyy_0, tpx_xzz_xx_0, \
                                     tpx_xzz_xxx_0, tpx_xzz_xxy_0, tpx_xzz_xxz_0, tpx_xzz_xy_0, tpx_xzz_xyy_0, \
                                     tpx_xzz_xyz_0, tpx_xzz_xz_0, tpx_xzz_xzz_0, tpx_xzz_yy_0, tpx_xzz_yyy_0, \
                                     tpx_xzz_yyz_0, tpx_xzz_yz_0, tpx_xzz_yzz_0, tpx_xzz_zz_0, tpx_xzz_zzz_0, \
                                     tpx_yyy_xx_0, tpx_yyy_xxx_0, tpx_yyy_xxy_0, tpx_yyy_xxz_0, tpx_yyy_xy_0, \
                                     tpx_yyy_xyy_0, tpx_yyy_xyz_0, tpx_yyy_xz_0, tpx_yyy_xzz_0, tpx_yyy_yy_0, \
                                     tpx_yyy_yyy_0, tpx_yyy_yz_0, tpx_yyy_zz_0, tpx_zz_xxx_0, tpx_zz_xxy_0, tpx_zz_xxz_0, \
                                     tpx_zz_xyy_0, tpx_zz_xyz_0, tpx_zz_xzz_0, tpx_zz_yyy_0, tpx_zz_yyz_0, tpx_zz_yzz_0, \
                                     tpx_zz_zzz_0, tpy_xxzz_xxx_0, tpy_xxzz_xxy_0, tpy_xxzz_xxz_0, tpy_xxzz_xyy_0, \
                                     tpy_xxzz_xyz_0, tpy_xxzz_xzz_0, tpy_xxzz_yyy_0, tpy_xxzz_yyz_0, tpy_xxzz_yzz_0, \
                                     tpy_xxzz_zzz_0, tpy_xyyy_xxx_0, tpy_xyyy_xxy_0, tpy_xyyy_xxz_0, tpy_xyyy_xyy_0, \
                                     tpy_xyyy_xyz_0, tpy_xyyy_xzz_0, tpy_xyyy_yyy_0, tpy_xzz_xx_0, tpy_xzz_xxx_0, \
                                     tpy_xzz_xxy_0, tpy_xzz_xxz_0, tpy_xzz_xy_0, tpy_xzz_xyy_0, tpy_xzz_xyz_0, \
                                     tpy_xzz_xz_0, tpy_xzz_xzz_0, tpy_xzz_yy_0, tpy_xzz_yyy_0, tpy_xzz_yyz_0, \
                                     tpy_xzz_yz_0, tpy_xzz_yzz_0, tpy_xzz_zz_0, tpy_xzz_zzz_0, tpy_yyy_xx_0, \
                                     tpy_yyy_xxx_0, tpy_yyy_xxy_0, tpy_yyy_xxz_0, tpy_yyy_xy_0, tpy_yyy_xyy_0, \
                                     tpy_yyy_xyz_0, tpy_yyy_xz_0, tpy_yyy_xzz_0, tpy_yyy_yy_0, tpy_yyy_yyy_0, \
                                     tpy_yyy_yz_0, tpy_yyy_zz_0, tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, tpy_zz_xyy_0, \
                                     tpy_zz_xyz_0, tpy_zz_xzz_0, tpy_zz_yyy_0, tpy_zz_yyz_0, tpy_zz_yzz_0, tpy_zz_zzz_0, \
                                     tpz_xxzz_xxx_0, tpz_xxzz_xxy_0, tpz_xxzz_xxz_0, tpz_xxzz_xyy_0, tpz_xxzz_xyz_0, \
                                     tpz_xxzz_xzz_0, tpz_xxzz_yyy_0, tpz_xxzz_yyz_0, tpz_xxzz_yzz_0, tpz_xxzz_zzz_0, \
                                     tpz_xyyy_xxx_0, tpz_xyyy_xxy_0, tpz_xyyy_xxz_0, tpz_xyyy_xyy_0, tpz_xyyy_xyz_0, \
                                     tpz_xyyy_xzz_0, tpz_xzz_xx_0, tpz_xzz_xxx_0, tpz_xzz_xxy_0, tpz_xzz_xxz_0, \
                                     tpz_xzz_xy_0, tpz_xzz_xyy_0, tpz_xzz_xyz_0, tpz_xzz_xz_0, tpz_xzz_xzz_0, \
                                     tpz_xzz_yy_0, tpz_xzz_yyy_0, tpz_xzz_yyz_0, tpz_xzz_yz_0, tpz_xzz_yzz_0, \
                                     tpz_xzz_zz_0, tpz_xzz_zzz_0, tpz_yyy_xx_0, tpz_yyy_xxx_0, tpz_yyy_xxy_0, \
                                     tpz_yyy_xxz_0, tpz_yyy_xy_0, tpz_yyy_xyy_0, tpz_yyy_xyz_0, tpz_yyy_xz_0, \
                                     tpz_yyy_xzz_0, tpz_yyy_yy_0, tpz_yyy_yz_0, tpz_yyy_zz_0, tpz_zz_xxx_0, tpz_zz_xxy_0, \
                                     tpz_zz_xxz_0, tpz_zz_xyy_0, tpz_zz_xyz_0, tpz_zz_xzz_0, tpz_zz_yyy_0, tpz_zz_yyz_0, \
                                     tpz_zz_yzz_0, tpz_zz_zzz_0, ts_xzz_xxx_0, ts_xzz_xxy_0, ts_xzz_xxz_0, ts_xzz_xyy_0, \
                                     ts_xzz_xyz_0, ts_xzz_xzz_0, ts_xzz_yyy_0, ts_xzz_yyz_0, ts_xzz_yzz_0, ts_xzz_zzz_0, \
                                     ts_yyy_xxx_0, ts_yyy_xxy_0, ts_yyy_xxz_0, ts_yyy_xyy_0, ts_yyy_xyz_0, ts_yyy_xzz_0, \
                                     ts_yyy_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxzz_xxx_0[j] =
                pa_x[j] * tpx_xzz_xxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxx_0[j] + 1.5 * fl1_fx * tpx_xzz_xx_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxx_0[j];

            tpy_xxzz_xxx_0[j] = pa_x[j] * tpy_xzz_xxx_0[j] + 0.5 * fl1_fx * tpy_zz_xxx_0[j] + 1.5 * fl1_fx * tpy_xzz_xx_0[j];

            tpz_xxzz_xxx_0[j] = pa_x[j] * tpz_xzz_xxx_0[j] + 0.5 * fl1_fx * tpz_zz_xxx_0[j] + 1.5 * fl1_fx * tpz_xzz_xx_0[j];

            tpx_xxzz_xxy_0[j] =
                pa_x[j] * tpx_xzz_xxy_0[j] + 0.5 * fl1_fx * tpx_zz_xxy_0[j] + fl1_fx * tpx_xzz_xy_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxy_0[j];

            tpy_xxzz_xxy_0[j] = pa_x[j] * tpy_xzz_xxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxy_0[j] + fl1_fx * tpy_xzz_xy_0[j];

            tpz_xxzz_xxy_0[j] = pa_x[j] * tpz_xzz_xxy_0[j] + 0.5 * fl1_fx * tpz_zz_xxy_0[j] + fl1_fx * tpz_xzz_xy_0[j];

            tpx_xxzz_xxz_0[j] =
                pa_x[j] * tpx_xzz_xxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxz_0[j] + fl1_fx * tpx_xzz_xz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xxz_0[j];

            tpy_xxzz_xxz_0[j] = pa_x[j] * tpy_xzz_xxz_0[j] + 0.5 * fl1_fx * tpy_zz_xxz_0[j] + fl1_fx * tpy_xzz_xz_0[j];

            tpz_xxzz_xxz_0[j] = pa_x[j] * tpz_xzz_xxz_0[j] + 0.5 * fl1_fx * tpz_zz_xxz_0[j] + fl1_fx * tpz_xzz_xz_0[j];

            tpx_xxzz_xyy_0[j] =
                pa_x[j] * tpx_xzz_xyy_0[j] + 0.5 * fl1_fx * tpx_zz_xyy_0[j] + 0.5 * fl1_fx * tpx_xzz_yy_0[j] - fl1_fgb * fl1_fx * ts_xzz_xyy_0[j];

            tpy_xxzz_xyy_0[j] = pa_x[j] * tpy_xzz_xyy_0[j] + 0.5 * fl1_fx * tpy_zz_xyy_0[j] + 0.5 * fl1_fx * tpy_xzz_yy_0[j];

            tpz_xxzz_xyy_0[j] = pa_x[j] * tpz_xzz_xyy_0[j] + 0.5 * fl1_fx * tpz_zz_xyy_0[j] + 0.5 * fl1_fx * tpz_xzz_yy_0[j];

            tpx_xxzz_xyz_0[j] =
                pa_x[j] * tpx_xzz_xyz_0[j] + 0.5 * fl1_fx * tpx_zz_xyz_0[j] + 0.5 * fl1_fx * tpx_xzz_yz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xyz_0[j];

            tpy_xxzz_xyz_0[j] = pa_x[j] * tpy_xzz_xyz_0[j] + 0.5 * fl1_fx * tpy_zz_xyz_0[j] + 0.5 * fl1_fx * tpy_xzz_yz_0[j];

            tpz_xxzz_xyz_0[j] = pa_x[j] * tpz_xzz_xyz_0[j] + 0.5 * fl1_fx * tpz_zz_xyz_0[j] + 0.5 * fl1_fx * tpz_xzz_yz_0[j];

            tpx_xxzz_xzz_0[j] =
                pa_x[j] * tpx_xzz_xzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzz_0[j] + 0.5 * fl1_fx * tpx_xzz_zz_0[j] - fl1_fgb * fl1_fx * ts_xzz_xzz_0[j];

            tpy_xxzz_xzz_0[j] = pa_x[j] * tpy_xzz_xzz_0[j] + 0.5 * fl1_fx * tpy_zz_xzz_0[j] + 0.5 * fl1_fx * tpy_xzz_zz_0[j];

            tpz_xxzz_xzz_0[j] = pa_x[j] * tpz_xzz_xzz_0[j] + 0.5 * fl1_fx * tpz_zz_xzz_0[j] + 0.5 * fl1_fx * tpz_xzz_zz_0[j];

            tpx_xxzz_yyy_0[j] = pa_x[j] * tpx_xzz_yyy_0[j] + 0.5 * fl1_fx * tpx_zz_yyy_0[j] - fl1_fgb * fl1_fx * ts_xzz_yyy_0[j];

            tpy_xxzz_yyy_0[j] = pa_x[j] * tpy_xzz_yyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyy_0[j];

            tpz_xxzz_yyy_0[j] = pa_x[j] * tpz_xzz_yyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyy_0[j];

            tpx_xxzz_yyz_0[j] = pa_x[j] * tpx_xzz_yyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyz_0[j] - fl1_fgb * fl1_fx * ts_xzz_yyz_0[j];

            tpy_xxzz_yyz_0[j] = pa_x[j] * tpy_xzz_yyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyz_0[j];

            tpz_xxzz_yyz_0[j] = pa_x[j] * tpz_xzz_yyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyz_0[j];

            tpx_xxzz_yzz_0[j] = pa_x[j] * tpx_xzz_yzz_0[j] + 0.5 * fl1_fx * tpx_zz_yzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_yzz_0[j];

            tpy_xxzz_yzz_0[j] = pa_x[j] * tpy_xzz_yzz_0[j] + 0.5 * fl1_fx * tpy_zz_yzz_0[j];

            tpz_xxzz_yzz_0[j] = pa_x[j] * tpz_xzz_yzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzz_0[j];

            tpx_xxzz_zzz_0[j] = pa_x[j] * tpx_xzz_zzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_xzz_zzz_0[j];

            tpy_xxzz_zzz_0[j] = pa_x[j] * tpy_xzz_zzz_0[j] + 0.5 * fl1_fx * tpy_zz_zzz_0[j];

            tpz_xxzz_zzz_0[j] = pa_x[j] * tpz_xzz_zzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzz_0[j];

            tpx_xyyy_xxx_0[j] = pa_x[j] * tpx_yyy_xxx_0[j] + 1.5 * fl1_fx * tpx_yyy_xx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxx_0[j];

            tpy_xyyy_xxx_0[j] = pa_x[j] * tpy_yyy_xxx_0[j] + 1.5 * fl1_fx * tpy_yyy_xx_0[j];

            tpz_xyyy_xxx_0[j] = pa_x[j] * tpz_yyy_xxx_0[j] + 1.5 * fl1_fx * tpz_yyy_xx_0[j];

            tpx_xyyy_xxy_0[j] = pa_x[j] * tpx_yyy_xxy_0[j] + fl1_fx * tpx_yyy_xy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxy_0[j];

            tpy_xyyy_xxy_0[j] = pa_x[j] * tpy_yyy_xxy_0[j] + fl1_fx * tpy_yyy_xy_0[j];

            tpz_xyyy_xxy_0[j] = pa_x[j] * tpz_yyy_xxy_0[j] + fl1_fx * tpz_yyy_xy_0[j];

            tpx_xyyy_xxz_0[j] = pa_x[j] * tpx_yyy_xxz_0[j] + fl1_fx * tpx_yyy_xz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxz_0[j];

            tpy_xyyy_xxz_0[j] = pa_x[j] * tpy_yyy_xxz_0[j] + fl1_fx * tpy_yyy_xz_0[j];

            tpz_xyyy_xxz_0[j] = pa_x[j] * tpz_yyy_xxz_0[j] + fl1_fx * tpz_yyy_xz_0[j];

            tpx_xyyy_xyy_0[j] = pa_x[j] * tpx_yyy_xyy_0[j] + 0.5 * fl1_fx * tpx_yyy_yy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyy_0[j];

            tpy_xyyy_xyy_0[j] = pa_x[j] * tpy_yyy_xyy_0[j] + 0.5 * fl1_fx * tpy_yyy_yy_0[j];

            tpz_xyyy_xyy_0[j] = pa_x[j] * tpz_yyy_xyy_0[j] + 0.5 * fl1_fx * tpz_yyy_yy_0[j];

            tpx_xyyy_xyz_0[j] = pa_x[j] * tpx_yyy_xyz_0[j] + 0.5 * fl1_fx * tpx_yyy_yz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyz_0[j];

            tpy_xyyy_xyz_0[j] = pa_x[j] * tpy_yyy_xyz_0[j] + 0.5 * fl1_fx * tpy_yyy_yz_0[j];

            tpz_xyyy_xyz_0[j] = pa_x[j] * tpz_yyy_xyz_0[j] + 0.5 * fl1_fx * tpz_yyy_yz_0[j];

            tpx_xyyy_xzz_0[j] = pa_x[j] * tpx_yyy_xzz_0[j] + 0.5 * fl1_fx * tpx_yyy_zz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xzz_0[j];

            tpy_xyyy_xzz_0[j] = pa_x[j] * tpy_yyy_xzz_0[j] + 0.5 * fl1_fx * tpy_yyy_zz_0[j];

            tpz_xyyy_xzz_0[j] = pa_x[j] * tpz_yyy_xzz_0[j] + 0.5 * fl1_fx * tpz_yyy_zz_0[j];

            tpx_xyyy_yyy_0[j] = pa_x[j] * tpx_yyy_yyy_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyy_0[j];

            tpy_xyyy_yyy_0[j] = pa_x[j] * tpy_yyy_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_200_250(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 42);

        auto tpy_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 43);

        auto tpy_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 44);

        auto tpy_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 45);

        auto tpy_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 46);

        auto tpy_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 47);

        auto tpy_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 48);

        auto tpy_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 49);

        auto tpy_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 50);

        auto tpy_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 51);

        auto ts_yyy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 67);

        auto ts_yyy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 68);

        auto ts_yyy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 69);

        auto ts_yyz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 70);

        auto ts_yyz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 71);

        auto ts_yyz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 72);

        auto ts_yyz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 73);

        auto ts_yyz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 74);

        auto ts_yyz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 75);

        auto ts_yyz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 76);

        auto ts_yyz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 77);

        auto ts_yyz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 78);

        auto ts_yyz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 79);

        auto ts_yzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 80);

        auto ts_yzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 81);

        auto ts_yzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 82);

        auto ts_yzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 83);

        // set up pointers to integrals

        auto tpz_xyyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 66);

        auto tpx_xyyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 67);

        auto tpy_xyyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 67);

        auto tpz_xyyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 67);

        auto tpx_xyyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 68);

        auto tpy_xyyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 68);

        auto tpz_xyyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 68);

        auto tpx_xyyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 69);

        auto tpy_xyyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 69);

        auto tpz_xyyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 69);

        auto tpx_xyyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 70);

        auto tpy_xyyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 70);

        auto tpz_xyyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 70);

        auto tpx_xyyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 71);

        auto tpy_xyyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 71);

        auto tpz_xyyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 71);

        auto tpx_xyyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 72);

        auto tpy_xyyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 72);

        auto tpz_xyyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 72);

        auto tpx_xyyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 73);

        auto tpy_xyyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 73);

        auto tpz_xyyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 73);

        auto tpx_xyyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 74);

        auto tpy_xyyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 74);

        auto tpz_xyyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 74);

        auto tpx_xyyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 75);

        auto tpy_xyyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 75);

        auto tpz_xyyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 75);

        auto tpx_xyyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 76);

        auto tpy_xyyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 76);

        auto tpz_xyyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 76);

        auto tpx_xyyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 77);

        auto tpy_xyyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 77);

        auto tpz_xyyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 77);

        auto tpx_xyyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 78);

        auto tpy_xyyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 78);

        auto tpz_xyyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 78);

        auto tpx_xyyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 79);

        auto tpy_xyyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 79);

        auto tpz_xyyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 79);

        auto tpx_xyzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 80);

        auto tpy_xyzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 80);

        auto tpz_xyzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 80);

        auto tpx_xyzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 81);

        auto tpy_xyzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 81);

        auto tpz_xyzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 81);

        auto tpx_xyzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 82);

        auto tpy_xyzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 82);

        auto tpz_xyzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 82);

        auto tpx_xyzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 83);

        // Batch of Integrals (200,250)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyyy_yyz_0, tpx_xyyy_yzz_0, tpx_xyyy_zzz_0, \
                                     tpx_xyyz_xxx_0, tpx_xyyz_xxy_0, tpx_xyyz_xxz_0, tpx_xyyz_xyy_0, tpx_xyyz_xyz_0, \
                                     tpx_xyyz_xzz_0, tpx_xyyz_yyy_0, tpx_xyyz_yyz_0, tpx_xyyz_yzz_0, tpx_xyyz_zzz_0, \
                                     tpx_xyzz_xxx_0, tpx_xyzz_xxy_0, tpx_xyzz_xxz_0, tpx_xyzz_xyy_0, tpx_yyy_yyz_0, \
                                     tpx_yyy_yzz_0, tpx_yyy_zzz_0, tpx_yyz_xx_0, tpx_yyz_xxx_0, tpx_yyz_xxy_0, \
                                     tpx_yyz_xxz_0, tpx_yyz_xy_0, tpx_yyz_xyy_0, tpx_yyz_xyz_0, tpx_yyz_xz_0, \
                                     tpx_yyz_xzz_0, tpx_yyz_yy_0, tpx_yyz_yyy_0, tpx_yyz_yyz_0, tpx_yyz_yz_0, \
                                     tpx_yyz_yzz_0, tpx_yyz_zz_0, tpx_yyz_zzz_0, tpx_yzz_xx_0, tpx_yzz_xxx_0, \
                                     tpx_yzz_xxy_0, tpx_yzz_xxz_0, tpx_yzz_xy_0, tpx_yzz_xyy_0, tpx_yzz_xz_0, \
                                     tpx_yzz_yy_0, tpy_xyyy_yyz_0, tpy_xyyy_yzz_0, tpy_xyyy_zzz_0, tpy_xyyz_xxx_0, \
                                     tpy_xyyz_xxy_0, tpy_xyyz_xxz_0, tpy_xyyz_xyy_0, tpy_xyyz_xyz_0, tpy_xyyz_xzz_0, \
                                     tpy_xyyz_yyy_0, tpy_xyyz_yyz_0, tpy_xyyz_yzz_0, tpy_xyyz_zzz_0, tpy_xyzz_xxx_0, \
                                     tpy_xyzz_xxy_0, tpy_xyzz_xxz_0, tpy_yyy_yyz_0, tpy_yyy_yzz_0, tpy_yyy_zzz_0, \
                                     tpy_yyz_xx_0, tpy_yyz_xxx_0, tpy_yyz_xxy_0, tpy_yyz_xxz_0, tpy_yyz_xy_0, \
                                     tpy_yyz_xyy_0, tpy_yyz_xyz_0, tpy_yyz_xz_0, tpy_yyz_xzz_0, tpy_yyz_yy_0, \
                                     tpy_yyz_yyy_0, tpy_yyz_yyz_0, tpy_yyz_yz_0, tpy_yyz_yzz_0, tpy_yyz_zz_0, \
                                     tpy_yyz_zzz_0, tpy_yzz_xx_0, tpy_yzz_xxx_0, tpy_yzz_xxy_0, tpy_yzz_xxz_0, \
                                     tpy_yzz_xy_0, tpy_yzz_xz_0, tpz_xyyy_yyy_0, tpz_xyyy_yyz_0, tpz_xyyy_yzz_0, \
                                     tpz_xyyy_zzz_0, tpz_xyyz_xxx_0, tpz_xyyz_xxy_0, tpz_xyyz_xxz_0, tpz_xyyz_xyy_0, \
                                     tpz_xyyz_xyz_0, tpz_xyyz_xzz_0, tpz_xyyz_yyy_0, tpz_xyyz_yyz_0, tpz_xyyz_yzz_0, \
                                     tpz_xyyz_zzz_0, tpz_xyzz_xxx_0, tpz_xyzz_xxy_0, tpz_xyzz_xxz_0, tpz_yyy_yyy_0, \
                                     tpz_yyy_yyz_0, tpz_yyy_yzz_0, tpz_yyy_zzz_0, tpz_yyz_xx_0, tpz_yyz_xxx_0, \
                                     tpz_yyz_xxy_0, tpz_yyz_xxz_0, tpz_yyz_xy_0, tpz_yyz_xyy_0, tpz_yyz_xyz_0, \
                                     tpz_yyz_xz_0, tpz_yyz_xzz_0, tpz_yyz_yy_0, tpz_yyz_yyy_0, tpz_yyz_yyz_0, \
                                     tpz_yyz_yz_0, tpz_yyz_yzz_0, tpz_yyz_zz_0, tpz_yyz_zzz_0, tpz_yzz_xx_0, \
                                     tpz_yzz_xxx_0, tpz_yzz_xxy_0, tpz_yzz_xxz_0, tpz_yzz_xy_0, tpz_yzz_xz_0, \
                                     ts_yyy_yyz_0, ts_yyy_yzz_0, ts_yyy_zzz_0, ts_yyz_xxx_0, ts_yyz_xxy_0, ts_yyz_xxz_0, \
                                     ts_yyz_xyy_0, ts_yyz_xyz_0, ts_yyz_xzz_0, ts_yyz_yyy_0, ts_yyz_yyz_0, ts_yyz_yzz_0, \
                                     ts_yyz_zzz_0, ts_yzz_xxx_0, ts_yzz_xxy_0, ts_yzz_xxz_0, ts_yzz_xyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_xyyy_yyy_0[j] = pa_x[j] * tpz_yyy_yyy_0[j];

            tpx_xyyy_yyz_0[j] = pa_x[j] * tpx_yyy_yyz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyz_0[j];

            tpy_xyyy_yyz_0[j] = pa_x[j] * tpy_yyy_yyz_0[j];

            tpz_xyyy_yyz_0[j] = pa_x[j] * tpz_yyy_yyz_0[j];

            tpx_xyyy_yzz_0[j] = pa_x[j] * tpx_yyy_yzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yzz_0[j];

            tpy_xyyy_yzz_0[j] = pa_x[j] * tpy_yyy_yzz_0[j];

            tpz_xyyy_yzz_0[j] = pa_x[j] * tpz_yyy_yzz_0[j];

            tpx_xyyy_zzz_0[j] = pa_x[j] * tpx_yyy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_zzz_0[j];

            tpy_xyyy_zzz_0[j] = pa_x[j] * tpy_yyy_zzz_0[j];

            tpz_xyyy_zzz_0[j] = pa_x[j] * tpz_yyy_zzz_0[j];

            tpx_xyyz_xxx_0[j] = pa_x[j] * tpx_yyz_xxx_0[j] + 1.5 * fl1_fx * tpx_yyz_xx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxx_0[j];

            tpy_xyyz_xxx_0[j] = pa_x[j] * tpy_yyz_xxx_0[j] + 1.5 * fl1_fx * tpy_yyz_xx_0[j];

            tpz_xyyz_xxx_0[j] = pa_x[j] * tpz_yyz_xxx_0[j] + 1.5 * fl1_fx * tpz_yyz_xx_0[j];

            tpx_xyyz_xxy_0[j] = pa_x[j] * tpx_yyz_xxy_0[j] + fl1_fx * tpx_yyz_xy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxy_0[j];

            tpy_xyyz_xxy_0[j] = pa_x[j] * tpy_yyz_xxy_0[j] + fl1_fx * tpy_yyz_xy_0[j];

            tpz_xyyz_xxy_0[j] = pa_x[j] * tpz_yyz_xxy_0[j] + fl1_fx * tpz_yyz_xy_0[j];

            tpx_xyyz_xxz_0[j] = pa_x[j] * tpx_yyz_xxz_0[j] + fl1_fx * tpx_yyz_xz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxz_0[j];

            tpy_xyyz_xxz_0[j] = pa_x[j] * tpy_yyz_xxz_0[j] + fl1_fx * tpy_yyz_xz_0[j];

            tpz_xyyz_xxz_0[j] = pa_x[j] * tpz_yyz_xxz_0[j] + fl1_fx * tpz_yyz_xz_0[j];

            tpx_xyyz_xyy_0[j] = pa_x[j] * tpx_yyz_xyy_0[j] + 0.5 * fl1_fx * tpx_yyz_yy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyy_0[j];

            tpy_xyyz_xyy_0[j] = pa_x[j] * tpy_yyz_xyy_0[j] + 0.5 * fl1_fx * tpy_yyz_yy_0[j];

            tpz_xyyz_xyy_0[j] = pa_x[j] * tpz_yyz_xyy_0[j] + 0.5 * fl1_fx * tpz_yyz_yy_0[j];

            tpx_xyyz_xyz_0[j] = pa_x[j] * tpx_yyz_xyz_0[j] + 0.5 * fl1_fx * tpx_yyz_yz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyz_0[j];

            tpy_xyyz_xyz_0[j] = pa_x[j] * tpy_yyz_xyz_0[j] + 0.5 * fl1_fx * tpy_yyz_yz_0[j];

            tpz_xyyz_xyz_0[j] = pa_x[j] * tpz_yyz_xyz_0[j] + 0.5 * fl1_fx * tpz_yyz_yz_0[j];

            tpx_xyyz_xzz_0[j] = pa_x[j] * tpx_yyz_xzz_0[j] + 0.5 * fl1_fx * tpx_yyz_zz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xzz_0[j];

            tpy_xyyz_xzz_0[j] = pa_x[j] * tpy_yyz_xzz_0[j] + 0.5 * fl1_fx * tpy_yyz_zz_0[j];

            tpz_xyyz_xzz_0[j] = pa_x[j] * tpz_yyz_xzz_0[j] + 0.5 * fl1_fx * tpz_yyz_zz_0[j];

            tpx_xyyz_yyy_0[j] = pa_x[j] * tpx_yyz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyy_0[j];

            tpy_xyyz_yyy_0[j] = pa_x[j] * tpy_yyz_yyy_0[j];

            tpz_xyyz_yyy_0[j] = pa_x[j] * tpz_yyz_yyy_0[j];

            tpx_xyyz_yyz_0[j] = pa_x[j] * tpx_yyz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyz_0[j];

            tpy_xyyz_yyz_0[j] = pa_x[j] * tpy_yyz_yyz_0[j];

            tpz_xyyz_yyz_0[j] = pa_x[j] * tpz_yyz_yyz_0[j];

            tpx_xyyz_yzz_0[j] = pa_x[j] * tpx_yyz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yzz_0[j];

            tpy_xyyz_yzz_0[j] = pa_x[j] * tpy_yyz_yzz_0[j];

            tpz_xyyz_yzz_0[j] = pa_x[j] * tpz_yyz_yzz_0[j];

            tpx_xyyz_zzz_0[j] = pa_x[j] * tpx_yyz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_zzz_0[j];

            tpy_xyyz_zzz_0[j] = pa_x[j] * tpy_yyz_zzz_0[j];

            tpz_xyyz_zzz_0[j] = pa_x[j] * tpz_yyz_zzz_0[j];

            tpx_xyzz_xxx_0[j] = pa_x[j] * tpx_yzz_xxx_0[j] + 1.5 * fl1_fx * tpx_yzz_xx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxx_0[j];

            tpy_xyzz_xxx_0[j] = pa_x[j] * tpy_yzz_xxx_0[j] + 1.5 * fl1_fx * tpy_yzz_xx_0[j];

            tpz_xyzz_xxx_0[j] = pa_x[j] * tpz_yzz_xxx_0[j] + 1.5 * fl1_fx * tpz_yzz_xx_0[j];

            tpx_xyzz_xxy_0[j] = pa_x[j] * tpx_yzz_xxy_0[j] + fl1_fx * tpx_yzz_xy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxy_0[j];

            tpy_xyzz_xxy_0[j] = pa_x[j] * tpy_yzz_xxy_0[j] + fl1_fx * tpy_yzz_xy_0[j];

            tpz_xyzz_xxy_0[j] = pa_x[j] * tpz_yzz_xxy_0[j] + fl1_fx * tpz_yzz_xy_0[j];

            tpx_xyzz_xxz_0[j] = pa_x[j] * tpx_yzz_xxz_0[j] + fl1_fx * tpx_yzz_xz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxz_0[j];

            tpy_xyzz_xxz_0[j] = pa_x[j] * tpy_yzz_xxz_0[j] + fl1_fx * tpy_yzz_xz_0[j];

            tpz_xyzz_xxz_0[j] = pa_x[j] * tpz_yzz_xxz_0[j] + fl1_fx * tpz_yzz_xz_0[j];

            tpx_xyzz_xyy_0[j] = pa_x[j] * tpx_yzz_xyy_0[j] + 0.5 * fl1_fx * tpx_yzz_yy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_250_300(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpy_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 52);

        auto tpy_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 53);

        auto tpy_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54);

        auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55);

        auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 56);

        auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 57);

        auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 58);

        auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 59);

        auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59);

        auto ts_yzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 84);

        auto ts_yzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 85);

        auto ts_yzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 86);

        auto ts_yzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 87);

        auto ts_yzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 88);

        auto ts_yzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 89);

        auto ts_zzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 90);

        auto ts_zzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 91);

        auto ts_zzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 92);

        auto ts_zzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 93);

        auto ts_zzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 94);

        auto ts_zzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 95);

        auto ts_zzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 96);

        auto ts_zzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 97);

        auto ts_zzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 98);

        auto ts_zzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 99);

        // set up pointers to integrals

        auto tpy_xyzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 83);

        auto tpz_xyzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 83);

        auto tpx_xyzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 84);

        auto tpy_xyzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 84);

        auto tpz_xyzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 84);

        auto tpx_xyzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 85);

        auto tpy_xyzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 85);

        auto tpz_xyzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 85);

        auto tpx_xyzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 86);

        auto tpy_xyzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 86);

        auto tpz_xyzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 86);

        auto tpx_xyzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 87);

        auto tpy_xyzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 87);

        auto tpz_xyzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 87);

        auto tpx_xyzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 88);

        auto tpy_xyzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 88);

        auto tpz_xyzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 88);

        auto tpx_xyzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 89);

        auto tpy_xyzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 89);

        auto tpz_xyzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 89);

        auto tpx_xzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 90);

        auto tpy_xzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 90);

        auto tpz_xzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 90);

        auto tpx_xzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 91);

        auto tpy_xzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 91);

        auto tpz_xzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 91);

        auto tpx_xzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 92);

        auto tpy_xzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 92);

        auto tpz_xzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 92);

        auto tpx_xzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 93);

        auto tpy_xzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 93);

        auto tpz_xzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 93);

        auto tpx_xzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 94);

        auto tpy_xzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 94);

        auto tpz_xzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 94);

        auto tpx_xzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 95);

        auto tpy_xzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 95);

        auto tpz_xzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 95);

        auto tpx_xzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 96);

        auto tpy_xzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 96);

        auto tpz_xzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 96);

        auto tpx_xzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 97);

        auto tpy_xzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 97);

        auto tpz_xzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 97);

        auto tpx_xzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 98);

        auto tpy_xzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 98);

        auto tpz_xzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 98);

        auto tpx_xzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 99);

        auto tpy_xzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 99);

        auto tpz_xzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 99);

        // Batch of Integrals (250,300)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xyzz_xyz_0, tpx_xyzz_xzz_0, tpx_xyzz_yyy_0, \
                                     tpx_xyzz_yyz_0, tpx_xyzz_yzz_0, tpx_xyzz_zzz_0, tpx_xzzz_xxx_0, tpx_xzzz_xxy_0, \
                                     tpx_xzzz_xxz_0, tpx_xzzz_xyy_0, tpx_xzzz_xyz_0, tpx_xzzz_xzz_0, tpx_xzzz_yyy_0, \
                                     tpx_xzzz_yyz_0, tpx_xzzz_yzz_0, tpx_xzzz_zzz_0, tpx_yzz_xyz_0, tpx_yzz_xzz_0, \
                                     tpx_yzz_yyy_0, tpx_yzz_yyz_0, tpx_yzz_yz_0, tpx_yzz_yzz_0, tpx_yzz_zz_0, \
                                     tpx_yzz_zzz_0, tpx_zzz_xx_0, tpx_zzz_xxx_0, tpx_zzz_xxy_0, tpx_zzz_xxz_0, \
                                     tpx_zzz_xy_0, tpx_zzz_xyy_0, tpx_zzz_xyz_0, tpx_zzz_xz_0, tpx_zzz_xzz_0, \
                                     tpx_zzz_yy_0, tpx_zzz_yyy_0, tpx_zzz_yyz_0, tpx_zzz_yz_0, tpx_zzz_yzz_0, \
                                     tpx_zzz_zz_0, tpx_zzz_zzz_0, tpy_xyzz_xyy_0, tpy_xyzz_xyz_0, tpy_xyzz_xzz_0, \
                                     tpy_xyzz_yyy_0, tpy_xyzz_yyz_0, tpy_xyzz_yzz_0, tpy_xyzz_zzz_0, tpy_xzzz_xxx_0, \
                                     tpy_xzzz_xxy_0, tpy_xzzz_xxz_0, tpy_xzzz_xyy_0, tpy_xzzz_xyz_0, tpy_xzzz_xzz_0, \
                                     tpy_xzzz_yyy_0, tpy_xzzz_yyz_0, tpy_xzzz_yzz_0, tpy_xzzz_zzz_0, tpy_yzz_xyy_0, \
                                     tpy_yzz_xyz_0, tpy_yzz_xzz_0, tpy_yzz_yy_0, tpy_yzz_yyy_0, tpy_yzz_yyz_0, \
                                     tpy_yzz_yz_0, tpy_yzz_yzz_0, tpy_yzz_zz_0, tpy_yzz_zzz_0, tpy_zzz_xx_0, \
                                     tpy_zzz_xxx_0, tpy_zzz_xxy_0, tpy_zzz_xxz_0, tpy_zzz_xy_0, tpy_zzz_xyy_0, \
                                     tpy_zzz_xyz_0, tpy_zzz_xz_0, tpy_zzz_xzz_0, tpy_zzz_yy_0, tpy_zzz_yyy_0, \
                                     tpy_zzz_yyz_0, tpy_zzz_yz_0, tpy_zzz_yzz_0, tpy_zzz_zz_0, tpy_zzz_zzz_0, \
                                     tpz_xyzz_xyy_0, tpz_xyzz_xyz_0, tpz_xyzz_xzz_0, tpz_xyzz_yyy_0, tpz_xyzz_yyz_0, \
                                     tpz_xyzz_yzz_0, tpz_xyzz_zzz_0, tpz_xzzz_xxx_0, tpz_xzzz_xxy_0, tpz_xzzz_xxz_0, \
                                     tpz_xzzz_xyy_0, tpz_xzzz_xyz_0, tpz_xzzz_xzz_0, tpz_xzzz_yyy_0, tpz_xzzz_yyz_0, \
                                     tpz_xzzz_yzz_0, tpz_xzzz_zzz_0, tpz_yzz_xyy_0, tpz_yzz_xyz_0, tpz_yzz_xzz_0, \
                                     tpz_yzz_yy_0, tpz_yzz_yyy_0, tpz_yzz_yyz_0, tpz_yzz_yz_0, tpz_yzz_yzz_0, \
                                     tpz_yzz_zz_0, tpz_yzz_zzz_0, tpz_zzz_xx_0, tpz_zzz_xxx_0, tpz_zzz_xxy_0, \
                                     tpz_zzz_xxz_0, tpz_zzz_xy_0, tpz_zzz_xyy_0, tpz_zzz_xyz_0, tpz_zzz_xz_0, \
                                     tpz_zzz_xzz_0, tpz_zzz_yy_0, tpz_zzz_yyy_0, tpz_zzz_yyz_0, tpz_zzz_yz_0, \
                                     tpz_zzz_yzz_0, tpz_zzz_zz_0, tpz_zzz_zzz_0, ts_yzz_xyz_0, ts_yzz_xzz_0, \
                                     ts_yzz_yyy_0, ts_yzz_yyz_0, ts_yzz_yzz_0, ts_yzz_zzz_0, ts_zzz_xxx_0, ts_zzz_xxy_0, \
                                     ts_zzz_xxz_0, ts_zzz_xyy_0, ts_zzz_xyz_0, ts_zzz_xzz_0, ts_zzz_yyy_0, ts_zzz_yyz_0, \
                                     ts_zzz_yzz_0, ts_zzz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_xyzz_xyy_0[j] = pa_x[j] * tpy_yzz_xyy_0[j] + 0.5 * fl1_fx * tpy_yzz_yy_0[j];

            tpz_xyzz_xyy_0[j] = pa_x[j] * tpz_yzz_xyy_0[j] + 0.5 * fl1_fx * tpz_yzz_yy_0[j];

            tpx_xyzz_xyz_0[j] = pa_x[j] * tpx_yzz_xyz_0[j] + 0.5 * fl1_fx * tpx_yzz_yz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyz_0[j];

            tpy_xyzz_xyz_0[j] = pa_x[j] * tpy_yzz_xyz_0[j] + 0.5 * fl1_fx * tpy_yzz_yz_0[j];

            tpz_xyzz_xyz_0[j] = pa_x[j] * tpz_yzz_xyz_0[j] + 0.5 * fl1_fx * tpz_yzz_yz_0[j];

            tpx_xyzz_xzz_0[j] = pa_x[j] * tpx_yzz_xzz_0[j] + 0.5 * fl1_fx * tpx_yzz_zz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xzz_0[j];

            tpy_xyzz_xzz_0[j] = pa_x[j] * tpy_yzz_xzz_0[j] + 0.5 * fl1_fx * tpy_yzz_zz_0[j];

            tpz_xyzz_xzz_0[j] = pa_x[j] * tpz_yzz_xzz_0[j] + 0.5 * fl1_fx * tpz_yzz_zz_0[j];

            tpx_xyzz_yyy_0[j] = pa_x[j] * tpx_yzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyy_0[j];

            tpy_xyzz_yyy_0[j] = pa_x[j] * tpy_yzz_yyy_0[j];

            tpz_xyzz_yyy_0[j] = pa_x[j] * tpz_yzz_yyy_0[j];

            tpx_xyzz_yyz_0[j] = pa_x[j] * tpx_yzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyz_0[j];

            tpy_xyzz_yyz_0[j] = pa_x[j] * tpy_yzz_yyz_0[j];

            tpz_xyzz_yyz_0[j] = pa_x[j] * tpz_yzz_yyz_0[j];

            tpx_xyzz_yzz_0[j] = pa_x[j] * tpx_yzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yzz_0[j];

            tpy_xyzz_yzz_0[j] = pa_x[j] * tpy_yzz_yzz_0[j];

            tpz_xyzz_yzz_0[j] = pa_x[j] * tpz_yzz_yzz_0[j];

            tpx_xyzz_zzz_0[j] = pa_x[j] * tpx_yzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_zzz_0[j];

            tpy_xyzz_zzz_0[j] = pa_x[j] * tpy_yzz_zzz_0[j];

            tpz_xyzz_zzz_0[j] = pa_x[j] * tpz_yzz_zzz_0[j];

            tpx_xzzz_xxx_0[j] = pa_x[j] * tpx_zzz_xxx_0[j] + 1.5 * fl1_fx * tpx_zzz_xx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxx_0[j];

            tpy_xzzz_xxx_0[j] = pa_x[j] * tpy_zzz_xxx_0[j] + 1.5 * fl1_fx * tpy_zzz_xx_0[j];

            tpz_xzzz_xxx_0[j] = pa_x[j] * tpz_zzz_xxx_0[j] + 1.5 * fl1_fx * tpz_zzz_xx_0[j];

            tpx_xzzz_xxy_0[j] = pa_x[j] * tpx_zzz_xxy_0[j] + fl1_fx * tpx_zzz_xy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxy_0[j];

            tpy_xzzz_xxy_0[j] = pa_x[j] * tpy_zzz_xxy_0[j] + fl1_fx * tpy_zzz_xy_0[j];

            tpz_xzzz_xxy_0[j] = pa_x[j] * tpz_zzz_xxy_0[j] + fl1_fx * tpz_zzz_xy_0[j];

            tpx_xzzz_xxz_0[j] = pa_x[j] * tpx_zzz_xxz_0[j] + fl1_fx * tpx_zzz_xz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxz_0[j];

            tpy_xzzz_xxz_0[j] = pa_x[j] * tpy_zzz_xxz_0[j] + fl1_fx * tpy_zzz_xz_0[j];

            tpz_xzzz_xxz_0[j] = pa_x[j] * tpz_zzz_xxz_0[j] + fl1_fx * tpz_zzz_xz_0[j];

            tpx_xzzz_xyy_0[j] = pa_x[j] * tpx_zzz_xyy_0[j] + 0.5 * fl1_fx * tpx_zzz_yy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyy_0[j];

            tpy_xzzz_xyy_0[j] = pa_x[j] * tpy_zzz_xyy_0[j] + 0.5 * fl1_fx * tpy_zzz_yy_0[j];

            tpz_xzzz_xyy_0[j] = pa_x[j] * tpz_zzz_xyy_0[j] + 0.5 * fl1_fx * tpz_zzz_yy_0[j];

            tpx_xzzz_xyz_0[j] = pa_x[j] * tpx_zzz_xyz_0[j] + 0.5 * fl1_fx * tpx_zzz_yz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyz_0[j];

            tpy_xzzz_xyz_0[j] = pa_x[j] * tpy_zzz_xyz_0[j] + 0.5 * fl1_fx * tpy_zzz_yz_0[j];

            tpz_xzzz_xyz_0[j] = pa_x[j] * tpz_zzz_xyz_0[j] + 0.5 * fl1_fx * tpz_zzz_yz_0[j];

            tpx_xzzz_xzz_0[j] = pa_x[j] * tpx_zzz_xzz_0[j] + 0.5 * fl1_fx * tpx_zzz_zz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xzz_0[j];

            tpy_xzzz_xzz_0[j] = pa_x[j] * tpy_zzz_xzz_0[j] + 0.5 * fl1_fx * tpy_zzz_zz_0[j];

            tpz_xzzz_xzz_0[j] = pa_x[j] * tpz_zzz_xzz_0[j] + 0.5 * fl1_fx * tpz_zzz_zz_0[j];

            tpx_xzzz_yyy_0[j] = pa_x[j] * tpx_zzz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyy_0[j];

            tpy_xzzz_yyy_0[j] = pa_x[j] * tpy_zzz_yyy_0[j];

            tpz_xzzz_yyy_0[j] = pa_x[j] * tpz_zzz_yyy_0[j];

            tpx_xzzz_yyz_0[j] = pa_x[j] * tpx_zzz_yyz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyz_0[j];

            tpy_xzzz_yyz_0[j] = pa_x[j] * tpy_zzz_yyz_0[j];

            tpz_xzzz_yyz_0[j] = pa_x[j] * tpz_zzz_yyz_0[j];

            tpx_xzzz_yzz_0[j] = pa_x[j] * tpx_zzz_yzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yzz_0[j];

            tpy_xzzz_yzz_0[j] = pa_x[j] * tpy_zzz_yzz_0[j];

            tpz_xzzz_yzz_0[j] = pa_x[j] * tpz_zzz_yzz_0[j];

            tpx_xzzz_zzz_0[j] = pa_x[j] * tpx_zzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zzz_0[j];

            tpy_xzzz_zzz_0[j] = pa_x[j] * tpy_zzz_zzz_0[j];

            tpz_xzzz_zzz_0[j] = pa_x[j] * tpz_zzz_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_300_350(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 36);

        auto tpy_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 36);

        auto tpz_yyy_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 36);

        auto tpx_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 37);

        auto tpy_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 37);

        auto tpz_yyy_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 37);

        auto tpx_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 38);

        auto tpy_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 38);

        auto tpz_yyy_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 38);

        auto tpx_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 39);

        auto tpy_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 39);

        auto tpz_yyy_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 39);

        auto tpx_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 40);

        auto tpy_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 40);

        auto tpz_yyy_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 40);

        auto tpx_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 41);

        auto tpy_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 41);

        auto tpz_yyy_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 41);

        auto tpx_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 42);

        auto tpy_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 42);

        auto tpz_yyz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 42);

        auto tpx_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 43);

        auto tpy_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 43);

        auto tpz_yyz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 43);

        auto tpx_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 44);

        auto tpy_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 44);

        auto tpz_yyz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 44);

        auto tpx_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 45);

        auto tpy_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 45);

        auto ts_yyy_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 60);

        auto ts_yyy_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 61);

        auto ts_yyy_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 62);

        auto ts_yyy_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 63);

        auto ts_yyy_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 64);

        auto ts_yyy_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 65);

        auto ts_yyy_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 66);

        auto ts_yyy_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 67);

        auto ts_yyy_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 68);

        auto ts_yyy_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 69);

        auto ts_yyz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 70);

        auto ts_yyz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 71);

        auto ts_yyz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 72);

        auto ts_yyz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 73);

        auto ts_yyz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 74);

        auto ts_yyz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 75);

        auto ts_yyz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 76);

        // set up pointers to integrals

        auto tpx_yyyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 100);

        auto tpy_yyyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 100);

        auto tpz_yyyy_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 100);

        auto tpx_yyyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 101);

        auto tpy_yyyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 101);

        auto tpz_yyyy_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 101);

        auto tpx_yyyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 102);

        auto tpy_yyyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 102);

        auto tpz_yyyy_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 102);

        auto tpx_yyyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 103);

        auto tpy_yyyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 103);

        auto tpz_yyyy_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 103);

        auto tpx_yyyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 104);

        auto tpy_yyyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 104);

        auto tpz_yyyy_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 104);

        auto tpx_yyyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 105);

        auto tpy_yyyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 105);

        auto tpz_yyyy_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 105);

        auto tpx_yyyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 106);

        auto tpy_yyyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 106);

        auto tpz_yyyy_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 106);

        auto tpx_yyyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 107);

        auto tpy_yyyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 107);

        auto tpz_yyyy_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 107);

        auto tpx_yyyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 108);

        auto tpy_yyyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 108);

        auto tpz_yyyy_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 108);

        auto tpx_yyyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 109);

        auto tpy_yyyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 109);

        auto tpz_yyyy_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 109);

        auto tpx_yyyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 110);

        auto tpy_yyyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 110);

        auto tpz_yyyz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 110);

        auto tpx_yyyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 111);

        auto tpy_yyyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 111);

        auto tpz_yyyz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 111);

        auto tpx_yyyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 112);

        auto tpy_yyyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 112);

        auto tpz_yyyz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 112);

        auto tpx_yyyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 113);

        auto tpy_yyyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 113);

        auto tpz_yyyz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 113);

        auto tpx_yyyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 114);

        auto tpy_yyyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 114);

        auto tpz_yyyz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 114);

        auto tpx_yyyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 115);

        auto tpy_yyyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 115);

        auto tpz_yyyz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 115);

        auto tpx_yyyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 116);

        auto tpy_yyyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 116);

        // Batch of Integrals (300,350)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yy_xxx_0, tpx_yy_xxy_0, tpx_yy_xxz_0, tpx_yy_xyy_0, \
                                     tpx_yy_xyz_0, tpx_yy_xzz_0, tpx_yy_yyy_0, tpx_yy_yyz_0, tpx_yy_yzz_0, tpx_yy_zzz_0, \
                                     tpx_yyy_xx_0, tpx_yyy_xxx_0, tpx_yyy_xxy_0, tpx_yyy_xxz_0, tpx_yyy_xy_0, \
                                     tpx_yyy_xyy_0, tpx_yyy_xyz_0, tpx_yyy_xz_0, tpx_yyy_xzz_0, tpx_yyy_yy_0, \
                                     tpx_yyy_yyy_0, tpx_yyy_yyz_0, tpx_yyy_yz_0, tpx_yyy_yzz_0, tpx_yyy_zz_0, \
                                     tpx_yyy_zzz_0, tpx_yyyy_xxx_0, tpx_yyyy_xxy_0, tpx_yyyy_xxz_0, tpx_yyyy_xyy_0, \
                                     tpx_yyyy_xyz_0, tpx_yyyy_xzz_0, tpx_yyyy_yyy_0, tpx_yyyy_yyz_0, tpx_yyyy_yzz_0, \
                                     tpx_yyyy_zzz_0, tpx_yyyz_xxx_0, tpx_yyyz_xxy_0, tpx_yyyz_xxz_0, tpx_yyyz_xyy_0, \
                                     tpx_yyyz_xyz_0, tpx_yyyz_xzz_0, tpx_yyyz_yyy_0, tpx_yyz_xx_0, tpx_yyz_xxx_0, \
                                     tpx_yyz_xxy_0, tpx_yyz_xxz_0, tpx_yyz_xy_0, tpx_yyz_xyy_0, tpx_yyz_xyz_0, \
                                     tpx_yyz_xz_0, tpx_yyz_xzz_0, tpx_yyz_yy_0, tpx_yyz_yyy_0, tpx_yz_xxx_0, \
                                     tpx_yz_xxy_0, tpx_yz_xxz_0, tpx_yz_xyy_0, tpx_yz_xyz_0, tpx_yz_xzz_0, tpx_yz_yyy_0, \
                                     tpy_yy_xxx_0, tpy_yy_xxy_0, tpy_yy_xxz_0, tpy_yy_xyy_0, tpy_yy_xyz_0, tpy_yy_xzz_0, \
                                     tpy_yy_yyy_0, tpy_yy_yyz_0, tpy_yy_yzz_0, tpy_yy_zzz_0, tpy_yyy_xx_0, \
                                     tpy_yyy_xxx_0, tpy_yyy_xxy_0, tpy_yyy_xxz_0, tpy_yyy_xy_0, tpy_yyy_xyy_0, \
                                     tpy_yyy_xyz_0, tpy_yyy_xz_0, tpy_yyy_xzz_0, tpy_yyy_yy_0, tpy_yyy_yyy_0, \
                                     tpy_yyy_yyz_0, tpy_yyy_yz_0, tpy_yyy_yzz_0, tpy_yyy_zz_0, tpy_yyy_zzz_0, \
                                     tpy_yyyy_xxx_0, tpy_yyyy_xxy_0, tpy_yyyy_xxz_0, tpy_yyyy_xyy_0, tpy_yyyy_xyz_0, \
                                     tpy_yyyy_xzz_0, tpy_yyyy_yyy_0, tpy_yyyy_yyz_0, tpy_yyyy_yzz_0, tpy_yyyy_zzz_0, \
                                     tpy_yyyz_xxx_0, tpy_yyyz_xxy_0, tpy_yyyz_xxz_0, tpy_yyyz_xyy_0, tpy_yyyz_xyz_0, \
                                     tpy_yyyz_xzz_0, tpy_yyyz_yyy_0, tpy_yyz_xx_0, tpy_yyz_xxx_0, tpy_yyz_xxy_0, \
                                     tpy_yyz_xxz_0, tpy_yyz_xy_0, tpy_yyz_xyy_0, tpy_yyz_xyz_0, tpy_yyz_xz_0, \
                                     tpy_yyz_xzz_0, tpy_yyz_yy_0, tpy_yyz_yyy_0, tpy_yz_xxx_0, tpy_yz_xxy_0, \
                                     tpy_yz_xxz_0, tpy_yz_xyy_0, tpy_yz_xyz_0, tpy_yz_xzz_0, tpy_yz_yyy_0, tpz_yy_xxx_0, \
                                     tpz_yy_xxy_0, tpz_yy_xxz_0, tpz_yy_xyy_0, tpz_yy_xyz_0, tpz_yy_xzz_0, tpz_yy_yyy_0, \
                                     tpz_yy_yyz_0, tpz_yy_yzz_0, tpz_yy_zzz_0, tpz_yyy_xx_0, tpz_yyy_xxx_0, \
                                     tpz_yyy_xxy_0, tpz_yyy_xxz_0, tpz_yyy_xy_0, tpz_yyy_xyy_0, tpz_yyy_xyz_0, \
                                     tpz_yyy_xz_0, tpz_yyy_xzz_0, tpz_yyy_yy_0, tpz_yyy_yyy_0, tpz_yyy_yyz_0, \
                                     tpz_yyy_yz_0, tpz_yyy_yzz_0, tpz_yyy_zz_0, tpz_yyy_zzz_0, tpz_yyyy_xxx_0, \
                                     tpz_yyyy_xxy_0, tpz_yyyy_xxz_0, tpz_yyyy_xyy_0, tpz_yyyy_xyz_0, tpz_yyyy_xzz_0, \
                                     tpz_yyyy_yyy_0, tpz_yyyy_yyz_0, tpz_yyyy_yzz_0, tpz_yyyy_zzz_0, tpz_yyyz_xxx_0, \
                                     tpz_yyyz_xxy_0, tpz_yyyz_xxz_0, tpz_yyyz_xyy_0, tpz_yyyz_xyz_0, tpz_yyyz_xzz_0, \
                                     tpz_yyz_xx_0, tpz_yyz_xxx_0, tpz_yyz_xxy_0, tpz_yyz_xxz_0, tpz_yyz_xy_0, \
                                     tpz_yyz_xyy_0, tpz_yyz_xyz_0, tpz_yyz_xz_0, tpz_yyz_xzz_0, tpz_yz_xxx_0, \
                                     tpz_yz_xxy_0, tpz_yz_xxz_0, tpz_yz_xyy_0, tpz_yz_xyz_0, tpz_yz_xzz_0, ts_yyy_xxx_0, \
                                     ts_yyy_xxy_0, ts_yyy_xxz_0, ts_yyy_xyy_0, ts_yyy_xyz_0, ts_yyy_xzz_0, ts_yyy_yyy_0, \
                                     ts_yyy_yyz_0, ts_yyy_yzz_0, ts_yyy_zzz_0, ts_yyz_xxx_0, ts_yyz_xxy_0, ts_yyz_xxz_0, \
                                     ts_yyz_xyy_0, ts_yyz_xyz_0, ts_yyz_xzz_0, ts_yyz_yyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyyy_xxx_0[j] = pa_y[j] * tpx_yyy_xxx_0[j] + 1.5 * fl1_fx * tpx_yy_xxx_0[j];

            tpy_yyyy_xxx_0[j] = pa_y[j] * tpy_yyy_xxx_0[j] + 1.5 * fl1_fx * tpy_yy_xxx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxx_0[j];

            tpz_yyyy_xxx_0[j] = pa_y[j] * tpz_yyy_xxx_0[j] + 1.5 * fl1_fx * tpz_yy_xxx_0[j];

            tpx_yyyy_xxy_0[j] = pa_y[j] * tpx_yyy_xxy_0[j] + 1.5 * fl1_fx * tpx_yy_xxy_0[j] + 0.5 * fl1_fx * tpx_yyy_xx_0[j];

            tpy_yyyy_xxy_0[j] =
                pa_y[j] * tpy_yyy_xxy_0[j] + 1.5 * fl1_fx * tpy_yy_xxy_0[j] + 0.5 * fl1_fx * tpy_yyy_xx_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxy_0[j];

            tpz_yyyy_xxy_0[j] = pa_y[j] * tpz_yyy_xxy_0[j] + 1.5 * fl1_fx * tpz_yy_xxy_0[j] + 0.5 * fl1_fx * tpz_yyy_xx_0[j];

            tpx_yyyy_xxz_0[j] = pa_y[j] * tpx_yyy_xxz_0[j] + 1.5 * fl1_fx * tpx_yy_xxz_0[j];

            tpy_yyyy_xxz_0[j] = pa_y[j] * tpy_yyy_xxz_0[j] + 1.5 * fl1_fx * tpy_yy_xxz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xxz_0[j];

            tpz_yyyy_xxz_0[j] = pa_y[j] * tpz_yyy_xxz_0[j] + 1.5 * fl1_fx * tpz_yy_xxz_0[j];

            tpx_yyyy_xyy_0[j] = pa_y[j] * tpx_yyy_xyy_0[j] + 1.5 * fl1_fx * tpx_yy_xyy_0[j] + fl1_fx * tpx_yyy_xy_0[j];

            tpy_yyyy_xyy_0[j] =
                pa_y[j] * tpy_yyy_xyy_0[j] + 1.5 * fl1_fx * tpy_yy_xyy_0[j] + fl1_fx * tpy_yyy_xy_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyy_0[j];

            tpz_yyyy_xyy_0[j] = pa_y[j] * tpz_yyy_xyy_0[j] + 1.5 * fl1_fx * tpz_yy_xyy_0[j] + fl1_fx * tpz_yyy_xy_0[j];

            tpx_yyyy_xyz_0[j] = pa_y[j] * tpx_yyy_xyz_0[j] + 1.5 * fl1_fx * tpx_yy_xyz_0[j] + 0.5 * fl1_fx * tpx_yyy_xz_0[j];

            tpy_yyyy_xyz_0[j] =
                pa_y[j] * tpy_yyy_xyz_0[j] + 1.5 * fl1_fx * tpy_yy_xyz_0[j] + 0.5 * fl1_fx * tpy_yyy_xz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xyz_0[j];

            tpz_yyyy_xyz_0[j] = pa_y[j] * tpz_yyy_xyz_0[j] + 1.5 * fl1_fx * tpz_yy_xyz_0[j] + 0.5 * fl1_fx * tpz_yyy_xz_0[j];

            tpx_yyyy_xzz_0[j] = pa_y[j] * tpx_yyy_xzz_0[j] + 1.5 * fl1_fx * tpx_yy_xzz_0[j];

            tpy_yyyy_xzz_0[j] = pa_y[j] * tpy_yyy_xzz_0[j] + 1.5 * fl1_fx * tpy_yy_xzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_xzz_0[j];

            tpz_yyyy_xzz_0[j] = pa_y[j] * tpz_yyy_xzz_0[j] + 1.5 * fl1_fx * tpz_yy_xzz_0[j];

            tpx_yyyy_yyy_0[j] = pa_y[j] * tpx_yyy_yyy_0[j] + 1.5 * fl1_fx * tpx_yy_yyy_0[j] + 1.5 * fl1_fx * tpx_yyy_yy_0[j];

            tpy_yyyy_yyy_0[j] =
                pa_y[j] * tpy_yyy_yyy_0[j] + 1.5 * fl1_fx * tpy_yy_yyy_0[j] + 1.5 * fl1_fx * tpy_yyy_yy_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyy_0[j];

            tpz_yyyy_yyy_0[j] = pa_y[j] * tpz_yyy_yyy_0[j] + 1.5 * fl1_fx * tpz_yy_yyy_0[j] + 1.5 * fl1_fx * tpz_yyy_yy_0[j];

            tpx_yyyy_yyz_0[j] = pa_y[j] * tpx_yyy_yyz_0[j] + 1.5 * fl1_fx * tpx_yy_yyz_0[j] + fl1_fx * tpx_yyy_yz_0[j];

            tpy_yyyy_yyz_0[j] =
                pa_y[j] * tpy_yyy_yyz_0[j] + 1.5 * fl1_fx * tpy_yy_yyz_0[j] + fl1_fx * tpy_yyy_yz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yyz_0[j];

            tpz_yyyy_yyz_0[j] = pa_y[j] * tpz_yyy_yyz_0[j] + 1.5 * fl1_fx * tpz_yy_yyz_0[j] + fl1_fx * tpz_yyy_yz_0[j];

            tpx_yyyy_yzz_0[j] = pa_y[j] * tpx_yyy_yzz_0[j] + 1.5 * fl1_fx * tpx_yy_yzz_0[j] + 0.5 * fl1_fx * tpx_yyy_zz_0[j];

            tpy_yyyy_yzz_0[j] =
                pa_y[j] * tpy_yyy_yzz_0[j] + 1.5 * fl1_fx * tpy_yy_yzz_0[j] + 0.5 * fl1_fx * tpy_yyy_zz_0[j] - fl1_fgb * fl1_fx * ts_yyy_yzz_0[j];

            tpz_yyyy_yzz_0[j] = pa_y[j] * tpz_yyy_yzz_0[j] + 1.5 * fl1_fx * tpz_yy_yzz_0[j] + 0.5 * fl1_fx * tpz_yyy_zz_0[j];

            tpx_yyyy_zzz_0[j] = pa_y[j] * tpx_yyy_zzz_0[j] + 1.5 * fl1_fx * tpx_yy_zzz_0[j];

            tpy_yyyy_zzz_0[j] = pa_y[j] * tpy_yyy_zzz_0[j] + 1.5 * fl1_fx * tpy_yy_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyy_zzz_0[j];

            tpz_yyyy_zzz_0[j] = pa_y[j] * tpz_yyy_zzz_0[j] + 1.5 * fl1_fx * tpz_yy_zzz_0[j];

            tpx_yyyz_xxx_0[j] = pa_y[j] * tpx_yyz_xxx_0[j] + fl1_fx * tpx_yz_xxx_0[j];

            tpy_yyyz_xxx_0[j] = pa_y[j] * tpy_yyz_xxx_0[j] + fl1_fx * tpy_yz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxx_0[j];

            tpz_yyyz_xxx_0[j] = pa_y[j] * tpz_yyz_xxx_0[j] + fl1_fx * tpz_yz_xxx_0[j];

            tpx_yyyz_xxy_0[j] = pa_y[j] * tpx_yyz_xxy_0[j] + fl1_fx * tpx_yz_xxy_0[j] + 0.5 * fl1_fx * tpx_yyz_xx_0[j];

            tpy_yyyz_xxy_0[j] =
                pa_y[j] * tpy_yyz_xxy_0[j] + fl1_fx * tpy_yz_xxy_0[j] + 0.5 * fl1_fx * tpy_yyz_xx_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxy_0[j];

            tpz_yyyz_xxy_0[j] = pa_y[j] * tpz_yyz_xxy_0[j] + fl1_fx * tpz_yz_xxy_0[j] + 0.5 * fl1_fx * tpz_yyz_xx_0[j];

            tpx_yyyz_xxz_0[j] = pa_y[j] * tpx_yyz_xxz_0[j] + fl1_fx * tpx_yz_xxz_0[j];

            tpy_yyyz_xxz_0[j] = pa_y[j] * tpy_yyz_xxz_0[j] + fl1_fx * tpy_yz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xxz_0[j];

            tpz_yyyz_xxz_0[j] = pa_y[j] * tpz_yyz_xxz_0[j] + fl1_fx * tpz_yz_xxz_0[j];

            tpx_yyyz_xyy_0[j] = pa_y[j] * tpx_yyz_xyy_0[j] + fl1_fx * tpx_yz_xyy_0[j] + fl1_fx * tpx_yyz_xy_0[j];

            tpy_yyyz_xyy_0[j] = pa_y[j] * tpy_yyz_xyy_0[j] + fl1_fx * tpy_yz_xyy_0[j] + fl1_fx * tpy_yyz_xy_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyy_0[j];

            tpz_yyyz_xyy_0[j] = pa_y[j] * tpz_yyz_xyy_0[j] + fl1_fx * tpz_yz_xyy_0[j] + fl1_fx * tpz_yyz_xy_0[j];

            tpx_yyyz_xyz_0[j] = pa_y[j] * tpx_yyz_xyz_0[j] + fl1_fx * tpx_yz_xyz_0[j] + 0.5 * fl1_fx * tpx_yyz_xz_0[j];

            tpy_yyyz_xyz_0[j] =
                pa_y[j] * tpy_yyz_xyz_0[j] + fl1_fx * tpy_yz_xyz_0[j] + 0.5 * fl1_fx * tpy_yyz_xz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xyz_0[j];

            tpz_yyyz_xyz_0[j] = pa_y[j] * tpz_yyz_xyz_0[j] + fl1_fx * tpz_yz_xyz_0[j] + 0.5 * fl1_fx * tpz_yyz_xz_0[j];

            tpx_yyyz_xzz_0[j] = pa_y[j] * tpx_yyz_xzz_0[j] + fl1_fx * tpx_yz_xzz_0[j];

            tpy_yyyz_xzz_0[j] = pa_y[j] * tpy_yyz_xzz_0[j] + fl1_fx * tpy_yz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_xzz_0[j];

            tpz_yyyz_xzz_0[j] = pa_y[j] * tpz_yyz_xzz_0[j] + fl1_fx * tpz_yz_xzz_0[j];

            tpx_yyyz_yyy_0[j] = pa_y[j] * tpx_yyz_yyy_0[j] + fl1_fx * tpx_yz_yyy_0[j] + 1.5 * fl1_fx * tpx_yyz_yy_0[j];

            tpy_yyyz_yyy_0[j] =
                pa_y[j] * tpy_yyz_yyy_0[j] + fl1_fx * tpy_yz_yyy_0[j] + 1.5 * fl1_fx * tpy_yyz_yy_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_350_400(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpy_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zz_yzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * idx + 59);

        auto tpy_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zz_zzz_0 = primBuffer.data(pidx_p_2_3_m0 + 120 * bdim + 60 * idx + 59);

        auto tpz_yyz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 45);

        auto tpx_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 46);

        auto tpy_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 46);

        auto tpz_yyz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 46);

        auto tpx_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 47);

        auto tpy_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 47);

        auto tpz_yyz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 47);

        auto tpx_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 48);

        auto tpy_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 48);

        auto tpz_yzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 48);

        auto tpx_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 49);

        auto tpy_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 49);

        auto tpz_yzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 49);

        auto tpx_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 50);

        auto tpy_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 50);

        auto tpz_yzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 50);

        auto tpx_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 51);

        auto tpy_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 51);

        auto tpz_yzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 51);

        auto tpx_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 52);

        auto tpy_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 52);

        auto tpz_yzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 52);

        auto tpx_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 53);

        auto tpy_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 53);

        auto tpz_yzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 53);

        auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54);

        auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55);

        auto ts_yyz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 77);

        auto ts_yyz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 78);

        auto ts_yyz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 79);

        auto ts_yzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 80);

        auto ts_yzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 81);

        auto ts_yzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 82);

        auto ts_yzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 83);

        auto ts_yzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 84);

        auto ts_yzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 85);

        auto ts_yzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 86);

        auto ts_yzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 87);

        auto ts_yzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 88);

        auto ts_yzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 89);

        auto ts_zzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 90);

        auto ts_zzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 91);

        auto ts_zzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 92);

        // set up pointers to integrals

        auto tpz_yyyz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 116);

        auto tpx_yyyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 117);

        auto tpy_yyyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 117);

        auto tpz_yyyz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 117);

        auto tpx_yyyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 118);

        auto tpy_yyyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 118);

        auto tpz_yyyz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 118);

        auto tpx_yyyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 119);

        auto tpy_yyyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 119);

        auto tpz_yyyz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 119);

        auto tpx_yyzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 120);

        auto tpy_yyzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 120);

        auto tpz_yyzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 120);

        auto tpx_yyzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 121);

        auto tpy_yyzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 121);

        auto tpz_yyzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 121);

        auto tpx_yyzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 122);

        auto tpy_yyzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 122);

        auto tpz_yyzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 122);

        auto tpx_yyzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 123);

        auto tpy_yyzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 123);

        auto tpz_yyzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 123);

        auto tpx_yyzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 124);

        auto tpy_yyzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 124);

        auto tpz_yyzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 124);

        auto tpx_yyzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 125);

        auto tpy_yyzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 125);

        auto tpz_yyzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 125);

        auto tpx_yyzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 126);

        auto tpy_yyzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 126);

        auto tpz_yyzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 126);

        auto tpx_yyzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 127);

        auto tpy_yyzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 127);

        auto tpz_yyzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 127);

        auto tpx_yyzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 128);

        auto tpy_yyzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 128);

        auto tpz_yyzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 128);

        auto tpx_yyzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 129);

        auto tpy_yyzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 129);

        auto tpz_yyzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 129);

        auto tpx_yzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 130);

        auto tpy_yzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 130);

        auto tpz_yzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 130);

        auto tpx_yzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 131);

        auto tpy_yzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 131);

        auto tpz_yzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 131);

        auto tpx_yzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 132);

        auto tpy_yzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 132);

        auto tpz_yzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 132);

        auto tpx_yzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 133);

        // Batch of Integrals (350,400)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_yyyz_yyz_0, tpx_yyyz_yzz_0, tpx_yyyz_zzz_0, \
                                     tpx_yyz_yyz_0, tpx_yyz_yz_0, tpx_yyz_yzz_0, tpx_yyz_zz_0, tpx_yyz_zzz_0, \
                                     tpx_yyzz_xxx_0, tpx_yyzz_xxy_0, tpx_yyzz_xxz_0, tpx_yyzz_xyy_0, tpx_yyzz_xyz_0, \
                                     tpx_yyzz_xzz_0, tpx_yyzz_yyy_0, tpx_yyzz_yyz_0, tpx_yyzz_yzz_0, tpx_yyzz_zzz_0, \
                                     tpx_yz_yyz_0, tpx_yz_yzz_0, tpx_yz_zzz_0, tpx_yzz_xx_0, tpx_yzz_xxx_0, \
                                     tpx_yzz_xxy_0, tpx_yzz_xxz_0, tpx_yzz_xy_0, tpx_yzz_xyy_0, tpx_yzz_xyz_0, \
                                     tpx_yzz_xz_0, tpx_yzz_xzz_0, tpx_yzz_yy_0, tpx_yzz_yyy_0, tpx_yzz_yyz_0, \
                                     tpx_yzz_yz_0, tpx_yzz_yzz_0, tpx_yzz_zz_0, tpx_yzz_zzz_0, tpx_yzzz_xxx_0, \
                                     tpx_yzzz_xxy_0, tpx_yzzz_xxz_0, tpx_yzzz_xyy_0, tpx_zz_xxx_0, tpx_zz_xxy_0, \
                                     tpx_zz_xxz_0, tpx_zz_xyy_0, tpx_zz_xyz_0, tpx_zz_xzz_0, tpx_zz_yyy_0, tpx_zz_yyz_0, \
                                     tpx_zz_yzz_0, tpx_zz_zzz_0, tpx_zzz_xx_0, tpx_zzz_xxx_0, tpx_zzz_xxy_0, \
                                     tpx_zzz_xxz_0, tpx_zzz_xy_0, tpx_zzz_xyy_0, tpy_yyyz_yyz_0, tpy_yyyz_yzz_0, \
                                     tpy_yyyz_zzz_0, tpy_yyz_yyz_0, tpy_yyz_yz_0, tpy_yyz_yzz_0, tpy_yyz_zz_0, \
                                     tpy_yyz_zzz_0, tpy_yyzz_xxx_0, tpy_yyzz_xxy_0, tpy_yyzz_xxz_0, tpy_yyzz_xyy_0, \
                                     tpy_yyzz_xyz_0, tpy_yyzz_xzz_0, tpy_yyzz_yyy_0, tpy_yyzz_yyz_0, tpy_yyzz_yzz_0, \
                                     tpy_yyzz_zzz_0, tpy_yz_yyz_0, tpy_yz_yzz_0, tpy_yz_zzz_0, tpy_yzz_xx_0, \
                                     tpy_yzz_xxx_0, tpy_yzz_xxy_0, tpy_yzz_xxz_0, tpy_yzz_xy_0, tpy_yzz_xyy_0, \
                                     tpy_yzz_xyz_0, tpy_yzz_xz_0, tpy_yzz_xzz_0, tpy_yzz_yy_0, tpy_yzz_yyy_0, \
                                     tpy_yzz_yyz_0, tpy_yzz_yz_0, tpy_yzz_yzz_0, tpy_yzz_zz_0, tpy_yzz_zzz_0, \
                                     tpy_yzzz_xxx_0, tpy_yzzz_xxy_0, tpy_yzzz_xxz_0, tpy_zz_xxx_0, tpy_zz_xxy_0, \
                                     tpy_zz_xxz_0, tpy_zz_xyy_0, tpy_zz_xyz_0, tpy_zz_xzz_0, tpy_zz_yyy_0, tpy_zz_yyz_0, \
                                     tpy_zz_yzz_0, tpy_zz_zzz_0, tpy_zzz_xx_0, tpy_zzz_xxx_0, tpy_zzz_xxy_0, \
                                     tpy_zzz_xxz_0, tpz_yyyz_yyy_0, tpz_yyyz_yyz_0, tpz_yyyz_yzz_0, tpz_yyyz_zzz_0, \
                                     tpz_yyz_yy_0, tpz_yyz_yyy_0, tpz_yyz_yyz_0, tpz_yyz_yz_0, tpz_yyz_yzz_0, \
                                     tpz_yyz_zz_0, tpz_yyz_zzz_0, tpz_yyzz_xxx_0, tpz_yyzz_xxy_0, tpz_yyzz_xxz_0, \
                                     tpz_yyzz_xyy_0, tpz_yyzz_xyz_0, tpz_yyzz_xzz_0, tpz_yyzz_yyy_0, tpz_yyzz_yyz_0, \
                                     tpz_yyzz_yzz_0, tpz_yyzz_zzz_0, tpz_yz_yyy_0, tpz_yz_yyz_0, tpz_yz_yzz_0, \
                                     tpz_yz_zzz_0, tpz_yzz_xx_0, tpz_yzz_xxx_0, tpz_yzz_xxy_0, tpz_yzz_xxz_0, \
                                     tpz_yzz_xy_0, tpz_yzz_xyy_0, tpz_yzz_xyz_0, tpz_yzz_xz_0, tpz_yzz_xzz_0, \
                                     tpz_yzz_yy_0, tpz_yzz_yyy_0, tpz_yzz_yyz_0, tpz_yzz_yz_0, tpz_yzz_yzz_0, \
                                     tpz_yzz_zz_0, tpz_yzz_zzz_0, tpz_yzzz_xxx_0, tpz_yzzz_xxy_0, tpz_yzzz_xxz_0, \
                                     tpz_zz_xxx_0, tpz_zz_xxy_0, tpz_zz_xxz_0, tpz_zz_xyy_0, tpz_zz_xyz_0, tpz_zz_xzz_0, \
                                     tpz_zz_yyy_0, tpz_zz_yyz_0, tpz_zz_yzz_0, tpz_zz_zzz_0, tpz_zzz_xx_0, \
                                     tpz_zzz_xxx_0, tpz_zzz_xxy_0, tpz_zzz_xxz_0, ts_yyz_yyz_0, ts_yyz_yzz_0, \
                                     ts_yyz_zzz_0, ts_yzz_xxx_0, ts_yzz_xxy_0, ts_yzz_xxz_0, ts_yzz_xyy_0, ts_yzz_xyz_0, \
                                     ts_yzz_xzz_0, ts_yzz_yyy_0, ts_yzz_yyz_0, ts_yzz_yzz_0, ts_yzz_zzz_0, ts_zzz_xxx_0, \
                                     ts_zzz_xxy_0, ts_zzz_xxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpz_yyyz_yyy_0[j] = pa_y[j] * tpz_yyz_yyy_0[j] + fl1_fx * tpz_yz_yyy_0[j] + 1.5 * fl1_fx * tpz_yyz_yy_0[j];

            tpx_yyyz_yyz_0[j] = pa_y[j] * tpx_yyz_yyz_0[j] + fl1_fx * tpx_yz_yyz_0[j] + fl1_fx * tpx_yyz_yz_0[j];

            tpy_yyyz_yyz_0[j] = pa_y[j] * tpy_yyz_yyz_0[j] + fl1_fx * tpy_yz_yyz_0[j] + fl1_fx * tpy_yyz_yz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yyz_0[j];

            tpz_yyyz_yyz_0[j] = pa_y[j] * tpz_yyz_yyz_0[j] + fl1_fx * tpz_yz_yyz_0[j] + fl1_fx * tpz_yyz_yz_0[j];

            tpx_yyyz_yzz_0[j] = pa_y[j] * tpx_yyz_yzz_0[j] + fl1_fx * tpx_yz_yzz_0[j] + 0.5 * fl1_fx * tpx_yyz_zz_0[j];

            tpy_yyyz_yzz_0[j] =
                pa_y[j] * tpy_yyz_yzz_0[j] + fl1_fx * tpy_yz_yzz_0[j] + 0.5 * fl1_fx * tpy_yyz_zz_0[j] - fl1_fgb * fl1_fx * ts_yyz_yzz_0[j];

            tpz_yyyz_yzz_0[j] = pa_y[j] * tpz_yyz_yzz_0[j] + fl1_fx * tpz_yz_yzz_0[j] + 0.5 * fl1_fx * tpz_yyz_zz_0[j];

            tpx_yyyz_zzz_0[j] = pa_y[j] * tpx_yyz_zzz_0[j] + fl1_fx * tpx_yz_zzz_0[j];

            tpy_yyyz_zzz_0[j] = pa_y[j] * tpy_yyz_zzz_0[j] + fl1_fx * tpy_yz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yyz_zzz_0[j];

            tpz_yyyz_zzz_0[j] = pa_y[j] * tpz_yyz_zzz_0[j] + fl1_fx * tpz_yz_zzz_0[j];

            tpx_yyzz_xxx_0[j] = pa_y[j] * tpx_yzz_xxx_0[j] + 0.5 * fl1_fx * tpx_zz_xxx_0[j];

            tpy_yyzz_xxx_0[j] = pa_y[j] * tpy_yzz_xxx_0[j] + 0.5 * fl1_fx * tpy_zz_xxx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxx_0[j];

            tpz_yyzz_xxx_0[j] = pa_y[j] * tpz_yzz_xxx_0[j] + 0.5 * fl1_fx * tpz_zz_xxx_0[j];

            tpx_yyzz_xxy_0[j] = pa_y[j] * tpx_yzz_xxy_0[j] + 0.5 * fl1_fx * tpx_zz_xxy_0[j] + 0.5 * fl1_fx * tpx_yzz_xx_0[j];

            tpy_yyzz_xxy_0[j] =
                pa_y[j] * tpy_yzz_xxy_0[j] + 0.5 * fl1_fx * tpy_zz_xxy_0[j] + 0.5 * fl1_fx * tpy_yzz_xx_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxy_0[j];

            tpz_yyzz_xxy_0[j] = pa_y[j] * tpz_yzz_xxy_0[j] + 0.5 * fl1_fx * tpz_zz_xxy_0[j] + 0.5 * fl1_fx * tpz_yzz_xx_0[j];

            tpx_yyzz_xxz_0[j] = pa_y[j] * tpx_yzz_xxz_0[j] + 0.5 * fl1_fx * tpx_zz_xxz_0[j];

            tpy_yyzz_xxz_0[j] = pa_y[j] * tpy_yzz_xxz_0[j] + 0.5 * fl1_fx * tpy_zz_xxz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xxz_0[j];

            tpz_yyzz_xxz_0[j] = pa_y[j] * tpz_yzz_xxz_0[j] + 0.5 * fl1_fx * tpz_zz_xxz_0[j];

            tpx_yyzz_xyy_0[j] = pa_y[j] * tpx_yzz_xyy_0[j] + 0.5 * fl1_fx * tpx_zz_xyy_0[j] + fl1_fx * tpx_yzz_xy_0[j];

            tpy_yyzz_xyy_0[j] =
                pa_y[j] * tpy_yzz_xyy_0[j] + 0.5 * fl1_fx * tpy_zz_xyy_0[j] + fl1_fx * tpy_yzz_xy_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyy_0[j];

            tpz_yyzz_xyy_0[j] = pa_y[j] * tpz_yzz_xyy_0[j] + 0.5 * fl1_fx * tpz_zz_xyy_0[j] + fl1_fx * tpz_yzz_xy_0[j];

            tpx_yyzz_xyz_0[j] = pa_y[j] * tpx_yzz_xyz_0[j] + 0.5 * fl1_fx * tpx_zz_xyz_0[j] + 0.5 * fl1_fx * tpx_yzz_xz_0[j];

            tpy_yyzz_xyz_0[j] =
                pa_y[j] * tpy_yzz_xyz_0[j] + 0.5 * fl1_fx * tpy_zz_xyz_0[j] + 0.5 * fl1_fx * tpy_yzz_xz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xyz_0[j];

            tpz_yyzz_xyz_0[j] = pa_y[j] * tpz_yzz_xyz_0[j] + 0.5 * fl1_fx * tpz_zz_xyz_0[j] + 0.5 * fl1_fx * tpz_yzz_xz_0[j];

            tpx_yyzz_xzz_0[j] = pa_y[j] * tpx_yzz_xzz_0[j] + 0.5 * fl1_fx * tpx_zz_xzz_0[j];

            tpy_yyzz_xzz_0[j] = pa_y[j] * tpy_yzz_xzz_0[j] + 0.5 * fl1_fx * tpy_zz_xzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_xzz_0[j];

            tpz_yyzz_xzz_0[j] = pa_y[j] * tpz_yzz_xzz_0[j] + 0.5 * fl1_fx * tpz_zz_xzz_0[j];

            tpx_yyzz_yyy_0[j] = pa_y[j] * tpx_yzz_yyy_0[j] + 0.5 * fl1_fx * tpx_zz_yyy_0[j] + 1.5 * fl1_fx * tpx_yzz_yy_0[j];

            tpy_yyzz_yyy_0[j] =
                pa_y[j] * tpy_yzz_yyy_0[j] + 0.5 * fl1_fx * tpy_zz_yyy_0[j] + 1.5 * fl1_fx * tpy_yzz_yy_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyy_0[j];

            tpz_yyzz_yyy_0[j] = pa_y[j] * tpz_yzz_yyy_0[j] + 0.5 * fl1_fx * tpz_zz_yyy_0[j] + 1.5 * fl1_fx * tpz_yzz_yy_0[j];

            tpx_yyzz_yyz_0[j] = pa_y[j] * tpx_yzz_yyz_0[j] + 0.5 * fl1_fx * tpx_zz_yyz_0[j] + fl1_fx * tpx_yzz_yz_0[j];

            tpy_yyzz_yyz_0[j] =
                pa_y[j] * tpy_yzz_yyz_0[j] + 0.5 * fl1_fx * tpy_zz_yyz_0[j] + fl1_fx * tpy_yzz_yz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yyz_0[j];

            tpz_yyzz_yyz_0[j] = pa_y[j] * tpz_yzz_yyz_0[j] + 0.5 * fl1_fx * tpz_zz_yyz_0[j] + fl1_fx * tpz_yzz_yz_0[j];

            tpx_yyzz_yzz_0[j] = pa_y[j] * tpx_yzz_yzz_0[j] + 0.5 * fl1_fx * tpx_zz_yzz_0[j] + 0.5 * fl1_fx * tpx_yzz_zz_0[j];

            tpy_yyzz_yzz_0[j] =
                pa_y[j] * tpy_yzz_yzz_0[j] + 0.5 * fl1_fx * tpy_zz_yzz_0[j] + 0.5 * fl1_fx * tpy_yzz_zz_0[j] - fl1_fgb * fl1_fx * ts_yzz_yzz_0[j];

            tpz_yyzz_yzz_0[j] = pa_y[j] * tpz_yzz_yzz_0[j] + 0.5 * fl1_fx * tpz_zz_yzz_0[j] + 0.5 * fl1_fx * tpz_yzz_zz_0[j];

            tpx_yyzz_zzz_0[j] = pa_y[j] * tpx_yzz_zzz_0[j] + 0.5 * fl1_fx * tpx_zz_zzz_0[j];

            tpy_yyzz_zzz_0[j] = pa_y[j] * tpy_yzz_zzz_0[j] + 0.5 * fl1_fx * tpy_zz_zzz_0[j] - fl1_fgb * fl1_fx * ts_yzz_zzz_0[j];

            tpz_yyzz_zzz_0[j] = pa_y[j] * tpz_yzz_zzz_0[j] + 0.5 * fl1_fx * tpz_zz_zzz_0[j];

            tpx_yzzz_xxx_0[j] = pa_y[j] * tpx_zzz_xxx_0[j];

            tpy_yzzz_xxx_0[j] = pa_y[j] * tpy_zzz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxx_0[j];

            tpz_yzzz_xxx_0[j] = pa_y[j] * tpz_zzz_xxx_0[j];

            tpx_yzzz_xxy_0[j] = pa_y[j] * tpx_zzz_xxy_0[j] + 0.5 * fl1_fx * tpx_zzz_xx_0[j];

            tpy_yzzz_xxy_0[j] = pa_y[j] * tpy_zzz_xxy_0[j] + 0.5 * fl1_fx * tpy_zzz_xx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxy_0[j];

            tpz_yzzz_xxy_0[j] = pa_y[j] * tpz_zzz_xxy_0[j] + 0.5 * fl1_fx * tpz_zzz_xx_0[j];

            tpx_yzzz_xxz_0[j] = pa_y[j] * tpx_zzz_xxz_0[j];

            tpy_yzzz_xxz_0[j] = pa_y[j] * tpy_zzz_xxz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxz_0[j];

            tpz_yzzz_xxz_0[j] = pa_y[j] * tpz_zzz_xxz_0[j];

            tpx_yzzz_xyy_0[j] = pa_y[j] * tpx_zzz_xyy_0[j] + fl1_fx * tpx_zzz_xy_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGF_400_450(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_4_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 54);

        auto tpy_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 54);

        auto tpz_zzz_xx_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 54);

        auto tpx_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 55);

        auto tpy_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 55);

        auto tpz_zzz_xy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 55);

        auto tpx_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 56);

        auto tpy_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 56);

        auto tpz_zzz_xz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 56);

        auto tpx_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 57);

        auto tpy_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 57);

        auto tpz_zzz_yy_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 57);

        auto tpx_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 58);

        auto tpy_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 58);

        auto tpz_zzz_yz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 58);

        auto tpx_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * idx + 59);

        auto tpy_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 60 * bdim + 60 * idx + 59);

        auto tpz_zzz_zz_0 = primBuffer.data(pidx_p_3_2_m0 + 120 * bdim + 60 * idx + 59);

        auto ts_zzz_xxx_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 90);

        auto ts_zzz_xxy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 91);

        auto ts_zzz_xxz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 92);

        auto ts_zzz_xyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 93);

        auto ts_zzz_xyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 94);

        auto ts_zzz_xzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 95);

        auto ts_zzz_yyy_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 96);

        auto ts_zzz_yyz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 97);

        auto ts_zzz_yzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 98);

        auto ts_zzz_zzz_0 = primBuffer.data(pidx_s_3_3_m0 + 100 * idx + 99);

        // set up pointers to integrals

        auto tpy_yzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 133);

        auto tpz_yzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 133);

        auto tpx_yzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 134);

        auto tpy_yzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 134);

        auto tpz_yzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 134);

        auto tpx_yzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 135);

        auto tpy_yzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 135);

        auto tpz_yzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 135);

        auto tpx_yzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 136);

        auto tpy_yzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 136);

        auto tpz_yzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 136);

        auto tpx_yzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 137);

        auto tpy_yzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 137);

        auto tpz_yzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 137);

        auto tpx_yzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 138);

        auto tpy_yzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 138);

        auto tpz_yzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 138);

        auto tpx_yzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 139);

        auto tpy_yzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 139);

        auto tpz_yzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 139);

        auto tpx_zzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 140);

        auto tpy_zzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 140);

        auto tpz_zzzz_xxx_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 140);

        auto tpx_zzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 141);

        auto tpy_zzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 141);

        auto tpz_zzzz_xxy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 141);

        auto tpx_zzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 142);

        auto tpy_zzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 142);

        auto tpz_zzzz_xxz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 142);

        auto tpx_zzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 143);

        auto tpy_zzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 143);

        auto tpz_zzzz_xyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 143);

        auto tpx_zzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 144);

        auto tpy_zzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 144);

        auto tpz_zzzz_xyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 144);

        auto tpx_zzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 145);

        auto tpy_zzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 145);

        auto tpz_zzzz_xzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 145);

        auto tpx_zzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 146);

        auto tpy_zzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 146);

        auto tpz_zzzz_yyy_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 146);

        auto tpx_zzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 147);

        auto tpy_zzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 147);

        auto tpz_zzzz_yyz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 147);

        auto tpx_zzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 148);

        auto tpy_zzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 148);

        auto tpz_zzzz_yzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 148);

        auto tpx_zzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * idx + 149);

        auto tpy_zzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 150 * bdim + 150 * idx + 149);

        auto tpz_zzzz_zzz_0 = primBuffer.data(pidx_p_4_3_m0 + 300 * bdim + 150 * idx + 149);

        // Batch of Integrals (400,450)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yzzz_xyz_0, tpx_yzzz_xzz_0, tpx_yzzz_yyy_0, \
                                     tpx_yzzz_yyz_0, tpx_yzzz_yzz_0, tpx_yzzz_zzz_0, tpx_zz_xxx_0, tpx_zz_xxy_0, \
                                     tpx_zz_xxz_0, tpx_zz_xyy_0, tpx_zz_xyz_0, tpx_zz_xzz_0, tpx_zz_yyy_0, tpx_zz_yyz_0, \
                                     tpx_zz_yzz_0, tpx_zz_zzz_0, tpx_zzz_xx_0, tpx_zzz_xxx_0, tpx_zzz_xxy_0, \
                                     tpx_zzz_xxz_0, tpx_zzz_xy_0, tpx_zzz_xyy_0, tpx_zzz_xyz_0, tpx_zzz_xz_0, \
                                     tpx_zzz_xzz_0, tpx_zzz_yy_0, tpx_zzz_yyy_0, tpx_zzz_yyz_0, tpx_zzz_yz_0, \
                                     tpx_zzz_yzz_0, tpx_zzz_zz_0, tpx_zzz_zzz_0, tpx_zzzz_xxx_0, tpx_zzzz_xxy_0, \
                                     tpx_zzzz_xxz_0, tpx_zzzz_xyy_0, tpx_zzzz_xyz_0, tpx_zzzz_xzz_0, tpx_zzzz_yyy_0, \
                                     tpx_zzzz_yyz_0, tpx_zzzz_yzz_0, tpx_zzzz_zzz_0, tpy_yzzz_xyy_0, tpy_yzzz_xyz_0, \
                                     tpy_yzzz_xzz_0, tpy_yzzz_yyy_0, tpy_yzzz_yyz_0, tpy_yzzz_yzz_0, tpy_yzzz_zzz_0, \
                                     tpy_zz_xxx_0, tpy_zz_xxy_0, tpy_zz_xxz_0, tpy_zz_xyy_0, tpy_zz_xyz_0, tpy_zz_xzz_0, \
                                     tpy_zz_yyy_0, tpy_zz_yyz_0, tpy_zz_yzz_0, tpy_zz_zzz_0, tpy_zzz_xx_0, \
                                     tpy_zzz_xxx_0, tpy_zzz_xxy_0, tpy_zzz_xxz_0, tpy_zzz_xy_0, tpy_zzz_xyy_0, \
                                     tpy_zzz_xyz_0, tpy_zzz_xz_0, tpy_zzz_xzz_0, tpy_zzz_yy_0, tpy_zzz_yyy_0, \
                                     tpy_zzz_yyz_0, tpy_zzz_yz_0, tpy_zzz_yzz_0, tpy_zzz_zz_0, tpy_zzz_zzz_0, \
                                     tpy_zzzz_xxx_0, tpy_zzzz_xxy_0, tpy_zzzz_xxz_0, tpy_zzzz_xyy_0, tpy_zzzz_xyz_0, \
                                     tpy_zzzz_xzz_0, tpy_zzzz_yyy_0, tpy_zzzz_yyz_0, tpy_zzzz_yzz_0, tpy_zzzz_zzz_0, \
                                     tpz_yzzz_xyy_0, tpz_yzzz_xyz_0, tpz_yzzz_xzz_0, tpz_yzzz_yyy_0, tpz_yzzz_yyz_0, \
                                     tpz_yzzz_yzz_0, tpz_yzzz_zzz_0, tpz_zz_xxx_0, tpz_zz_xxy_0, tpz_zz_xxz_0, \
                                     tpz_zz_xyy_0, tpz_zz_xyz_0, tpz_zz_xzz_0, tpz_zz_yyy_0, tpz_zz_yyz_0, tpz_zz_yzz_0, \
                                     tpz_zz_zzz_0, tpz_zzz_xx_0, tpz_zzz_xxx_0, tpz_zzz_xxy_0, tpz_zzz_xxz_0, \
                                     tpz_zzz_xy_0, tpz_zzz_xyy_0, tpz_zzz_xyz_0, tpz_zzz_xz_0, tpz_zzz_xzz_0, \
                                     tpz_zzz_yy_0, tpz_zzz_yyy_0, tpz_zzz_yyz_0, tpz_zzz_yz_0, tpz_zzz_yzz_0, \
                                     tpz_zzz_zz_0, tpz_zzz_zzz_0, tpz_zzzz_xxx_0, tpz_zzzz_xxy_0, tpz_zzzz_xxz_0, \
                                     tpz_zzzz_xyy_0, tpz_zzzz_xyz_0, tpz_zzzz_xzz_0, tpz_zzzz_yyy_0, tpz_zzzz_yyz_0, \
                                     tpz_zzzz_yzz_0, tpz_zzzz_zzz_0, ts_zzz_xxx_0, ts_zzz_xxy_0, ts_zzz_xxz_0, \
                                     ts_zzz_xyy_0, ts_zzz_xyz_0, ts_zzz_xzz_0, ts_zzz_yyy_0, ts_zzz_yyz_0, ts_zzz_yzz_0, \
                                     ts_zzz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpy_yzzz_xyy_0[j] = pa_y[j] * tpy_zzz_xyy_0[j] + fl1_fx * tpy_zzz_xy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyy_0[j];

            tpz_yzzz_xyy_0[j] = pa_y[j] * tpz_zzz_xyy_0[j] + fl1_fx * tpz_zzz_xy_0[j];

            tpx_yzzz_xyz_0[j] = pa_y[j] * tpx_zzz_xyz_0[j] + 0.5 * fl1_fx * tpx_zzz_xz_0[j];

            tpy_yzzz_xyz_0[j] = pa_y[j] * tpy_zzz_xyz_0[j] + 0.5 * fl1_fx * tpy_zzz_xz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyz_0[j];

            tpz_yzzz_xyz_0[j] = pa_y[j] * tpz_zzz_xyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xz_0[j];

            tpx_yzzz_xzz_0[j] = pa_y[j] * tpx_zzz_xzz_0[j];

            tpy_yzzz_xzz_0[j] = pa_y[j] * tpy_zzz_xzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xzz_0[j];

            tpz_yzzz_xzz_0[j] = pa_y[j] * tpz_zzz_xzz_0[j];

            tpx_yzzz_yyy_0[j] = pa_y[j] * tpx_zzz_yyy_0[j] + 1.5 * fl1_fx * tpx_zzz_yy_0[j];

            tpy_yzzz_yyy_0[j] = pa_y[j] * tpy_zzz_yyy_0[j] + 1.5 * fl1_fx * tpy_zzz_yy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyy_0[j];

            tpz_yzzz_yyy_0[j] = pa_y[j] * tpz_zzz_yyy_0[j] + 1.5 * fl1_fx * tpz_zzz_yy_0[j];

            tpx_yzzz_yyz_0[j] = pa_y[j] * tpx_zzz_yyz_0[j] + fl1_fx * tpx_zzz_yz_0[j];

            tpy_yzzz_yyz_0[j] = pa_y[j] * tpy_zzz_yyz_0[j] + fl1_fx * tpy_zzz_yz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyz_0[j];

            tpz_yzzz_yyz_0[j] = pa_y[j] * tpz_zzz_yyz_0[j] + fl1_fx * tpz_zzz_yz_0[j];

            tpx_yzzz_yzz_0[j] = pa_y[j] * tpx_zzz_yzz_0[j] + 0.5 * fl1_fx * tpx_zzz_zz_0[j];

            tpy_yzzz_yzz_0[j] = pa_y[j] * tpy_zzz_yzz_0[j] + 0.5 * fl1_fx * tpy_zzz_zz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yzz_0[j];

            tpz_yzzz_yzz_0[j] = pa_y[j] * tpz_zzz_yzz_0[j] + 0.5 * fl1_fx * tpz_zzz_zz_0[j];

            tpx_yzzz_zzz_0[j] = pa_y[j] * tpx_zzz_zzz_0[j];

            tpy_yzzz_zzz_0[j] = pa_y[j] * tpy_zzz_zzz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zzz_0[j];

            tpz_yzzz_zzz_0[j] = pa_y[j] * tpz_zzz_zzz_0[j];

            tpx_zzzz_xxx_0[j] = pa_z[j] * tpx_zzz_xxx_0[j] + 1.5 * fl1_fx * tpx_zz_xxx_0[j];

            tpy_zzzz_xxx_0[j] = pa_z[j] * tpy_zzz_xxx_0[j] + 1.5 * fl1_fx * tpy_zz_xxx_0[j];

            tpz_zzzz_xxx_0[j] = pa_z[j] * tpz_zzz_xxx_0[j] + 1.5 * fl1_fx * tpz_zz_xxx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxx_0[j];

            tpx_zzzz_xxy_0[j] = pa_z[j] * tpx_zzz_xxy_0[j] + 1.5 * fl1_fx * tpx_zz_xxy_0[j];

            tpy_zzzz_xxy_0[j] = pa_z[j] * tpy_zzz_xxy_0[j] + 1.5 * fl1_fx * tpy_zz_xxy_0[j];

            tpz_zzzz_xxy_0[j] = pa_z[j] * tpz_zzz_xxy_0[j] + 1.5 * fl1_fx * tpz_zz_xxy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxy_0[j];

            tpx_zzzz_xxz_0[j] = pa_z[j] * tpx_zzz_xxz_0[j] + 1.5 * fl1_fx * tpx_zz_xxz_0[j] + 0.5 * fl1_fx * tpx_zzz_xx_0[j];

            tpy_zzzz_xxz_0[j] = pa_z[j] * tpy_zzz_xxz_0[j] + 1.5 * fl1_fx * tpy_zz_xxz_0[j] + 0.5 * fl1_fx * tpy_zzz_xx_0[j];

            tpz_zzzz_xxz_0[j] =
                pa_z[j] * tpz_zzz_xxz_0[j] + 1.5 * fl1_fx * tpz_zz_xxz_0[j] + 0.5 * fl1_fx * tpz_zzz_xx_0[j] - fl1_fgb * fl1_fx * ts_zzz_xxz_0[j];

            tpx_zzzz_xyy_0[j] = pa_z[j] * tpx_zzz_xyy_0[j] + 1.5 * fl1_fx * tpx_zz_xyy_0[j];

            tpy_zzzz_xyy_0[j] = pa_z[j] * tpy_zzz_xyy_0[j] + 1.5 * fl1_fx * tpy_zz_xyy_0[j];

            tpz_zzzz_xyy_0[j] = pa_z[j] * tpz_zzz_xyy_0[j] + 1.5 * fl1_fx * tpz_zz_xyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyy_0[j];

            tpx_zzzz_xyz_0[j] = pa_z[j] * tpx_zzz_xyz_0[j] + 1.5 * fl1_fx * tpx_zz_xyz_0[j] + 0.5 * fl1_fx * tpx_zzz_xy_0[j];

            tpy_zzzz_xyz_0[j] = pa_z[j] * tpy_zzz_xyz_0[j] + 1.5 * fl1_fx * tpy_zz_xyz_0[j] + 0.5 * fl1_fx * tpy_zzz_xy_0[j];

            tpz_zzzz_xyz_0[j] =
                pa_z[j] * tpz_zzz_xyz_0[j] + 1.5 * fl1_fx * tpz_zz_xyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xy_0[j] - fl1_fgb * fl1_fx * ts_zzz_xyz_0[j];

            tpx_zzzz_xzz_0[j] = pa_z[j] * tpx_zzz_xzz_0[j] + 1.5 * fl1_fx * tpx_zz_xzz_0[j] + fl1_fx * tpx_zzz_xz_0[j];

            tpy_zzzz_xzz_0[j] = pa_z[j] * tpy_zzz_xzz_0[j] + 1.5 * fl1_fx * tpy_zz_xzz_0[j] + fl1_fx * tpy_zzz_xz_0[j];

            tpz_zzzz_xzz_0[j] =
                pa_z[j] * tpz_zzz_xzz_0[j] + 1.5 * fl1_fx * tpz_zz_xzz_0[j] + fl1_fx * tpz_zzz_xz_0[j] - fl1_fgb * fl1_fx * ts_zzz_xzz_0[j];

            tpx_zzzz_yyy_0[j] = pa_z[j] * tpx_zzz_yyy_0[j] + 1.5 * fl1_fx * tpx_zz_yyy_0[j];

            tpy_zzzz_yyy_0[j] = pa_z[j] * tpy_zzz_yyy_0[j] + 1.5 * fl1_fx * tpy_zz_yyy_0[j];

            tpz_zzzz_yyy_0[j] = pa_z[j] * tpz_zzz_yyy_0[j] + 1.5 * fl1_fx * tpz_zz_yyy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyy_0[j];

            tpx_zzzz_yyz_0[j] = pa_z[j] * tpx_zzz_yyz_0[j] + 1.5 * fl1_fx * tpx_zz_yyz_0[j] + 0.5 * fl1_fx * tpx_zzz_yy_0[j];

            tpy_zzzz_yyz_0[j] = pa_z[j] * tpy_zzz_yyz_0[j] + 1.5 * fl1_fx * tpy_zz_yyz_0[j] + 0.5 * fl1_fx * tpy_zzz_yy_0[j];

            tpz_zzzz_yyz_0[j] =
                pa_z[j] * tpz_zzz_yyz_0[j] + 1.5 * fl1_fx * tpz_zz_yyz_0[j] + 0.5 * fl1_fx * tpz_zzz_yy_0[j] - fl1_fgb * fl1_fx * ts_zzz_yyz_0[j];

            tpx_zzzz_yzz_0[j] = pa_z[j] * tpx_zzz_yzz_0[j] + 1.5 * fl1_fx * tpx_zz_yzz_0[j] + fl1_fx * tpx_zzz_yz_0[j];

            tpy_zzzz_yzz_0[j] = pa_z[j] * tpy_zzz_yzz_0[j] + 1.5 * fl1_fx * tpy_zz_yzz_0[j] + fl1_fx * tpy_zzz_yz_0[j];

            tpz_zzzz_yzz_0[j] =
                pa_z[j] * tpz_zzz_yzz_0[j] + 1.5 * fl1_fx * tpz_zz_yzz_0[j] + fl1_fx * tpz_zzz_yz_0[j] - fl1_fgb * fl1_fx * ts_zzz_yzz_0[j];

            tpx_zzzz_zzz_0[j] = pa_z[j] * tpx_zzz_zzz_0[j] + 1.5 * fl1_fx * tpx_zz_zzz_0[j] + 1.5 * fl1_fx * tpx_zzz_zz_0[j];

            tpy_zzzz_zzz_0[j] = pa_z[j] * tpy_zzz_zzz_0[j] + 1.5 * fl1_fx * tpy_zz_zzz_0[j] + 1.5 * fl1_fx * tpy_zzz_zz_0[j];

            tpz_zzzz_zzz_0[j] =
                pa_z[j] * tpz_zzz_zzz_0[j] + 1.5 * fl1_fx * tpz_zz_zzz_0[j] + 1.5 * fl1_fx * tpz_zzz_zz_0[j] - fl1_fgb * fl1_fx * ts_zzz_zzz_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
