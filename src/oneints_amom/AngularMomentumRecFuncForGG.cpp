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

#include "AngularMomentumRecFuncForGG.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

void
compAngularMomentumForGG(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
                         const CGtoBlock&           braGtoBlock,
                         const CGtoBlock&           ketGtoBlock,
                         const int32_t              iContrGto)
{
    amomrecfunc::compAngularMomentumForGG_0_49(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_49_98(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_98_147(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_147_195(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_195_243(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_243_291(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_291_339(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_339_387(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_387_435(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_435_483(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_483_531(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_531_579(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_579_627(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    amomrecfunc::compAngularMomentumForGG_627_675(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compAngularMomentumForGG_0_49(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xxx_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx);

        auto tly_xxx_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx);

        auto tlz_xxx_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx);

        auto tlx_xxx_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 1);

        auto tly_xxx_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 1);

        auto tlz_xxx_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 1);

        auto tlx_xxx_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 2);

        auto tly_xxx_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 2);

        auto tlz_xxx_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 2);

        auto tlx_xxx_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 3);

        auto tly_xxx_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 3);

        auto tlz_xxx_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 3);

        auto tlx_xxx_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 4);

        auto tly_xxx_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 4);

        auto tlz_xxx_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 4);

        auto tlx_xxx_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 5);

        auto tly_xxx_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 5);

        auto tlz_xxx_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 5);

        auto tlx_xxx_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 6);

        auto tly_xxx_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 6);

        auto tlz_xxx_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 6);

        auto tlx_xxx_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 7);

        auto tly_xxx_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 7);

        auto tlz_xxx_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 7);

        auto tlx_xxx_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 8);

        auto tly_xxx_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 8);

        auto tlz_xxx_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 8);

        auto tlx_xxx_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 9);

        auto tly_xxx_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 9);

        auto tlz_xxx_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 9);

        auto tlx_xxx_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 10);

        auto tly_xxx_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 10);

        auto tlz_xxx_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 10);

        auto tlx_xxx_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 11);

        auto tly_xxx_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 11);

        auto tlz_xxx_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 11);

        auto tlx_xxx_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 12);

        auto tly_xxx_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 12);

        auto tlz_xxx_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 12);

        auto tlx_xxx_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 13);

        auto tly_xxx_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 13);

        auto tlz_xxx_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 13);

        auto tlx_xxx_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 14);

        auto tly_xxx_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 14);

        auto tlz_xxx_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 14);

        auto tlx_xxy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 15);

        auto tly_xxy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 15);

        auto tlz_xxy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 15);

        auto tlx_xxy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 16);

        auto tlx_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx);

        auto tly_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx);

        auto tlz_xx_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx);

        auto tlx_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 1);

        auto tly_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 1);

        auto tlz_xx_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 1);

        auto tlx_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 2);

        auto tly_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 2);

        auto tlz_xx_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 2);

        auto tlx_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 3);

        auto tly_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 3);

        auto tlz_xx_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 3);

        auto tlx_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 4);

        auto tly_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 4);

        auto tlz_xx_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 4);

        auto tlx_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 5);

        auto tly_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 5);

        auto tlz_xx_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 5);

        auto tlx_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 6);

        auto tly_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 6);

        auto tlz_xx_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 6);

        auto tlx_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 7);

        auto tly_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 7);

        auto tlz_xx_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 7);

        auto tlx_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 8);

        auto tly_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 8);

        auto tlz_xx_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 8);

        auto tlx_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 9);

        auto tly_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 9);

        auto tlz_xx_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 9);

        auto tlx_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 10);

        auto tly_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 10);

        auto tlz_xx_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 10);

        auto tlx_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 11);

        auto tly_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 11);

        auto tlz_xx_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 11);

        auto tlx_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 12);

        auto tly_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 12);

        auto tlz_xx_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 12);

        auto tlx_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 13);

        auto tly_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 13);

        auto tlz_xx_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 13);

        auto tlx_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 14);

        auto tly_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 14);

        auto tlz_xx_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 14);

        auto tlx_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 15);

        auto tly_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 15);

        auto tlz_xy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 15);

        auto tlx_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 16);

        auto tlx_xxx_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx);

        auto tly_xxx_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx);

        auto tlz_xxx_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx);

        auto tlx_xxx_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 1);

        auto tly_xxx_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 1);

        auto tlz_xxx_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 1);

        auto tlx_xxx_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 2);

        auto tly_xxx_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 2);

        auto tlz_xxx_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 2);

        auto tlx_xxx_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 3);

        auto tly_xxx_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 3);

        auto tlz_xxx_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 3);

        auto tlx_xxx_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 4);

        auto tly_xxx_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 4);

        auto tlz_xxx_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 4);

        auto tlx_xxx_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 5);

        auto tly_xxx_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 5);

        auto tlz_xxx_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 5);

        auto tlx_xxx_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 6);

        auto tly_xxx_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 6);

        auto tlz_xxx_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 6);

        auto tlx_xxx_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 7);

        auto tly_xxx_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 7);

        auto tlz_xxx_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 7);

        auto tlx_xxx_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 8);

        auto tly_xxx_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 8);

        auto tlz_xxx_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 8);

        auto tlx_xxx_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 9);

        auto tly_xxx_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 9);

        auto tlz_xxx_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 9);

        auto tlx_xxy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 10);

        auto tly_xxy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 10);

        auto tlz_xxy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 10);

        auto tlx_xxy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 11);

        auto tpy_xxx_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx);

        auto tpz_xxx_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx);

        auto tpy_xxx_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 1);

        auto tpz_xxx_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 1);

        auto tpy_xxx_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 2);

        auto tpz_xxx_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 2);

        auto tpy_xxx_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 3);

        auto tpz_xxx_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 3);

        auto tpy_xxx_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 4);

        auto tpz_xxx_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 4);

        auto tpy_xxx_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 5);

        auto tpz_xxx_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 5);

        auto tpy_xxx_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 6);

        auto tpz_xxx_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 6);

        auto tpy_xxx_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 7);

        auto tpz_xxx_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 7);

        auto tpy_xxx_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 8);

        auto tpz_xxx_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 8);

        auto tpy_xxx_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 9);

        auto tpz_xxx_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 9);

        auto tpy_xxx_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 10);

        auto tpz_xxx_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 10);

        auto tpy_xxx_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 11);

        auto tpz_xxx_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 11);

        auto tpy_xxx_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 12);

        auto tpz_xxx_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 12);

        auto tpy_xxx_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 13);

        auto tpz_xxx_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 13);

        auto tpy_xxx_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 14);

        auto tpz_xxx_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 14);

        auto tpy_xxy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 15);

        auto tpz_xxy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 15);

        auto tdy_xxx_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx);

        auto tdz_xxx_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx);

        auto tdy_xxx_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 1);

        auto tdz_xxx_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 1);

        auto tdy_xxx_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 2);

        auto tdz_xxx_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 2);

        auto tdy_xxx_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 3);

        auto tdz_xxx_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 3);

        auto tdy_xxx_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 4);

        auto tdz_xxx_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 4);

        auto tdy_xxx_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 5);

        auto tdz_xxx_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 5);

        auto tdy_xxx_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 6);

        auto tdz_xxx_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 6);

        auto tdy_xxx_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 7);

        auto tdz_xxx_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 7);

        auto tdy_xxx_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 8);

        auto tdz_xxx_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 8);

        auto tdy_xxx_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 9);

        auto tdz_xxx_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 9);

        auto tdy_xxx_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 10);

        auto tdz_xxx_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 10);

        auto tdy_xxx_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 11);

        auto tdz_xxx_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 11);

        auto tdy_xxx_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 12);

        auto tdz_xxx_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 12);

        auto tdy_xxx_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 13);

        auto tdz_xxx_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 13);

        auto tdy_xxx_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 14);

        auto tdz_xxx_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 14);

        auto tdy_xxy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 15);

        auto tdz_xxy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 15);

        // set up pointers to integrals

        auto tlx_xxxx_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx);

        auto tly_xxxx_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx);

        auto tlz_xxxx_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx);

        auto tlx_xxxx_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 1);

        auto tly_xxxx_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 1);

        auto tlz_xxxx_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 1);

        auto tlx_xxxx_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 2);

        auto tly_xxxx_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 2);

        auto tlz_xxxx_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 2);

        auto tlx_xxxx_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 3);

        auto tly_xxxx_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 3);

        auto tlz_xxxx_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 3);

        auto tlx_xxxx_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 4);

        auto tly_xxxx_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 4);

        auto tlz_xxxx_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 4);

        auto tlx_xxxx_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 5);

        auto tly_xxxx_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 5);

        auto tlz_xxxx_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 5);

        auto tlx_xxxx_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 6);

        auto tly_xxxx_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 6);

        auto tlz_xxxx_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 6);

        auto tlx_xxxx_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 7);

        auto tly_xxxx_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 7);

        auto tlz_xxxx_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 7);

        auto tlx_xxxx_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 8);

        auto tly_xxxx_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 8);

        auto tlz_xxxx_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 8);

        auto tlx_xxxx_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 9);

        auto tly_xxxx_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 9);

        auto tlz_xxxx_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 9);

        auto tlx_xxxx_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 10);

        auto tly_xxxx_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 10);

        auto tlz_xxxx_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 10);

        auto tlx_xxxx_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 11);

        auto tly_xxxx_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 11);

        auto tlz_xxxx_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 11);

        auto tlx_xxxx_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 12);

        auto tly_xxxx_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 12);

        auto tlz_xxxx_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 12);

        auto tlx_xxxx_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 13);

        auto tly_xxxx_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 13);

        auto tlz_xxxx_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 13);

        auto tlx_xxxx_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 14);

        auto tly_xxxx_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 14);

        auto tlz_xxxx_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 14);

        auto tlx_xxxy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 15);

        auto tly_xxxy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 15);

        auto tlz_xxxy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 15);

        auto tlx_xxxy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 16);

        // Batch of Integrals (0,49)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxx_xxxx_0, tdy_xxx_xxxy_0, tdy_xxx_xxxz_0, \
                                     tdy_xxx_xxyy_0, tdy_xxx_xxyz_0, tdy_xxx_xxzz_0, tdy_xxx_xyyy_0, tdy_xxx_xyyz_0, \
                                     tdy_xxx_xyzz_0, tdy_xxx_xzzz_0, tdy_xxx_yyyy_0, tdy_xxx_yyyz_0, tdy_xxx_yyzz_0, \
                                     tdy_xxx_yzzz_0, tdy_xxx_zzzz_0, tdy_xxy_xxxx_0, tdz_xxx_xxxx_0, tdz_xxx_xxxy_0, \
                                     tdz_xxx_xxxz_0, tdz_xxx_xxyy_0, tdz_xxx_xxyz_0, tdz_xxx_xxzz_0, tdz_xxx_xyyy_0, \
                                     tdz_xxx_xyyz_0, tdz_xxx_xyzz_0, tdz_xxx_xzzz_0, tdz_xxx_yyyy_0, tdz_xxx_yyyz_0, \
                                     tdz_xxx_yyzz_0, tdz_xxx_yzzz_0, tdz_xxx_zzzz_0, tdz_xxy_xxxx_0, tlx_xx_xxxx_0, \
                                     tlx_xx_xxxy_0, tlx_xx_xxxz_0, tlx_xx_xxyy_0, tlx_xx_xxyz_0, tlx_xx_xxzz_0, \
                                     tlx_xx_xyyy_0, tlx_xx_xyyz_0, tlx_xx_xyzz_0, tlx_xx_xzzz_0, tlx_xx_yyyy_0, \
                                     tlx_xx_yyyz_0, tlx_xx_yyzz_0, tlx_xx_yzzz_0, tlx_xx_zzzz_0, tlx_xxx_xxx_0, \
                                     tlx_xxx_xxxx_0, tlx_xxx_xxxy_0, tlx_xxx_xxxz_0, tlx_xxx_xxy_0, tlx_xxx_xxyy_0, \
                                     tlx_xxx_xxyz_0, tlx_xxx_xxz_0, tlx_xxx_xxzz_0, tlx_xxx_xyy_0, tlx_xxx_xyyy_0, \
                                     tlx_xxx_xyyz_0, tlx_xxx_xyz_0, tlx_xxx_xyzz_0, tlx_xxx_xzz_0, tlx_xxx_xzzz_0, \
                                     tlx_xxx_yyy_0, tlx_xxx_yyyy_0, tlx_xxx_yyyz_0, tlx_xxx_yyz_0, tlx_xxx_yyzz_0, \
                                     tlx_xxx_yzz_0, tlx_xxx_yzzz_0, tlx_xxx_zzz_0, tlx_xxx_zzzz_0, tlx_xxxx_xxxx_0, \
                                     tlx_xxxx_xxxy_0, tlx_xxxx_xxxz_0, tlx_xxxx_xxyy_0, tlx_xxxx_xxyz_0, tlx_xxxx_xxzz_0, \
                                     tlx_xxxx_xyyy_0, tlx_xxxx_xyyz_0, tlx_xxxx_xyzz_0, tlx_xxxx_xzzz_0, tlx_xxxx_yyyy_0, \
                                     tlx_xxxx_yyyz_0, tlx_xxxx_yyzz_0, tlx_xxxx_yzzz_0, tlx_xxxx_zzzz_0, tlx_xxxy_xxxx_0, \
                                     tlx_xxxy_xxxy_0, tlx_xxy_xxx_0, tlx_xxy_xxxx_0, tlx_xxy_xxxy_0, tlx_xxy_xxy_0, \
                                     tlx_xy_xxxx_0, tlx_xy_xxxy_0, tly_xx_xxxx_0, tly_xx_xxxy_0, tly_xx_xxxz_0, \
                                     tly_xx_xxyy_0, tly_xx_xxyz_0, tly_xx_xxzz_0, tly_xx_xyyy_0, tly_xx_xyyz_0, \
                                     tly_xx_xyzz_0, tly_xx_xzzz_0, tly_xx_yyyy_0, tly_xx_yyyz_0, tly_xx_yyzz_0, \
                                     tly_xx_yzzz_0, tly_xx_zzzz_0, tly_xxx_xxx_0, tly_xxx_xxxx_0, tly_xxx_xxxy_0, \
                                     tly_xxx_xxxz_0, tly_xxx_xxy_0, tly_xxx_xxyy_0, tly_xxx_xxyz_0, tly_xxx_xxz_0, \
                                     tly_xxx_xxzz_0, tly_xxx_xyy_0, tly_xxx_xyyy_0, tly_xxx_xyyz_0, tly_xxx_xyz_0, \
                                     tly_xxx_xyzz_0, tly_xxx_xzz_0, tly_xxx_xzzz_0, tly_xxx_yyy_0, tly_xxx_yyyy_0, \
                                     tly_xxx_yyyz_0, tly_xxx_yyz_0, tly_xxx_yyzz_0, tly_xxx_yzz_0, tly_xxx_yzzz_0, \
                                     tly_xxx_zzz_0, tly_xxx_zzzz_0, tly_xxxx_xxxx_0, tly_xxxx_xxxy_0, tly_xxxx_xxxz_0, \
                                     tly_xxxx_xxyy_0, tly_xxxx_xxyz_0, tly_xxxx_xxzz_0, tly_xxxx_xyyy_0, tly_xxxx_xyyz_0, \
                                     tly_xxxx_xyzz_0, tly_xxxx_xzzz_0, tly_xxxx_yyyy_0, tly_xxxx_yyyz_0, tly_xxxx_yyzz_0, \
                                     tly_xxxx_yzzz_0, tly_xxxx_zzzz_0, tly_xxxy_xxxx_0, tly_xxy_xxx_0, tly_xxy_xxxx_0, \
                                     tly_xy_xxxx_0, tlz_xx_xxxx_0, tlz_xx_xxxy_0, tlz_xx_xxxz_0, tlz_xx_xxyy_0, \
                                     tlz_xx_xxyz_0, tlz_xx_xxzz_0, tlz_xx_xyyy_0, tlz_xx_xyyz_0, tlz_xx_xyzz_0, \
                                     tlz_xx_xzzz_0, tlz_xx_yyyy_0, tlz_xx_yyyz_0, tlz_xx_yyzz_0, tlz_xx_yzzz_0, \
                                     tlz_xx_zzzz_0, tlz_xxx_xxx_0, tlz_xxx_xxxx_0, tlz_xxx_xxxy_0, tlz_xxx_xxxz_0, \
                                     tlz_xxx_xxy_0, tlz_xxx_xxyy_0, tlz_xxx_xxyz_0, tlz_xxx_xxz_0, tlz_xxx_xxzz_0, \
                                     tlz_xxx_xyy_0, tlz_xxx_xyyy_0, tlz_xxx_xyyz_0, tlz_xxx_xyz_0, tlz_xxx_xyzz_0, \
                                     tlz_xxx_xzz_0, tlz_xxx_xzzz_0, tlz_xxx_yyy_0, tlz_xxx_yyyy_0, tlz_xxx_yyyz_0, \
                                     tlz_xxx_yyz_0, tlz_xxx_yyzz_0, tlz_xxx_yzz_0, tlz_xxx_yzzz_0, tlz_xxx_zzz_0, \
                                     tlz_xxx_zzzz_0, tlz_xxxx_xxxx_0, tlz_xxxx_xxxy_0, tlz_xxxx_xxxz_0, tlz_xxxx_xxyy_0, \
                                     tlz_xxxx_xxyz_0, tlz_xxxx_xxzz_0, tlz_xxxx_xyyy_0, tlz_xxxx_xyyz_0, tlz_xxxx_xyzz_0, \
                                     tlz_xxxx_xzzz_0, tlz_xxxx_yyyy_0, tlz_xxxx_yyyz_0, tlz_xxxx_yyzz_0, tlz_xxxx_yzzz_0, \
                                     tlz_xxxx_zzzz_0, tlz_xxxy_xxxx_0, tlz_xxy_xxx_0, tlz_xxy_xxxx_0, tlz_xy_xxxx_0, \
                                     tpy_xxx_xxxx_0, tpy_xxx_xxxy_0, tpy_xxx_xxxz_0, tpy_xxx_xxyy_0, tpy_xxx_xxyz_0, \
                                     tpy_xxx_xxzz_0, tpy_xxx_xyyy_0, tpy_xxx_xyyz_0, tpy_xxx_xyzz_0, tpy_xxx_xzzz_0, \
                                     tpy_xxx_yyyy_0, tpy_xxx_yyyz_0, tpy_xxx_yyzz_0, tpy_xxx_yzzz_0, tpy_xxx_zzzz_0, \
                                     tpy_xxy_xxxx_0, tpz_xxx_xxxx_0, tpz_xxx_xxxy_0, tpz_xxx_xxxz_0, tpz_xxx_xxyy_0, \
                                     tpz_xxx_xxyz_0, tpz_xxx_xxzz_0, tpz_xxx_xyyy_0, tpz_xxx_xyyz_0, tpz_xxx_xyzz_0, \
                                     tpz_xxx_xzzz_0, tpz_xxx_yyyy_0, tpz_xxx_yyyz_0, tpz_xxx_yyzz_0, tpz_xxx_yzzz_0, \
                                     tpz_xxx_zzzz_0, tpz_xxy_xxxx_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxxx_xxxx_0[j] = pa_x[j] * tlx_xxx_xxxx_0[j] + 1.5 * fl1_fx * tlx_xx_xxxx_0[j] + 2.0 * fl1_fx * tlx_xxx_xxx_0[j];

            tly_xxxx_xxxx_0[j] = pa_x[j] * tly_xxx_xxxx_0[j] + 1.5 * fl1_fx * tly_xx_xxxx_0[j] + 2.0 * fl1_fx * tly_xxx_xxx_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxxx_0[j];

            tlz_xxxx_xxxx_0[j] = pa_x[j] * tlz_xxx_xxxx_0[j] + 1.5 * fl1_fx * tlz_xx_xxxx_0[j] + 2.0 * fl1_fx * tlz_xxx_xxx_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxxx_0[j];

            tlx_xxxx_xxxy_0[j] = pa_x[j] * tlx_xxx_xxxy_0[j] + 1.5 * fl1_fx * tlx_xx_xxxy_0[j] + 1.5 * fl1_fx * tlx_xxx_xxy_0[j];

            tly_xxxx_xxxy_0[j] = pa_x[j] * tly_xxx_xxxy_0[j] + 1.5 * fl1_fx * tly_xx_xxxy_0[j] + 1.5 * fl1_fx * tly_xxx_xxy_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxxy_0[j];

            tlz_xxxx_xxxy_0[j] = pa_x[j] * tlz_xxx_xxxy_0[j] + 1.5 * fl1_fx * tlz_xx_xxxy_0[j] + 1.5 * fl1_fx * tlz_xxx_xxy_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxxy_0[j];

            tlx_xxxx_xxxz_0[j] = pa_x[j] * tlx_xxx_xxxz_0[j] + 1.5 * fl1_fx * tlx_xx_xxxz_0[j] + 1.5 * fl1_fx * tlx_xxx_xxz_0[j];

            tly_xxxx_xxxz_0[j] = pa_x[j] * tly_xxx_xxxz_0[j] + 1.5 * fl1_fx * tly_xx_xxxz_0[j] + 1.5 * fl1_fx * tly_xxx_xxz_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxxz_0[j];

            tlz_xxxx_xxxz_0[j] = pa_x[j] * tlz_xxx_xxxz_0[j] + 1.5 * fl1_fx * tlz_xx_xxxz_0[j] + 1.5 * fl1_fx * tlz_xxx_xxz_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxxz_0[j];

            tlx_xxxx_xxyy_0[j] = pa_x[j] * tlx_xxx_xxyy_0[j] + 1.5 * fl1_fx * tlx_xx_xxyy_0[j] + fl1_fx * tlx_xxx_xyy_0[j];

            tly_xxxx_xxyy_0[j] = pa_x[j] * tly_xxx_xxyy_0[j] + 1.5 * fl1_fx * tly_xx_xxyy_0[j] + fl1_fx * tly_xxx_xyy_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxyy_0[j];

            tlz_xxxx_xxyy_0[j] = pa_x[j] * tlz_xxx_xxyy_0[j] + 1.5 * fl1_fx * tlz_xx_xxyy_0[j] + fl1_fx * tlz_xxx_xyy_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxyy_0[j];

            tlx_xxxx_xxyz_0[j] = pa_x[j] * tlx_xxx_xxyz_0[j] + 1.5 * fl1_fx * tlx_xx_xxyz_0[j] + fl1_fx * tlx_xxx_xyz_0[j];

            tly_xxxx_xxyz_0[j] = pa_x[j] * tly_xxx_xxyz_0[j] + 1.5 * fl1_fx * tly_xx_xxyz_0[j] + fl1_fx * tly_xxx_xyz_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxyz_0[j];

            tlz_xxxx_xxyz_0[j] = pa_x[j] * tlz_xxx_xxyz_0[j] + 1.5 * fl1_fx * tlz_xx_xxyz_0[j] + fl1_fx * tlz_xxx_xyz_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxyz_0[j];

            tlx_xxxx_xxzz_0[j] = pa_x[j] * tlx_xxx_xxzz_0[j] + 1.5 * fl1_fx * tlx_xx_xxzz_0[j] + fl1_fx * tlx_xxx_xzz_0[j];

            tly_xxxx_xxzz_0[j] = pa_x[j] * tly_xxx_xxzz_0[j] + 1.5 * fl1_fx * tly_xx_xxzz_0[j] + fl1_fx * tly_xxx_xzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xxzz_0[j];

            tlz_xxxx_xxzz_0[j] = pa_x[j] * tlz_xxx_xxzz_0[j] + 1.5 * fl1_fx * tlz_xx_xxzz_0[j] + fl1_fx * tlz_xxx_xzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xxzz_0[j];

            tlx_xxxx_xyyy_0[j] = pa_x[j] * tlx_xxx_xyyy_0[j] + 1.5 * fl1_fx * tlx_xx_xyyy_0[j] + 0.5 * fl1_fx * tlx_xxx_yyy_0[j];

            tly_xxxx_xyyy_0[j] = pa_x[j] * tly_xxx_xyyy_0[j] + 1.5 * fl1_fx * tly_xx_xyyy_0[j] + 0.5 * fl1_fx * tly_xxx_yyy_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xyyy_0[j];

            tlz_xxxx_xyyy_0[j] = pa_x[j] * tlz_xxx_xyyy_0[j] + 1.5 * fl1_fx * tlz_xx_xyyy_0[j] + 0.5 * fl1_fx * tlz_xxx_yyy_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xyyy_0[j];

            tlx_xxxx_xyyz_0[j] = pa_x[j] * tlx_xxx_xyyz_0[j] + 1.5 * fl1_fx * tlx_xx_xyyz_0[j] + 0.5 * fl1_fx * tlx_xxx_yyz_0[j];

            tly_xxxx_xyyz_0[j] = pa_x[j] * tly_xxx_xyyz_0[j] + 1.5 * fl1_fx * tly_xx_xyyz_0[j] + 0.5 * fl1_fx * tly_xxx_yyz_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xyyz_0[j];

            tlz_xxxx_xyyz_0[j] = pa_x[j] * tlz_xxx_xyyz_0[j] + 1.5 * fl1_fx * tlz_xx_xyyz_0[j] + 0.5 * fl1_fx * tlz_xxx_yyz_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xyyz_0[j];

            tlx_xxxx_xyzz_0[j] = pa_x[j] * tlx_xxx_xyzz_0[j] + 1.5 * fl1_fx * tlx_xx_xyzz_0[j] + 0.5 * fl1_fx * tlx_xxx_yzz_0[j];

            tly_xxxx_xyzz_0[j] = pa_x[j] * tly_xxx_xyzz_0[j] + 1.5 * fl1_fx * tly_xx_xyzz_0[j] + 0.5 * fl1_fx * tly_xxx_yzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xyzz_0[j];

            tlz_xxxx_xyzz_0[j] = pa_x[j] * tlz_xxx_xyzz_0[j] + 1.5 * fl1_fx * tlz_xx_xyzz_0[j] + 0.5 * fl1_fx * tlz_xxx_yzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xyzz_0[j];

            tlx_xxxx_xzzz_0[j] = pa_x[j] * tlx_xxx_xzzz_0[j] + 1.5 * fl1_fx * tlx_xx_xzzz_0[j] + 0.5 * fl1_fx * tlx_xxx_zzz_0[j];

            tly_xxxx_xzzz_0[j] = pa_x[j] * tly_xxx_xzzz_0[j] + 1.5 * fl1_fx * tly_xx_xzzz_0[j] + 0.5 * fl1_fx * tly_xxx_zzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxx_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxx_xzzz_0[j];

            tlz_xxxx_xzzz_0[j] = pa_x[j] * tlz_xxx_xzzz_0[j] + 1.5 * fl1_fx * tlz_xx_xzzz_0[j] + 0.5 * fl1_fx * tlz_xxx_zzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxx_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxx_xzzz_0[j];

            tlx_xxxx_yyyy_0[j] = pa_x[j] * tlx_xxx_yyyy_0[j] + 1.5 * fl1_fx * tlx_xx_yyyy_0[j];

            tly_xxxx_yyyy_0[j] = pa_x[j] * tly_xxx_yyyy_0[j] + 1.5 * fl1_fx * tly_xx_yyyy_0[j] + 0.5 * fl1_fx * tpz_xxx_yyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xxx_yyyy_0[j];

            tlz_xxxx_yyyy_0[j] = pa_x[j] * tlz_xxx_yyyy_0[j] + 1.5 * fl1_fx * tlz_xx_yyyy_0[j] - 0.5 * fl1_fx * tpy_xxx_yyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xxx_yyyy_0[j];

            tlx_xxxx_yyyz_0[j] = pa_x[j] * tlx_xxx_yyyz_0[j] + 1.5 * fl1_fx * tlx_xx_yyyz_0[j];

            tly_xxxx_yyyz_0[j] = pa_x[j] * tly_xxx_yyyz_0[j] + 1.5 * fl1_fx * tly_xx_yyyz_0[j] + 0.5 * fl1_fx * tpz_xxx_yyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xxx_yyyz_0[j];

            tlz_xxxx_yyyz_0[j] = pa_x[j] * tlz_xxx_yyyz_0[j] + 1.5 * fl1_fx * tlz_xx_yyyz_0[j] - 0.5 * fl1_fx * tpy_xxx_yyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xxx_yyyz_0[j];

            tlx_xxxx_yyzz_0[j] = pa_x[j] * tlx_xxx_yyzz_0[j] + 1.5 * fl1_fx * tlx_xx_yyzz_0[j];

            tly_xxxx_yyzz_0[j] = pa_x[j] * tly_xxx_yyzz_0[j] + 1.5 * fl1_fx * tly_xx_yyzz_0[j] + 0.5 * fl1_fx * tpz_xxx_yyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xxx_yyzz_0[j];

            tlz_xxxx_yyzz_0[j] = pa_x[j] * tlz_xxx_yyzz_0[j] + 1.5 * fl1_fx * tlz_xx_yyzz_0[j] - 0.5 * fl1_fx * tpy_xxx_yyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xxx_yyzz_0[j];

            tlx_xxxx_yzzz_0[j] = pa_x[j] * tlx_xxx_yzzz_0[j] + 1.5 * fl1_fx * tlx_xx_yzzz_0[j];

            tly_xxxx_yzzz_0[j] = pa_x[j] * tly_xxx_yzzz_0[j] + 1.5 * fl1_fx * tly_xx_yzzz_0[j] + 0.5 * fl1_fx * tpz_xxx_yzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xxx_yzzz_0[j];

            tlz_xxxx_yzzz_0[j] = pa_x[j] * tlz_xxx_yzzz_0[j] + 1.5 * fl1_fx * tlz_xx_yzzz_0[j] - 0.5 * fl1_fx * tpy_xxx_yzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xxx_yzzz_0[j];

            tlx_xxxx_zzzz_0[j] = pa_x[j] * tlx_xxx_zzzz_0[j] + 1.5 * fl1_fx * tlx_xx_zzzz_0[j];

            tly_xxxx_zzzz_0[j] = pa_x[j] * tly_xxx_zzzz_0[j] + 1.5 * fl1_fx * tly_xx_zzzz_0[j] + 0.5 * fl1_fx * tpz_xxx_zzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xxx_zzzz_0[j];

            tlz_xxxx_zzzz_0[j] = pa_x[j] * tlz_xxx_zzzz_0[j] + 1.5 * fl1_fx * tlz_xx_zzzz_0[j] - 0.5 * fl1_fx * tpy_xxx_zzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xxx_zzzz_0[j];

            tlx_xxxy_xxxx_0[j] = pa_x[j] * tlx_xxy_xxxx_0[j] + fl1_fx * tlx_xy_xxxx_0[j] + 2.0 * fl1_fx * tlx_xxy_xxx_0[j];

            tly_xxxy_xxxx_0[j] = pa_x[j] * tly_xxy_xxxx_0[j] + fl1_fx * tly_xy_xxxx_0[j] + 2.0 * fl1_fx * tly_xxy_xxx_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxxx_0[j];

            tlz_xxxy_xxxx_0[j] = pa_x[j] * tlz_xxy_xxxx_0[j] + fl1_fx * tlz_xy_xxxx_0[j] + 2.0 * fl1_fx * tlz_xxy_xxx_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxxx_0[j];

            tlx_xxxy_xxxy_0[j] = pa_x[j] * tlx_xxy_xxxy_0[j] + fl1_fx * tlx_xy_xxxy_0[j] + 1.5 * fl1_fx * tlx_xxy_xxy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_49_98(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tly_xxy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 16);

        auto tlz_xxy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 16);

        auto tlx_xxy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 17);

        auto tly_xxy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 17);

        auto tlz_xxy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 17);

        auto tlx_xxy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 18);

        auto tly_xxy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 18);

        auto tlz_xxy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 18);

        auto tlx_xxy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 19);

        auto tly_xxy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 19);

        auto tlz_xxy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 19);

        auto tlx_xxy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 20);

        auto tly_xxy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 20);

        auto tlz_xxy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 20);

        auto tlx_xxy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 21);

        auto tly_xxy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 21);

        auto tlz_xxy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 21);

        auto tlx_xxy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 22);

        auto tly_xxy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 22);

        auto tlz_xxy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 22);

        auto tlx_xxy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 23);

        auto tly_xxy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 23);

        auto tlz_xxy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 23);

        auto tlx_xxy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 24);

        auto tly_xxy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 24);

        auto tlz_xxy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 24);

        auto tlx_xxy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 25);

        auto tly_xxy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 25);

        auto tlz_xxy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 25);

        auto tlx_xxy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 26);

        auto tly_xxy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 26);

        auto tlz_xxy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 26);

        auto tlx_xxy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 27);

        auto tly_xxy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 27);

        auto tlz_xxy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 27);

        auto tlx_xxy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 28);

        auto tly_xxy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 28);

        auto tlz_xxy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 28);

        auto tlx_xxy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 29);

        auto tly_xxy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 29);

        auto tlz_xxy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 29);

        auto tlx_xxz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 30);

        auto tly_xxz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 30);

        auto tlz_xxz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 30);

        auto tlx_xxz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 31);

        auto tly_xxz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 31);

        auto tlz_xxz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 31);

        auto tlx_xxz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 32);

        auto tly_xxz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 32);

        auto tly_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 16);

        auto tlz_xy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 16);

        auto tlx_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 17);

        auto tly_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 17);

        auto tlz_xy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 17);

        auto tlx_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 18);

        auto tly_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 18);

        auto tlz_xy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 18);

        auto tlx_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 19);

        auto tly_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 19);

        auto tlz_xy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 19);

        auto tlx_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 20);

        auto tly_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 20);

        auto tlz_xy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 20);

        auto tlx_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 21);

        auto tly_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 21);

        auto tlz_xy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 21);

        auto tlx_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 22);

        auto tly_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 22);

        auto tlz_xy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 22);

        auto tlx_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 23);

        auto tly_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 23);

        auto tlz_xy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 23);

        auto tlx_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 24);

        auto tly_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 24);

        auto tlz_xy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 24);

        auto tlx_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 25);

        auto tly_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 25);

        auto tlz_xy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 25);

        auto tlx_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 26);

        auto tly_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 26);

        auto tlz_xy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 26);

        auto tlx_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 27);

        auto tly_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 27);

        auto tlz_xy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 27);

        auto tlx_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 28);

        auto tly_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 28);

        auto tlz_xy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 28);

        auto tlx_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 29);

        auto tly_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 29);

        auto tlz_xy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 29);

        auto tlx_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 30);

        auto tly_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 30);

        auto tlz_xz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 30);

        auto tlx_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 31);

        auto tly_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 31);

        auto tlz_xz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 31);

        auto tlx_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 32);

        auto tly_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 32);

        auto tly_xxy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 11);

        auto tlz_xxy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 11);

        auto tlx_xxy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 12);

        auto tly_xxy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 12);

        auto tlz_xxy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 12);

        auto tlx_xxy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 13);

        auto tly_xxy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 13);

        auto tlz_xxy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 13);

        auto tlx_xxy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 14);

        auto tly_xxy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 14);

        auto tlz_xxy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 14);

        auto tlx_xxy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 15);

        auto tly_xxy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 15);

        auto tlz_xxy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 15);

        auto tlx_xxy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 16);

        auto tly_xxy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 16);

        auto tlz_xxy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 16);

        auto tlx_xxy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 17);

        auto tly_xxy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 17);

        auto tlz_xxy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 17);

        auto tlx_xxy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 18);

        auto tly_xxy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 18);

        auto tlz_xxy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 18);

        auto tlx_xxy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 19);

        auto tly_xxy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 19);

        auto tlz_xxy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 19);

        auto tlx_xxz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 20);

        auto tly_xxz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 20);

        auto tlz_xxz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 20);

        auto tlx_xxz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 21);

        auto tly_xxz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 21);

        auto tlz_xxz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 21);

        auto tlx_xxz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 22);

        auto tly_xxz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 22);

        auto tpy_xxy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 16);

        auto tpz_xxy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 16);

        auto tpy_xxy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 17);

        auto tpz_xxy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 17);

        auto tpy_xxy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 18);

        auto tpz_xxy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 18);

        auto tpy_xxy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 19);

        auto tpz_xxy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 19);

        auto tpy_xxy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 20);

        auto tpz_xxy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 20);

        auto tpy_xxy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 21);

        auto tpz_xxy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 21);

        auto tpy_xxy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 22);

        auto tpz_xxy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 22);

        auto tpy_xxy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 23);

        auto tpz_xxy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 23);

        auto tpy_xxy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 24);

        auto tpz_xxy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 24);

        auto tpy_xxy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 25);

        auto tpz_xxy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 25);

        auto tpy_xxy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 26);

        auto tpz_xxy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 26);

        auto tpy_xxy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 27);

        auto tpz_xxy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 27);

        auto tpy_xxy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 28);

        auto tpz_xxy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 28);

        auto tpy_xxy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 29);

        auto tpz_xxy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 29);

        auto tpy_xxz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 30);

        auto tpz_xxz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 30);

        auto tpy_xxz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 31);

        auto tpz_xxz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 31);

        auto tpz_xxz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tdy_xxy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 16);

        auto tdz_xxy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 16);

        auto tdy_xxy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 17);

        auto tdz_xxy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 17);

        auto tdy_xxy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 18);

        auto tdz_xxy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 18);

        auto tdy_xxy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 19);

        auto tdz_xxy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 19);

        auto tdy_xxy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 20);

        auto tdz_xxy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 20);

        auto tdy_xxy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 21);

        auto tdz_xxy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 21);

        auto tdy_xxy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 22);

        auto tdz_xxy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 22);

        auto tdy_xxy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 23);

        auto tdz_xxy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 23);

        auto tdy_xxy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 24);

        auto tdz_xxy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 24);

        auto tdy_xxy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 25);

        auto tdz_xxy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 25);

        auto tdy_xxy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 26);

        auto tdz_xxy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 26);

        auto tdy_xxy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 27);

        auto tdz_xxy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 27);

        auto tdy_xxy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 28);

        auto tdz_xxy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 28);

        auto tdy_xxy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 29);

        auto tdz_xxy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 29);

        auto tdy_xxz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 30);

        auto tdz_xxz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 30);

        auto tdy_xxz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 31);

        auto tdz_xxz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 31);

        auto tdz_xxz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 32);

        // set up pointers to integrals

        auto tly_xxxy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 16);

        auto tlz_xxxy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 16);

        auto tlx_xxxy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 17);

        auto tly_xxxy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 17);

        auto tlz_xxxy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 17);

        auto tlx_xxxy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 18);

        auto tly_xxxy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 18);

        auto tlz_xxxy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 18);

        auto tlx_xxxy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 19);

        auto tly_xxxy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 19);

        auto tlz_xxxy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 19);

        auto tlx_xxxy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 20);

        auto tly_xxxy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 20);

        auto tlz_xxxy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 20);

        auto tlx_xxxy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 21);

        auto tly_xxxy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 21);

        auto tlz_xxxy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 21);

        auto tlx_xxxy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 22);

        auto tly_xxxy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 22);

        auto tlz_xxxy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 22);

        auto tlx_xxxy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 23);

        auto tly_xxxy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 23);

        auto tlz_xxxy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 23);

        auto tlx_xxxy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 24);

        auto tly_xxxy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 24);

        auto tlz_xxxy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 24);

        auto tlx_xxxy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 25);

        auto tly_xxxy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 25);

        auto tlz_xxxy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 25);

        auto tlx_xxxy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 26);

        auto tly_xxxy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 26);

        auto tlz_xxxy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 26);

        auto tlx_xxxy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 27);

        auto tly_xxxy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 27);

        auto tlz_xxxy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 27);

        auto tlx_xxxy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 28);

        auto tly_xxxy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 28);

        auto tlz_xxxy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 28);

        auto tlx_xxxy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 29);

        auto tly_xxxy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 29);

        auto tlz_xxxy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 29);

        auto tlx_xxxz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 30);

        auto tly_xxxz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 30);

        auto tlz_xxxz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 30);

        auto tlx_xxxz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 31);

        auto tly_xxxz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 31);

        auto tlz_xxxz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 31);

        auto tlx_xxxz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 32);

        auto tly_xxxz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 32);

        // Batch of Integrals (49,98)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxy_xxxy_0, tdy_xxy_xxxz_0, tdy_xxy_xxyy_0, \
                                     tdy_xxy_xxyz_0, tdy_xxy_xxzz_0, tdy_xxy_xyyy_0, tdy_xxy_xyyz_0, tdy_xxy_xyzz_0, \
                                     tdy_xxy_xzzz_0, tdy_xxy_yyyy_0, tdy_xxy_yyyz_0, tdy_xxy_yyzz_0, tdy_xxy_yzzz_0, \
                                     tdy_xxy_zzzz_0, tdy_xxz_xxxx_0, tdy_xxz_xxxy_0, tdz_xxy_xxxy_0, tdz_xxy_xxxz_0, \
                                     tdz_xxy_xxyy_0, tdz_xxy_xxyz_0, tdz_xxy_xxzz_0, tdz_xxy_xyyy_0, tdz_xxy_xyyz_0, \
                                     tdz_xxy_xyzz_0, tdz_xxy_xzzz_0, tdz_xxy_yyyy_0, tdz_xxy_yyyz_0, tdz_xxy_yyzz_0, \
                                     tdz_xxy_yzzz_0, tdz_xxy_zzzz_0, tdz_xxz_xxxx_0, tdz_xxz_xxxy_0, tdz_xxz_xxxz_0, \
                                     tlx_xxxy_xxxz_0, tlx_xxxy_xxyy_0, tlx_xxxy_xxyz_0, tlx_xxxy_xxzz_0, tlx_xxxy_xyyy_0, \
                                     tlx_xxxy_xyyz_0, tlx_xxxy_xyzz_0, tlx_xxxy_xzzz_0, tlx_xxxy_yyyy_0, tlx_xxxy_yyyz_0, \
                                     tlx_xxxy_yyzz_0, tlx_xxxy_yzzz_0, tlx_xxxy_zzzz_0, tlx_xxxz_xxxx_0, tlx_xxxz_xxxy_0, \
                                     tlx_xxxz_xxxz_0, tlx_xxy_xxxz_0, tlx_xxy_xxyy_0, tlx_xxy_xxyz_0, tlx_xxy_xxz_0, \
                                     tlx_xxy_xxzz_0, tlx_xxy_xyy_0, tlx_xxy_xyyy_0, tlx_xxy_xyyz_0, tlx_xxy_xyz_0, \
                                     tlx_xxy_xyzz_0, tlx_xxy_xzz_0, tlx_xxy_xzzz_0, tlx_xxy_yyy_0, tlx_xxy_yyyy_0, \
                                     tlx_xxy_yyyz_0, tlx_xxy_yyz_0, tlx_xxy_yyzz_0, tlx_xxy_yzz_0, tlx_xxy_yzzz_0, \
                                     tlx_xxy_zzz_0, tlx_xxy_zzzz_0, tlx_xxz_xxx_0, tlx_xxz_xxxx_0, tlx_xxz_xxxy_0, \
                                     tlx_xxz_xxxz_0, tlx_xxz_xxy_0, tlx_xxz_xxz_0, tlx_xy_xxxz_0, tlx_xy_xxyy_0, \
                                     tlx_xy_xxyz_0, tlx_xy_xxzz_0, tlx_xy_xyyy_0, tlx_xy_xyyz_0, tlx_xy_xyzz_0, \
                                     tlx_xy_xzzz_0, tlx_xy_yyyy_0, tlx_xy_yyyz_0, tlx_xy_yyzz_0, tlx_xy_yzzz_0, \
                                     tlx_xy_zzzz_0, tlx_xz_xxxx_0, tlx_xz_xxxy_0, tlx_xz_xxxz_0, tly_xxxy_xxxy_0, \
                                     tly_xxxy_xxxz_0, tly_xxxy_xxyy_0, tly_xxxy_xxyz_0, tly_xxxy_xxzz_0, tly_xxxy_xyyy_0, \
                                     tly_xxxy_xyyz_0, tly_xxxy_xyzz_0, tly_xxxy_xzzz_0, tly_xxxy_yyyy_0, tly_xxxy_yyyz_0, \
                                     tly_xxxy_yyzz_0, tly_xxxy_yzzz_0, tly_xxxy_zzzz_0, tly_xxxz_xxxx_0, tly_xxxz_xxxy_0, \
                                     tly_xxxz_xxxz_0, tly_xxy_xxxy_0, tly_xxy_xxxz_0, tly_xxy_xxy_0, tly_xxy_xxyy_0, \
                                     tly_xxy_xxyz_0, tly_xxy_xxz_0, tly_xxy_xxzz_0, tly_xxy_xyy_0, tly_xxy_xyyy_0, \
                                     tly_xxy_xyyz_0, tly_xxy_xyz_0, tly_xxy_xyzz_0, tly_xxy_xzz_0, tly_xxy_xzzz_0, \
                                     tly_xxy_yyy_0, tly_xxy_yyyy_0, tly_xxy_yyyz_0, tly_xxy_yyz_0, tly_xxy_yyzz_0, \
                                     tly_xxy_yzz_0, tly_xxy_yzzz_0, tly_xxy_zzz_0, tly_xxy_zzzz_0, tly_xxz_xxx_0, \
                                     tly_xxz_xxxx_0, tly_xxz_xxxy_0, tly_xxz_xxxz_0, tly_xxz_xxy_0, tly_xxz_xxz_0, \
                                     tly_xy_xxxy_0, tly_xy_xxxz_0, tly_xy_xxyy_0, tly_xy_xxyz_0, tly_xy_xxzz_0, \
                                     tly_xy_xyyy_0, tly_xy_xyyz_0, tly_xy_xyzz_0, tly_xy_xzzz_0, tly_xy_yyyy_0, \
                                     tly_xy_yyyz_0, tly_xy_yyzz_0, tly_xy_yzzz_0, tly_xy_zzzz_0, tly_xz_xxxx_0, \
                                     tly_xz_xxxy_0, tly_xz_xxxz_0, tlz_xxxy_xxxy_0, tlz_xxxy_xxxz_0, tlz_xxxy_xxyy_0, \
                                     tlz_xxxy_xxyz_0, tlz_xxxy_xxzz_0, tlz_xxxy_xyyy_0, tlz_xxxy_xyyz_0, tlz_xxxy_xyzz_0, \
                                     tlz_xxxy_xzzz_0, tlz_xxxy_yyyy_0, tlz_xxxy_yyyz_0, tlz_xxxy_yyzz_0, tlz_xxxy_yzzz_0, \
                                     tlz_xxxy_zzzz_0, tlz_xxxz_xxxx_0, tlz_xxxz_xxxy_0, tlz_xxy_xxxy_0, tlz_xxy_xxxz_0, \
                                     tlz_xxy_xxy_0, tlz_xxy_xxyy_0, tlz_xxy_xxyz_0, tlz_xxy_xxz_0, tlz_xxy_xxzz_0, \
                                     tlz_xxy_xyy_0, tlz_xxy_xyyy_0, tlz_xxy_xyyz_0, tlz_xxy_xyz_0, tlz_xxy_xyzz_0, \
                                     tlz_xxy_xzz_0, tlz_xxy_xzzz_0, tlz_xxy_yyy_0, tlz_xxy_yyyy_0, tlz_xxy_yyyz_0, \
                                     tlz_xxy_yyz_0, tlz_xxy_yyzz_0, tlz_xxy_yzz_0, tlz_xxy_yzzz_0, tlz_xxy_zzz_0, \
                                     tlz_xxy_zzzz_0, tlz_xxz_xxx_0, tlz_xxz_xxxx_0, tlz_xxz_xxxy_0, tlz_xxz_xxy_0, \
                                     tlz_xy_xxxy_0, tlz_xy_xxxz_0, tlz_xy_xxyy_0, tlz_xy_xxyz_0, tlz_xy_xxzz_0, \
                                     tlz_xy_xyyy_0, tlz_xy_xyyz_0, tlz_xy_xyzz_0, tlz_xy_xzzz_0, tlz_xy_yyyy_0, \
                                     tlz_xy_yyyz_0, tlz_xy_yyzz_0, tlz_xy_yzzz_0, tlz_xy_zzzz_0, tlz_xz_xxxx_0, \
                                     tlz_xz_xxxy_0, tpy_xxy_xxxy_0, tpy_xxy_xxxz_0, tpy_xxy_xxyy_0, tpy_xxy_xxyz_0, \
                                     tpy_xxy_xxzz_0, tpy_xxy_xyyy_0, tpy_xxy_xyyz_0, tpy_xxy_xyzz_0, tpy_xxy_xzzz_0, \
                                     tpy_xxy_yyyy_0, tpy_xxy_yyyz_0, tpy_xxy_yyzz_0, tpy_xxy_yzzz_0, tpy_xxy_zzzz_0, \
                                     tpy_xxz_xxxx_0, tpy_xxz_xxxy_0, tpz_xxy_xxxy_0, tpz_xxy_xxxz_0, tpz_xxy_xxyy_0, \
                                     tpz_xxy_xxyz_0, tpz_xxy_xxzz_0, tpz_xxy_xyyy_0, tpz_xxy_xyyz_0, tpz_xxy_xyzz_0, \
                                     tpz_xxy_xzzz_0, tpz_xxy_yyyy_0, tpz_xxy_yyyz_0, tpz_xxy_yyzz_0, tpz_xxy_yzzz_0, \
                                     tpz_xxy_zzzz_0, tpz_xxz_xxxx_0, tpz_xxz_xxxy_0, tpz_xxz_xxxz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tly_xxxy_xxxy_0[j] = pa_x[j] * tly_xxy_xxxy_0[j] + fl1_fx * tly_xy_xxxy_0[j] + 1.5 * fl1_fx * tly_xxy_xxy_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxxy_0[j];

            tlz_xxxy_xxxy_0[j] = pa_x[j] * tlz_xxy_xxxy_0[j] + fl1_fx * tlz_xy_xxxy_0[j] + 1.5 * fl1_fx * tlz_xxy_xxy_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxxy_0[j];

            tlx_xxxy_xxxz_0[j] = pa_x[j] * tlx_xxy_xxxz_0[j] + fl1_fx * tlx_xy_xxxz_0[j] + 1.5 * fl1_fx * tlx_xxy_xxz_0[j];

            tly_xxxy_xxxz_0[j] = pa_x[j] * tly_xxy_xxxz_0[j] + fl1_fx * tly_xy_xxxz_0[j] + 1.5 * fl1_fx * tly_xxy_xxz_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxxz_0[j];

            tlz_xxxy_xxxz_0[j] = pa_x[j] * tlz_xxy_xxxz_0[j] + fl1_fx * tlz_xy_xxxz_0[j] + 1.5 * fl1_fx * tlz_xxy_xxz_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxxz_0[j];

            tlx_xxxy_xxyy_0[j] = pa_x[j] * tlx_xxy_xxyy_0[j] + fl1_fx * tlx_xy_xxyy_0[j] + fl1_fx * tlx_xxy_xyy_0[j];

            tly_xxxy_xxyy_0[j] = pa_x[j] * tly_xxy_xxyy_0[j] + fl1_fx * tly_xy_xxyy_0[j] + fl1_fx * tly_xxy_xyy_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxyy_0[j];

            tlz_xxxy_xxyy_0[j] = pa_x[j] * tlz_xxy_xxyy_0[j] + fl1_fx * tlz_xy_xxyy_0[j] + fl1_fx * tlz_xxy_xyy_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxyy_0[j];

            tlx_xxxy_xxyz_0[j] = pa_x[j] * tlx_xxy_xxyz_0[j] + fl1_fx * tlx_xy_xxyz_0[j] + fl1_fx * tlx_xxy_xyz_0[j];

            tly_xxxy_xxyz_0[j] = pa_x[j] * tly_xxy_xxyz_0[j] + fl1_fx * tly_xy_xxyz_0[j] + fl1_fx * tly_xxy_xyz_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxyz_0[j];

            tlz_xxxy_xxyz_0[j] = pa_x[j] * tlz_xxy_xxyz_0[j] + fl1_fx * tlz_xy_xxyz_0[j] + fl1_fx * tlz_xxy_xyz_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxyz_0[j];

            tlx_xxxy_xxzz_0[j] = pa_x[j] * tlx_xxy_xxzz_0[j] + fl1_fx * tlx_xy_xxzz_0[j] + fl1_fx * tlx_xxy_xzz_0[j];

            tly_xxxy_xxzz_0[j] = pa_x[j] * tly_xxy_xxzz_0[j] + fl1_fx * tly_xy_xxzz_0[j] + fl1_fx * tly_xxy_xzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xxzz_0[j];

            tlz_xxxy_xxzz_0[j] = pa_x[j] * tlz_xxy_xxzz_0[j] + fl1_fx * tlz_xy_xxzz_0[j] + fl1_fx * tlz_xxy_xzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xxzz_0[j];

            tlx_xxxy_xyyy_0[j] = pa_x[j] * tlx_xxy_xyyy_0[j] + fl1_fx * tlx_xy_xyyy_0[j] + 0.5 * fl1_fx * tlx_xxy_yyy_0[j];

            tly_xxxy_xyyy_0[j] = pa_x[j] * tly_xxy_xyyy_0[j] + fl1_fx * tly_xy_xyyy_0[j] + 0.5 * fl1_fx * tly_xxy_yyy_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xyyy_0[j];

            tlz_xxxy_xyyy_0[j] = pa_x[j] * tlz_xxy_xyyy_0[j] + fl1_fx * tlz_xy_xyyy_0[j] + 0.5 * fl1_fx * tlz_xxy_yyy_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xyyy_0[j];

            tlx_xxxy_xyyz_0[j] = pa_x[j] * tlx_xxy_xyyz_0[j] + fl1_fx * tlx_xy_xyyz_0[j] + 0.5 * fl1_fx * tlx_xxy_yyz_0[j];

            tly_xxxy_xyyz_0[j] = pa_x[j] * tly_xxy_xyyz_0[j] + fl1_fx * tly_xy_xyyz_0[j] + 0.5 * fl1_fx * tly_xxy_yyz_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xyyz_0[j];

            tlz_xxxy_xyyz_0[j] = pa_x[j] * tlz_xxy_xyyz_0[j] + fl1_fx * tlz_xy_xyyz_0[j] + 0.5 * fl1_fx * tlz_xxy_yyz_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xyyz_0[j];

            tlx_xxxy_xyzz_0[j] = pa_x[j] * tlx_xxy_xyzz_0[j] + fl1_fx * tlx_xy_xyzz_0[j] + 0.5 * fl1_fx * tlx_xxy_yzz_0[j];

            tly_xxxy_xyzz_0[j] = pa_x[j] * tly_xxy_xyzz_0[j] + fl1_fx * tly_xy_xyzz_0[j] + 0.5 * fl1_fx * tly_xxy_yzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xyzz_0[j];

            tlz_xxxy_xyzz_0[j] = pa_x[j] * tlz_xxy_xyzz_0[j] + fl1_fx * tlz_xy_xyzz_0[j] + 0.5 * fl1_fx * tlz_xxy_yzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xyzz_0[j];

            tlx_xxxy_xzzz_0[j] = pa_x[j] * tlx_xxy_xzzz_0[j] + fl1_fx * tlx_xy_xzzz_0[j] + 0.5 * fl1_fx * tlx_xxy_zzz_0[j];

            tly_xxxy_xzzz_0[j] = pa_x[j] * tly_xxy_xzzz_0[j] + fl1_fx * tly_xy_xzzz_0[j] + 0.5 * fl1_fx * tly_xxy_zzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxy_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_xzzz_0[j];

            tlz_xxxy_xzzz_0[j] = pa_x[j] * tlz_xxy_xzzz_0[j] + fl1_fx * tlz_xy_xzzz_0[j] + 0.5 * fl1_fx * tlz_xxy_zzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxy_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_xzzz_0[j];

            tlx_xxxy_yyyy_0[j] = pa_x[j] * tlx_xxy_yyyy_0[j] + fl1_fx * tlx_xy_yyyy_0[j];

            tly_xxxy_yyyy_0[j] =
                pa_x[j] * tly_xxy_yyyy_0[j] + fl1_fx * tly_xy_yyyy_0[j] + 0.5 * fl1_fx * tpz_xxy_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yyyy_0[j];

            tlz_xxxy_yyyy_0[j] =
                pa_x[j] * tlz_xxy_yyyy_0[j] + fl1_fx * tlz_xy_yyyy_0[j] - 0.5 * fl1_fx * tpy_xxy_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yyyy_0[j];

            tlx_xxxy_yyyz_0[j] = pa_x[j] * tlx_xxy_yyyz_0[j] + fl1_fx * tlx_xy_yyyz_0[j];

            tly_xxxy_yyyz_0[j] =
                pa_x[j] * tly_xxy_yyyz_0[j] + fl1_fx * tly_xy_yyyz_0[j] + 0.5 * fl1_fx * tpz_xxy_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yyyz_0[j];

            tlz_xxxy_yyyz_0[j] =
                pa_x[j] * tlz_xxy_yyyz_0[j] + fl1_fx * tlz_xy_yyyz_0[j] - 0.5 * fl1_fx * tpy_xxy_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yyyz_0[j];

            tlx_xxxy_yyzz_0[j] = pa_x[j] * tlx_xxy_yyzz_0[j] + fl1_fx * tlx_xy_yyzz_0[j];

            tly_xxxy_yyzz_0[j] =
                pa_x[j] * tly_xxy_yyzz_0[j] + fl1_fx * tly_xy_yyzz_0[j] + 0.5 * fl1_fx * tpz_xxy_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yyzz_0[j];

            tlz_xxxy_yyzz_0[j] =
                pa_x[j] * tlz_xxy_yyzz_0[j] + fl1_fx * tlz_xy_yyzz_0[j] - 0.5 * fl1_fx * tpy_xxy_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yyzz_0[j];

            tlx_xxxy_yzzz_0[j] = pa_x[j] * tlx_xxy_yzzz_0[j] + fl1_fx * tlx_xy_yzzz_0[j];

            tly_xxxy_yzzz_0[j] =
                pa_x[j] * tly_xxy_yzzz_0[j] + fl1_fx * tly_xy_yzzz_0[j] + 0.5 * fl1_fx * tpz_xxy_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_yzzz_0[j];

            tlz_xxxy_yzzz_0[j] =
                pa_x[j] * tlz_xxy_yzzz_0[j] + fl1_fx * tlz_xy_yzzz_0[j] - 0.5 * fl1_fx * tpy_xxy_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_yzzz_0[j];

            tlx_xxxy_zzzz_0[j] = pa_x[j] * tlx_xxy_zzzz_0[j] + fl1_fx * tlx_xy_zzzz_0[j];

            tly_xxxy_zzzz_0[j] =
                pa_x[j] * tly_xxy_zzzz_0[j] + fl1_fx * tly_xy_zzzz_0[j] + 0.5 * fl1_fx * tpz_xxy_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxy_zzzz_0[j];

            tlz_xxxy_zzzz_0[j] =
                pa_x[j] * tlz_xxy_zzzz_0[j] + fl1_fx * tlz_xy_zzzz_0[j] - 0.5 * fl1_fx * tpy_xxy_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxy_zzzz_0[j];

            tlx_xxxz_xxxx_0[j] = pa_x[j] * tlx_xxz_xxxx_0[j] + fl1_fx * tlx_xz_xxxx_0[j] + 2.0 * fl1_fx * tlx_xxz_xxx_0[j];

            tly_xxxz_xxxx_0[j] = pa_x[j] * tly_xxz_xxxx_0[j] + fl1_fx * tly_xz_xxxx_0[j] + 2.0 * fl1_fx * tly_xxz_xxx_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxxx_0[j];

            tlz_xxxz_xxxx_0[j] = pa_x[j] * tlz_xxz_xxxx_0[j] + fl1_fx * tlz_xz_xxxx_0[j] + 2.0 * fl1_fx * tlz_xxz_xxx_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxxx_0[j];

            tlx_xxxz_xxxy_0[j] = pa_x[j] * tlx_xxz_xxxy_0[j] + fl1_fx * tlx_xz_xxxy_0[j] + 1.5 * fl1_fx * tlx_xxz_xxy_0[j];

            tly_xxxz_xxxy_0[j] = pa_x[j] * tly_xxz_xxxy_0[j] + fl1_fx * tly_xz_xxxy_0[j] + 1.5 * fl1_fx * tly_xxz_xxy_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxxy_0[j];

            tlz_xxxz_xxxy_0[j] = pa_x[j] * tlz_xxz_xxxy_0[j] + fl1_fx * tlz_xz_xxxy_0[j] + 1.5 * fl1_fx * tlz_xxz_xxy_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxxy_0[j];

            tlx_xxxz_xxxz_0[j] = pa_x[j] * tlx_xxz_xxxz_0[j] + fl1_fx * tlx_xz_xxxz_0[j] + 1.5 * fl1_fx * tlx_xxz_xxz_0[j];

            tly_xxxz_xxxz_0[j] = pa_x[j] * tly_xxz_xxxz_0[j] + fl1_fx * tly_xz_xxxz_0[j] + 1.5 * fl1_fx * tly_xxz_xxz_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxxz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_98_147(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlz_xxz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 32);

        auto tlx_xxz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 33);

        auto tly_xxz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 33);

        auto tlz_xxz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 33);

        auto tlx_xxz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 34);

        auto tly_xxz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 34);

        auto tlz_xxz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 34);

        auto tlx_xxz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 35);

        auto tly_xxz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 35);

        auto tlz_xxz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 35);

        auto tlx_xxz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 36);

        auto tly_xxz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 36);

        auto tlz_xxz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 36);

        auto tlx_xxz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 37);

        auto tly_xxz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 37);

        auto tlz_xxz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 37);

        auto tlx_xxz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 38);

        auto tly_xxz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 38);

        auto tlz_xxz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 38);

        auto tlx_xxz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 39);

        auto tly_xxz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 39);

        auto tlz_xxz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 39);

        auto tlx_xxz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 40);

        auto tly_xxz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 40);

        auto tlz_xxz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 40);

        auto tlx_xxz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 41);

        auto tly_xxz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 41);

        auto tlz_xxz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 41);

        auto tlx_xxz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 42);

        auto tly_xxz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 42);

        auto tlz_xxz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 42);

        auto tlx_xxz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 43);

        auto tly_xxz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 43);

        auto tlz_xxz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 43);

        auto tlx_xxz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 44);

        auto tly_xxz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 44);

        auto tlz_xxz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 44);

        auto tlx_xyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 45);

        auto tly_xyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 45);

        auto tlz_xyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 45);

        auto tlx_xyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 46);

        auto tly_xyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 46);

        auto tlz_xyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 46);

        auto tlx_xyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 47);

        auto tly_xyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 47);

        auto tlz_xyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 47);

        auto tlx_xyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 48);

        auto tly_xyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 48);

        auto tlz_xyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 48);

        auto tlz_xz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 32);

        auto tlx_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 33);

        auto tly_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 33);

        auto tlz_xz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 33);

        auto tlx_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 34);

        auto tly_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 34);

        auto tlz_xz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 34);

        auto tlx_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 35);

        auto tly_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 35);

        auto tlz_xz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 35);

        auto tlx_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 36);

        auto tly_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 36);

        auto tlz_xz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 36);

        auto tlx_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 37);

        auto tly_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 37);

        auto tlz_xz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 37);

        auto tlx_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 38);

        auto tly_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 38);

        auto tlz_xz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 38);

        auto tlx_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 39);

        auto tly_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 39);

        auto tlz_xz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 39);

        auto tlx_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 40);

        auto tly_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 40);

        auto tlz_xz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 40);

        auto tlx_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 41);

        auto tly_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 41);

        auto tlz_xz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 41);

        auto tlx_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 42);

        auto tly_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 42);

        auto tlz_xz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 42);

        auto tlx_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 43);

        auto tly_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 43);

        auto tlz_xz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 43);

        auto tlx_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 44);

        auto tly_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 44);

        auto tlz_xz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 44);

        auto tlx_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 45);

        auto tly_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tlz_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tlx_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 46);

        auto tly_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tlz_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tlx_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 47);

        auto tly_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tlz_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tlx_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 48);

        auto tly_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tlz_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tlz_xxz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 22);

        auto tlx_xxz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 23);

        auto tly_xxz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 23);

        auto tlz_xxz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 23);

        auto tlx_xxz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 24);

        auto tly_xxz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 24);

        auto tlz_xxz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 24);

        auto tlx_xxz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 25);

        auto tly_xxz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 25);

        auto tlz_xxz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 25);

        auto tlx_xxz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 26);

        auto tly_xxz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 26);

        auto tlz_xxz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 26);

        auto tlx_xxz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 27);

        auto tly_xxz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 27);

        auto tlz_xxz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 27);

        auto tlx_xxz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 28);

        auto tly_xxz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 28);

        auto tlz_xxz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 28);

        auto tlx_xxz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 29);

        auto tly_xxz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 29);

        auto tlz_xxz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 29);

        auto tlx_xyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 30);

        auto tly_xyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 30);

        auto tlz_xyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 30);

        auto tlx_xyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 31);

        auto tly_xyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 31);

        auto tlz_xyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 31);

        auto tlx_xyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 32);

        auto tly_xyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 32);

        auto tlz_xyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 32);

        auto tlx_xyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 33);

        auto tly_xyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 33);

        auto tlz_xyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 33);

        auto tpy_xxz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 32);

        auto tpy_xxz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 33);

        auto tpz_xxz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 33);

        auto tpy_xxz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 34);

        auto tpz_xxz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 34);

        auto tpy_xxz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 35);

        auto tpz_xxz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 35);

        auto tpy_xxz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 36);

        auto tpz_xxz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 36);

        auto tpy_xxz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 37);

        auto tpz_xxz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 37);

        auto tpy_xxz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 38);

        auto tpz_xxz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 38);

        auto tpy_xxz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 39);

        auto tpz_xxz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 39);

        auto tpy_xxz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 40);

        auto tpz_xxz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 40);

        auto tpy_xxz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 41);

        auto tpz_xxz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 41);

        auto tpy_xxz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 42);

        auto tpz_xxz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 42);

        auto tpy_xxz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 43);

        auto tpz_xxz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 43);

        auto tpy_xxz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 44);

        auto tpz_xxz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 44);

        auto tpy_xyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 45);

        auto tpz_xyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 45);

        auto tpy_xyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 46);

        auto tpz_xyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 46);

        auto tpy_xyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 47);

        auto tpz_xyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 47);

        auto tpy_xyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 48);

        auto tpz_xyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 48);

        auto tdy_xxz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 32);

        auto tdy_xxz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 33);

        auto tdz_xxz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 33);

        auto tdy_xxz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 34);

        auto tdz_xxz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 34);

        auto tdy_xxz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 35);

        auto tdz_xxz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 35);

        auto tdy_xxz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 36);

        auto tdz_xxz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 36);

        auto tdy_xxz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 37);

        auto tdz_xxz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 37);

        auto tdy_xxz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 38);

        auto tdz_xxz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 38);

        auto tdy_xxz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 39);

        auto tdz_xxz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 39);

        auto tdy_xxz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 40);

        auto tdz_xxz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 40);

        auto tdy_xxz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 41);

        auto tdz_xxz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 41);

        auto tdy_xxz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 42);

        auto tdz_xxz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 42);

        auto tdy_xxz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 43);

        auto tdz_xxz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 43);

        auto tdy_xxz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 44);

        auto tdz_xxz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 44);

        auto tdy_xyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 45);

        auto tdz_xyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 45);

        auto tdy_xyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 46);

        auto tdz_xyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 46);

        auto tdy_xyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 47);

        auto tdz_xyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 47);

        auto tdy_xyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 48);

        auto tdz_xyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 48);

        // set up pointers to integrals

        auto tlz_xxxz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 32);

        auto tlx_xxxz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 33);

        auto tly_xxxz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 33);

        auto tlz_xxxz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 33);

        auto tlx_xxxz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 34);

        auto tly_xxxz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 34);

        auto tlz_xxxz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 34);

        auto tlx_xxxz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 35);

        auto tly_xxxz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 35);

        auto tlz_xxxz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 35);

        auto tlx_xxxz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 36);

        auto tly_xxxz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 36);

        auto tlz_xxxz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 36);

        auto tlx_xxxz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 37);

        auto tly_xxxz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 37);

        auto tlz_xxxz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 37);

        auto tlx_xxxz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 38);

        auto tly_xxxz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 38);

        auto tlz_xxxz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 38);

        auto tlx_xxxz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 39);

        auto tly_xxxz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 39);

        auto tlz_xxxz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 39);

        auto tlx_xxxz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 40);

        auto tly_xxxz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 40);

        auto tlz_xxxz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 40);

        auto tlx_xxxz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 41);

        auto tly_xxxz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 41);

        auto tlz_xxxz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 41);

        auto tlx_xxxz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 42);

        auto tly_xxxz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 42);

        auto tlz_xxxz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 42);

        auto tlx_xxxz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 43);

        auto tly_xxxz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 43);

        auto tlz_xxxz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 43);

        auto tlx_xxxz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 44);

        auto tly_xxxz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 44);

        auto tlz_xxxz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 44);

        auto tlx_xxyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 45);

        auto tly_xxyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 45);

        auto tlz_xxyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 45);

        auto tlx_xxyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 46);

        auto tly_xxyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 46);

        auto tlz_xxyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 46);

        auto tlx_xxyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 47);

        auto tly_xxyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 47);

        auto tlz_xxyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 47);

        auto tlx_xxyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 48);

        auto tly_xxyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 48);

        auto tlz_xxyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 48);

        // Batch of Integrals (98,147)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xxz_xxxz_0, tdy_xxz_xxyy_0, tdy_xxz_xxyz_0, \
                                     tdy_xxz_xxzz_0, tdy_xxz_xyyy_0, tdy_xxz_xyyz_0, tdy_xxz_xyzz_0, tdy_xxz_xzzz_0, \
                                     tdy_xxz_yyyy_0, tdy_xxz_yyyz_0, tdy_xxz_yyzz_0, tdy_xxz_yzzz_0, tdy_xxz_zzzz_0, \
                                     tdy_xyy_xxxx_0, tdy_xyy_xxxy_0, tdy_xyy_xxxz_0, tdy_xyy_xxyy_0, tdz_xxz_xxyy_0, \
                                     tdz_xxz_xxyz_0, tdz_xxz_xxzz_0, tdz_xxz_xyyy_0, tdz_xxz_xyyz_0, tdz_xxz_xyzz_0, \
                                     tdz_xxz_xzzz_0, tdz_xxz_yyyy_0, tdz_xxz_yyyz_0, tdz_xxz_yyzz_0, tdz_xxz_yzzz_0, \
                                     tdz_xxz_zzzz_0, tdz_xyy_xxxx_0, tdz_xyy_xxxy_0, tdz_xyy_xxxz_0, tdz_xyy_xxyy_0, \
                                     tlx_xxxz_xxyy_0, tlx_xxxz_xxyz_0, tlx_xxxz_xxzz_0, tlx_xxxz_xyyy_0, tlx_xxxz_xyyz_0, \
                                     tlx_xxxz_xyzz_0, tlx_xxxz_xzzz_0, tlx_xxxz_yyyy_0, tlx_xxxz_yyyz_0, tlx_xxxz_yyzz_0, \
                                     tlx_xxxz_yzzz_0, tlx_xxxz_zzzz_0, tlx_xxyy_xxxx_0, tlx_xxyy_xxxy_0, tlx_xxyy_xxxz_0, \
                                     tlx_xxyy_xxyy_0, tlx_xxz_xxyy_0, tlx_xxz_xxyz_0, tlx_xxz_xxzz_0, tlx_xxz_xyy_0, \
                                     tlx_xxz_xyyy_0, tlx_xxz_xyyz_0, tlx_xxz_xyz_0, tlx_xxz_xyzz_0, tlx_xxz_xzz_0, \
                                     tlx_xxz_xzzz_0, tlx_xxz_yyy_0, tlx_xxz_yyyy_0, tlx_xxz_yyyz_0, tlx_xxz_yyz_0, \
                                     tlx_xxz_yyzz_0, tlx_xxz_yzz_0, tlx_xxz_yzzz_0, tlx_xxz_zzz_0, tlx_xxz_zzzz_0, \
                                     tlx_xyy_xxx_0, tlx_xyy_xxxx_0, tlx_xyy_xxxy_0, tlx_xyy_xxxz_0, tlx_xyy_xxy_0, \
                                     tlx_xyy_xxyy_0, tlx_xyy_xxz_0, tlx_xyy_xyy_0, tlx_xz_xxyy_0, tlx_xz_xxyz_0, \
                                     tlx_xz_xxzz_0, tlx_xz_xyyy_0, tlx_xz_xyyz_0, tlx_xz_xyzz_0, tlx_xz_xzzz_0, \
                                     tlx_xz_yyyy_0, tlx_xz_yyyz_0, tlx_xz_yyzz_0, tlx_xz_yzzz_0, tlx_xz_zzzz_0, \
                                     tlx_yy_xxxx_0, tlx_yy_xxxy_0, tlx_yy_xxxz_0, tlx_yy_xxyy_0, tly_xxxz_xxyy_0, \
                                     tly_xxxz_xxyz_0, tly_xxxz_xxzz_0, tly_xxxz_xyyy_0, tly_xxxz_xyyz_0, tly_xxxz_xyzz_0, \
                                     tly_xxxz_xzzz_0, tly_xxxz_yyyy_0, tly_xxxz_yyyz_0, tly_xxxz_yyzz_0, tly_xxxz_yzzz_0, \
                                     tly_xxxz_zzzz_0, tly_xxyy_xxxx_0, tly_xxyy_xxxy_0, tly_xxyy_xxxz_0, tly_xxyy_xxyy_0, \
                                     tly_xxz_xxyy_0, tly_xxz_xxyz_0, tly_xxz_xxzz_0, tly_xxz_xyy_0, tly_xxz_xyyy_0, \
                                     tly_xxz_xyyz_0, tly_xxz_xyz_0, tly_xxz_xyzz_0, tly_xxz_xzz_0, tly_xxz_xzzz_0, \
                                     tly_xxz_yyy_0, tly_xxz_yyyy_0, tly_xxz_yyyz_0, tly_xxz_yyz_0, tly_xxz_yyzz_0, \
                                     tly_xxz_yzz_0, tly_xxz_yzzz_0, tly_xxz_zzz_0, tly_xxz_zzzz_0, tly_xyy_xxx_0, \
                                     tly_xyy_xxxx_0, tly_xyy_xxxy_0, tly_xyy_xxxz_0, tly_xyy_xxy_0, tly_xyy_xxyy_0, \
                                     tly_xyy_xxz_0, tly_xyy_xyy_0, tly_xz_xxyy_0, tly_xz_xxyz_0, tly_xz_xxzz_0, \
                                     tly_xz_xyyy_0, tly_xz_xyyz_0, tly_xz_xyzz_0, tly_xz_xzzz_0, tly_xz_yyyy_0, \
                                     tly_xz_yyyz_0, tly_xz_yyzz_0, tly_xz_yzzz_0, tly_xz_zzzz_0, tly_yy_xxxx_0, \
                                     tly_yy_xxxy_0, tly_yy_xxxz_0, tly_yy_xxyy_0, tlz_xxxz_xxxz_0, tlz_xxxz_xxyy_0, \
                                     tlz_xxxz_xxyz_0, tlz_xxxz_xxzz_0, tlz_xxxz_xyyy_0, tlz_xxxz_xyyz_0, tlz_xxxz_xyzz_0, \
                                     tlz_xxxz_xzzz_0, tlz_xxxz_yyyy_0, tlz_xxxz_yyyz_0, tlz_xxxz_yyzz_0, tlz_xxxz_yzzz_0, \
                                     tlz_xxxz_zzzz_0, tlz_xxyy_xxxx_0, tlz_xxyy_xxxy_0, tlz_xxyy_xxxz_0, tlz_xxyy_xxyy_0, \
                                     tlz_xxz_xxxz_0, tlz_xxz_xxyy_0, tlz_xxz_xxyz_0, tlz_xxz_xxz_0, tlz_xxz_xxzz_0, \
                                     tlz_xxz_xyy_0, tlz_xxz_xyyy_0, tlz_xxz_xyyz_0, tlz_xxz_xyz_0, tlz_xxz_xyzz_0, \
                                     tlz_xxz_xzz_0, tlz_xxz_xzzz_0, tlz_xxz_yyy_0, tlz_xxz_yyyy_0, tlz_xxz_yyyz_0, \
                                     tlz_xxz_yyz_0, tlz_xxz_yyzz_0, tlz_xxz_yzz_0, tlz_xxz_yzzz_0, tlz_xxz_zzz_0, \
                                     tlz_xxz_zzzz_0, tlz_xyy_xxx_0, tlz_xyy_xxxx_0, tlz_xyy_xxxy_0, tlz_xyy_xxxz_0, \
                                     tlz_xyy_xxy_0, tlz_xyy_xxyy_0, tlz_xyy_xxz_0, tlz_xyy_xyy_0, tlz_xz_xxxz_0, \
                                     tlz_xz_xxyy_0, tlz_xz_xxyz_0, tlz_xz_xxzz_0, tlz_xz_xyyy_0, tlz_xz_xyyz_0, \
                                     tlz_xz_xyzz_0, tlz_xz_xzzz_0, tlz_xz_yyyy_0, tlz_xz_yyyz_0, tlz_xz_yyzz_0, \
                                     tlz_xz_yzzz_0, tlz_xz_zzzz_0, tlz_yy_xxxx_0, tlz_yy_xxxy_0, tlz_yy_xxxz_0, \
                                     tlz_yy_xxyy_0, tpy_xxz_xxxz_0, tpy_xxz_xxyy_0, tpy_xxz_xxyz_0, tpy_xxz_xxzz_0, \
                                     tpy_xxz_xyyy_0, tpy_xxz_xyyz_0, tpy_xxz_xyzz_0, tpy_xxz_xzzz_0, tpy_xxz_yyyy_0, \
                                     tpy_xxz_yyyz_0, tpy_xxz_yyzz_0, tpy_xxz_yzzz_0, tpy_xxz_zzzz_0, tpy_xyy_xxxx_0, \
                                     tpy_xyy_xxxy_0, tpy_xyy_xxxz_0, tpy_xyy_xxyy_0, tpz_xxz_xxyy_0, tpz_xxz_xxyz_0, \
                                     tpz_xxz_xxzz_0, tpz_xxz_xyyy_0, tpz_xxz_xyyz_0, tpz_xxz_xyzz_0, tpz_xxz_xzzz_0, \
                                     tpz_xxz_yyyy_0, tpz_xxz_yyyz_0, tpz_xxz_yyzz_0, tpz_xxz_yzzz_0, tpz_xxz_zzzz_0, \
                                     tpz_xyy_xxxx_0, tpz_xyy_xxxy_0, tpz_xyy_xxxz_0, tpz_xyy_xxyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlz_xxxz_xxxz_0[j] = pa_x[j] * tlz_xxz_xxxz_0[j] + fl1_fx * tlz_xz_xxxz_0[j] + 1.5 * fl1_fx * tlz_xxz_xxz_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxxz_0[j];

            tlx_xxxz_xxyy_0[j] = pa_x[j] * tlx_xxz_xxyy_0[j] + fl1_fx * tlx_xz_xxyy_0[j] + fl1_fx * tlx_xxz_xyy_0[j];

            tly_xxxz_xxyy_0[j] = pa_x[j] * tly_xxz_xxyy_0[j] + fl1_fx * tly_xz_xxyy_0[j] + fl1_fx * tly_xxz_xyy_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxyy_0[j];

            tlz_xxxz_xxyy_0[j] = pa_x[j] * tlz_xxz_xxyy_0[j] + fl1_fx * tlz_xz_xxyy_0[j] + fl1_fx * tlz_xxz_xyy_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxyy_0[j];

            tlx_xxxz_xxyz_0[j] = pa_x[j] * tlx_xxz_xxyz_0[j] + fl1_fx * tlx_xz_xxyz_0[j] + fl1_fx * tlx_xxz_xyz_0[j];

            tly_xxxz_xxyz_0[j] = pa_x[j] * tly_xxz_xxyz_0[j] + fl1_fx * tly_xz_xxyz_0[j] + fl1_fx * tly_xxz_xyz_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxyz_0[j];

            tlz_xxxz_xxyz_0[j] = pa_x[j] * tlz_xxz_xxyz_0[j] + fl1_fx * tlz_xz_xxyz_0[j] + fl1_fx * tlz_xxz_xyz_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxyz_0[j];

            tlx_xxxz_xxzz_0[j] = pa_x[j] * tlx_xxz_xxzz_0[j] + fl1_fx * tlx_xz_xxzz_0[j] + fl1_fx * tlx_xxz_xzz_0[j];

            tly_xxxz_xxzz_0[j] = pa_x[j] * tly_xxz_xxzz_0[j] + fl1_fx * tly_xz_xxzz_0[j] + fl1_fx * tly_xxz_xzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xxzz_0[j];

            tlz_xxxz_xxzz_0[j] = pa_x[j] * tlz_xxz_xxzz_0[j] + fl1_fx * tlz_xz_xxzz_0[j] + fl1_fx * tlz_xxz_xzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xxzz_0[j];

            tlx_xxxz_xyyy_0[j] = pa_x[j] * tlx_xxz_xyyy_0[j] + fl1_fx * tlx_xz_xyyy_0[j] + 0.5 * fl1_fx * tlx_xxz_yyy_0[j];

            tly_xxxz_xyyy_0[j] = pa_x[j] * tly_xxz_xyyy_0[j] + fl1_fx * tly_xz_xyyy_0[j] + 0.5 * fl1_fx * tly_xxz_yyy_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xyyy_0[j];

            tlz_xxxz_xyyy_0[j] = pa_x[j] * tlz_xxz_xyyy_0[j] + fl1_fx * tlz_xz_xyyy_0[j] + 0.5 * fl1_fx * tlz_xxz_yyy_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xyyy_0[j];

            tlx_xxxz_xyyz_0[j] = pa_x[j] * tlx_xxz_xyyz_0[j] + fl1_fx * tlx_xz_xyyz_0[j] + 0.5 * fl1_fx * tlx_xxz_yyz_0[j];

            tly_xxxz_xyyz_0[j] = pa_x[j] * tly_xxz_xyyz_0[j] + fl1_fx * tly_xz_xyyz_0[j] + 0.5 * fl1_fx * tly_xxz_yyz_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xyyz_0[j];

            tlz_xxxz_xyyz_0[j] = pa_x[j] * tlz_xxz_xyyz_0[j] + fl1_fx * tlz_xz_xyyz_0[j] + 0.5 * fl1_fx * tlz_xxz_yyz_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xyyz_0[j];

            tlx_xxxz_xyzz_0[j] = pa_x[j] * tlx_xxz_xyzz_0[j] + fl1_fx * tlx_xz_xyzz_0[j] + 0.5 * fl1_fx * tlx_xxz_yzz_0[j];

            tly_xxxz_xyzz_0[j] = pa_x[j] * tly_xxz_xyzz_0[j] + fl1_fx * tly_xz_xyzz_0[j] + 0.5 * fl1_fx * tly_xxz_yzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xyzz_0[j];

            tlz_xxxz_xyzz_0[j] = pa_x[j] * tlz_xxz_xyzz_0[j] + fl1_fx * tlz_xz_xyzz_0[j] + 0.5 * fl1_fx * tlz_xxz_yzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xyzz_0[j];

            tlx_xxxz_xzzz_0[j] = pa_x[j] * tlx_xxz_xzzz_0[j] + fl1_fx * tlx_xz_xzzz_0[j] + 0.5 * fl1_fx * tlx_xxz_zzz_0[j];

            tly_xxxz_xzzz_0[j] = pa_x[j] * tly_xxz_xzzz_0[j] + fl1_fx * tly_xz_xzzz_0[j] + 0.5 * fl1_fx * tly_xxz_zzz_0[j] +
                                 0.5 * fl1_fx * tpz_xxz_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_xzzz_0[j];

            tlz_xxxz_xzzz_0[j] = pa_x[j] * tlz_xxz_xzzz_0[j] + fl1_fx * tlz_xz_xzzz_0[j] + 0.5 * fl1_fx * tlz_xxz_zzz_0[j] -
                                 0.5 * fl1_fx * tpy_xxz_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_xzzz_0[j];

            tlx_xxxz_yyyy_0[j] = pa_x[j] * tlx_xxz_yyyy_0[j] + fl1_fx * tlx_xz_yyyy_0[j];

            tly_xxxz_yyyy_0[j] =
                pa_x[j] * tly_xxz_yyyy_0[j] + fl1_fx * tly_xz_yyyy_0[j] + 0.5 * fl1_fx * tpz_xxz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yyyy_0[j];

            tlz_xxxz_yyyy_0[j] =
                pa_x[j] * tlz_xxz_yyyy_0[j] + fl1_fx * tlz_xz_yyyy_0[j] - 0.5 * fl1_fx * tpy_xxz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yyyy_0[j];

            tlx_xxxz_yyyz_0[j] = pa_x[j] * tlx_xxz_yyyz_0[j] + fl1_fx * tlx_xz_yyyz_0[j];

            tly_xxxz_yyyz_0[j] =
                pa_x[j] * tly_xxz_yyyz_0[j] + fl1_fx * tly_xz_yyyz_0[j] + 0.5 * fl1_fx * tpz_xxz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yyyz_0[j];

            tlz_xxxz_yyyz_0[j] =
                pa_x[j] * tlz_xxz_yyyz_0[j] + fl1_fx * tlz_xz_yyyz_0[j] - 0.5 * fl1_fx * tpy_xxz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yyyz_0[j];

            tlx_xxxz_yyzz_0[j] = pa_x[j] * tlx_xxz_yyzz_0[j] + fl1_fx * tlx_xz_yyzz_0[j];

            tly_xxxz_yyzz_0[j] =
                pa_x[j] * tly_xxz_yyzz_0[j] + fl1_fx * tly_xz_yyzz_0[j] + 0.5 * fl1_fx * tpz_xxz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yyzz_0[j];

            tlz_xxxz_yyzz_0[j] =
                pa_x[j] * tlz_xxz_yyzz_0[j] + fl1_fx * tlz_xz_yyzz_0[j] - 0.5 * fl1_fx * tpy_xxz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yyzz_0[j];

            tlx_xxxz_yzzz_0[j] = pa_x[j] * tlx_xxz_yzzz_0[j] + fl1_fx * tlx_xz_yzzz_0[j];

            tly_xxxz_yzzz_0[j] =
                pa_x[j] * tly_xxz_yzzz_0[j] + fl1_fx * tly_xz_yzzz_0[j] + 0.5 * fl1_fx * tpz_xxz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_yzzz_0[j];

            tlz_xxxz_yzzz_0[j] =
                pa_x[j] * tlz_xxz_yzzz_0[j] + fl1_fx * tlz_xz_yzzz_0[j] - 0.5 * fl1_fx * tpy_xxz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_yzzz_0[j];

            tlx_xxxz_zzzz_0[j] = pa_x[j] * tlx_xxz_zzzz_0[j] + fl1_fx * tlx_xz_zzzz_0[j];

            tly_xxxz_zzzz_0[j] =
                pa_x[j] * tly_xxz_zzzz_0[j] + fl1_fx * tly_xz_zzzz_0[j] + 0.5 * fl1_fx * tpz_xxz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_xxz_zzzz_0[j];

            tlz_xxxz_zzzz_0[j] =
                pa_x[j] * tlz_xxz_zzzz_0[j] + fl1_fx * tlz_xz_zzzz_0[j] - 0.5 * fl1_fx * tpy_xxz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_xxz_zzzz_0[j];

            tlx_xxyy_xxxx_0[j] = pa_x[j] * tlx_xyy_xxxx_0[j] + 0.5 * fl1_fx * tlx_yy_xxxx_0[j] + 2.0 * fl1_fx * tlx_xyy_xxx_0[j];

            tly_xxyy_xxxx_0[j] = pa_x[j] * tly_xyy_xxxx_0[j] + 0.5 * fl1_fx * tly_yy_xxxx_0[j] + 2.0 * fl1_fx * tly_xyy_xxx_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxxx_0[j];

            tlz_xxyy_xxxx_0[j] = pa_x[j] * tlz_xyy_xxxx_0[j] + 0.5 * fl1_fx * tlz_yy_xxxx_0[j] + 2.0 * fl1_fx * tlz_xyy_xxx_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxxx_0[j];

            tlx_xxyy_xxxy_0[j] = pa_x[j] * tlx_xyy_xxxy_0[j] + 0.5 * fl1_fx * tlx_yy_xxxy_0[j] + 1.5 * fl1_fx * tlx_xyy_xxy_0[j];

            tly_xxyy_xxxy_0[j] = pa_x[j] * tly_xyy_xxxy_0[j] + 0.5 * fl1_fx * tly_yy_xxxy_0[j] + 1.5 * fl1_fx * tly_xyy_xxy_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxxy_0[j];

            tlz_xxyy_xxxy_0[j] = pa_x[j] * tlz_xyy_xxxy_0[j] + 0.5 * fl1_fx * tlz_yy_xxxy_0[j] + 1.5 * fl1_fx * tlz_xyy_xxy_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxxy_0[j];

            tlx_xxyy_xxxz_0[j] = pa_x[j] * tlx_xyy_xxxz_0[j] + 0.5 * fl1_fx * tlx_yy_xxxz_0[j] + 1.5 * fl1_fx * tlx_xyy_xxz_0[j];

            tly_xxyy_xxxz_0[j] = pa_x[j] * tly_xyy_xxxz_0[j] + 0.5 * fl1_fx * tly_yy_xxxz_0[j] + 1.5 * fl1_fx * tly_xyy_xxz_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxxz_0[j];

            tlz_xxyy_xxxz_0[j] = pa_x[j] * tlz_xyy_xxxz_0[j] + 0.5 * fl1_fx * tlz_yy_xxxz_0[j] + 1.5 * fl1_fx * tlz_xyy_xxz_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxxz_0[j];

            tlx_xxyy_xxyy_0[j] = pa_x[j] * tlx_xyy_xxyy_0[j] + 0.5 * fl1_fx * tlx_yy_xxyy_0[j] + fl1_fx * tlx_xyy_xyy_0[j];

            tly_xxyy_xxyy_0[j] = pa_x[j] * tly_xyy_xxyy_0[j] + 0.5 * fl1_fx * tly_yy_xxyy_0[j] + fl1_fx * tly_xyy_xyy_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxyy_0[j];

            tlz_xxyy_xxyy_0[j] = pa_x[j] * tlz_xyy_xxyy_0[j] + 0.5 * fl1_fx * tlz_yy_xxyy_0[j] + fl1_fx * tlz_xyy_xyy_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_147_195(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 49);

        auto tly_xyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tlz_xyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 49);

        auto tlx_xyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 50);

        auto tly_xyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 50);

        auto tlz_xyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 50);

        auto tlx_xyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 51);

        auto tly_xyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 51);

        auto tlz_xyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 51);

        auto tlx_xyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 52);

        auto tly_xyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 52);

        auto tlz_xyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 52);

        auto tlx_xyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 53);

        auto tly_xyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 53);

        auto tlz_xyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 53);

        auto tlx_xyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 54);

        auto tly_xyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 54);

        auto tlz_xyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 54);

        auto tlx_xyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 55);

        auto tly_xyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 55);

        auto tlz_xyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 55);

        auto tlx_xyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 56);

        auto tly_xyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 56);

        auto tlz_xyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 56);

        auto tlx_xyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 57);

        auto tly_xyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 57);

        auto tlz_xyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 57);

        auto tlx_xyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 58);

        auto tly_xyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 58);

        auto tlz_xyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 58);

        auto tlx_xyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 59);

        auto tly_xyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 59);

        auto tlz_xyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 59);

        auto tlx_xyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 60);

        auto tly_xyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 60);

        auto tlz_xyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 60);

        auto tlx_xyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 61);

        auto tly_xyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 61);

        auto tlz_xyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 61);

        auto tlx_xyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 62);

        auto tly_xyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 62);

        auto tlz_xyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 62);

        auto tlx_xyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 63);

        auto tly_xyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 63);

        auto tlz_xyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 63);

        auto tlx_xyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 64);

        auto tly_xyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 64);

        auto tlz_xyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 64);

        auto tlx_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 49);

        auto tly_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tlz_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tlx_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 50);

        auto tly_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tlz_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tlx_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 51);

        auto tly_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tlz_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tlx_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 52);

        auto tly_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tlz_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tlx_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 53);

        auto tly_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tlz_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tlx_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 54);

        auto tly_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tlz_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tlx_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 55);

        auto tly_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tlz_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tlx_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 56);

        auto tly_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tlz_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tlx_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 57);

        auto tly_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tlz_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tlx_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 58);

        auto tly_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tlz_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tlx_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 59);

        auto tly_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tlz_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tlx_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 60);

        auto tly_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tlz_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tlx_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 61);

        auto tly_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tlz_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tlx_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 62);

        auto tly_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tlz_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tlx_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 63);

        auto tly_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tlz_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tlx_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 64);

        auto tly_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tlz_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tlx_xyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 34);

        auto tly_xyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 34);

        auto tlz_xyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 34);

        auto tlx_xyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 35);

        auto tly_xyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 35);

        auto tlz_xyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 35);

        auto tlx_xyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 36);

        auto tly_xyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 36);

        auto tlz_xyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 36);

        auto tlx_xyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 37);

        auto tly_xyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 37);

        auto tlz_xyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 37);

        auto tlx_xyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 38);

        auto tly_xyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 38);

        auto tlz_xyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 38);

        auto tlx_xyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 39);

        auto tly_xyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 39);

        auto tlz_xyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 39);

        auto tlx_xyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 40);

        auto tly_xyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 40);

        auto tlz_xyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 40);

        auto tlx_xyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 41);

        auto tly_xyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 41);

        auto tlz_xyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 41);

        auto tlx_xyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 42);

        auto tly_xyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 42);

        auto tlz_xyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 42);

        auto tlx_xyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 43);

        auto tly_xyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 43);

        auto tlz_xyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 43);

        auto tlx_xyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 44);

        auto tly_xyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 44);

        auto tlz_xyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 44);

        auto tpy_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tpz_xyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 49);

        auto tpy_xyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 50);

        auto tpz_xyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 50);

        auto tpy_xyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 51);

        auto tpz_xyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 51);

        auto tpy_xyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 52);

        auto tpz_xyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 52);

        auto tpy_xyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 53);

        auto tpz_xyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 53);

        auto tpy_xyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 54);

        auto tpz_xyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 54);

        auto tpy_xyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 55);

        auto tpz_xyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 55);

        auto tpy_xyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 56);

        auto tpz_xyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 56);

        auto tpy_xyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 57);

        auto tpz_xyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 57);

        auto tpy_xyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 58);

        auto tpz_xyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 58);

        auto tpy_xyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 59);

        auto tpz_xyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 59);

        auto tpy_xyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 60);

        auto tpz_xyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 60);

        auto tpy_xyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 61);

        auto tpz_xyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 61);

        auto tpy_xyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 62);

        auto tpz_xyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 62);

        auto tpy_xyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 63);

        auto tpz_xyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 63);

        auto tpy_xyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 64);

        auto tpz_xyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 64);

        auto tdy_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 49);

        auto tdz_xyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 49);

        auto tdy_xyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 50);

        auto tdz_xyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 50);

        auto tdy_xyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 51);

        auto tdz_xyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 51);

        auto tdy_xyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 52);

        auto tdz_xyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 52);

        auto tdy_xyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 53);

        auto tdz_xyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 53);

        auto tdy_xyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 54);

        auto tdz_xyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 54);

        auto tdy_xyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 55);

        auto tdz_xyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 55);

        auto tdy_xyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 56);

        auto tdz_xyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 56);

        auto tdy_xyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 57);

        auto tdz_xyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 57);

        auto tdy_xyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 58);

        auto tdz_xyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 58);

        auto tdy_xyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 59);

        auto tdz_xyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 59);

        auto tdy_xyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 60);

        auto tdz_xyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 60);

        auto tdy_xyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 61);

        auto tdz_xyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 61);

        auto tdy_xyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 62);

        auto tdz_xyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 62);

        auto tdy_xyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 63);

        auto tdz_xyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 63);

        auto tdy_xyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 64);

        auto tdz_xyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 64);

        // set up pointers to integrals

        auto tlx_xxyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 49);

        auto tly_xxyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 49);

        auto tlz_xxyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 49);

        auto tlx_xxyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 50);

        auto tly_xxyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 50);

        auto tlz_xxyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 50);

        auto tlx_xxyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 51);

        auto tly_xxyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 51);

        auto tlz_xxyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 51);

        auto tlx_xxyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 52);

        auto tly_xxyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 52);

        auto tlz_xxyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 52);

        auto tlx_xxyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 53);

        auto tly_xxyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 53);

        auto tlz_xxyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 53);

        auto tlx_xxyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 54);

        auto tly_xxyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 54);

        auto tlz_xxyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 54);

        auto tlx_xxyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 55);

        auto tly_xxyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 55);

        auto tlz_xxyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 55);

        auto tlx_xxyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 56);

        auto tly_xxyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 56);

        auto tlz_xxyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 56);

        auto tlx_xxyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 57);

        auto tly_xxyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 57);

        auto tlz_xxyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 57);

        auto tlx_xxyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 58);

        auto tly_xxyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 58);

        auto tlz_xxyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 58);

        auto tlx_xxyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 59);

        auto tly_xxyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 59);

        auto tlz_xxyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 59);

        auto tlx_xxyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 60);

        auto tly_xxyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 60);

        auto tlz_xxyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 60);

        auto tlx_xxyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 61);

        auto tly_xxyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 61);

        auto tlz_xxyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 61);

        auto tlx_xxyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 62);

        auto tly_xxyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 62);

        auto tlz_xxyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 62);

        auto tlx_xxyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 63);

        auto tly_xxyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 63);

        auto tlz_xxyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 63);

        auto tlx_xxyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 64);

        auto tly_xxyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 64);

        auto tlz_xxyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 64);

        // Batch of Integrals (147,195)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xyy_xxyz_0, tdy_xyy_xxzz_0, tdy_xyy_xyyy_0, \
                                     tdy_xyy_xyyz_0, tdy_xyy_xyzz_0, tdy_xyy_xzzz_0, tdy_xyy_yyyy_0, tdy_xyy_yyyz_0, \
                                     tdy_xyy_yyzz_0, tdy_xyy_yzzz_0, tdy_xyy_zzzz_0, tdy_xyz_xxxx_0, tdy_xyz_xxxy_0, \
                                     tdy_xyz_xxxz_0, tdy_xyz_xxyy_0, tdy_xyz_xxyz_0, tdz_xyy_xxyz_0, tdz_xyy_xxzz_0, \
                                     tdz_xyy_xyyy_0, tdz_xyy_xyyz_0, tdz_xyy_xyzz_0, tdz_xyy_xzzz_0, tdz_xyy_yyyy_0, \
                                     tdz_xyy_yyyz_0, tdz_xyy_yyzz_0, tdz_xyy_yzzz_0, tdz_xyy_zzzz_0, tdz_xyz_xxxx_0, \
                                     tdz_xyz_xxxy_0, tdz_xyz_xxxz_0, tdz_xyz_xxyy_0, tdz_xyz_xxyz_0, tlx_xxyy_xxyz_0, \
                                     tlx_xxyy_xxzz_0, tlx_xxyy_xyyy_0, tlx_xxyy_xyyz_0, tlx_xxyy_xyzz_0, tlx_xxyy_xzzz_0, \
                                     tlx_xxyy_yyyy_0, tlx_xxyy_yyyz_0, tlx_xxyy_yyzz_0, tlx_xxyy_yzzz_0, tlx_xxyy_zzzz_0, \
                                     tlx_xxyz_xxxx_0, tlx_xxyz_xxxy_0, tlx_xxyz_xxxz_0, tlx_xxyz_xxyy_0, tlx_xxyz_xxyz_0, \
                                     tlx_xyy_xxyz_0, tlx_xyy_xxzz_0, tlx_xyy_xyyy_0, tlx_xyy_xyyz_0, tlx_xyy_xyz_0, \
                                     tlx_xyy_xyzz_0, tlx_xyy_xzz_0, tlx_xyy_xzzz_0, tlx_xyy_yyy_0, tlx_xyy_yyyy_0, \
                                     tlx_xyy_yyyz_0, tlx_xyy_yyz_0, tlx_xyy_yyzz_0, tlx_xyy_yzz_0, tlx_xyy_yzzz_0, \
                                     tlx_xyy_zzz_0, tlx_xyy_zzzz_0, tlx_xyz_xxx_0, tlx_xyz_xxxx_0, tlx_xyz_xxxy_0, \
                                     tlx_xyz_xxxz_0, tlx_xyz_xxy_0, tlx_xyz_xxyy_0, tlx_xyz_xxyz_0, tlx_xyz_xxz_0, \
                                     tlx_xyz_xyy_0, tlx_xyz_xyz_0, tlx_yy_xxyz_0, tlx_yy_xxzz_0, tlx_yy_xyyy_0, \
                                     tlx_yy_xyyz_0, tlx_yy_xyzz_0, tlx_yy_xzzz_0, tlx_yy_yyyy_0, tlx_yy_yyyz_0, \
                                     tlx_yy_yyzz_0, tlx_yy_yzzz_0, tlx_yy_zzzz_0, tlx_yz_xxxx_0, tlx_yz_xxxy_0, \
                                     tlx_yz_xxxz_0, tlx_yz_xxyy_0, tlx_yz_xxyz_0, tly_xxyy_xxyz_0, tly_xxyy_xxzz_0, \
                                     tly_xxyy_xyyy_0, tly_xxyy_xyyz_0, tly_xxyy_xyzz_0, tly_xxyy_xzzz_0, tly_xxyy_yyyy_0, \
                                     tly_xxyy_yyyz_0, tly_xxyy_yyzz_0, tly_xxyy_yzzz_0, tly_xxyy_zzzz_0, tly_xxyz_xxxx_0, \
                                     tly_xxyz_xxxy_0, tly_xxyz_xxxz_0, tly_xxyz_xxyy_0, tly_xxyz_xxyz_0, tly_xyy_xxyz_0, \
                                     tly_xyy_xxzz_0, tly_xyy_xyyy_0, tly_xyy_xyyz_0, tly_xyy_xyz_0, tly_xyy_xyzz_0, \
                                     tly_xyy_xzz_0, tly_xyy_xzzz_0, tly_xyy_yyy_0, tly_xyy_yyyy_0, tly_xyy_yyyz_0, \
                                     tly_xyy_yyz_0, tly_xyy_yyzz_0, tly_xyy_yzz_0, tly_xyy_yzzz_0, tly_xyy_zzz_0, \
                                     tly_xyy_zzzz_0, tly_xyz_xxx_0, tly_xyz_xxxx_0, tly_xyz_xxxy_0, tly_xyz_xxxz_0, \
                                     tly_xyz_xxy_0, tly_xyz_xxyy_0, tly_xyz_xxyz_0, tly_xyz_xxz_0, tly_xyz_xyy_0, \
                                     tly_xyz_xyz_0, tly_yy_xxyz_0, tly_yy_xxzz_0, tly_yy_xyyy_0, tly_yy_xyyz_0, \
                                     tly_yy_xyzz_0, tly_yy_xzzz_0, tly_yy_yyyy_0, tly_yy_yyyz_0, tly_yy_yyzz_0, \
                                     tly_yy_yzzz_0, tly_yy_zzzz_0, tly_yz_xxxx_0, tly_yz_xxxy_0, tly_yz_xxxz_0, \
                                     tly_yz_xxyy_0, tly_yz_xxyz_0, tlz_xxyy_xxyz_0, tlz_xxyy_xxzz_0, tlz_xxyy_xyyy_0, \
                                     tlz_xxyy_xyyz_0, tlz_xxyy_xyzz_0, tlz_xxyy_xzzz_0, tlz_xxyy_yyyy_0, tlz_xxyy_yyyz_0, \
                                     tlz_xxyy_yyzz_0, tlz_xxyy_yzzz_0, tlz_xxyy_zzzz_0, tlz_xxyz_xxxx_0, tlz_xxyz_xxxy_0, \
                                     tlz_xxyz_xxxz_0, tlz_xxyz_xxyy_0, tlz_xxyz_xxyz_0, tlz_xyy_xxyz_0, tlz_xyy_xxzz_0, \
                                     tlz_xyy_xyyy_0, tlz_xyy_xyyz_0, tlz_xyy_xyz_0, tlz_xyy_xyzz_0, tlz_xyy_xzz_0, \
                                     tlz_xyy_xzzz_0, tlz_xyy_yyy_0, tlz_xyy_yyyy_0, tlz_xyy_yyyz_0, tlz_xyy_yyz_0, \
                                     tlz_xyy_yyzz_0, tlz_xyy_yzz_0, tlz_xyy_yzzz_0, tlz_xyy_zzz_0, tlz_xyy_zzzz_0, \
                                     tlz_xyz_xxx_0, tlz_xyz_xxxx_0, tlz_xyz_xxxy_0, tlz_xyz_xxxz_0, tlz_xyz_xxy_0, \
                                     tlz_xyz_xxyy_0, tlz_xyz_xxyz_0, tlz_xyz_xxz_0, tlz_xyz_xyy_0, tlz_xyz_xyz_0, \
                                     tlz_yy_xxyz_0, tlz_yy_xxzz_0, tlz_yy_xyyy_0, tlz_yy_xyyz_0, tlz_yy_xyzz_0, \
                                     tlz_yy_xzzz_0, tlz_yy_yyyy_0, tlz_yy_yyyz_0, tlz_yy_yyzz_0, tlz_yy_yzzz_0, \
                                     tlz_yy_zzzz_0, tlz_yz_xxxx_0, tlz_yz_xxxy_0, tlz_yz_xxxz_0, tlz_yz_xxyy_0, \
                                     tlz_yz_xxyz_0, tpy_xyy_xxyz_0, tpy_xyy_xxzz_0, tpy_xyy_xyyy_0, tpy_xyy_xyyz_0, \
                                     tpy_xyy_xyzz_0, tpy_xyy_xzzz_0, tpy_xyy_yyyy_0, tpy_xyy_yyyz_0, tpy_xyy_yyzz_0, \
                                     tpy_xyy_yzzz_0, tpy_xyy_zzzz_0, tpy_xyz_xxxx_0, tpy_xyz_xxxy_0, tpy_xyz_xxxz_0, \
                                     tpy_xyz_xxyy_0, tpy_xyz_xxyz_0, tpz_xyy_xxyz_0, tpz_xyy_xxzz_0, tpz_xyy_xyyy_0, \
                                     tpz_xyy_xyyz_0, tpz_xyy_xyzz_0, tpz_xyy_xzzz_0, tpz_xyy_yyyy_0, tpz_xyy_yyyz_0, \
                                     tpz_xyy_yyzz_0, tpz_xyy_yzzz_0, tpz_xyy_zzzz_0, tpz_xyz_xxxx_0, tpz_xyz_xxxy_0, \
                                     tpz_xyz_xxxz_0, tpz_xyz_xxyy_0, tpz_xyz_xxyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxyy_xxyz_0[j] = pa_x[j] * tlx_xyy_xxyz_0[j] + 0.5 * fl1_fx * tlx_yy_xxyz_0[j] + fl1_fx * tlx_xyy_xyz_0[j];

            tly_xxyy_xxyz_0[j] = pa_x[j] * tly_xyy_xxyz_0[j] + 0.5 * fl1_fx * tly_yy_xxyz_0[j] + fl1_fx * tly_xyy_xyz_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxyz_0[j];

            tlz_xxyy_xxyz_0[j] = pa_x[j] * tlz_xyy_xxyz_0[j] + 0.5 * fl1_fx * tlz_yy_xxyz_0[j] + fl1_fx * tlz_xyy_xyz_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxyz_0[j];

            tlx_xxyy_xxzz_0[j] = pa_x[j] * tlx_xyy_xxzz_0[j] + 0.5 * fl1_fx * tlx_yy_xxzz_0[j] + fl1_fx * tlx_xyy_xzz_0[j];

            tly_xxyy_xxzz_0[j] = pa_x[j] * tly_xyy_xxzz_0[j] + 0.5 * fl1_fx * tly_yy_xxzz_0[j] + fl1_fx * tly_xyy_xzz_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xxzz_0[j];

            tlz_xxyy_xxzz_0[j] = pa_x[j] * tlz_xyy_xxzz_0[j] + 0.5 * fl1_fx * tlz_yy_xxzz_0[j] + fl1_fx * tlz_xyy_xzz_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xxzz_0[j];

            tlx_xxyy_xyyy_0[j] = pa_x[j] * tlx_xyy_xyyy_0[j] + 0.5 * fl1_fx * tlx_yy_xyyy_0[j] + 0.5 * fl1_fx * tlx_xyy_yyy_0[j];

            tly_xxyy_xyyy_0[j] = pa_x[j] * tly_xyy_xyyy_0[j] + 0.5 * fl1_fx * tly_yy_xyyy_0[j] + 0.5 * fl1_fx * tly_xyy_yyy_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xyyy_0[j];

            tlz_xxyy_xyyy_0[j] = pa_x[j] * tlz_xyy_xyyy_0[j] + 0.5 * fl1_fx * tlz_yy_xyyy_0[j] + 0.5 * fl1_fx * tlz_xyy_yyy_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xyyy_0[j];

            tlx_xxyy_xyyz_0[j] = pa_x[j] * tlx_xyy_xyyz_0[j] + 0.5 * fl1_fx * tlx_yy_xyyz_0[j] + 0.5 * fl1_fx * tlx_xyy_yyz_0[j];

            tly_xxyy_xyyz_0[j] = pa_x[j] * tly_xyy_xyyz_0[j] + 0.5 * fl1_fx * tly_yy_xyyz_0[j] + 0.5 * fl1_fx * tly_xyy_yyz_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xyyz_0[j];

            tlz_xxyy_xyyz_0[j] = pa_x[j] * tlz_xyy_xyyz_0[j] + 0.5 * fl1_fx * tlz_yy_xyyz_0[j] + 0.5 * fl1_fx * tlz_xyy_yyz_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xyyz_0[j];

            tlx_xxyy_xyzz_0[j] = pa_x[j] * tlx_xyy_xyzz_0[j] + 0.5 * fl1_fx * tlx_yy_xyzz_0[j] + 0.5 * fl1_fx * tlx_xyy_yzz_0[j];

            tly_xxyy_xyzz_0[j] = pa_x[j] * tly_xyy_xyzz_0[j] + 0.5 * fl1_fx * tly_yy_xyzz_0[j] + 0.5 * fl1_fx * tly_xyy_yzz_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xyzz_0[j];

            tlz_xxyy_xyzz_0[j] = pa_x[j] * tlz_xyy_xyzz_0[j] + 0.5 * fl1_fx * tlz_yy_xyzz_0[j] + 0.5 * fl1_fx * tlz_xyy_yzz_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xyzz_0[j];

            tlx_xxyy_xzzz_0[j] = pa_x[j] * tlx_xyy_xzzz_0[j] + 0.5 * fl1_fx * tlx_yy_xzzz_0[j] + 0.5 * fl1_fx * tlx_xyy_zzz_0[j];

            tly_xxyy_xzzz_0[j] = pa_x[j] * tly_xyy_xzzz_0[j] + 0.5 * fl1_fx * tly_yy_xzzz_0[j] + 0.5 * fl1_fx * tly_xyy_zzz_0[j] +
                                 0.5 * fl1_fx * tpz_xyy_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xyy_xzzz_0[j];

            tlz_xxyy_xzzz_0[j] = pa_x[j] * tlz_xyy_xzzz_0[j] + 0.5 * fl1_fx * tlz_yy_xzzz_0[j] + 0.5 * fl1_fx * tlz_xyy_zzz_0[j] -
                                 0.5 * fl1_fx * tpy_xyy_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xyy_xzzz_0[j];

            tlx_xxyy_yyyy_0[j] = pa_x[j] * tlx_xyy_yyyy_0[j] + 0.5 * fl1_fx * tlx_yy_yyyy_0[j];

            tly_xxyy_yyyy_0[j] = pa_x[j] * tly_xyy_yyyy_0[j] + 0.5 * fl1_fx * tly_yy_yyyy_0[j] + 0.5 * fl1_fx * tpz_xyy_yyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyy_yyyy_0[j];

            tlz_xxyy_yyyy_0[j] = pa_x[j] * tlz_xyy_yyyy_0[j] + 0.5 * fl1_fx * tlz_yy_yyyy_0[j] - 0.5 * fl1_fx * tpy_xyy_yyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyy_yyyy_0[j];

            tlx_xxyy_yyyz_0[j] = pa_x[j] * tlx_xyy_yyyz_0[j] + 0.5 * fl1_fx * tlx_yy_yyyz_0[j];

            tly_xxyy_yyyz_0[j] = pa_x[j] * tly_xyy_yyyz_0[j] + 0.5 * fl1_fx * tly_yy_yyyz_0[j] + 0.5 * fl1_fx * tpz_xyy_yyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyy_yyyz_0[j];

            tlz_xxyy_yyyz_0[j] = pa_x[j] * tlz_xyy_yyyz_0[j] + 0.5 * fl1_fx * tlz_yy_yyyz_0[j] - 0.5 * fl1_fx * tpy_xyy_yyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyy_yyyz_0[j];

            tlx_xxyy_yyzz_0[j] = pa_x[j] * tlx_xyy_yyzz_0[j] + 0.5 * fl1_fx * tlx_yy_yyzz_0[j];

            tly_xxyy_yyzz_0[j] = pa_x[j] * tly_xyy_yyzz_0[j] + 0.5 * fl1_fx * tly_yy_yyzz_0[j] + 0.5 * fl1_fx * tpz_xyy_yyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyy_yyzz_0[j];

            tlz_xxyy_yyzz_0[j] = pa_x[j] * tlz_xyy_yyzz_0[j] + 0.5 * fl1_fx * tlz_yy_yyzz_0[j] - 0.5 * fl1_fx * tpy_xyy_yyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyy_yyzz_0[j];

            tlx_xxyy_yzzz_0[j] = pa_x[j] * tlx_xyy_yzzz_0[j] + 0.5 * fl1_fx * tlx_yy_yzzz_0[j];

            tly_xxyy_yzzz_0[j] = pa_x[j] * tly_xyy_yzzz_0[j] + 0.5 * fl1_fx * tly_yy_yzzz_0[j] + 0.5 * fl1_fx * tpz_xyy_yzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyy_yzzz_0[j];

            tlz_xxyy_yzzz_0[j] = pa_x[j] * tlz_xyy_yzzz_0[j] + 0.5 * fl1_fx * tlz_yy_yzzz_0[j] - 0.5 * fl1_fx * tpy_xyy_yzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyy_yzzz_0[j];

            tlx_xxyy_zzzz_0[j] = pa_x[j] * tlx_xyy_zzzz_0[j] + 0.5 * fl1_fx * tlx_yy_zzzz_0[j];

            tly_xxyy_zzzz_0[j] = pa_x[j] * tly_xyy_zzzz_0[j] + 0.5 * fl1_fx * tly_yy_zzzz_0[j] + 0.5 * fl1_fx * tpz_xyy_zzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyy_zzzz_0[j];

            tlz_xxyy_zzzz_0[j] = pa_x[j] * tlz_xyy_zzzz_0[j] + 0.5 * fl1_fx * tlz_yy_zzzz_0[j] - 0.5 * fl1_fx * tpy_xyy_zzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyy_zzzz_0[j];

            tlx_xxyz_xxxx_0[j] = pa_x[j] * tlx_xyz_xxxx_0[j] + 0.5 * fl1_fx * tlx_yz_xxxx_0[j] + 2.0 * fl1_fx * tlx_xyz_xxx_0[j];

            tly_xxyz_xxxx_0[j] = pa_x[j] * tly_xyz_xxxx_0[j] + 0.5 * fl1_fx * tly_yz_xxxx_0[j] + 2.0 * fl1_fx * tly_xyz_xxx_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxxx_0[j];

            tlz_xxyz_xxxx_0[j] = pa_x[j] * tlz_xyz_xxxx_0[j] + 0.5 * fl1_fx * tlz_yz_xxxx_0[j] + 2.0 * fl1_fx * tlz_xyz_xxx_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxxx_0[j];

            tlx_xxyz_xxxy_0[j] = pa_x[j] * tlx_xyz_xxxy_0[j] + 0.5 * fl1_fx * tlx_yz_xxxy_0[j] + 1.5 * fl1_fx * tlx_xyz_xxy_0[j];

            tly_xxyz_xxxy_0[j] = pa_x[j] * tly_xyz_xxxy_0[j] + 0.5 * fl1_fx * tly_yz_xxxy_0[j] + 1.5 * fl1_fx * tly_xyz_xxy_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxxy_0[j];

            tlz_xxyz_xxxy_0[j] = pa_x[j] * tlz_xyz_xxxy_0[j] + 0.5 * fl1_fx * tlz_yz_xxxy_0[j] + 1.5 * fl1_fx * tlz_xyz_xxy_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxxy_0[j];

            tlx_xxyz_xxxz_0[j] = pa_x[j] * tlx_xyz_xxxz_0[j] + 0.5 * fl1_fx * tlx_yz_xxxz_0[j] + 1.5 * fl1_fx * tlx_xyz_xxz_0[j];

            tly_xxyz_xxxz_0[j] = pa_x[j] * tly_xyz_xxxz_0[j] + 0.5 * fl1_fx * tly_yz_xxxz_0[j] + 1.5 * fl1_fx * tly_xyz_xxz_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxxz_0[j];

            tlz_xxyz_xxxz_0[j] = pa_x[j] * tlz_xyz_xxxz_0[j] + 0.5 * fl1_fx * tlz_yz_xxxz_0[j] + 1.5 * fl1_fx * tlz_xyz_xxz_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxxz_0[j];

            tlx_xxyz_xxyy_0[j] = pa_x[j] * tlx_xyz_xxyy_0[j] + 0.5 * fl1_fx * tlx_yz_xxyy_0[j] + fl1_fx * tlx_xyz_xyy_0[j];

            tly_xxyz_xxyy_0[j] = pa_x[j] * tly_xyz_xxyy_0[j] + 0.5 * fl1_fx * tly_yz_xxyy_0[j] + fl1_fx * tly_xyz_xyy_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxyy_0[j];

            tlz_xxyz_xxyy_0[j] = pa_x[j] * tlz_xyz_xxyy_0[j] + 0.5 * fl1_fx * tlz_yz_xxyy_0[j] + fl1_fx * tlz_xyz_xyy_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxyy_0[j];

            tlx_xxyz_xxyz_0[j] = pa_x[j] * tlx_xyz_xxyz_0[j] + 0.5 * fl1_fx * tlx_yz_xxyz_0[j] + fl1_fx * tlx_xyz_xyz_0[j];

            tly_xxyz_xxyz_0[j] = pa_x[j] * tly_xyz_xxyz_0[j] + 0.5 * fl1_fx * tly_yz_xxyz_0[j] + fl1_fx * tly_xyz_xyz_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxyz_0[j];

            tlz_xxyz_xxyz_0[j] = pa_x[j] * tlz_xyz_xxyz_0[j] + 0.5 * fl1_fx * tlz_yz_xxyz_0[j] + fl1_fx * tlz_xyz_xyz_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxyz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_195_243(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 65);

        auto tly_xyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tlz_xyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tlx_xyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 66);

        auto tly_xyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 66);

        auto tlz_xyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 66);

        auto tlx_xyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 67);

        auto tly_xyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 67);

        auto tlz_xyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 67);

        auto tlx_xyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 68);

        auto tly_xyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 68);

        auto tlz_xyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 68);

        auto tlx_xyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 69);

        auto tly_xyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 69);

        auto tlz_xyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 69);

        auto tlx_xyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 70);

        auto tly_xyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 70);

        auto tlz_xyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 70);

        auto tlx_xyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 71);

        auto tly_xyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 71);

        auto tlz_xyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 71);

        auto tlx_xyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 72);

        auto tly_xyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 72);

        auto tlz_xyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 72);

        auto tlx_xyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 73);

        auto tly_xyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 73);

        auto tlz_xyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 73);

        auto tlx_xyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 74);

        auto tly_xyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 74);

        auto tlz_xyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 74);

        auto tlx_xzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 75);

        auto tly_xzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 75);

        auto tlz_xzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 75);

        auto tlx_xzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 76);

        auto tly_xzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 76);

        auto tlz_xzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 76);

        auto tlx_xzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 77);

        auto tly_xzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 77);

        auto tlz_xzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 77);

        auto tlx_xzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 78);

        auto tly_xzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 78);

        auto tlz_xzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 78);

        auto tlx_xzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 79);

        auto tly_xzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 79);

        auto tlz_xzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 79);

        auto tlx_xzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 80);

        auto tly_xzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 80);

        auto tlz_xzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 80);

        auto tlx_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 65);

        auto tly_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tlz_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tlx_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 66);

        auto tly_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tlz_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tlx_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 67);

        auto tly_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tlz_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tlx_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 68);

        auto tly_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tlz_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tlx_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 69);

        auto tly_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tlz_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tlx_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 70);

        auto tly_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tlz_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tlx_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 71);

        auto tly_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tlz_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tlx_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 72);

        auto tly_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tlz_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tlx_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 73);

        auto tly_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tlz_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tlx_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 74);

        auto tly_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tlz_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75);

        auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76);

        auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77);

        auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78);

        auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79);

        auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80);

        auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tlx_xyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 45);

        auto tly_xyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 45);

        auto tlz_xyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 45);

        auto tlx_xyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 46);

        auto tly_xyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 46);

        auto tlz_xyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 46);

        auto tlx_xyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 47);

        auto tly_xyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 47);

        auto tlz_xyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 47);

        auto tlx_xyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 48);

        auto tly_xyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 48);

        auto tlz_xyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 48);

        auto tlx_xyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 49);

        auto tly_xyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 49);

        auto tlz_xyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 49);

        auto tlx_xzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 50);

        auto tly_xzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 50);

        auto tlz_xzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 50);

        auto tlx_xzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 51);

        auto tly_xzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 51);

        auto tlz_xzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 51);

        auto tlx_xzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 52);

        auto tly_xzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 52);

        auto tlz_xzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 52);

        auto tlx_xzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 53);

        auto tly_xzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 53);

        auto tlz_xzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 53);

        auto tlx_xzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 54);

        auto tly_xzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 54);

        auto tlz_xzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 54);

        auto tlx_xzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 55);

        auto tly_xzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 55);

        auto tlz_xzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 55);

        auto tpy_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tpz_xyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tpy_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 66);

        auto tpz_xyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 66);

        auto tpy_xyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 67);

        auto tpz_xyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 67);

        auto tpy_xyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 68);

        auto tpz_xyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 68);

        auto tpy_xyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 69);

        auto tpz_xyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 69);

        auto tpy_xyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 70);

        auto tpz_xyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 70);

        auto tpy_xyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 71);

        auto tpz_xyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 71);

        auto tpy_xyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 72);

        auto tpz_xyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 72);

        auto tpy_xyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 73);

        auto tpz_xyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 73);

        auto tpy_xyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 74);

        auto tpz_xyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 74);

        auto tpy_xzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 75);

        auto tpz_xzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 75);

        auto tpy_xzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 76);

        auto tpz_xzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 76);

        auto tpy_xzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 77);

        auto tpz_xzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 77);

        auto tpy_xzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 78);

        auto tpz_xzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 78);

        auto tpy_xzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 79);

        auto tpz_xzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 79);

        auto tpy_xzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 80);

        auto tpz_xzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 80);

        auto tdy_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 65);

        auto tdz_xyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 65);

        auto tdy_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 66);

        auto tdz_xyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 66);

        auto tdy_xyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 67);

        auto tdz_xyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 67);

        auto tdy_xyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 68);

        auto tdz_xyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 68);

        auto tdy_xyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 69);

        auto tdz_xyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 69);

        auto tdy_xyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 70);

        auto tdz_xyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 70);

        auto tdy_xyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 71);

        auto tdz_xyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 71);

        auto tdy_xyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 72);

        auto tdz_xyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 72);

        auto tdy_xyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 73);

        auto tdz_xyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 73);

        auto tdy_xyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 74);

        auto tdz_xyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 74);

        auto tdy_xzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 75);

        auto tdz_xzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 75);

        auto tdy_xzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 76);

        auto tdz_xzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 76);

        auto tdy_xzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 77);

        auto tdz_xzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 77);

        auto tdy_xzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 78);

        auto tdz_xzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 78);

        auto tdy_xzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 79);

        auto tdz_xzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 79);

        auto tdy_xzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 80);

        auto tdz_xzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 80);

        // set up pointers to integrals

        auto tlx_xxyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 65);

        auto tly_xxyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 65);

        auto tlz_xxyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 65);

        auto tlx_xxyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 66);

        auto tly_xxyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 66);

        auto tlz_xxyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 66);

        auto tlx_xxyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 67);

        auto tly_xxyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 67);

        auto tlz_xxyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 67);

        auto tlx_xxyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 68);

        auto tly_xxyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 68);

        auto tlz_xxyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 68);

        auto tlx_xxyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 69);

        auto tly_xxyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 69);

        auto tlz_xxyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 69);

        auto tlx_xxyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 70);

        auto tly_xxyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 70);

        auto tlz_xxyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 70);

        auto tlx_xxyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 71);

        auto tly_xxyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 71);

        auto tlz_xxyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 71);

        auto tlx_xxyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 72);

        auto tly_xxyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 72);

        auto tlz_xxyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 72);

        auto tlx_xxyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 73);

        auto tly_xxyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 73);

        auto tlz_xxyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 73);

        auto tlx_xxyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 74);

        auto tly_xxyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 74);

        auto tlz_xxyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 74);

        auto tlx_xxzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 75);

        auto tly_xxzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 75);

        auto tlz_xxzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 75);

        auto tlx_xxzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 76);

        auto tly_xxzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 76);

        auto tlz_xxzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 76);

        auto tlx_xxzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 77);

        auto tly_xxzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 77);

        auto tlz_xxzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 77);

        auto tlx_xxzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 78);

        auto tly_xxzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 78);

        auto tlz_xxzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 78);

        auto tlx_xxzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 79);

        auto tly_xxzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 79);

        auto tlz_xxzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 79);

        auto tlx_xxzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 80);

        auto tly_xxzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 80);

        auto tlz_xxzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 80);

        // Batch of Integrals (195,243)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xyz_xxzz_0, tdy_xyz_xyyy_0, tdy_xyz_xyyz_0, \
                                     tdy_xyz_xyzz_0, tdy_xyz_xzzz_0, tdy_xyz_yyyy_0, tdy_xyz_yyyz_0, tdy_xyz_yyzz_0, \
                                     tdy_xyz_yzzz_0, tdy_xyz_zzzz_0, tdy_xzz_xxxx_0, tdy_xzz_xxxy_0, tdy_xzz_xxxz_0, \
                                     tdy_xzz_xxyy_0, tdy_xzz_xxyz_0, tdy_xzz_xxzz_0, tdz_xyz_xxzz_0, tdz_xyz_xyyy_0, \
                                     tdz_xyz_xyyz_0, tdz_xyz_xyzz_0, tdz_xyz_xzzz_0, tdz_xyz_yyyy_0, tdz_xyz_yyyz_0, \
                                     tdz_xyz_yyzz_0, tdz_xyz_yzzz_0, tdz_xyz_zzzz_0, tdz_xzz_xxxx_0, tdz_xzz_xxxy_0, \
                                     tdz_xzz_xxxz_0, tdz_xzz_xxyy_0, tdz_xzz_xxyz_0, tdz_xzz_xxzz_0, tlx_xxyz_xxzz_0, \
                                     tlx_xxyz_xyyy_0, tlx_xxyz_xyyz_0, tlx_xxyz_xyzz_0, tlx_xxyz_xzzz_0, tlx_xxyz_yyyy_0, \
                                     tlx_xxyz_yyyz_0, tlx_xxyz_yyzz_0, tlx_xxyz_yzzz_0, tlx_xxyz_zzzz_0, tlx_xxzz_xxxx_0, \
                                     tlx_xxzz_xxxy_0, tlx_xxzz_xxxz_0, tlx_xxzz_xxyy_0, tlx_xxzz_xxyz_0, tlx_xxzz_xxzz_0, \
                                     tlx_xyz_xxzz_0, tlx_xyz_xyyy_0, tlx_xyz_xyyz_0, tlx_xyz_xyzz_0, tlx_xyz_xzz_0, \
                                     tlx_xyz_xzzz_0, tlx_xyz_yyy_0, tlx_xyz_yyyy_0, tlx_xyz_yyyz_0, tlx_xyz_yyz_0, \
                                     tlx_xyz_yyzz_0, tlx_xyz_yzz_0, tlx_xyz_yzzz_0, tlx_xyz_zzz_0, tlx_xyz_zzzz_0, \
                                     tlx_xzz_xxx_0, tlx_xzz_xxxx_0, tlx_xzz_xxxy_0, tlx_xzz_xxxz_0, tlx_xzz_xxy_0, \
                                     tlx_xzz_xxyy_0, tlx_xzz_xxyz_0, tlx_xzz_xxz_0, tlx_xzz_xxzz_0, tlx_xzz_xyy_0, \
                                     tlx_xzz_xyz_0, tlx_xzz_xzz_0, tlx_yz_xxzz_0, tlx_yz_xyyy_0, tlx_yz_xyyz_0, \
                                     tlx_yz_xyzz_0, tlx_yz_xzzz_0, tlx_yz_yyyy_0, tlx_yz_yyyz_0, tlx_yz_yyzz_0, \
                                     tlx_yz_yzzz_0, tlx_yz_zzzz_0, tlx_zz_xxxx_0, tlx_zz_xxxy_0, tlx_zz_xxxz_0, \
                                     tlx_zz_xxyy_0, tlx_zz_xxyz_0, tlx_zz_xxzz_0, tly_xxyz_xxzz_0, tly_xxyz_xyyy_0, \
                                     tly_xxyz_xyyz_0, tly_xxyz_xyzz_0, tly_xxyz_xzzz_0, tly_xxyz_yyyy_0, tly_xxyz_yyyz_0, \
                                     tly_xxyz_yyzz_0, tly_xxyz_yzzz_0, tly_xxyz_zzzz_0, tly_xxzz_xxxx_0, tly_xxzz_xxxy_0, \
                                     tly_xxzz_xxxz_0, tly_xxzz_xxyy_0, tly_xxzz_xxyz_0, tly_xxzz_xxzz_0, tly_xyz_xxzz_0, \
                                     tly_xyz_xyyy_0, tly_xyz_xyyz_0, tly_xyz_xyzz_0, tly_xyz_xzz_0, tly_xyz_xzzz_0, \
                                     tly_xyz_yyy_0, tly_xyz_yyyy_0, tly_xyz_yyyz_0, tly_xyz_yyz_0, tly_xyz_yyzz_0, \
                                     tly_xyz_yzz_0, tly_xyz_yzzz_0, tly_xyz_zzz_0, tly_xyz_zzzz_0, tly_xzz_xxx_0, \
                                     tly_xzz_xxxx_0, tly_xzz_xxxy_0, tly_xzz_xxxz_0, tly_xzz_xxy_0, tly_xzz_xxyy_0, \
                                     tly_xzz_xxyz_0, tly_xzz_xxz_0, tly_xzz_xxzz_0, tly_xzz_xyy_0, tly_xzz_xyz_0, \
                                     tly_xzz_xzz_0, tly_yz_xxzz_0, tly_yz_xyyy_0, tly_yz_xyyz_0, tly_yz_xyzz_0, \
                                     tly_yz_xzzz_0, tly_yz_yyyy_0, tly_yz_yyyz_0, tly_yz_yyzz_0, tly_yz_yzzz_0, \
                                     tly_yz_zzzz_0, tly_zz_xxxx_0, tly_zz_xxxy_0, tly_zz_xxxz_0, tly_zz_xxyy_0, \
                                     tly_zz_xxyz_0, tly_zz_xxzz_0, tlz_xxyz_xxzz_0, tlz_xxyz_xyyy_0, tlz_xxyz_xyyz_0, \
                                     tlz_xxyz_xyzz_0, tlz_xxyz_xzzz_0, tlz_xxyz_yyyy_0, tlz_xxyz_yyyz_0, tlz_xxyz_yyzz_0, \
                                     tlz_xxyz_yzzz_0, tlz_xxyz_zzzz_0, tlz_xxzz_xxxx_0, tlz_xxzz_xxxy_0, tlz_xxzz_xxxz_0, \
                                     tlz_xxzz_xxyy_0, tlz_xxzz_xxyz_0, tlz_xxzz_xxzz_0, tlz_xyz_xxzz_0, tlz_xyz_xyyy_0, \
                                     tlz_xyz_xyyz_0, tlz_xyz_xyzz_0, tlz_xyz_xzz_0, tlz_xyz_xzzz_0, tlz_xyz_yyy_0, \
                                     tlz_xyz_yyyy_0, tlz_xyz_yyyz_0, tlz_xyz_yyz_0, tlz_xyz_yyzz_0, tlz_xyz_yzz_0, \
                                     tlz_xyz_yzzz_0, tlz_xyz_zzz_0, tlz_xyz_zzzz_0, tlz_xzz_xxx_0, tlz_xzz_xxxx_0, \
                                     tlz_xzz_xxxy_0, tlz_xzz_xxxz_0, tlz_xzz_xxy_0, tlz_xzz_xxyy_0, tlz_xzz_xxyz_0, \
                                     tlz_xzz_xxz_0, tlz_xzz_xxzz_0, tlz_xzz_xyy_0, tlz_xzz_xyz_0, tlz_xzz_xzz_0, \
                                     tlz_yz_xxzz_0, tlz_yz_xyyy_0, tlz_yz_xyyz_0, tlz_yz_xyzz_0, tlz_yz_xzzz_0, \
                                     tlz_yz_yyyy_0, tlz_yz_yyyz_0, tlz_yz_yyzz_0, tlz_yz_yzzz_0, tlz_yz_zzzz_0, \
                                     tlz_zz_xxxx_0, tlz_zz_xxxy_0, tlz_zz_xxxz_0, tlz_zz_xxyy_0, tlz_zz_xxyz_0, \
                                     tlz_zz_xxzz_0, tpy_xyz_xxzz_0, tpy_xyz_xyyy_0, tpy_xyz_xyyz_0, tpy_xyz_xyzz_0, \
                                     tpy_xyz_xzzz_0, tpy_xyz_yyyy_0, tpy_xyz_yyyz_0, tpy_xyz_yyzz_0, tpy_xyz_yzzz_0, \
                                     tpy_xyz_zzzz_0, tpy_xzz_xxxx_0, tpy_xzz_xxxy_0, tpy_xzz_xxxz_0, tpy_xzz_xxyy_0, \
                                     tpy_xzz_xxyz_0, tpy_xzz_xxzz_0, tpz_xyz_xxzz_0, tpz_xyz_xyyy_0, tpz_xyz_xyyz_0, \
                                     tpz_xyz_xyzz_0, tpz_xyz_xzzz_0, tpz_xyz_yyyy_0, tpz_xyz_yyyz_0, tpz_xyz_yyzz_0, \
                                     tpz_xyz_yzzz_0, tpz_xyz_zzzz_0, tpz_xzz_xxxx_0, tpz_xzz_xxxy_0, tpz_xzz_xxxz_0, \
                                     tpz_xzz_xxyy_0, tpz_xzz_xxyz_0, tpz_xzz_xxzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxyz_xxzz_0[j] = pa_x[j] * tlx_xyz_xxzz_0[j] + 0.5 * fl1_fx * tlx_yz_xxzz_0[j] + fl1_fx * tlx_xyz_xzz_0[j];

            tly_xxyz_xxzz_0[j] = pa_x[j] * tly_xyz_xxzz_0[j] + 0.5 * fl1_fx * tly_yz_xxzz_0[j] + fl1_fx * tly_xyz_xzz_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xxzz_0[j];

            tlz_xxyz_xxzz_0[j] = pa_x[j] * tlz_xyz_xxzz_0[j] + 0.5 * fl1_fx * tlz_yz_xxzz_0[j] + fl1_fx * tlz_xyz_xzz_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xxzz_0[j];

            tlx_xxyz_xyyy_0[j] = pa_x[j] * tlx_xyz_xyyy_0[j] + 0.5 * fl1_fx * tlx_yz_xyyy_0[j] + 0.5 * fl1_fx * tlx_xyz_yyy_0[j];

            tly_xxyz_xyyy_0[j] = pa_x[j] * tly_xyz_xyyy_0[j] + 0.5 * fl1_fx * tly_yz_xyyy_0[j] + 0.5 * fl1_fx * tly_xyz_yyy_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xyyy_0[j];

            tlz_xxyz_xyyy_0[j] = pa_x[j] * tlz_xyz_xyyy_0[j] + 0.5 * fl1_fx * tlz_yz_xyyy_0[j] + 0.5 * fl1_fx * tlz_xyz_yyy_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xyyy_0[j];

            tlx_xxyz_xyyz_0[j] = pa_x[j] * tlx_xyz_xyyz_0[j] + 0.5 * fl1_fx * tlx_yz_xyyz_0[j] + 0.5 * fl1_fx * tlx_xyz_yyz_0[j];

            tly_xxyz_xyyz_0[j] = pa_x[j] * tly_xyz_xyyz_0[j] + 0.5 * fl1_fx * tly_yz_xyyz_0[j] + 0.5 * fl1_fx * tly_xyz_yyz_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xyyz_0[j];

            tlz_xxyz_xyyz_0[j] = pa_x[j] * tlz_xyz_xyyz_0[j] + 0.5 * fl1_fx * tlz_yz_xyyz_0[j] + 0.5 * fl1_fx * tlz_xyz_yyz_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xyyz_0[j];

            tlx_xxyz_xyzz_0[j] = pa_x[j] * tlx_xyz_xyzz_0[j] + 0.5 * fl1_fx * tlx_yz_xyzz_0[j] + 0.5 * fl1_fx * tlx_xyz_yzz_0[j];

            tly_xxyz_xyzz_0[j] = pa_x[j] * tly_xyz_xyzz_0[j] + 0.5 * fl1_fx * tly_yz_xyzz_0[j] + 0.5 * fl1_fx * tly_xyz_yzz_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xyzz_0[j];

            tlz_xxyz_xyzz_0[j] = pa_x[j] * tlz_xyz_xyzz_0[j] + 0.5 * fl1_fx * tlz_yz_xyzz_0[j] + 0.5 * fl1_fx * tlz_xyz_yzz_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xyzz_0[j];

            tlx_xxyz_xzzz_0[j] = pa_x[j] * tlx_xyz_xzzz_0[j] + 0.5 * fl1_fx * tlx_yz_xzzz_0[j] + 0.5 * fl1_fx * tlx_xyz_zzz_0[j];

            tly_xxyz_xzzz_0[j] = pa_x[j] * tly_xyz_xzzz_0[j] + 0.5 * fl1_fx * tly_yz_xzzz_0[j] + 0.5 * fl1_fx * tly_xyz_zzz_0[j] +
                                 0.5 * fl1_fx * tpz_xyz_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xyz_xzzz_0[j];

            tlz_xxyz_xzzz_0[j] = pa_x[j] * tlz_xyz_xzzz_0[j] + 0.5 * fl1_fx * tlz_yz_xzzz_0[j] + 0.5 * fl1_fx * tlz_xyz_zzz_0[j] -
                                 0.5 * fl1_fx * tpy_xyz_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xyz_xzzz_0[j];

            tlx_xxyz_yyyy_0[j] = pa_x[j] * tlx_xyz_yyyy_0[j] + 0.5 * fl1_fx * tlx_yz_yyyy_0[j];

            tly_xxyz_yyyy_0[j] = pa_x[j] * tly_xyz_yyyy_0[j] + 0.5 * fl1_fx * tly_yz_yyyy_0[j] + 0.5 * fl1_fx * tpz_xyz_yyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyz_yyyy_0[j];

            tlz_xxyz_yyyy_0[j] = pa_x[j] * tlz_xyz_yyyy_0[j] + 0.5 * fl1_fx * tlz_yz_yyyy_0[j] - 0.5 * fl1_fx * tpy_xyz_yyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyz_yyyy_0[j];

            tlx_xxyz_yyyz_0[j] = pa_x[j] * tlx_xyz_yyyz_0[j] + 0.5 * fl1_fx * tlx_yz_yyyz_0[j];

            tly_xxyz_yyyz_0[j] = pa_x[j] * tly_xyz_yyyz_0[j] + 0.5 * fl1_fx * tly_yz_yyyz_0[j] + 0.5 * fl1_fx * tpz_xyz_yyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyz_yyyz_0[j];

            tlz_xxyz_yyyz_0[j] = pa_x[j] * tlz_xyz_yyyz_0[j] + 0.5 * fl1_fx * tlz_yz_yyyz_0[j] - 0.5 * fl1_fx * tpy_xyz_yyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyz_yyyz_0[j];

            tlx_xxyz_yyzz_0[j] = pa_x[j] * tlx_xyz_yyzz_0[j] + 0.5 * fl1_fx * tlx_yz_yyzz_0[j];

            tly_xxyz_yyzz_0[j] = pa_x[j] * tly_xyz_yyzz_0[j] + 0.5 * fl1_fx * tly_yz_yyzz_0[j] + 0.5 * fl1_fx * tpz_xyz_yyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyz_yyzz_0[j];

            tlz_xxyz_yyzz_0[j] = pa_x[j] * tlz_xyz_yyzz_0[j] + 0.5 * fl1_fx * tlz_yz_yyzz_0[j] - 0.5 * fl1_fx * tpy_xyz_yyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyz_yyzz_0[j];

            tlx_xxyz_yzzz_0[j] = pa_x[j] * tlx_xyz_yzzz_0[j] + 0.5 * fl1_fx * tlx_yz_yzzz_0[j];

            tly_xxyz_yzzz_0[j] = pa_x[j] * tly_xyz_yzzz_0[j] + 0.5 * fl1_fx * tly_yz_yzzz_0[j] + 0.5 * fl1_fx * tpz_xyz_yzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyz_yzzz_0[j];

            tlz_xxyz_yzzz_0[j] = pa_x[j] * tlz_xyz_yzzz_0[j] + 0.5 * fl1_fx * tlz_yz_yzzz_0[j] - 0.5 * fl1_fx * tpy_xyz_yzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyz_yzzz_0[j];

            tlx_xxyz_zzzz_0[j] = pa_x[j] * tlx_xyz_zzzz_0[j] + 0.5 * fl1_fx * tlx_yz_zzzz_0[j];

            tly_xxyz_zzzz_0[j] = pa_x[j] * tly_xyz_zzzz_0[j] + 0.5 * fl1_fx * tly_yz_zzzz_0[j] + 0.5 * fl1_fx * tpz_xyz_zzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xyz_zzzz_0[j];

            tlz_xxyz_zzzz_0[j] = pa_x[j] * tlz_xyz_zzzz_0[j] + 0.5 * fl1_fx * tlz_yz_zzzz_0[j] - 0.5 * fl1_fx * tpy_xyz_zzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xyz_zzzz_0[j];

            tlx_xxzz_xxxx_0[j] = pa_x[j] * tlx_xzz_xxxx_0[j] + 0.5 * fl1_fx * tlx_zz_xxxx_0[j] + 2.0 * fl1_fx * tlx_xzz_xxx_0[j];

            tly_xxzz_xxxx_0[j] = pa_x[j] * tly_xzz_xxxx_0[j] + 0.5 * fl1_fx * tly_zz_xxxx_0[j] + 2.0 * fl1_fx * tly_xzz_xxx_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xxxx_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxxx_0[j];

            tlz_xxzz_xxxx_0[j] = pa_x[j] * tlz_xzz_xxxx_0[j] + 0.5 * fl1_fx * tlz_zz_xxxx_0[j] + 2.0 * fl1_fx * tlz_xzz_xxx_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xxxx_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxxx_0[j];

            tlx_xxzz_xxxy_0[j] = pa_x[j] * tlx_xzz_xxxy_0[j] + 0.5 * fl1_fx * tlx_zz_xxxy_0[j] + 1.5 * fl1_fx * tlx_xzz_xxy_0[j];

            tly_xxzz_xxxy_0[j] = pa_x[j] * tly_xzz_xxxy_0[j] + 0.5 * fl1_fx * tly_zz_xxxy_0[j] + 1.5 * fl1_fx * tly_xzz_xxy_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xxxy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxxy_0[j];

            tlz_xxzz_xxxy_0[j] = pa_x[j] * tlz_xzz_xxxy_0[j] + 0.5 * fl1_fx * tlz_zz_xxxy_0[j] + 1.5 * fl1_fx * tlz_xzz_xxy_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xxxy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxxy_0[j];

            tlx_xxzz_xxxz_0[j] = pa_x[j] * tlx_xzz_xxxz_0[j] + 0.5 * fl1_fx * tlx_zz_xxxz_0[j] + 1.5 * fl1_fx * tlx_xzz_xxz_0[j];

            tly_xxzz_xxxz_0[j] = pa_x[j] * tly_xzz_xxxz_0[j] + 0.5 * fl1_fx * tly_zz_xxxz_0[j] + 1.5 * fl1_fx * tly_xzz_xxz_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xxxz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxxz_0[j];

            tlz_xxzz_xxxz_0[j] = pa_x[j] * tlz_xzz_xxxz_0[j] + 0.5 * fl1_fx * tlz_zz_xxxz_0[j] + 1.5 * fl1_fx * tlz_xzz_xxz_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xxxz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxxz_0[j];

            tlx_xxzz_xxyy_0[j] = pa_x[j] * tlx_xzz_xxyy_0[j] + 0.5 * fl1_fx * tlx_zz_xxyy_0[j] + fl1_fx * tlx_xzz_xyy_0[j];

            tly_xxzz_xxyy_0[j] = pa_x[j] * tly_xzz_xxyy_0[j] + 0.5 * fl1_fx * tly_zz_xxyy_0[j] + fl1_fx * tly_xzz_xyy_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxyy_0[j];

            tlz_xxzz_xxyy_0[j] = pa_x[j] * tlz_xzz_xxyy_0[j] + 0.5 * fl1_fx * tlz_zz_xxyy_0[j] + fl1_fx * tlz_xzz_xyy_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxyy_0[j];

            tlx_xxzz_xxyz_0[j] = pa_x[j] * tlx_xzz_xxyz_0[j] + 0.5 * fl1_fx * tlx_zz_xxyz_0[j] + fl1_fx * tlx_xzz_xyz_0[j];

            tly_xxzz_xxyz_0[j] = pa_x[j] * tly_xzz_xxyz_0[j] + 0.5 * fl1_fx * tly_zz_xxyz_0[j] + fl1_fx * tly_xzz_xyz_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxyz_0[j];

            tlz_xxzz_xxyz_0[j] = pa_x[j] * tlz_xzz_xxyz_0[j] + 0.5 * fl1_fx * tlz_zz_xxyz_0[j] + fl1_fx * tlz_xzz_xyz_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxyz_0[j];

            tlx_xxzz_xxzz_0[j] = pa_x[j] * tlx_xzz_xxzz_0[j] + 0.5 * fl1_fx * tlx_zz_xxzz_0[j] + fl1_fx * tlx_xzz_xzz_0[j];

            tly_xxzz_xxzz_0[j] = pa_x[j] * tly_xzz_xxzz_0[j] + 0.5 * fl1_fx * tly_zz_xxzz_0[j] + fl1_fx * tly_xzz_xzz_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xxzz_0[j];

            tlz_xxzz_xxzz_0[j] = pa_x[j] * tlz_xzz_xxzz_0[j] + 0.5 * fl1_fx * tlz_zz_xxzz_0[j] + fl1_fx * tlz_xzz_xzz_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xxzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_243_291(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_xzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 81);

        auto tly_xzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tlz_xzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tlx_xzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 82);

        auto tly_xzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tlz_xzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tlx_xzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 83);

        auto tly_xzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 83);

        auto tlz_xzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 83);

        auto tlx_xzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 84);

        auto tly_xzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 84);

        auto tlz_xzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 84);

        auto tlx_xzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 85);

        auto tly_xzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 85);

        auto tlz_xzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 85);

        auto tlx_xzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 86);

        auto tly_xzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 86);

        auto tlz_xzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 86);

        auto tlx_xzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 87);

        auto tly_xzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 87);

        auto tlz_xzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 87);

        auto tlx_xzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 88);

        auto tly_xzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 88);

        auto tlz_xzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 88);

        auto tlx_xzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 89);

        auto tly_xzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 89);

        auto tlz_xzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 89);

        auto tlx_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 90);

        auto tly_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tlz_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tlx_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 91);

        auto tly_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tlz_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tlx_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 92);

        auto tly_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tlz_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tlx_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 93);

        auto tly_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tlz_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tlx_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 94);

        auto tly_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tlz_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tlx_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 95);

        auto tly_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tlz_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tlx_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 96);

        auto tly_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tlz_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81);

        auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82);

        auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83);

        auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84);

        auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85);

        auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86);

        auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87);

        auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88);

        auto tly_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tlz_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tlx_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 89);

        auto tly_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tlz_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tlx_xzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 56);

        auto tly_xzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 56);

        auto tlz_xzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 56);

        auto tlx_xzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 57);

        auto tly_xzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 57);

        auto tlz_xzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 57);

        auto tlx_xzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 58);

        auto tly_xzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 58);

        auto tlz_xzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 58);

        auto tlx_xzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 59);

        auto tly_xzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 59);

        auto tlz_xzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 59);

        auto tlx_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 60);

        auto tly_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tlz_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tlx_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 61);

        auto tly_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tlz_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tlx_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 62);

        auto tly_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tlz_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tlx_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 63);

        auto tly_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tlz_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tlx_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 64);

        auto tly_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tlz_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tlx_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 65);

        auto tly_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tlz_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tlx_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 66);

        auto tly_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tlz_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto tpy_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tpz_xzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tpy_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tpz_xzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tpy_xzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 83);

        auto tpz_xzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 83);

        auto tpy_xzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 84);

        auto tpz_xzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 84);

        auto tpy_xzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 85);

        auto tpz_xzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 85);

        auto tpy_xzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 86);

        auto tpz_xzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 86);

        auto tpy_xzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 87);

        auto tpz_xzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 87);

        auto tpy_xzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 88);

        auto tpz_xzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 88);

        auto tpy_xzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 89);

        auto tpz_xzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 89);

        auto tpy_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tpz_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tpy_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tpz_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tpy_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tpz_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tpy_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tpz_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tpy_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tpz_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tpy_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tpz_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tpy_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tpz_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tdy_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 81);

        auto tdz_xzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 81);

        auto tdy_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 82);

        auto tdz_xzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 82);

        auto tdy_xzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 83);

        auto tdz_xzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 83);

        auto tdy_xzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 84);

        auto tdz_xzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 84);

        auto tdy_xzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 85);

        auto tdz_xzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 85);

        auto tdy_xzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 86);

        auto tdz_xzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 86);

        auto tdy_xzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 87);

        auto tdz_xzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 87);

        auto tdy_xzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 88);

        auto tdz_xzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 88);

        auto tdy_xzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 89);

        auto tdz_xzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 89);

        auto tdy_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tdz_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tdy_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tdz_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tdy_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tdz_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tdy_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tdz_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tdy_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tdz_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tdy_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tdz_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tdy_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tdz_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 96);

        // set up pointers to integrals

        auto tlx_xxzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 81);

        auto tly_xxzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 81);

        auto tlz_xxzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 81);

        auto tlx_xxzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 82);

        auto tly_xxzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 82);

        auto tlz_xxzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 82);

        auto tlx_xxzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 83);

        auto tly_xxzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 83);

        auto tlz_xxzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 83);

        auto tlx_xxzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 84);

        auto tly_xxzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 84);

        auto tlz_xxzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 84);

        auto tlx_xxzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 85);

        auto tly_xxzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 85);

        auto tlz_xxzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 85);

        auto tlx_xxzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 86);

        auto tly_xxzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 86);

        auto tlz_xxzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 86);

        auto tlx_xxzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 87);

        auto tly_xxzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 87);

        auto tlz_xxzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 87);

        auto tlx_xxzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 88);

        auto tly_xxzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 88);

        auto tlz_xxzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 88);

        auto tlx_xxzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 89);

        auto tly_xxzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 89);

        auto tlz_xxzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 89);

        auto tlx_xyyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 90);

        auto tly_xyyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 90);

        auto tlz_xyyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 90);

        auto tlx_xyyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 91);

        auto tly_xyyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 91);

        auto tlz_xyyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 91);

        auto tlx_xyyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 92);

        auto tly_xyyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 92);

        auto tlz_xyyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 92);

        auto tlx_xyyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 93);

        auto tly_xyyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 93);

        auto tlz_xyyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 93);

        auto tlx_xyyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 94);

        auto tly_xyyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 94);

        auto tlz_xyyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 94);

        auto tlx_xyyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 95);

        auto tly_xyyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 95);

        auto tlz_xyyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 95);

        auto tlx_xyyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 96);

        auto tly_xyyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 96);

        auto tlz_xyyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 96);

        // Batch of Integrals (243,291)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_xzz_xyyy_0, tdy_xzz_xyyz_0, tdy_xzz_xyzz_0, \
                                     tdy_xzz_xzzz_0, tdy_xzz_yyyy_0, tdy_xzz_yyyz_0, tdy_xzz_yyzz_0, tdy_xzz_yzzz_0, \
                                     tdy_xzz_zzzz_0, tdy_yyy_xxxx_0, tdy_yyy_xxxy_0, tdy_yyy_xxxz_0, tdy_yyy_xxyy_0, \
                                     tdy_yyy_xxyz_0, tdy_yyy_xxzz_0, tdy_yyy_xyyy_0, tdz_xzz_xyyy_0, tdz_xzz_xyyz_0, \
                                     tdz_xzz_xyzz_0, tdz_xzz_xzzz_0, tdz_xzz_yyyy_0, tdz_xzz_yyyz_0, tdz_xzz_yyzz_0, \
                                     tdz_xzz_yzzz_0, tdz_xzz_zzzz_0, tdz_yyy_xxxx_0, tdz_yyy_xxxy_0, tdz_yyy_xxxz_0, \
                                     tdz_yyy_xxyy_0, tdz_yyy_xxyz_0, tdz_yyy_xxzz_0, tdz_yyy_xyyy_0, tlx_xxzz_xyyy_0, \
                                     tlx_xxzz_xyyz_0, tlx_xxzz_xyzz_0, tlx_xxzz_xzzz_0, tlx_xxzz_yyyy_0, tlx_xxzz_yyyz_0, \
                                     tlx_xxzz_yyzz_0, tlx_xxzz_yzzz_0, tlx_xxzz_zzzz_0, tlx_xyyy_xxxx_0, tlx_xyyy_xxxy_0, \
                                     tlx_xyyy_xxxz_0, tlx_xyyy_xxyy_0, tlx_xyyy_xxyz_0, tlx_xyyy_xxzz_0, tlx_xyyy_xyyy_0, \
                                     tlx_xzz_xyyy_0, tlx_xzz_xyyz_0, tlx_xzz_xyzz_0, tlx_xzz_xzzz_0, tlx_xzz_yyy_0, \
                                     tlx_xzz_yyyy_0, tlx_xzz_yyyz_0, tlx_xzz_yyz_0, tlx_xzz_yyzz_0, tlx_xzz_yzz_0, \
                                     tlx_xzz_yzzz_0, tlx_xzz_zzz_0, tlx_xzz_zzzz_0, tlx_yyy_xxx_0, tlx_yyy_xxxx_0, \
                                     tlx_yyy_xxxy_0, tlx_yyy_xxxz_0, tlx_yyy_xxy_0, tlx_yyy_xxyy_0, tlx_yyy_xxyz_0, \
                                     tlx_yyy_xxz_0, tlx_yyy_xxzz_0, tlx_yyy_xyy_0, tlx_yyy_xyyy_0, tlx_yyy_xyz_0, \
                                     tlx_yyy_xzz_0, tlx_yyy_yyy_0, tlx_zz_xyyy_0, tlx_zz_xyyz_0, tlx_zz_xyzz_0, \
                                     tlx_zz_xzzz_0, tlx_zz_yyyy_0, tlx_zz_yyyz_0, tlx_zz_yyzz_0, tlx_zz_yzzz_0, \
                                     tlx_zz_zzzz_0, tly_xxzz_xyyy_0, tly_xxzz_xyyz_0, tly_xxzz_xyzz_0, tly_xxzz_xzzz_0, \
                                     tly_xxzz_yyyy_0, tly_xxzz_yyyz_0, tly_xxzz_yyzz_0, tly_xxzz_yzzz_0, tly_xxzz_zzzz_0, \
                                     tly_xyyy_xxxx_0, tly_xyyy_xxxy_0, tly_xyyy_xxxz_0, tly_xyyy_xxyy_0, tly_xyyy_xxyz_0, \
                                     tly_xyyy_xxzz_0, tly_xyyy_xyyy_0, tly_xzz_xyyy_0, tly_xzz_xyyz_0, tly_xzz_xyzz_0, \
                                     tly_xzz_xzzz_0, tly_xzz_yyy_0, tly_xzz_yyyy_0, tly_xzz_yyyz_0, tly_xzz_yyz_0, \
                                     tly_xzz_yyzz_0, tly_xzz_yzz_0, tly_xzz_yzzz_0, tly_xzz_zzz_0, tly_xzz_zzzz_0, \
                                     tly_yyy_xxx_0, tly_yyy_xxxx_0, tly_yyy_xxxy_0, tly_yyy_xxxz_0, tly_yyy_xxy_0, \
                                     tly_yyy_xxyy_0, tly_yyy_xxyz_0, tly_yyy_xxz_0, tly_yyy_xxzz_0, tly_yyy_xyy_0, \
                                     tly_yyy_xyyy_0, tly_yyy_xyz_0, tly_yyy_xzz_0, tly_yyy_yyy_0, tly_zz_xyyy_0, \
                                     tly_zz_xyyz_0, tly_zz_xyzz_0, tly_zz_xzzz_0, tly_zz_yyyy_0, tly_zz_yyyz_0, \
                                     tly_zz_yyzz_0, tly_zz_yzzz_0, tly_zz_zzzz_0, tlz_xxzz_xyyy_0, tlz_xxzz_xyyz_0, \
                                     tlz_xxzz_xyzz_0, tlz_xxzz_xzzz_0, tlz_xxzz_yyyy_0, tlz_xxzz_yyyz_0, tlz_xxzz_yyzz_0, \
                                     tlz_xxzz_yzzz_0, tlz_xxzz_zzzz_0, tlz_xyyy_xxxx_0, tlz_xyyy_xxxy_0, tlz_xyyy_xxxz_0, \
                                     tlz_xyyy_xxyy_0, tlz_xyyy_xxyz_0, tlz_xyyy_xxzz_0, tlz_xyyy_xyyy_0, tlz_xzz_xyyy_0, \
                                     tlz_xzz_xyyz_0, tlz_xzz_xyzz_0, tlz_xzz_xzzz_0, tlz_xzz_yyy_0, tlz_xzz_yyyy_0, \
                                     tlz_xzz_yyyz_0, tlz_xzz_yyz_0, tlz_xzz_yyzz_0, tlz_xzz_yzz_0, tlz_xzz_yzzz_0, \
                                     tlz_xzz_zzz_0, tlz_xzz_zzzz_0, tlz_yyy_xxx_0, tlz_yyy_xxxx_0, tlz_yyy_xxxy_0, \
                                     tlz_yyy_xxxz_0, tlz_yyy_xxy_0, tlz_yyy_xxyy_0, tlz_yyy_xxyz_0, tlz_yyy_xxz_0, \
                                     tlz_yyy_xxzz_0, tlz_yyy_xyy_0, tlz_yyy_xyyy_0, tlz_yyy_xyz_0, tlz_yyy_xzz_0, \
                                     tlz_yyy_yyy_0, tlz_zz_xyyy_0, tlz_zz_xyyz_0, tlz_zz_xyzz_0, tlz_zz_xzzz_0, \
                                     tlz_zz_yyyy_0, tlz_zz_yyyz_0, tlz_zz_yyzz_0, tlz_zz_yzzz_0, tlz_zz_zzzz_0, \
                                     tpy_xzz_xyyy_0, tpy_xzz_xyyz_0, tpy_xzz_xyzz_0, tpy_xzz_xzzz_0, tpy_xzz_yyyy_0, \
                                     tpy_xzz_yyyz_0, tpy_xzz_yyzz_0, tpy_xzz_yzzz_0, tpy_xzz_zzzz_0, tpy_yyy_xxxx_0, \
                                     tpy_yyy_xxxy_0, tpy_yyy_xxxz_0, tpy_yyy_xxyy_0, tpy_yyy_xxyz_0, tpy_yyy_xxzz_0, \
                                     tpy_yyy_xyyy_0, tpz_xzz_xyyy_0, tpz_xzz_xyyz_0, tpz_xzz_xyzz_0, tpz_xzz_xzzz_0, \
                                     tpz_xzz_yyyy_0, tpz_xzz_yyyz_0, tpz_xzz_yyzz_0, tpz_xzz_yzzz_0, tpz_xzz_zzzz_0, \
                                     tpz_yyy_xxxx_0, tpz_yyy_xxxy_0, tpz_yyy_xxxz_0, tpz_yyy_xxyy_0, tpz_yyy_xxyz_0, \
                                     tpz_yyy_xxzz_0, tpz_yyy_xyyy_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxzz_xyyy_0[j] = pa_x[j] * tlx_xzz_xyyy_0[j] + 0.5 * fl1_fx * tlx_zz_xyyy_0[j] + 0.5 * fl1_fx * tlx_xzz_yyy_0[j];

            tly_xxzz_xyyy_0[j] = pa_x[j] * tly_xzz_xyyy_0[j] + 0.5 * fl1_fx * tly_zz_xyyy_0[j] + 0.5 * fl1_fx * tly_xzz_yyy_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xyyy_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xyyy_0[j];

            tlz_xxzz_xyyy_0[j] = pa_x[j] * tlz_xzz_xyyy_0[j] + 0.5 * fl1_fx * tlz_zz_xyyy_0[j] + 0.5 * fl1_fx * tlz_xzz_yyy_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xyyy_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xyyy_0[j];

            tlx_xxzz_xyyz_0[j] = pa_x[j] * tlx_xzz_xyyz_0[j] + 0.5 * fl1_fx * tlx_zz_xyyz_0[j] + 0.5 * fl1_fx * tlx_xzz_yyz_0[j];

            tly_xxzz_xyyz_0[j] = pa_x[j] * tly_xzz_xyyz_0[j] + 0.5 * fl1_fx * tly_zz_xyyz_0[j] + 0.5 * fl1_fx * tly_xzz_yyz_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xyyz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xyyz_0[j];

            tlz_xxzz_xyyz_0[j] = pa_x[j] * tlz_xzz_xyyz_0[j] + 0.5 * fl1_fx * tlz_zz_xyyz_0[j] + 0.5 * fl1_fx * tlz_xzz_yyz_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xyyz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xyyz_0[j];

            tlx_xxzz_xyzz_0[j] = pa_x[j] * tlx_xzz_xyzz_0[j] + 0.5 * fl1_fx * tlx_zz_xyzz_0[j] + 0.5 * fl1_fx * tlx_xzz_yzz_0[j];

            tly_xxzz_xyzz_0[j] = pa_x[j] * tly_xzz_xyzz_0[j] + 0.5 * fl1_fx * tly_zz_xyzz_0[j] + 0.5 * fl1_fx * tly_xzz_yzz_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xyzz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xyzz_0[j];

            tlz_xxzz_xyzz_0[j] = pa_x[j] * tlz_xzz_xyzz_0[j] + 0.5 * fl1_fx * tlz_zz_xyzz_0[j] + 0.5 * fl1_fx * tlz_xzz_yzz_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xyzz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xyzz_0[j];

            tlx_xxzz_xzzz_0[j] = pa_x[j] * tlx_xzz_xzzz_0[j] + 0.5 * fl1_fx * tlx_zz_xzzz_0[j] + 0.5 * fl1_fx * tlx_xzz_zzz_0[j];

            tly_xxzz_xzzz_0[j] = pa_x[j] * tly_xzz_xzzz_0[j] + 0.5 * fl1_fx * tly_zz_xzzz_0[j] + 0.5 * fl1_fx * tly_xzz_zzz_0[j] +
                                 0.5 * fl1_fx * tpz_xzz_xzzz_0[j] + fl1_fx * fl1_fgb * tdz_xzz_xzzz_0[j];

            tlz_xxzz_xzzz_0[j] = pa_x[j] * tlz_xzz_xzzz_0[j] + 0.5 * fl1_fx * tlz_zz_xzzz_0[j] + 0.5 * fl1_fx * tlz_xzz_zzz_0[j] -
                                 0.5 * fl1_fx * tpy_xzz_xzzz_0[j] - fl1_fx * fl1_fgb * tdy_xzz_xzzz_0[j];

            tlx_xxzz_yyyy_0[j] = pa_x[j] * tlx_xzz_yyyy_0[j] + 0.5 * fl1_fx * tlx_zz_yyyy_0[j];

            tly_xxzz_yyyy_0[j] = pa_x[j] * tly_xzz_yyyy_0[j] + 0.5 * fl1_fx * tly_zz_yyyy_0[j] + 0.5 * fl1_fx * tpz_xzz_yyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xzz_yyyy_0[j];

            tlz_xxzz_yyyy_0[j] = pa_x[j] * tlz_xzz_yyyy_0[j] + 0.5 * fl1_fx * tlz_zz_yyyy_0[j] - 0.5 * fl1_fx * tpy_xzz_yyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xzz_yyyy_0[j];

            tlx_xxzz_yyyz_0[j] = pa_x[j] * tlx_xzz_yyyz_0[j] + 0.5 * fl1_fx * tlx_zz_yyyz_0[j];

            tly_xxzz_yyyz_0[j] = pa_x[j] * tly_xzz_yyyz_0[j] + 0.5 * fl1_fx * tly_zz_yyyz_0[j] + 0.5 * fl1_fx * tpz_xzz_yyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xzz_yyyz_0[j];

            tlz_xxzz_yyyz_0[j] = pa_x[j] * tlz_xzz_yyyz_0[j] + 0.5 * fl1_fx * tlz_zz_yyyz_0[j] - 0.5 * fl1_fx * tpy_xzz_yyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xzz_yyyz_0[j];

            tlx_xxzz_yyzz_0[j] = pa_x[j] * tlx_xzz_yyzz_0[j] + 0.5 * fl1_fx * tlx_zz_yyzz_0[j];

            tly_xxzz_yyzz_0[j] = pa_x[j] * tly_xzz_yyzz_0[j] + 0.5 * fl1_fx * tly_zz_yyzz_0[j] + 0.5 * fl1_fx * tpz_xzz_yyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xzz_yyzz_0[j];

            tlz_xxzz_yyzz_0[j] = pa_x[j] * tlz_xzz_yyzz_0[j] + 0.5 * fl1_fx * tlz_zz_yyzz_0[j] - 0.5 * fl1_fx * tpy_xzz_yyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xzz_yyzz_0[j];

            tlx_xxzz_yzzz_0[j] = pa_x[j] * tlx_xzz_yzzz_0[j] + 0.5 * fl1_fx * tlx_zz_yzzz_0[j];

            tly_xxzz_yzzz_0[j] = pa_x[j] * tly_xzz_yzzz_0[j] + 0.5 * fl1_fx * tly_zz_yzzz_0[j] + 0.5 * fl1_fx * tpz_xzz_yzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xzz_yzzz_0[j];

            tlz_xxzz_yzzz_0[j] = pa_x[j] * tlz_xzz_yzzz_0[j] + 0.5 * fl1_fx * tlz_zz_yzzz_0[j] - 0.5 * fl1_fx * tpy_xzz_yzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xzz_yzzz_0[j];

            tlx_xxzz_zzzz_0[j] = pa_x[j] * tlx_xzz_zzzz_0[j] + 0.5 * fl1_fx * tlx_zz_zzzz_0[j];

            tly_xxzz_zzzz_0[j] = pa_x[j] * tly_xzz_zzzz_0[j] + 0.5 * fl1_fx * tly_zz_zzzz_0[j] + 0.5 * fl1_fx * tpz_xzz_zzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_xzz_zzzz_0[j];

            tlz_xxzz_zzzz_0[j] = pa_x[j] * tlz_xzz_zzzz_0[j] + 0.5 * fl1_fx * tlz_zz_zzzz_0[j] - 0.5 * fl1_fx * tpy_xzz_zzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_xzz_zzzz_0[j];

            tlx_xyyy_xxxx_0[j] = pa_x[j] * tlx_yyy_xxxx_0[j] + 2.0 * fl1_fx * tlx_yyy_xxx_0[j];

            tly_xyyy_xxxx_0[j] = pa_x[j] * tly_yyy_xxxx_0[j] + 2.0 * fl1_fx * tly_yyy_xxx_0[j] + 0.5 * fl1_fx * tpz_yyy_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xxxx_0[j];

            tlz_xyyy_xxxx_0[j] = pa_x[j] * tlz_yyy_xxxx_0[j] + 2.0 * fl1_fx * tlz_yyy_xxx_0[j] - 0.5 * fl1_fx * tpy_yyy_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xxxx_0[j];

            tlx_xyyy_xxxy_0[j] = pa_x[j] * tlx_yyy_xxxy_0[j] + 1.5 * fl1_fx * tlx_yyy_xxy_0[j];

            tly_xyyy_xxxy_0[j] = pa_x[j] * tly_yyy_xxxy_0[j] + 1.5 * fl1_fx * tly_yyy_xxy_0[j] + 0.5 * fl1_fx * tpz_yyy_xxxy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xxxy_0[j];

            tlz_xyyy_xxxy_0[j] = pa_x[j] * tlz_yyy_xxxy_0[j] + 1.5 * fl1_fx * tlz_yyy_xxy_0[j] - 0.5 * fl1_fx * tpy_yyy_xxxy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xxxy_0[j];

            tlx_xyyy_xxxz_0[j] = pa_x[j] * tlx_yyy_xxxz_0[j] + 1.5 * fl1_fx * tlx_yyy_xxz_0[j];

            tly_xyyy_xxxz_0[j] = pa_x[j] * tly_yyy_xxxz_0[j] + 1.5 * fl1_fx * tly_yyy_xxz_0[j] + 0.5 * fl1_fx * tpz_yyy_xxxz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xxxz_0[j];

            tlz_xyyy_xxxz_0[j] = pa_x[j] * tlz_yyy_xxxz_0[j] + 1.5 * fl1_fx * tlz_yyy_xxz_0[j] - 0.5 * fl1_fx * tpy_yyy_xxxz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xxxz_0[j];

            tlx_xyyy_xxyy_0[j] = pa_x[j] * tlx_yyy_xxyy_0[j] + fl1_fx * tlx_yyy_xyy_0[j];

            tly_xyyy_xxyy_0[j] =
                pa_x[j] * tly_yyy_xxyy_0[j] + fl1_fx * tly_yyy_xyy_0[j] + 0.5 * fl1_fx * tpz_yyy_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xxyy_0[j];

            tlz_xyyy_xxyy_0[j] =
                pa_x[j] * tlz_yyy_xxyy_0[j] + fl1_fx * tlz_yyy_xyy_0[j] - 0.5 * fl1_fx * tpy_yyy_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xxyy_0[j];

            tlx_xyyy_xxyz_0[j] = pa_x[j] * tlx_yyy_xxyz_0[j] + fl1_fx * tlx_yyy_xyz_0[j];

            tly_xyyy_xxyz_0[j] =
                pa_x[j] * tly_yyy_xxyz_0[j] + fl1_fx * tly_yyy_xyz_0[j] + 0.5 * fl1_fx * tpz_yyy_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xxyz_0[j];

            tlz_xyyy_xxyz_0[j] =
                pa_x[j] * tlz_yyy_xxyz_0[j] + fl1_fx * tlz_yyy_xyz_0[j] - 0.5 * fl1_fx * tpy_yyy_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xxyz_0[j];

            tlx_xyyy_xxzz_0[j] = pa_x[j] * tlx_yyy_xxzz_0[j] + fl1_fx * tlx_yyy_xzz_0[j];

            tly_xyyy_xxzz_0[j] =
                pa_x[j] * tly_yyy_xxzz_0[j] + fl1_fx * tly_yyy_xzz_0[j] + 0.5 * fl1_fx * tpz_yyy_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_xxzz_0[j];

            tlz_xyyy_xxzz_0[j] =
                pa_x[j] * tlz_yyy_xxzz_0[j] + fl1_fx * tlz_yyy_xzz_0[j] - 0.5 * fl1_fx * tpy_yyy_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_xxzz_0[j];

            tlx_xyyy_xyyy_0[j] = pa_x[j] * tlx_yyy_xyyy_0[j] + 0.5 * fl1_fx * tlx_yyy_yyy_0[j];

            tly_xyyy_xyyy_0[j] = pa_x[j] * tly_yyy_xyyy_0[j] + 0.5 * fl1_fx * tly_yyy_yyy_0[j] + 0.5 * fl1_fx * tpz_yyy_xyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xyyy_0[j];

            tlz_xyyy_xyyy_0[j] = pa_x[j] * tlz_yyy_xyyy_0[j] + 0.5 * fl1_fx * tlz_yyy_yyy_0[j] - 0.5 * fl1_fx * tpy_yyy_xyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xyyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_291_339(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 97);

        auto tly_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tlz_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tlx_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 98);

        auto tly_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tlz_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tlx_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 99);

        auto tly_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tlz_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tlx_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 100);

        auto tly_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tlz_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tlx_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 101);

        auto tly_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tlz_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tlx_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 102);

        auto tly_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tlz_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tlx_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 103);

        auto tly_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tlz_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tlx_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 104);

        auto tly_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tlz_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tlx_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 105);

        auto tly_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tlz_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tlx_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 106);

        auto tly_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tlz_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tlx_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 107);

        auto tly_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tlz_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tlx_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 108);

        auto tly_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tlz_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tlx_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 109);

        auto tly_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tlz_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tlx_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 110);

        auto tly_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tlz_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tlx_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 111);

        auto tly_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tlz_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tlx_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 112);

        auto tly_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tlz_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tlx_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 67);

        auto tly_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tlz_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tlx_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 68);

        auto tly_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tlz_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tlx_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 69);

        auto tly_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tlz_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tlx_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 70);

        auto tly_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tlz_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tlx_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 71);

        auto tly_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tlz_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tlx_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 72);

        auto tly_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tlz_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tlx_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 73);

        auto tly_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tlz_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tlx_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 74);

        auto tly_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tlz_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tlx_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 75);

        auto tly_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tlz_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tlx_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 76);

        auto tly_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tlz_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tlx_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 77);

        auto tly_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tlz_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tpy_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tpz_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tpy_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tpz_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tpy_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tpz_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tpy_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tpz_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tpy_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tpz_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tpy_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tpz_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tpy_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tpz_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tpy_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tpz_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tpy_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tpz_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tpy_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tpz_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tpy_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tpz_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tpy_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tpz_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tpy_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tpz_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tpy_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tpz_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tpy_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tpz_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tpy_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tpz_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tdy_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tdz_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tdy_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tdz_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tdy_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tdz_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tdy_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tdz_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tdy_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tdz_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tdy_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tdz_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tdy_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tdz_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tdy_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tdz_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tdy_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tdz_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tdy_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tdz_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tdy_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tdz_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tdy_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tdz_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tdy_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tdz_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tdy_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tdz_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tdy_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tdz_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tdy_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tdz_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 112);

        // set up pointers to integrals

        auto tlx_xyyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 97);

        auto tly_xyyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 97);

        auto tlz_xyyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 97);

        auto tlx_xyyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 98);

        auto tly_xyyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 98);

        auto tlz_xyyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 98);

        auto tlx_xyyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 99);

        auto tly_xyyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 99);

        auto tlz_xyyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 99);

        auto tlx_xyyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 100);

        auto tly_xyyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 100);

        auto tlz_xyyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 100);

        auto tlx_xyyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 101);

        auto tly_xyyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 101);

        auto tlz_xyyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 101);

        auto tlx_xyyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 102);

        auto tly_xyyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 102);

        auto tlz_xyyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 102);

        auto tlx_xyyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 103);

        auto tly_xyyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 103);

        auto tlz_xyyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 103);

        auto tlx_xyyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 104);

        auto tly_xyyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 104);

        auto tlz_xyyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 104);

        auto tlx_xyyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 105);

        auto tly_xyyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 105);

        auto tlz_xyyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 105);

        auto tlx_xyyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 106);

        auto tly_xyyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 106);

        auto tlz_xyyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 106);

        auto tlx_xyyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 107);

        auto tly_xyyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 107);

        auto tlz_xyyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 107);

        auto tlx_xyyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 108);

        auto tly_xyyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 108);

        auto tlz_xyyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 108);

        auto tlx_xyyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 109);

        auto tly_xyyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 109);

        auto tlz_xyyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 109);

        auto tlx_xyyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 110);

        auto tly_xyyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 110);

        auto tlz_xyyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 110);

        auto tlx_xyyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 111);

        auto tly_xyyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 111);

        auto tlz_xyyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 111);

        auto tlx_xyyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 112);

        auto tly_xyyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 112);

        auto tlz_xyyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 112);

        // Batch of Integrals (291,339)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yyy_xyyz_0, tdy_yyy_xyzz_0, tdy_yyy_xzzz_0, \
                                     tdy_yyy_yyyy_0, tdy_yyy_yyyz_0, tdy_yyy_yyzz_0, tdy_yyy_yzzz_0, tdy_yyy_zzzz_0, \
                                     tdy_yyz_xxxx_0, tdy_yyz_xxxy_0, tdy_yyz_xxxz_0, tdy_yyz_xxyy_0, tdy_yyz_xxyz_0, \
                                     tdy_yyz_xxzz_0, tdy_yyz_xyyy_0, tdy_yyz_xyyz_0, tdz_yyy_xyyz_0, tdz_yyy_xyzz_0, \
                                     tdz_yyy_xzzz_0, tdz_yyy_yyyy_0, tdz_yyy_yyyz_0, tdz_yyy_yyzz_0, tdz_yyy_yzzz_0, \
                                     tdz_yyy_zzzz_0, tdz_yyz_xxxx_0, tdz_yyz_xxxy_0, tdz_yyz_xxxz_0, tdz_yyz_xxyy_0, \
                                     tdz_yyz_xxyz_0, tdz_yyz_xxzz_0, tdz_yyz_xyyy_0, tdz_yyz_xyyz_0, tlx_xyyy_xyyz_0, \
                                     tlx_xyyy_xyzz_0, tlx_xyyy_xzzz_0, tlx_xyyy_yyyy_0, tlx_xyyy_yyyz_0, tlx_xyyy_yyzz_0, \
                                     tlx_xyyy_yzzz_0, tlx_xyyy_zzzz_0, tlx_xyyz_xxxx_0, tlx_xyyz_xxxy_0, tlx_xyyz_xxxz_0, \
                                     tlx_xyyz_xxyy_0, tlx_xyyz_xxyz_0, tlx_xyyz_xxzz_0, tlx_xyyz_xyyy_0, tlx_xyyz_xyyz_0, \
                                     tlx_yyy_xyyz_0, tlx_yyy_xyzz_0, tlx_yyy_xzzz_0, tlx_yyy_yyyy_0, tlx_yyy_yyyz_0, \
                                     tlx_yyy_yyz_0, tlx_yyy_yyzz_0, tlx_yyy_yzz_0, tlx_yyy_yzzz_0, tlx_yyy_zzz_0, \
                                     tlx_yyy_zzzz_0, tlx_yyz_xxx_0, tlx_yyz_xxxx_0, tlx_yyz_xxxy_0, tlx_yyz_xxxz_0, \
                                     tlx_yyz_xxy_0, tlx_yyz_xxyy_0, tlx_yyz_xxyz_0, tlx_yyz_xxz_0, tlx_yyz_xxzz_0, \
                                     tlx_yyz_xyy_0, tlx_yyz_xyyy_0, tlx_yyz_xyyz_0, tlx_yyz_xyz_0, tlx_yyz_xzz_0, \
                                     tlx_yyz_yyy_0, tlx_yyz_yyz_0, tly_xyyy_xyyz_0, tly_xyyy_xyzz_0, tly_xyyy_xzzz_0, \
                                     tly_xyyy_yyyy_0, tly_xyyy_yyyz_0, tly_xyyy_yyzz_0, tly_xyyy_yzzz_0, tly_xyyy_zzzz_0, \
                                     tly_xyyz_xxxx_0, tly_xyyz_xxxy_0, tly_xyyz_xxxz_0, tly_xyyz_xxyy_0, tly_xyyz_xxyz_0, \
                                     tly_xyyz_xxzz_0, tly_xyyz_xyyy_0, tly_xyyz_xyyz_0, tly_yyy_xyyz_0, tly_yyy_xyzz_0, \
                                     tly_yyy_xzzz_0, tly_yyy_yyyy_0, tly_yyy_yyyz_0, tly_yyy_yyz_0, tly_yyy_yyzz_0, \
                                     tly_yyy_yzz_0, tly_yyy_yzzz_0, tly_yyy_zzz_0, tly_yyy_zzzz_0, tly_yyz_xxx_0, \
                                     tly_yyz_xxxx_0, tly_yyz_xxxy_0, tly_yyz_xxxz_0, tly_yyz_xxy_0, tly_yyz_xxyy_0, \
                                     tly_yyz_xxyz_0, tly_yyz_xxz_0, tly_yyz_xxzz_0, tly_yyz_xyy_0, tly_yyz_xyyy_0, \
                                     tly_yyz_xyyz_0, tly_yyz_xyz_0, tly_yyz_xzz_0, tly_yyz_yyy_0, tly_yyz_yyz_0, \
                                     tlz_xyyy_xyyz_0, tlz_xyyy_xyzz_0, tlz_xyyy_xzzz_0, tlz_xyyy_yyyy_0, tlz_xyyy_yyyz_0, \
                                     tlz_xyyy_yyzz_0, tlz_xyyy_yzzz_0, tlz_xyyy_zzzz_0, tlz_xyyz_xxxx_0, tlz_xyyz_xxxy_0, \
                                     tlz_xyyz_xxxz_0, tlz_xyyz_xxyy_0, tlz_xyyz_xxyz_0, tlz_xyyz_xxzz_0, tlz_xyyz_xyyy_0, \
                                     tlz_xyyz_xyyz_0, tlz_yyy_xyyz_0, tlz_yyy_xyzz_0, tlz_yyy_xzzz_0, tlz_yyy_yyyy_0, \
                                     tlz_yyy_yyyz_0, tlz_yyy_yyz_0, tlz_yyy_yyzz_0, tlz_yyy_yzz_0, tlz_yyy_yzzz_0, \
                                     tlz_yyy_zzz_0, tlz_yyy_zzzz_0, tlz_yyz_xxx_0, tlz_yyz_xxxx_0, tlz_yyz_xxxy_0, \
                                     tlz_yyz_xxxz_0, tlz_yyz_xxy_0, tlz_yyz_xxyy_0, tlz_yyz_xxyz_0, tlz_yyz_xxz_0, \
                                     tlz_yyz_xxzz_0, tlz_yyz_xyy_0, tlz_yyz_xyyy_0, tlz_yyz_xyyz_0, tlz_yyz_xyz_0, \
                                     tlz_yyz_xzz_0, tlz_yyz_yyy_0, tlz_yyz_yyz_0, tpy_yyy_xyyz_0, tpy_yyy_xyzz_0, \
                                     tpy_yyy_xzzz_0, tpy_yyy_yyyy_0, tpy_yyy_yyyz_0, tpy_yyy_yyzz_0, tpy_yyy_yzzz_0, \
                                     tpy_yyy_zzzz_0, tpy_yyz_xxxx_0, tpy_yyz_xxxy_0, tpy_yyz_xxxz_0, tpy_yyz_xxyy_0, \
                                     tpy_yyz_xxyz_0, tpy_yyz_xxzz_0, tpy_yyz_xyyy_0, tpy_yyz_xyyz_0, tpz_yyy_xyyz_0, \
                                     tpz_yyy_xyzz_0, tpz_yyy_xzzz_0, tpz_yyy_yyyy_0, tpz_yyy_yyyz_0, tpz_yyy_yyzz_0, \
                                     tpz_yyy_yzzz_0, tpz_yyy_zzzz_0, tpz_yyz_xxxx_0, tpz_yyz_xxxy_0, tpz_yyz_xxxz_0, \
                                     tpz_yyz_xxyy_0, tpz_yyz_xxyz_0, tpz_yyz_xxzz_0, tpz_yyz_xyyy_0, tpz_yyz_xyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xyyy_xyyz_0[j] = pa_x[j] * tlx_yyy_xyyz_0[j] + 0.5 * fl1_fx * tlx_yyy_yyz_0[j];

            tly_xyyy_xyyz_0[j] = pa_x[j] * tly_yyy_xyyz_0[j] + 0.5 * fl1_fx * tly_yyy_yyz_0[j] + 0.5 * fl1_fx * tpz_yyy_xyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xyyz_0[j];

            tlz_xyyy_xyyz_0[j] = pa_x[j] * tlz_yyy_xyyz_0[j] + 0.5 * fl1_fx * tlz_yyy_yyz_0[j] - 0.5 * fl1_fx * tpy_yyy_xyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xyyz_0[j];

            tlx_xyyy_xyzz_0[j] = pa_x[j] * tlx_yyy_xyzz_0[j] + 0.5 * fl1_fx * tlx_yyy_yzz_0[j];

            tly_xyyy_xyzz_0[j] = pa_x[j] * tly_yyy_xyzz_0[j] + 0.5 * fl1_fx * tly_yyy_yzz_0[j] + 0.5 * fl1_fx * tpz_yyy_xyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xyzz_0[j];

            tlz_xyyy_xyzz_0[j] = pa_x[j] * tlz_yyy_xyzz_0[j] + 0.5 * fl1_fx * tlz_yyy_yzz_0[j] - 0.5 * fl1_fx * tpy_yyy_xyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xyzz_0[j];

            tlx_xyyy_xzzz_0[j] = pa_x[j] * tlx_yyy_xzzz_0[j] + 0.5 * fl1_fx * tlx_yyy_zzz_0[j];

            tly_xyyy_xzzz_0[j] = pa_x[j] * tly_yyy_xzzz_0[j] + 0.5 * fl1_fx * tly_yyy_zzz_0[j] + 0.5 * fl1_fx * tpz_yyy_xzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyy_xzzz_0[j];

            tlz_xyyy_xzzz_0[j] = pa_x[j] * tlz_yyy_xzzz_0[j] + 0.5 * fl1_fx * tlz_yyy_zzz_0[j] - 0.5 * fl1_fx * tpy_yyy_xzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyy_xzzz_0[j];

            tlx_xyyy_yyyy_0[j] = pa_x[j] * tlx_yyy_yyyy_0[j];

            tly_xyyy_yyyy_0[j] = pa_x[j] * tly_yyy_yyyy_0[j] + 0.5 * fl1_fx * tpz_yyy_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yyyy_0[j];

            tlz_xyyy_yyyy_0[j] = pa_x[j] * tlz_yyy_yyyy_0[j] - 0.5 * fl1_fx * tpy_yyy_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yyyy_0[j];

            tlx_xyyy_yyyz_0[j] = pa_x[j] * tlx_yyy_yyyz_0[j];

            tly_xyyy_yyyz_0[j] = pa_x[j] * tly_yyy_yyyz_0[j] + 0.5 * fl1_fx * tpz_yyy_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yyyz_0[j];

            tlz_xyyy_yyyz_0[j] = pa_x[j] * tlz_yyy_yyyz_0[j] - 0.5 * fl1_fx * tpy_yyy_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yyyz_0[j];

            tlx_xyyy_yyzz_0[j] = pa_x[j] * tlx_yyy_yyzz_0[j];

            tly_xyyy_yyzz_0[j] = pa_x[j] * tly_yyy_yyzz_0[j] + 0.5 * fl1_fx * tpz_yyy_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yyzz_0[j];

            tlz_xyyy_yyzz_0[j] = pa_x[j] * tlz_yyy_yyzz_0[j] - 0.5 * fl1_fx * tpy_yyy_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yyzz_0[j];

            tlx_xyyy_yzzz_0[j] = pa_x[j] * tlx_yyy_yzzz_0[j];

            tly_xyyy_yzzz_0[j] = pa_x[j] * tly_yyy_yzzz_0[j] + 0.5 * fl1_fx * tpz_yyy_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_yzzz_0[j];

            tlz_xyyy_yzzz_0[j] = pa_x[j] * tlz_yyy_yzzz_0[j] - 0.5 * fl1_fx * tpy_yyy_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_yzzz_0[j];

            tlx_xyyy_zzzz_0[j] = pa_x[j] * tlx_yyy_zzzz_0[j];

            tly_xyyy_zzzz_0[j] = pa_x[j] * tly_yyy_zzzz_0[j] + 0.5 * fl1_fx * tpz_yyy_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_yyy_zzzz_0[j];

            tlz_xyyy_zzzz_0[j] = pa_x[j] * tlz_yyy_zzzz_0[j] - 0.5 * fl1_fx * tpy_yyy_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_yyy_zzzz_0[j];

            tlx_xyyz_xxxx_0[j] = pa_x[j] * tlx_yyz_xxxx_0[j] + 2.0 * fl1_fx * tlx_yyz_xxx_0[j];

            tly_xyyz_xxxx_0[j] = pa_x[j] * tly_yyz_xxxx_0[j] + 2.0 * fl1_fx * tly_yyz_xxx_0[j] + 0.5 * fl1_fx * tpz_yyz_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xxxx_0[j];

            tlz_xyyz_xxxx_0[j] = pa_x[j] * tlz_yyz_xxxx_0[j] + 2.0 * fl1_fx * tlz_yyz_xxx_0[j] - 0.5 * fl1_fx * tpy_yyz_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xxxx_0[j];

            tlx_xyyz_xxxy_0[j] = pa_x[j] * tlx_yyz_xxxy_0[j] + 1.5 * fl1_fx * tlx_yyz_xxy_0[j];

            tly_xyyz_xxxy_0[j] = pa_x[j] * tly_yyz_xxxy_0[j] + 1.5 * fl1_fx * tly_yyz_xxy_0[j] + 0.5 * fl1_fx * tpz_yyz_xxxy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xxxy_0[j];

            tlz_xyyz_xxxy_0[j] = pa_x[j] * tlz_yyz_xxxy_0[j] + 1.5 * fl1_fx * tlz_yyz_xxy_0[j] - 0.5 * fl1_fx * tpy_yyz_xxxy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xxxy_0[j];

            tlx_xyyz_xxxz_0[j] = pa_x[j] * tlx_yyz_xxxz_0[j] + 1.5 * fl1_fx * tlx_yyz_xxz_0[j];

            tly_xyyz_xxxz_0[j] = pa_x[j] * tly_yyz_xxxz_0[j] + 1.5 * fl1_fx * tly_yyz_xxz_0[j] + 0.5 * fl1_fx * tpz_yyz_xxxz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xxxz_0[j];

            tlz_xyyz_xxxz_0[j] = pa_x[j] * tlz_yyz_xxxz_0[j] + 1.5 * fl1_fx * tlz_yyz_xxz_0[j] - 0.5 * fl1_fx * tpy_yyz_xxxz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xxxz_0[j];

            tlx_xyyz_xxyy_0[j] = pa_x[j] * tlx_yyz_xxyy_0[j] + fl1_fx * tlx_yyz_xyy_0[j];

            tly_xyyz_xxyy_0[j] =
                pa_x[j] * tly_yyz_xxyy_0[j] + fl1_fx * tly_yyz_xyy_0[j] + 0.5 * fl1_fx * tpz_yyz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xxyy_0[j];

            tlz_xyyz_xxyy_0[j] =
                pa_x[j] * tlz_yyz_xxyy_0[j] + fl1_fx * tlz_yyz_xyy_0[j] - 0.5 * fl1_fx * tpy_yyz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xxyy_0[j];

            tlx_xyyz_xxyz_0[j] = pa_x[j] * tlx_yyz_xxyz_0[j] + fl1_fx * tlx_yyz_xyz_0[j];

            tly_xyyz_xxyz_0[j] =
                pa_x[j] * tly_yyz_xxyz_0[j] + fl1_fx * tly_yyz_xyz_0[j] + 0.5 * fl1_fx * tpz_yyz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xxyz_0[j];

            tlz_xyyz_xxyz_0[j] =
                pa_x[j] * tlz_yyz_xxyz_0[j] + fl1_fx * tlz_yyz_xyz_0[j] - 0.5 * fl1_fx * tpy_yyz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xxyz_0[j];

            tlx_xyyz_xxzz_0[j] = pa_x[j] * tlx_yyz_xxzz_0[j] + fl1_fx * tlx_yyz_xzz_0[j];

            tly_xyyz_xxzz_0[j] =
                pa_x[j] * tly_yyz_xxzz_0[j] + fl1_fx * tly_yyz_xzz_0[j] + 0.5 * fl1_fx * tpz_yyz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_xxzz_0[j];

            tlz_xyyz_xxzz_0[j] =
                pa_x[j] * tlz_yyz_xxzz_0[j] + fl1_fx * tlz_yyz_xzz_0[j] - 0.5 * fl1_fx * tpy_yyz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_xxzz_0[j];

            tlx_xyyz_xyyy_0[j] = pa_x[j] * tlx_yyz_xyyy_0[j] + 0.5 * fl1_fx * tlx_yyz_yyy_0[j];

            tly_xyyz_xyyy_0[j] = pa_x[j] * tly_yyz_xyyy_0[j] + 0.5 * fl1_fx * tly_yyz_yyy_0[j] + 0.5 * fl1_fx * tpz_yyz_xyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xyyy_0[j];

            tlz_xyyz_xyyy_0[j] = pa_x[j] * tlz_yyz_xyyy_0[j] + 0.5 * fl1_fx * tlz_yyz_yyy_0[j] - 0.5 * fl1_fx * tpy_yyz_xyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xyyy_0[j];

            tlx_xyyz_xyyz_0[j] = pa_x[j] * tlx_yyz_xyyz_0[j] + 0.5 * fl1_fx * tlx_yyz_yyz_0[j];

            tly_xyyz_xyyz_0[j] = pa_x[j] * tly_yyz_xyyz_0[j] + 0.5 * fl1_fx * tly_yyz_yyz_0[j] + 0.5 * fl1_fx * tpz_yyz_xyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xyyz_0[j];

            tlz_xyyz_xyyz_0[j] = pa_x[j] * tlz_yyz_xyyz_0[j] + 0.5 * fl1_fx * tlz_yyz_yyz_0[j] - 0.5 * fl1_fx * tpy_yyz_xyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xyyz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_339_387(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 113);

        auto tly_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tlz_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tlx_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 114);

        auto tly_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tlz_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tlx_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 115);

        auto tly_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tlz_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tlx_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 116);

        auto tly_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 116);

        auto tlz_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tlx_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 117);

        auto tly_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tlz_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tlx_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 118);

        auto tly_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tlz_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tlx_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 119);

        auto tly_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tlz_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tlx_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 120);

        auto tly_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tlz_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tlx_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 121);

        auto tly_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tlz_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tlx_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 122);

        auto tly_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tlz_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tlx_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 123);

        auto tly_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tlz_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tlx_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 124);

        auto tly_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tlz_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tlx_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 125);

        auto tly_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tlz_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tlx_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 126);

        auto tly_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tlz_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tlx_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 127);

        auto tly_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tlz_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tlx_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 128);

        auto tly_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tlz_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tlx_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 78);

        auto tly_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tlz_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tlx_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 79);

        auto tly_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tlz_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tlx_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 80);

        auto tly_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tlz_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tlx_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 81);

        auto tly_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tlz_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tlx_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 82);

        auto tly_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tlz_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tlx_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 83);

        auto tly_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tlz_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tlx_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 84);

        auto tly_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tlz_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tlx_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 85);

        auto tly_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tlz_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tlx_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 86);

        auto tly_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tlz_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tlx_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 87);

        auto tly_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tlz_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tlx_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 88);

        auto tly_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tlz_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tpy_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tpz_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tpy_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tpz_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tpy_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tpz_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tpy_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 116);

        auto tpz_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tpy_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tpz_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tpy_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tpz_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tpy_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tpz_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tpy_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tpz_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tpy_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tpz_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tpy_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tpz_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tpy_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tpz_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tpy_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tpz_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tpy_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tpz_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tpy_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tpz_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tpy_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tpz_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tpy_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tpz_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tdy_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tdz_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tdy_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tdz_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tdy_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tdz_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tdy_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 116);

        auto tdz_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tdy_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tdz_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tdy_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tdz_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tdy_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tdz_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tdy_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tdz_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tdy_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tdz_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tdy_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tdz_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tdy_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tdz_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tdy_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tdz_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tdy_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tdz_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tdy_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tdz_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tdy_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tdz_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tdy_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tdz_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 128);

        // set up pointers to integrals

        auto tlx_xyyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 113);

        auto tly_xyyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 113);

        auto tlz_xyyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 113);

        auto tlx_xyyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 114);

        auto tly_xyyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 114);

        auto tlz_xyyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 114);

        auto tlx_xyyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 115);

        auto tly_xyyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 115);

        auto tlz_xyyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 115);

        auto tlx_xyyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 116);

        auto tly_xyyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 116);

        auto tlz_xyyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 116);

        auto tlx_xyyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 117);

        auto tly_xyyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 117);

        auto tlz_xyyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 117);

        auto tlx_xyyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 118);

        auto tly_xyyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 118);

        auto tlz_xyyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 118);

        auto tlx_xyyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 119);

        auto tly_xyyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 119);

        auto tlz_xyyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 119);

        auto tlx_xyzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 120);

        auto tly_xyzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 120);

        auto tlz_xyzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 120);

        auto tlx_xyzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 121);

        auto tly_xyzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 121);

        auto tlz_xyzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 121);

        auto tlx_xyzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 122);

        auto tly_xyzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 122);

        auto tlz_xyzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 122);

        auto tlx_xyzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 123);

        auto tly_xyzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 123);

        auto tlz_xyzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 123);

        auto tlx_xyzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 124);

        auto tly_xyzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 124);

        auto tlz_xyzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 124);

        auto tlx_xyzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 125);

        auto tly_xyzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 125);

        auto tlz_xyzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 125);

        auto tlx_xyzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 126);

        auto tly_xyzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 126);

        auto tlz_xyzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 126);

        auto tlx_xyzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 127);

        auto tly_xyzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 127);

        auto tlz_xyzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 127);

        auto tlx_xyzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 128);

        auto tly_xyzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 128);

        auto tlz_xyzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 128);

        // Batch of Integrals (339,387)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yyz_xyzz_0, tdy_yyz_xzzz_0, tdy_yyz_yyyy_0, \
                                     tdy_yyz_yyyz_0, tdy_yyz_yyzz_0, tdy_yyz_yzzz_0, tdy_yyz_zzzz_0, tdy_yzz_xxxx_0, \
                                     tdy_yzz_xxxy_0, tdy_yzz_xxxz_0, tdy_yzz_xxyy_0, tdy_yzz_xxyz_0, tdy_yzz_xxzz_0, \
                                     tdy_yzz_xyyy_0, tdy_yzz_xyyz_0, tdy_yzz_xyzz_0, tdz_yyz_xyzz_0, tdz_yyz_xzzz_0, \
                                     tdz_yyz_yyyy_0, tdz_yyz_yyyz_0, tdz_yyz_yyzz_0, tdz_yyz_yzzz_0, tdz_yyz_zzzz_0, \
                                     tdz_yzz_xxxx_0, tdz_yzz_xxxy_0, tdz_yzz_xxxz_0, tdz_yzz_xxyy_0, tdz_yzz_xxyz_0, \
                                     tdz_yzz_xxzz_0, tdz_yzz_xyyy_0, tdz_yzz_xyyz_0, tdz_yzz_xyzz_0, tlx_xyyz_xyzz_0, \
                                     tlx_xyyz_xzzz_0, tlx_xyyz_yyyy_0, tlx_xyyz_yyyz_0, tlx_xyyz_yyzz_0, tlx_xyyz_yzzz_0, \
                                     tlx_xyyz_zzzz_0, tlx_xyzz_xxxx_0, tlx_xyzz_xxxy_0, tlx_xyzz_xxxz_0, tlx_xyzz_xxyy_0, \
                                     tlx_xyzz_xxyz_0, tlx_xyzz_xxzz_0, tlx_xyzz_xyyy_0, tlx_xyzz_xyyz_0, tlx_xyzz_xyzz_0, \
                                     tlx_yyz_xyzz_0, tlx_yyz_xzzz_0, tlx_yyz_yyyy_0, tlx_yyz_yyyz_0, tlx_yyz_yyzz_0, \
                                     tlx_yyz_yzz_0, tlx_yyz_yzzz_0, tlx_yyz_zzz_0, tlx_yyz_zzzz_0, tlx_yzz_xxx_0, \
                                     tlx_yzz_xxxx_0, tlx_yzz_xxxy_0, tlx_yzz_xxxz_0, tlx_yzz_xxy_0, tlx_yzz_xxyy_0, \
                                     tlx_yzz_xxyz_0, tlx_yzz_xxz_0, tlx_yzz_xxzz_0, tlx_yzz_xyy_0, tlx_yzz_xyyy_0, \
                                     tlx_yzz_xyyz_0, tlx_yzz_xyz_0, tlx_yzz_xyzz_0, tlx_yzz_xzz_0, tlx_yzz_yyy_0, \
                                     tlx_yzz_yyz_0, tlx_yzz_yzz_0, tly_xyyz_xyzz_0, tly_xyyz_xzzz_0, tly_xyyz_yyyy_0, \
                                     tly_xyyz_yyyz_0, tly_xyyz_yyzz_0, tly_xyyz_yzzz_0, tly_xyyz_zzzz_0, tly_xyzz_xxxx_0, \
                                     tly_xyzz_xxxy_0, tly_xyzz_xxxz_0, tly_xyzz_xxyy_0, tly_xyzz_xxyz_0, tly_xyzz_xxzz_0, \
                                     tly_xyzz_xyyy_0, tly_xyzz_xyyz_0, tly_xyzz_xyzz_0, tly_yyz_xyzz_0, tly_yyz_xzzz_0, \
                                     tly_yyz_yyyy_0, tly_yyz_yyyz_0, tly_yyz_yyzz_0, tly_yyz_yzz_0, tly_yyz_yzzz_0, \
                                     tly_yyz_zzz_0, tly_yyz_zzzz_0, tly_yzz_xxx_0, tly_yzz_xxxx_0, tly_yzz_xxxy_0, \
                                     tly_yzz_xxxz_0, tly_yzz_xxy_0, tly_yzz_xxyy_0, tly_yzz_xxyz_0, tly_yzz_xxz_0, \
                                     tly_yzz_xxzz_0, tly_yzz_xyy_0, tly_yzz_xyyy_0, tly_yzz_xyyz_0, tly_yzz_xyz_0, \
                                     tly_yzz_xyzz_0, tly_yzz_xzz_0, tly_yzz_yyy_0, tly_yzz_yyz_0, tly_yzz_yzz_0, \
                                     tlz_xyyz_xyzz_0, tlz_xyyz_xzzz_0, tlz_xyyz_yyyy_0, tlz_xyyz_yyyz_0, tlz_xyyz_yyzz_0, \
                                     tlz_xyyz_yzzz_0, tlz_xyyz_zzzz_0, tlz_xyzz_xxxx_0, tlz_xyzz_xxxy_0, tlz_xyzz_xxxz_0, \
                                     tlz_xyzz_xxyy_0, tlz_xyzz_xxyz_0, tlz_xyzz_xxzz_0, tlz_xyzz_xyyy_0, tlz_xyzz_xyyz_0, \
                                     tlz_xyzz_xyzz_0, tlz_yyz_xyzz_0, tlz_yyz_xzzz_0, tlz_yyz_yyyy_0, tlz_yyz_yyyz_0, \
                                     tlz_yyz_yyzz_0, tlz_yyz_yzz_0, tlz_yyz_yzzz_0, tlz_yyz_zzz_0, tlz_yyz_zzzz_0, \
                                     tlz_yzz_xxx_0, tlz_yzz_xxxx_0, tlz_yzz_xxxy_0, tlz_yzz_xxxz_0, tlz_yzz_xxy_0, \
                                     tlz_yzz_xxyy_0, tlz_yzz_xxyz_0, tlz_yzz_xxz_0, tlz_yzz_xxzz_0, tlz_yzz_xyy_0, \
                                     tlz_yzz_xyyy_0, tlz_yzz_xyyz_0, tlz_yzz_xyz_0, tlz_yzz_xyzz_0, tlz_yzz_xzz_0, \
                                     tlz_yzz_yyy_0, tlz_yzz_yyz_0, tlz_yzz_yzz_0, tpy_yyz_xyzz_0, tpy_yyz_xzzz_0, \
                                     tpy_yyz_yyyy_0, tpy_yyz_yyyz_0, tpy_yyz_yyzz_0, tpy_yyz_yzzz_0, tpy_yyz_zzzz_0, \
                                     tpy_yzz_xxxx_0, tpy_yzz_xxxy_0, tpy_yzz_xxxz_0, tpy_yzz_xxyy_0, tpy_yzz_xxyz_0, \
                                     tpy_yzz_xxzz_0, tpy_yzz_xyyy_0, tpy_yzz_xyyz_0, tpy_yzz_xyzz_0, tpz_yyz_xyzz_0, \
                                     tpz_yyz_xzzz_0, tpz_yyz_yyyy_0, tpz_yyz_yyyz_0, tpz_yyz_yyzz_0, tpz_yyz_yzzz_0, \
                                     tpz_yyz_zzzz_0, tpz_yzz_xxxx_0, tpz_yzz_xxxy_0, tpz_yzz_xxxz_0, tpz_yzz_xxyy_0, \
                                     tpz_yzz_xxyz_0, tpz_yzz_xxzz_0, tpz_yzz_xyyy_0, tpz_yzz_xyyz_0, tpz_yzz_xyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xyyz_xyzz_0[j] = pa_x[j] * tlx_yyz_xyzz_0[j] + 0.5 * fl1_fx * tlx_yyz_yzz_0[j];

            tly_xyyz_xyzz_0[j] = pa_x[j] * tly_yyz_xyzz_0[j] + 0.5 * fl1_fx * tly_yyz_yzz_0[j] + 0.5 * fl1_fx * tpz_yyz_xyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xyzz_0[j];

            tlz_xyyz_xyzz_0[j] = pa_x[j] * tlz_yyz_xyzz_0[j] + 0.5 * fl1_fx * tlz_yyz_yzz_0[j] - 0.5 * fl1_fx * tpy_yyz_xyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xyzz_0[j];

            tlx_xyyz_xzzz_0[j] = pa_x[j] * tlx_yyz_xzzz_0[j] + 0.5 * fl1_fx * tlx_yyz_zzz_0[j];

            tly_xyyz_xzzz_0[j] = pa_x[j] * tly_yyz_xzzz_0[j] + 0.5 * fl1_fx * tly_yyz_zzz_0[j] + 0.5 * fl1_fx * tpz_yyz_xzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yyz_xzzz_0[j];

            tlz_xyyz_xzzz_0[j] = pa_x[j] * tlz_yyz_xzzz_0[j] + 0.5 * fl1_fx * tlz_yyz_zzz_0[j] - 0.5 * fl1_fx * tpy_yyz_xzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yyz_xzzz_0[j];

            tlx_xyyz_yyyy_0[j] = pa_x[j] * tlx_yyz_yyyy_0[j];

            tly_xyyz_yyyy_0[j] = pa_x[j] * tly_yyz_yyyy_0[j] + 0.5 * fl1_fx * tpz_yyz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yyyy_0[j];

            tlz_xyyz_yyyy_0[j] = pa_x[j] * tlz_yyz_yyyy_0[j] - 0.5 * fl1_fx * tpy_yyz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yyyy_0[j];

            tlx_xyyz_yyyz_0[j] = pa_x[j] * tlx_yyz_yyyz_0[j];

            tly_xyyz_yyyz_0[j] = pa_x[j] * tly_yyz_yyyz_0[j] + 0.5 * fl1_fx * tpz_yyz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yyyz_0[j];

            tlz_xyyz_yyyz_0[j] = pa_x[j] * tlz_yyz_yyyz_0[j] - 0.5 * fl1_fx * tpy_yyz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yyyz_0[j];

            tlx_xyyz_yyzz_0[j] = pa_x[j] * tlx_yyz_yyzz_0[j];

            tly_xyyz_yyzz_0[j] = pa_x[j] * tly_yyz_yyzz_0[j] + 0.5 * fl1_fx * tpz_yyz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yyzz_0[j];

            tlz_xyyz_yyzz_0[j] = pa_x[j] * tlz_yyz_yyzz_0[j] - 0.5 * fl1_fx * tpy_yyz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yyzz_0[j];

            tlx_xyyz_yzzz_0[j] = pa_x[j] * tlx_yyz_yzzz_0[j];

            tly_xyyz_yzzz_0[j] = pa_x[j] * tly_yyz_yzzz_0[j] + 0.5 * fl1_fx * tpz_yyz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_yzzz_0[j];

            tlz_xyyz_yzzz_0[j] = pa_x[j] * tlz_yyz_yzzz_0[j] - 0.5 * fl1_fx * tpy_yyz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_yzzz_0[j];

            tlx_xyyz_zzzz_0[j] = pa_x[j] * tlx_yyz_zzzz_0[j];

            tly_xyyz_zzzz_0[j] = pa_x[j] * tly_yyz_zzzz_0[j] + 0.5 * fl1_fx * tpz_yyz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_yyz_zzzz_0[j];

            tlz_xyyz_zzzz_0[j] = pa_x[j] * tlz_yyz_zzzz_0[j] - 0.5 * fl1_fx * tpy_yyz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_yyz_zzzz_0[j];

            tlx_xyzz_xxxx_0[j] = pa_x[j] * tlx_yzz_xxxx_0[j] + 2.0 * fl1_fx * tlx_yzz_xxx_0[j];

            tly_xyzz_xxxx_0[j] = pa_x[j] * tly_yzz_xxxx_0[j] + 2.0 * fl1_fx * tly_yzz_xxx_0[j] + 0.5 * fl1_fx * tpz_yzz_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xxxx_0[j];

            tlz_xyzz_xxxx_0[j] = pa_x[j] * tlz_yzz_xxxx_0[j] + 2.0 * fl1_fx * tlz_yzz_xxx_0[j] - 0.5 * fl1_fx * tpy_yzz_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xxxx_0[j];

            tlx_xyzz_xxxy_0[j] = pa_x[j] * tlx_yzz_xxxy_0[j] + 1.5 * fl1_fx * tlx_yzz_xxy_0[j];

            tly_xyzz_xxxy_0[j] = pa_x[j] * tly_yzz_xxxy_0[j] + 1.5 * fl1_fx * tly_yzz_xxy_0[j] + 0.5 * fl1_fx * tpz_yzz_xxxy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xxxy_0[j];

            tlz_xyzz_xxxy_0[j] = pa_x[j] * tlz_yzz_xxxy_0[j] + 1.5 * fl1_fx * tlz_yzz_xxy_0[j] - 0.5 * fl1_fx * tpy_yzz_xxxy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xxxy_0[j];

            tlx_xyzz_xxxz_0[j] = pa_x[j] * tlx_yzz_xxxz_0[j] + 1.5 * fl1_fx * tlx_yzz_xxz_0[j];

            tly_xyzz_xxxz_0[j] = pa_x[j] * tly_yzz_xxxz_0[j] + 1.5 * fl1_fx * tly_yzz_xxz_0[j] + 0.5 * fl1_fx * tpz_yzz_xxxz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xxxz_0[j];

            tlz_xyzz_xxxz_0[j] = pa_x[j] * tlz_yzz_xxxz_0[j] + 1.5 * fl1_fx * tlz_yzz_xxz_0[j] - 0.5 * fl1_fx * tpy_yzz_xxxz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xxxz_0[j];

            tlx_xyzz_xxyy_0[j] = pa_x[j] * tlx_yzz_xxyy_0[j] + fl1_fx * tlx_yzz_xyy_0[j];

            tly_xyzz_xxyy_0[j] =
                pa_x[j] * tly_yzz_xxyy_0[j] + fl1_fx * tly_yzz_xyy_0[j] + 0.5 * fl1_fx * tpz_yzz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xxyy_0[j];

            tlz_xyzz_xxyy_0[j] =
                pa_x[j] * tlz_yzz_xxyy_0[j] + fl1_fx * tlz_yzz_xyy_0[j] - 0.5 * fl1_fx * tpy_yzz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xxyy_0[j];

            tlx_xyzz_xxyz_0[j] = pa_x[j] * tlx_yzz_xxyz_0[j] + fl1_fx * tlx_yzz_xyz_0[j];

            tly_xyzz_xxyz_0[j] =
                pa_x[j] * tly_yzz_xxyz_0[j] + fl1_fx * tly_yzz_xyz_0[j] + 0.5 * fl1_fx * tpz_yzz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xxyz_0[j];

            tlz_xyzz_xxyz_0[j] =
                pa_x[j] * tlz_yzz_xxyz_0[j] + fl1_fx * tlz_yzz_xyz_0[j] - 0.5 * fl1_fx * tpy_yzz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xxyz_0[j];

            tlx_xyzz_xxzz_0[j] = pa_x[j] * tlx_yzz_xxzz_0[j] + fl1_fx * tlx_yzz_xzz_0[j];

            tly_xyzz_xxzz_0[j] =
                pa_x[j] * tly_yzz_xxzz_0[j] + fl1_fx * tly_yzz_xzz_0[j] + 0.5 * fl1_fx * tpz_yzz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_xxzz_0[j];

            tlz_xyzz_xxzz_0[j] =
                pa_x[j] * tlz_yzz_xxzz_0[j] + fl1_fx * tlz_yzz_xzz_0[j] - 0.5 * fl1_fx * tpy_yzz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_xxzz_0[j];

            tlx_xyzz_xyyy_0[j] = pa_x[j] * tlx_yzz_xyyy_0[j] + 0.5 * fl1_fx * tlx_yzz_yyy_0[j];

            tly_xyzz_xyyy_0[j] = pa_x[j] * tly_yzz_xyyy_0[j] + 0.5 * fl1_fx * tly_yzz_yyy_0[j] + 0.5 * fl1_fx * tpz_yzz_xyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xyyy_0[j];

            tlz_xyzz_xyyy_0[j] = pa_x[j] * tlz_yzz_xyyy_0[j] + 0.5 * fl1_fx * tlz_yzz_yyy_0[j] - 0.5 * fl1_fx * tpy_yzz_xyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xyyy_0[j];

            tlx_xyzz_xyyz_0[j] = pa_x[j] * tlx_yzz_xyyz_0[j] + 0.5 * fl1_fx * tlx_yzz_yyz_0[j];

            tly_xyzz_xyyz_0[j] = pa_x[j] * tly_yzz_xyyz_0[j] + 0.5 * fl1_fx * tly_yzz_yyz_0[j] + 0.5 * fl1_fx * tpz_yzz_xyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xyyz_0[j];

            tlz_xyzz_xyyz_0[j] = pa_x[j] * tlz_yzz_xyyz_0[j] + 0.5 * fl1_fx * tlz_yzz_yyz_0[j] - 0.5 * fl1_fx * tpy_yzz_xyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xyyz_0[j];

            tlx_xyzz_xyzz_0[j] = pa_x[j] * tlx_yzz_xyzz_0[j] + 0.5 * fl1_fx * tlx_yzz_yzz_0[j];

            tly_xyzz_xyzz_0[j] = pa_x[j] * tly_yzz_xyzz_0[j] + 0.5 * fl1_fx * tly_yzz_yzz_0[j] + 0.5 * fl1_fx * tpz_yzz_xyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xyzz_0[j];

            tlz_xyzz_xyzz_0[j] = pa_x[j] * tlz_yzz_xyzz_0[j] + 0.5 * fl1_fx * tlz_yzz_yzz_0[j] - 0.5 * fl1_fx * tpy_yzz_xyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xyzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_387_435(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 129);

        auto tly_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tlz_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tlx_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 130);

        auto tly_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tlz_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tlx_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 131);

        auto tly_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tlz_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tlx_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 132);

        auto tly_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tlz_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tlx_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 133);

        auto tly_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tlz_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tlx_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 134);

        auto tly_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tlz_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tlx_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 135);

        auto tly_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tlz_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tlx_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 136);

        auto tly_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tlz_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tlx_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 137);

        auto tly_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tlz_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tlx_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 138);

        auto tly_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tlz_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tlx_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 139);

        auto tly_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tlz_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tlx_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 140);

        auto tly_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tlz_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tlx_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 141);

        auto tly_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tlz_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tlx_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 142);

        auto tly_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tlz_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tlx_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 143);

        auto tly_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tlz_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tlx_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 144);

        auto tly_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tlz_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tlx_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 89);

        auto tly_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tlz_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tlx_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 90);

        auto tly_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tlz_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tlx_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 91);

        auto tly_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tlz_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tlx_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 92);

        auto tly_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tlz_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tlx_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 93);

        auto tly_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tlz_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tlx_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 94);

        auto tly_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tlz_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tlx_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 95);

        auto tly_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tlz_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tlx_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 96);

        auto tly_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tlz_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tlx_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 97);

        auto tly_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tlz_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tlx_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 98);

        auto tly_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tlz_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tlx_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 99);

        auto tly_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tlz_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto tpy_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tpz_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tpy_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tpz_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tpy_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tpz_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tpy_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tpz_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tpy_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tpz_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tpy_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tpz_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tpy_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tpz_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tpy_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tpz_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tpy_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tpz_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tpy_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tpz_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tpy_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tpz_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tpy_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tpz_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tpy_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tpz_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tpy_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tpz_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tpy_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tpz_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tpy_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tpz_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tdy_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tdz_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tdy_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tdz_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tdy_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tdz_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tdy_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tdz_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tdy_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tdz_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tdy_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tdz_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tdy_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tdz_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tdy_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tdz_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tdy_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tdz_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tdy_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tdz_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tdy_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tdz_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tdy_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tdz_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tdy_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tdz_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tdy_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tdz_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tdy_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tdz_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tdy_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tdz_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 144);

        // set up pointers to integrals

        auto tlx_xyzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 129);

        auto tly_xyzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 129);

        auto tlz_xyzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 129);

        auto tlx_xyzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 130);

        auto tly_xyzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 130);

        auto tlz_xyzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 130);

        auto tlx_xyzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 131);

        auto tly_xyzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 131);

        auto tlz_xyzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 131);

        auto tlx_xyzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 132);

        auto tly_xyzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 132);

        auto tlz_xyzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 132);

        auto tlx_xyzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 133);

        auto tly_xyzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 133);

        auto tlz_xyzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 133);

        auto tlx_xyzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 134);

        auto tly_xyzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 134);

        auto tlz_xyzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 134);

        auto tlx_xzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 135);

        auto tly_xzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 135);

        auto tlz_xzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 135);

        auto tlx_xzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 136);

        auto tly_xzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 136);

        auto tlz_xzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 136);

        auto tlx_xzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 137);

        auto tly_xzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 137);

        auto tlz_xzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 137);

        auto tlx_xzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 138);

        auto tly_xzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 138);

        auto tlz_xzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 138);

        auto tlx_xzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 139);

        auto tly_xzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 139);

        auto tlz_xzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 139);

        auto tlx_xzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 140);

        auto tly_xzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 140);

        auto tlz_xzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 140);

        auto tlx_xzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 141);

        auto tly_xzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 141);

        auto tlz_xzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 141);

        auto tlx_xzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 142);

        auto tly_xzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 142);

        auto tlz_xzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 142);

        auto tlx_xzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 143);

        auto tly_xzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 143);

        auto tlz_xzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 143);

        auto tlx_xzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 144);

        auto tly_xzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 144);

        auto tlz_xzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 144);

        // Batch of Integrals (387,435)

        #pragma omp simd aligned(fgb, fx, pa_x, tdy_yzz_xzzz_0, tdy_yzz_yyyy_0, tdy_yzz_yyyz_0, \
                                     tdy_yzz_yyzz_0, tdy_yzz_yzzz_0, tdy_yzz_zzzz_0, tdy_zzz_xxxx_0, tdy_zzz_xxxy_0, \
                                     tdy_zzz_xxxz_0, tdy_zzz_xxyy_0, tdy_zzz_xxyz_0, tdy_zzz_xxzz_0, tdy_zzz_xyyy_0, \
                                     tdy_zzz_xyyz_0, tdy_zzz_xyzz_0, tdy_zzz_xzzz_0, tdz_yzz_xzzz_0, tdz_yzz_yyyy_0, \
                                     tdz_yzz_yyyz_0, tdz_yzz_yyzz_0, tdz_yzz_yzzz_0, tdz_yzz_zzzz_0, tdz_zzz_xxxx_0, \
                                     tdz_zzz_xxxy_0, tdz_zzz_xxxz_0, tdz_zzz_xxyy_0, tdz_zzz_xxyz_0, tdz_zzz_xxzz_0, \
                                     tdz_zzz_xyyy_0, tdz_zzz_xyyz_0, tdz_zzz_xyzz_0, tdz_zzz_xzzz_0, tlx_xyzz_xzzz_0, \
                                     tlx_xyzz_yyyy_0, tlx_xyzz_yyyz_0, tlx_xyzz_yyzz_0, tlx_xyzz_yzzz_0, tlx_xyzz_zzzz_0, \
                                     tlx_xzzz_xxxx_0, tlx_xzzz_xxxy_0, tlx_xzzz_xxxz_0, tlx_xzzz_xxyy_0, tlx_xzzz_xxyz_0, \
                                     tlx_xzzz_xxzz_0, tlx_xzzz_xyyy_0, tlx_xzzz_xyyz_0, tlx_xzzz_xyzz_0, tlx_xzzz_xzzz_0, \
                                     tlx_yzz_xzzz_0, tlx_yzz_yyyy_0, tlx_yzz_yyyz_0, tlx_yzz_yyzz_0, tlx_yzz_yzzz_0, \
                                     tlx_yzz_zzz_0, tlx_yzz_zzzz_0, tlx_zzz_xxx_0, tlx_zzz_xxxx_0, tlx_zzz_xxxy_0, \
                                     tlx_zzz_xxxz_0, tlx_zzz_xxy_0, tlx_zzz_xxyy_0, tlx_zzz_xxyz_0, tlx_zzz_xxz_0, \
                                     tlx_zzz_xxzz_0, tlx_zzz_xyy_0, tlx_zzz_xyyy_0, tlx_zzz_xyyz_0, tlx_zzz_xyz_0, \
                                     tlx_zzz_xyzz_0, tlx_zzz_xzz_0, tlx_zzz_xzzz_0, tlx_zzz_yyy_0, tlx_zzz_yyz_0, \
                                     tlx_zzz_yzz_0, tlx_zzz_zzz_0, tly_xyzz_xzzz_0, tly_xyzz_yyyy_0, tly_xyzz_yyyz_0, \
                                     tly_xyzz_yyzz_0, tly_xyzz_yzzz_0, tly_xyzz_zzzz_0, tly_xzzz_xxxx_0, tly_xzzz_xxxy_0, \
                                     tly_xzzz_xxxz_0, tly_xzzz_xxyy_0, tly_xzzz_xxyz_0, tly_xzzz_xxzz_0, tly_xzzz_xyyy_0, \
                                     tly_xzzz_xyyz_0, tly_xzzz_xyzz_0, tly_xzzz_xzzz_0, tly_yzz_xzzz_0, tly_yzz_yyyy_0, \
                                     tly_yzz_yyyz_0, tly_yzz_yyzz_0, tly_yzz_yzzz_0, tly_yzz_zzz_0, tly_yzz_zzzz_0, \
                                     tly_zzz_xxx_0, tly_zzz_xxxx_0, tly_zzz_xxxy_0, tly_zzz_xxxz_0, tly_zzz_xxy_0, \
                                     tly_zzz_xxyy_0, tly_zzz_xxyz_0, tly_zzz_xxz_0, tly_zzz_xxzz_0, tly_zzz_xyy_0, \
                                     tly_zzz_xyyy_0, tly_zzz_xyyz_0, tly_zzz_xyz_0, tly_zzz_xyzz_0, tly_zzz_xzz_0, \
                                     tly_zzz_xzzz_0, tly_zzz_yyy_0, tly_zzz_yyz_0, tly_zzz_yzz_0, tly_zzz_zzz_0, \
                                     tlz_xyzz_xzzz_0, tlz_xyzz_yyyy_0, tlz_xyzz_yyyz_0, tlz_xyzz_yyzz_0, tlz_xyzz_yzzz_0, \
                                     tlz_xyzz_zzzz_0, tlz_xzzz_xxxx_0, tlz_xzzz_xxxy_0, tlz_xzzz_xxxz_0, tlz_xzzz_xxyy_0, \
                                     tlz_xzzz_xxyz_0, tlz_xzzz_xxzz_0, tlz_xzzz_xyyy_0, tlz_xzzz_xyyz_0, tlz_xzzz_xyzz_0, \
                                     tlz_xzzz_xzzz_0, tlz_yzz_xzzz_0, tlz_yzz_yyyy_0, tlz_yzz_yyyz_0, tlz_yzz_yyzz_0, \
                                     tlz_yzz_yzzz_0, tlz_yzz_zzz_0, tlz_yzz_zzzz_0, tlz_zzz_xxx_0, tlz_zzz_xxxx_0, \
                                     tlz_zzz_xxxy_0, tlz_zzz_xxxz_0, tlz_zzz_xxy_0, tlz_zzz_xxyy_0, tlz_zzz_xxyz_0, \
                                     tlz_zzz_xxz_0, tlz_zzz_xxzz_0, tlz_zzz_xyy_0, tlz_zzz_xyyy_0, tlz_zzz_xyyz_0, \
                                     tlz_zzz_xyz_0, tlz_zzz_xyzz_0, tlz_zzz_xzz_0, tlz_zzz_xzzz_0, tlz_zzz_yyy_0, \
                                     tlz_zzz_yyz_0, tlz_zzz_yzz_0, tlz_zzz_zzz_0, tpy_yzz_xzzz_0, tpy_yzz_yyyy_0, \
                                     tpy_yzz_yyyz_0, tpy_yzz_yyzz_0, tpy_yzz_yzzz_0, tpy_yzz_zzzz_0, tpy_zzz_xxxx_0, \
                                     tpy_zzz_xxxy_0, tpy_zzz_xxxz_0, tpy_zzz_xxyy_0, tpy_zzz_xxyz_0, tpy_zzz_xxzz_0, \
                                     tpy_zzz_xyyy_0, tpy_zzz_xyyz_0, tpy_zzz_xyzz_0, tpy_zzz_xzzz_0, tpz_yzz_xzzz_0, \
                                     tpz_yzz_yyyy_0, tpz_yzz_yyyz_0, tpz_yzz_yyzz_0, tpz_yzz_yzzz_0, tpz_yzz_zzzz_0, \
                                     tpz_zzz_xxxx_0, tpz_zzz_xxxy_0, tpz_zzz_xxxz_0, tpz_zzz_xxyy_0, tpz_zzz_xxyz_0, \
                                     tpz_zzz_xxzz_0, tpz_zzz_xyyy_0, tpz_zzz_xyyz_0, tpz_zzz_xyzz_0, tpz_zzz_xzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xyzz_xzzz_0[j] = pa_x[j] * tlx_yzz_xzzz_0[j] + 0.5 * fl1_fx * tlx_yzz_zzz_0[j];

            tly_xyzz_xzzz_0[j] = pa_x[j] * tly_yzz_xzzz_0[j] + 0.5 * fl1_fx * tly_yzz_zzz_0[j] + 0.5 * fl1_fx * tpz_yzz_xzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_yzz_xzzz_0[j];

            tlz_xyzz_xzzz_0[j] = pa_x[j] * tlz_yzz_xzzz_0[j] + 0.5 * fl1_fx * tlz_yzz_zzz_0[j] - 0.5 * fl1_fx * tpy_yzz_xzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_yzz_xzzz_0[j];

            tlx_xyzz_yyyy_0[j] = pa_x[j] * tlx_yzz_yyyy_0[j];

            tly_xyzz_yyyy_0[j] = pa_x[j] * tly_yzz_yyyy_0[j] + 0.5 * fl1_fx * tpz_yzz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yyyy_0[j];

            tlz_xyzz_yyyy_0[j] = pa_x[j] * tlz_yzz_yyyy_0[j] - 0.5 * fl1_fx * tpy_yzz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yyyy_0[j];

            tlx_xyzz_yyyz_0[j] = pa_x[j] * tlx_yzz_yyyz_0[j];

            tly_xyzz_yyyz_0[j] = pa_x[j] * tly_yzz_yyyz_0[j] + 0.5 * fl1_fx * tpz_yzz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yyyz_0[j];

            tlz_xyzz_yyyz_0[j] = pa_x[j] * tlz_yzz_yyyz_0[j] - 0.5 * fl1_fx * tpy_yzz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yyyz_0[j];

            tlx_xyzz_yyzz_0[j] = pa_x[j] * tlx_yzz_yyzz_0[j];

            tly_xyzz_yyzz_0[j] = pa_x[j] * tly_yzz_yyzz_0[j] + 0.5 * fl1_fx * tpz_yzz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yyzz_0[j];

            tlz_xyzz_yyzz_0[j] = pa_x[j] * tlz_yzz_yyzz_0[j] - 0.5 * fl1_fx * tpy_yzz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yyzz_0[j];

            tlx_xyzz_yzzz_0[j] = pa_x[j] * tlx_yzz_yzzz_0[j];

            tly_xyzz_yzzz_0[j] = pa_x[j] * tly_yzz_yzzz_0[j] + 0.5 * fl1_fx * tpz_yzz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_yzzz_0[j];

            tlz_xyzz_yzzz_0[j] = pa_x[j] * tlz_yzz_yzzz_0[j] - 0.5 * fl1_fx * tpy_yzz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_yzzz_0[j];

            tlx_xyzz_zzzz_0[j] = pa_x[j] * tlx_yzz_zzzz_0[j];

            tly_xyzz_zzzz_0[j] = pa_x[j] * tly_yzz_zzzz_0[j] + 0.5 * fl1_fx * tpz_yzz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_yzz_zzzz_0[j];

            tlz_xyzz_zzzz_0[j] = pa_x[j] * tlz_yzz_zzzz_0[j] - 0.5 * fl1_fx * tpy_yzz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_yzz_zzzz_0[j];

            tlx_xzzz_xxxx_0[j] = pa_x[j] * tlx_zzz_xxxx_0[j] + 2.0 * fl1_fx * tlx_zzz_xxx_0[j];

            tly_xzzz_xxxx_0[j] = pa_x[j] * tly_zzz_xxxx_0[j] + 2.0 * fl1_fx * tly_zzz_xxx_0[j] + 0.5 * fl1_fx * tpz_zzz_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xxxx_0[j];

            tlz_xzzz_xxxx_0[j] = pa_x[j] * tlz_zzz_xxxx_0[j] + 2.0 * fl1_fx * tlz_zzz_xxx_0[j] - 0.5 * fl1_fx * tpy_zzz_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xxxx_0[j];

            tlx_xzzz_xxxy_0[j] = pa_x[j] * tlx_zzz_xxxy_0[j] + 1.5 * fl1_fx * tlx_zzz_xxy_0[j];

            tly_xzzz_xxxy_0[j] = pa_x[j] * tly_zzz_xxxy_0[j] + 1.5 * fl1_fx * tly_zzz_xxy_0[j] + 0.5 * fl1_fx * tpz_zzz_xxxy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xxxy_0[j];

            tlz_xzzz_xxxy_0[j] = pa_x[j] * tlz_zzz_xxxy_0[j] + 1.5 * fl1_fx * tlz_zzz_xxy_0[j] - 0.5 * fl1_fx * tpy_zzz_xxxy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xxxy_0[j];

            tlx_xzzz_xxxz_0[j] = pa_x[j] * tlx_zzz_xxxz_0[j] + 1.5 * fl1_fx * tlx_zzz_xxz_0[j];

            tly_xzzz_xxxz_0[j] = pa_x[j] * tly_zzz_xxxz_0[j] + 1.5 * fl1_fx * tly_zzz_xxz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxxz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xxxz_0[j];

            tlz_xzzz_xxxz_0[j] = pa_x[j] * tlz_zzz_xxxz_0[j] + 1.5 * fl1_fx * tlz_zzz_xxz_0[j] - 0.5 * fl1_fx * tpy_zzz_xxxz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xxxz_0[j];

            tlx_xzzz_xxyy_0[j] = pa_x[j] * tlx_zzz_xxyy_0[j] + fl1_fx * tlx_zzz_xyy_0[j];

            tly_xzzz_xxyy_0[j] =
                pa_x[j] * tly_zzz_xxyy_0[j] + fl1_fx * tly_zzz_xyy_0[j] + 0.5 * fl1_fx * tpz_zzz_xxyy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xxyy_0[j];

            tlz_xzzz_xxyy_0[j] =
                pa_x[j] * tlz_zzz_xxyy_0[j] + fl1_fx * tlz_zzz_xyy_0[j] - 0.5 * fl1_fx * tpy_zzz_xxyy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xxyy_0[j];

            tlx_xzzz_xxyz_0[j] = pa_x[j] * tlx_zzz_xxyz_0[j] + fl1_fx * tlx_zzz_xyz_0[j];

            tly_xzzz_xxyz_0[j] =
                pa_x[j] * tly_zzz_xxyz_0[j] + fl1_fx * tly_zzz_xyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxyz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xxyz_0[j];

            tlz_xzzz_xxyz_0[j] =
                pa_x[j] * tlz_zzz_xxyz_0[j] + fl1_fx * tlz_zzz_xyz_0[j] - 0.5 * fl1_fx * tpy_zzz_xxyz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xxyz_0[j];

            tlx_xzzz_xxzz_0[j] = pa_x[j] * tlx_zzz_xxzz_0[j] + fl1_fx * tlx_zzz_xzz_0[j];

            tly_xzzz_xxzz_0[j] =
                pa_x[j] * tly_zzz_xxzz_0[j] + fl1_fx * tly_zzz_xzz_0[j] + 0.5 * fl1_fx * tpz_zzz_xxzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_xxzz_0[j];

            tlz_xzzz_xxzz_0[j] =
                pa_x[j] * tlz_zzz_xxzz_0[j] + fl1_fx * tlz_zzz_xzz_0[j] - 0.5 * fl1_fx * tpy_zzz_xxzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_xxzz_0[j];

            tlx_xzzz_xyyy_0[j] = pa_x[j] * tlx_zzz_xyyy_0[j] + 0.5 * fl1_fx * tlx_zzz_yyy_0[j];

            tly_xzzz_xyyy_0[j] = pa_x[j] * tly_zzz_xyyy_0[j] + 0.5 * fl1_fx * tly_zzz_yyy_0[j] + 0.5 * fl1_fx * tpz_zzz_xyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xyyy_0[j];

            tlz_xzzz_xyyy_0[j] = pa_x[j] * tlz_zzz_xyyy_0[j] + 0.5 * fl1_fx * tlz_zzz_yyy_0[j] - 0.5 * fl1_fx * tpy_zzz_xyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xyyy_0[j];

            tlx_xzzz_xyyz_0[j] = pa_x[j] * tlx_zzz_xyyz_0[j] + 0.5 * fl1_fx * tlx_zzz_yyz_0[j];

            tly_xzzz_xyyz_0[j] = pa_x[j] * tly_zzz_xyyz_0[j] + 0.5 * fl1_fx * tly_zzz_yyz_0[j] + 0.5 * fl1_fx * tpz_zzz_xyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xyyz_0[j];

            tlz_xzzz_xyyz_0[j] = pa_x[j] * tlz_zzz_xyyz_0[j] + 0.5 * fl1_fx * tlz_zzz_yyz_0[j] - 0.5 * fl1_fx * tpy_zzz_xyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xyyz_0[j];

            tlx_xzzz_xyzz_0[j] = pa_x[j] * tlx_zzz_xyzz_0[j] + 0.5 * fl1_fx * tlx_zzz_yzz_0[j];

            tly_xzzz_xyzz_0[j] = pa_x[j] * tly_zzz_xyzz_0[j] + 0.5 * fl1_fx * tly_zzz_yzz_0[j] + 0.5 * fl1_fx * tpz_zzz_xyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xyzz_0[j];

            tlz_xzzz_xyzz_0[j] = pa_x[j] * tlz_zzz_xyzz_0[j] + 0.5 * fl1_fx * tlz_zzz_yzz_0[j] - 0.5 * fl1_fx * tpy_zzz_xyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xyzz_0[j];

            tlx_xzzz_xzzz_0[j] = pa_x[j] * tlx_zzz_xzzz_0[j] + 0.5 * fl1_fx * tlx_zzz_zzz_0[j];

            tly_xzzz_xzzz_0[j] = pa_x[j] * tly_zzz_xzzz_0[j] + 0.5 * fl1_fx * tly_zzz_zzz_0[j] + 0.5 * fl1_fx * tpz_zzz_xzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdz_zzz_xzzz_0[j];

            tlz_xzzz_xzzz_0[j] = pa_x[j] * tlz_zzz_xzzz_0[j] + 0.5 * fl1_fx * tlz_zzz_zzz_0[j] - 0.5 * fl1_fx * tpy_zzz_xzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdy_zzz_xzzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_435_483(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 90);

        auto tly_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 90);

        auto tlz_yyy_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tlx_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 91);

        auto tly_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 91);

        auto tlz_yyy_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tlx_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 92);

        auto tly_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 92);

        auto tlz_yyy_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tlx_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 93);

        auto tly_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 93);

        auto tlz_yyy_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tlx_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 94);

        auto tly_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 94);

        auto tlz_yyy_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tlx_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 95);

        auto tly_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 95);

        auto tlz_yyy_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tlx_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 96);

        auto tly_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 96);

        auto tlz_yyy_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tlx_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 97);

        auto tly_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 97);

        auto tlz_yyy_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tlx_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 98);

        auto tly_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 98);

        auto tlz_yyy_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tlx_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 99);

        auto tly_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 99);

        auto tlz_yyy_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tlx_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 100);

        auto tly_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 100);

        auto tlz_yyy_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tlx_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 145);

        auto tly_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tlz_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tlx_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 146);

        auto tly_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tlz_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tlx_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 147);

        auto tly_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tlz_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tlx_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 148);

        auto tly_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tlz_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tlx_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 149);

        auto tly_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tlz_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 149);

        auto tlx_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 45);

        auto tly_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 45);

        auto tlz_yy_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 45);

        auto tlx_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 46);

        auto tly_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 46);

        auto tlz_yy_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 46);

        auto tlx_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 47);

        auto tly_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 47);

        auto tlz_yy_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 47);

        auto tlx_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 48);

        auto tly_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 48);

        auto tlz_yy_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 48);

        auto tlx_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 49);

        auto tly_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 49);

        auto tlz_yy_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 49);

        auto tlx_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 50);

        auto tly_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 50);

        auto tlz_yy_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 50);

        auto tlx_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 51);

        auto tly_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 51);

        auto tlz_yy_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 51);

        auto tlx_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 52);

        auto tly_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 52);

        auto tlz_yy_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 52);

        auto tlx_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 53);

        auto tly_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 53);

        auto tlz_yy_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 53);

        auto tlx_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 54);

        auto tly_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 54);

        auto tlz_yy_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 54);

        auto tlx_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 55);

        auto tly_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 55);

        auto tlz_yy_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 55);

        auto tlx_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 60);

        auto tly_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 60);

        auto tlz_yyy_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 60);

        auto tlx_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 61);

        auto tly_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 61);

        auto tlz_yyy_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 61);

        auto tlx_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 62);

        auto tly_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 62);

        auto tlz_yyy_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 62);

        auto tlx_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 63);

        auto tly_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 63);

        auto tlz_yyy_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 63);

        auto tlx_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 64);

        auto tly_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 64);

        auto tlz_yyy_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 64);

        auto tlx_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 65);

        auto tly_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 65);

        auto tlz_yyy_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 65);

        auto tlx_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 66);

        auto tly_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 66);

        auto tlz_yyy_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 66);

        auto tpx_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 90);

        auto tpz_yyy_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tpx_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 91);

        auto tpz_yyy_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tpx_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 92);

        auto tpz_yyy_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tpx_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 93);

        auto tpz_yyy_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tpx_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 94);

        auto tpz_yyy_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tpx_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 95);

        auto tpz_yyy_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tpx_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 96);

        auto tpz_yyy_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tpx_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 97);

        auto tpz_yyy_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tpx_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 98);

        auto tpz_yyy_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tpx_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 99);

        auto tpz_yyy_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tpx_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 100);

        auto tpz_yyy_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tpy_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tpz_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tpy_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tpz_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tpy_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tpz_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tpy_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tpz_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tpy_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tpz_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 149);

        auto tdx_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 90);

        auto tdz_yyy_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 90);

        auto tdx_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 91);

        auto tdz_yyy_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 91);

        auto tdx_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 92);

        auto tdz_yyy_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 92);

        auto tdx_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 93);

        auto tdz_yyy_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 93);

        auto tdx_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 94);

        auto tdz_yyy_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 94);

        auto tdx_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 95);

        auto tdz_yyy_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 95);

        auto tdx_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 96);

        auto tdz_yyy_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 96);

        auto tdx_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 97);

        auto tdz_yyy_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 97);

        auto tdx_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 98);

        auto tdz_yyy_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 98);

        auto tdx_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 99);

        auto tdz_yyy_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 99);

        auto tdx_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 100);

        auto tdz_yyy_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 100);

        auto tdy_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tdz_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tdy_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tdz_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tdy_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tdz_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tdy_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tdz_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tdy_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tdz_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 149);

        // set up pointers to integrals

        auto tlx_xzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 145);

        auto tly_xzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 145);

        auto tlz_xzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 145);

        auto tlx_xzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 146);

        auto tly_xzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 146);

        auto tlz_xzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 146);

        auto tlx_xzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 147);

        auto tly_xzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 147);

        auto tlz_xzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 147);

        auto tlx_xzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 148);

        auto tly_xzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 148);

        auto tlz_xzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 148);

        auto tlx_xzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 149);

        auto tly_xzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 149);

        auto tlz_xzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 149);

        auto tlx_yyyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 150);

        auto tly_yyyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 150);

        auto tlz_yyyy_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 150);

        auto tlx_yyyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 151);

        auto tly_yyyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 151);

        auto tlz_yyyy_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 151);

        auto tlx_yyyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 152);

        auto tly_yyyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 152);

        auto tlz_yyyy_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 152);

        auto tlx_yyyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 153);

        auto tly_yyyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 153);

        auto tlz_yyyy_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 153);

        auto tlx_yyyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 154);

        auto tly_yyyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 154);

        auto tlz_yyyy_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 154);

        auto tlx_yyyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 155);

        auto tly_yyyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 155);

        auto tlz_yyyy_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 155);

        auto tlx_yyyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 156);

        auto tly_yyyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 156);

        auto tlz_yyyy_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 156);

        auto tlx_yyyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 157);

        auto tly_yyyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 157);

        auto tlz_yyyy_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 157);

        auto tlx_yyyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 158);

        auto tly_yyyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 158);

        auto tlz_yyyy_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 158);

        auto tlx_yyyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 159);

        auto tly_yyyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 159);

        auto tlz_yyyy_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 159);

        auto tlx_yyyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 160);

        auto tly_yyyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 160);

        auto tlz_yyyy_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 160);

        // Batch of Integrals (435,483)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tdx_yyy_xxxx_0, tdx_yyy_xxxy_0, tdx_yyy_xxxz_0, \
                                     tdx_yyy_xxyy_0, tdx_yyy_xxyz_0, tdx_yyy_xxzz_0, tdx_yyy_xyyy_0, tdx_yyy_xyyz_0, \
                                     tdx_yyy_xyzz_0, tdx_yyy_xzzz_0, tdx_yyy_yyyy_0, tdy_zzz_yyyy_0, tdy_zzz_yyyz_0, \
                                     tdy_zzz_yyzz_0, tdy_zzz_yzzz_0, tdy_zzz_zzzz_0, tdz_yyy_xxxx_0, tdz_yyy_xxxy_0, \
                                     tdz_yyy_xxxz_0, tdz_yyy_xxyy_0, tdz_yyy_xxyz_0, tdz_yyy_xxzz_0, tdz_yyy_xyyy_0, \
                                     tdz_yyy_xyyz_0, tdz_yyy_xyzz_0, tdz_yyy_xzzz_0, tdz_yyy_yyyy_0, tdz_zzz_yyyy_0, \
                                     tdz_zzz_yyyz_0, tdz_zzz_yyzz_0, tdz_zzz_yzzz_0, tdz_zzz_zzzz_0, tlx_xzzz_yyyy_0, \
                                     tlx_xzzz_yyyz_0, tlx_xzzz_yyzz_0, tlx_xzzz_yzzz_0, tlx_xzzz_zzzz_0, tlx_yy_xxxx_0, \
                                     tlx_yy_xxxy_0, tlx_yy_xxxz_0, tlx_yy_xxyy_0, tlx_yy_xxyz_0, tlx_yy_xxzz_0, \
                                     tlx_yy_xyyy_0, tlx_yy_xyyz_0, tlx_yy_xyzz_0, tlx_yy_xzzz_0, tlx_yy_yyyy_0, \
                                     tlx_yyy_xxx_0, tlx_yyy_xxxx_0, tlx_yyy_xxxy_0, tlx_yyy_xxxz_0, tlx_yyy_xxy_0, \
                                     tlx_yyy_xxyy_0, tlx_yyy_xxyz_0, tlx_yyy_xxz_0, tlx_yyy_xxzz_0, tlx_yyy_xyy_0, \
                                     tlx_yyy_xyyy_0, tlx_yyy_xyyz_0, tlx_yyy_xyz_0, tlx_yyy_xyzz_0, tlx_yyy_xzz_0, \
                                     tlx_yyy_xzzz_0, tlx_yyy_yyy_0, tlx_yyy_yyyy_0, tlx_yyyy_xxxx_0, tlx_yyyy_xxxy_0, \
                                     tlx_yyyy_xxxz_0, tlx_yyyy_xxyy_0, tlx_yyyy_xxyz_0, tlx_yyyy_xxzz_0, tlx_yyyy_xyyy_0, \
                                     tlx_yyyy_xyyz_0, tlx_yyyy_xyzz_0, tlx_yyyy_xzzz_0, tlx_yyyy_yyyy_0, tlx_zzz_yyyy_0, \
                                     tlx_zzz_yyyz_0, tlx_zzz_yyzz_0, tlx_zzz_yzzz_0, tlx_zzz_zzzz_0, tly_xzzz_yyyy_0, \
                                     tly_xzzz_yyyz_0, tly_xzzz_yyzz_0, tly_xzzz_yzzz_0, tly_xzzz_zzzz_0, tly_yy_xxxx_0, \
                                     tly_yy_xxxy_0, tly_yy_xxxz_0, tly_yy_xxyy_0, tly_yy_xxyz_0, tly_yy_xxzz_0, \
                                     tly_yy_xyyy_0, tly_yy_xyyz_0, tly_yy_xyzz_0, tly_yy_xzzz_0, tly_yy_yyyy_0, \
                                     tly_yyy_xxx_0, tly_yyy_xxxx_0, tly_yyy_xxxy_0, tly_yyy_xxxz_0, tly_yyy_xxy_0, \
                                     tly_yyy_xxyy_0, tly_yyy_xxyz_0, tly_yyy_xxz_0, tly_yyy_xxzz_0, tly_yyy_xyy_0, \
                                     tly_yyy_xyyy_0, tly_yyy_xyyz_0, tly_yyy_xyz_0, tly_yyy_xyzz_0, tly_yyy_xzz_0, \
                                     tly_yyy_xzzz_0, tly_yyy_yyy_0, tly_yyy_yyyy_0, tly_yyyy_xxxx_0, tly_yyyy_xxxy_0, \
                                     tly_yyyy_xxxz_0, tly_yyyy_xxyy_0, tly_yyyy_xxyz_0, tly_yyyy_xxzz_0, tly_yyyy_xyyy_0, \
                                     tly_yyyy_xyyz_0, tly_yyyy_xyzz_0, tly_yyyy_xzzz_0, tly_yyyy_yyyy_0, tly_zzz_yyyy_0, \
                                     tly_zzz_yyyz_0, tly_zzz_yyzz_0, tly_zzz_yzzz_0, tly_zzz_zzzz_0, tlz_xzzz_yyyy_0, \
                                     tlz_xzzz_yyyz_0, tlz_xzzz_yyzz_0, tlz_xzzz_yzzz_0, tlz_xzzz_zzzz_0, tlz_yy_xxxx_0, \
                                     tlz_yy_xxxy_0, tlz_yy_xxxz_0, tlz_yy_xxyy_0, tlz_yy_xxyz_0, tlz_yy_xxzz_0, \
                                     tlz_yy_xyyy_0, tlz_yy_xyyz_0, tlz_yy_xyzz_0, tlz_yy_xzzz_0, tlz_yy_yyyy_0, \
                                     tlz_yyy_xxx_0, tlz_yyy_xxxx_0, tlz_yyy_xxxy_0, tlz_yyy_xxxz_0, tlz_yyy_xxy_0, \
                                     tlz_yyy_xxyy_0, tlz_yyy_xxyz_0, tlz_yyy_xxz_0, tlz_yyy_xxzz_0, tlz_yyy_xyy_0, \
                                     tlz_yyy_xyyy_0, tlz_yyy_xyyz_0, tlz_yyy_xyz_0, tlz_yyy_xyzz_0, tlz_yyy_xzz_0, \
                                     tlz_yyy_xzzz_0, tlz_yyy_yyy_0, tlz_yyy_yyyy_0, tlz_yyyy_xxxx_0, tlz_yyyy_xxxy_0, \
                                     tlz_yyyy_xxxz_0, tlz_yyyy_xxyy_0, tlz_yyyy_xxyz_0, tlz_yyyy_xxzz_0, tlz_yyyy_xyyy_0, \
                                     tlz_yyyy_xyyz_0, tlz_yyyy_xyzz_0, tlz_yyyy_xzzz_0, tlz_yyyy_yyyy_0, tlz_zzz_yyyy_0, \
                                     tlz_zzz_yyyz_0, tlz_zzz_yyzz_0, tlz_zzz_yzzz_0, tlz_zzz_zzzz_0, tpx_yyy_xxxx_0, \
                                     tpx_yyy_xxxy_0, tpx_yyy_xxxz_0, tpx_yyy_xxyy_0, tpx_yyy_xxyz_0, tpx_yyy_xxzz_0, \
                                     tpx_yyy_xyyy_0, tpx_yyy_xyyz_0, tpx_yyy_xyzz_0, tpx_yyy_xzzz_0, tpx_yyy_yyyy_0, \
                                     tpy_zzz_yyyy_0, tpy_zzz_yyyz_0, tpy_zzz_yyzz_0, tpy_zzz_yzzz_0, tpy_zzz_zzzz_0, \
                                     tpz_yyy_xxxx_0, tpz_yyy_xxxy_0, tpz_yyy_xxxz_0, tpz_yyy_xxyy_0, tpz_yyy_xxyz_0, \
                                     tpz_yyy_xxzz_0, tpz_yyy_xyyy_0, tpz_yyy_xyyz_0, tpz_yyy_xyzz_0, tpz_yyy_xzzz_0, \
                                     tpz_yyy_yyyy_0, tpz_zzz_yyyy_0, tpz_zzz_yyyz_0, tpz_zzz_yyzz_0, tpz_zzz_yzzz_0, \
                                     tpz_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xzzz_yyyy_0[j] = pa_x[j] * tlx_zzz_yyyy_0[j];

            tly_xzzz_yyyy_0[j] = pa_x[j] * tly_zzz_yyyy_0[j] + 0.5 * fl1_fx * tpz_zzz_yyyy_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yyyy_0[j];

            tlz_xzzz_yyyy_0[j] = pa_x[j] * tlz_zzz_yyyy_0[j] - 0.5 * fl1_fx * tpy_zzz_yyyy_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yyyy_0[j];

            tlx_xzzz_yyyz_0[j] = pa_x[j] * tlx_zzz_yyyz_0[j];

            tly_xzzz_yyyz_0[j] = pa_x[j] * tly_zzz_yyyz_0[j] + 0.5 * fl1_fx * tpz_zzz_yyyz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yyyz_0[j];

            tlz_xzzz_yyyz_0[j] = pa_x[j] * tlz_zzz_yyyz_0[j] - 0.5 * fl1_fx * tpy_zzz_yyyz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yyyz_0[j];

            tlx_xzzz_yyzz_0[j] = pa_x[j] * tlx_zzz_yyzz_0[j];

            tly_xzzz_yyzz_0[j] = pa_x[j] * tly_zzz_yyzz_0[j] + 0.5 * fl1_fx * tpz_zzz_yyzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yyzz_0[j];

            tlz_xzzz_yyzz_0[j] = pa_x[j] * tlz_zzz_yyzz_0[j] - 0.5 * fl1_fx * tpy_zzz_yyzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yyzz_0[j];

            tlx_xzzz_yzzz_0[j] = pa_x[j] * tlx_zzz_yzzz_0[j];

            tly_xzzz_yzzz_0[j] = pa_x[j] * tly_zzz_yzzz_0[j] + 0.5 * fl1_fx * tpz_zzz_yzzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_yzzz_0[j];

            tlz_xzzz_yzzz_0[j] = pa_x[j] * tlz_zzz_yzzz_0[j] - 0.5 * fl1_fx * tpy_zzz_yzzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_yzzz_0[j];

            tlx_xzzz_zzzz_0[j] = pa_x[j] * tlx_zzz_zzzz_0[j];

            tly_xzzz_zzzz_0[j] = pa_x[j] * tly_zzz_zzzz_0[j] + 0.5 * fl1_fx * tpz_zzz_zzzz_0[j] + fl1_fx * fl1_fgb * tdz_zzz_zzzz_0[j];

            tlz_xzzz_zzzz_0[j] = pa_x[j] * tlz_zzz_zzzz_0[j] - 0.5 * fl1_fx * tpy_zzz_zzzz_0[j] - fl1_fx * fl1_fgb * tdy_zzz_zzzz_0[j];

            tlx_yyyy_xxxx_0[j] = pa_y[j] * tlx_yyy_xxxx_0[j] + 1.5 * fl1_fx * tlx_yy_xxxx_0[j] - 0.5 * fl1_fx * tpz_yyy_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yyy_xxxx_0[j];

            tly_yyyy_xxxx_0[j] = pa_y[j] * tly_yyy_xxxx_0[j] + 1.5 * fl1_fx * tly_yy_xxxx_0[j];

            tlz_yyyy_xxxx_0[j] = pa_y[j] * tlz_yyy_xxxx_0[j] + 1.5 * fl1_fx * tlz_yy_xxxx_0[j] + 0.5 * fl1_fx * tpx_yyy_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yyy_xxxx_0[j];

            tlx_yyyy_xxxy_0[j] = pa_y[j] * tlx_yyy_xxxy_0[j] + 1.5 * fl1_fx * tlx_yy_xxxy_0[j] + 0.5 * fl1_fx * tlx_yyy_xxx_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xxxy_0[j];

            tly_yyyy_xxxy_0[j] = pa_y[j] * tly_yyy_xxxy_0[j] + 1.5 * fl1_fx * tly_yy_xxxy_0[j] + 0.5 * fl1_fx * tly_yyy_xxx_0[j];

            tlz_yyyy_xxxy_0[j] = pa_y[j] * tlz_yyy_xxxy_0[j] + 1.5 * fl1_fx * tlz_yy_xxxy_0[j] + 0.5 * fl1_fx * tlz_yyy_xxx_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xxxy_0[j];

            tlx_yyyy_xxxz_0[j] = pa_y[j] * tlx_yyy_xxxz_0[j] + 1.5 * fl1_fx * tlx_yy_xxxz_0[j] - 0.5 * fl1_fx * tpz_yyy_xxxz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yyy_xxxz_0[j];

            tly_yyyy_xxxz_0[j] = pa_y[j] * tly_yyy_xxxz_0[j] + 1.5 * fl1_fx * tly_yy_xxxz_0[j];

            tlz_yyyy_xxxz_0[j] = pa_y[j] * tlz_yyy_xxxz_0[j] + 1.5 * fl1_fx * tlz_yy_xxxz_0[j] + 0.5 * fl1_fx * tpx_yyy_xxxz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yyy_xxxz_0[j];

            tlx_yyyy_xxyy_0[j] = pa_y[j] * tlx_yyy_xxyy_0[j] + 1.5 * fl1_fx * tlx_yy_xxyy_0[j] + fl1_fx * tlx_yyy_xxy_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xxyy_0[j];

            tly_yyyy_xxyy_0[j] = pa_y[j] * tly_yyy_xxyy_0[j] + 1.5 * fl1_fx * tly_yy_xxyy_0[j] + fl1_fx * tly_yyy_xxy_0[j];

            tlz_yyyy_xxyy_0[j] = pa_y[j] * tlz_yyy_xxyy_0[j] + 1.5 * fl1_fx * tlz_yy_xxyy_0[j] + fl1_fx * tlz_yyy_xxy_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xxyy_0[j];

            tlx_yyyy_xxyz_0[j] = pa_y[j] * tlx_yyy_xxyz_0[j] + 1.5 * fl1_fx * tlx_yy_xxyz_0[j] + 0.5 * fl1_fx * tlx_yyy_xxz_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xxyz_0[j];

            tly_yyyy_xxyz_0[j] = pa_y[j] * tly_yyy_xxyz_0[j] + 1.5 * fl1_fx * tly_yy_xxyz_0[j] + 0.5 * fl1_fx * tly_yyy_xxz_0[j];

            tlz_yyyy_xxyz_0[j] = pa_y[j] * tlz_yyy_xxyz_0[j] + 1.5 * fl1_fx * tlz_yy_xxyz_0[j] + 0.5 * fl1_fx * tlz_yyy_xxz_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xxyz_0[j];

            tlx_yyyy_xxzz_0[j] = pa_y[j] * tlx_yyy_xxzz_0[j] + 1.5 * fl1_fx * tlx_yy_xxzz_0[j] - 0.5 * fl1_fx * tpz_yyy_xxzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yyy_xxzz_0[j];

            tly_yyyy_xxzz_0[j] = pa_y[j] * tly_yyy_xxzz_0[j] + 1.5 * fl1_fx * tly_yy_xxzz_0[j];

            tlz_yyyy_xxzz_0[j] = pa_y[j] * tlz_yyy_xxzz_0[j] + 1.5 * fl1_fx * tlz_yy_xxzz_0[j] + 0.5 * fl1_fx * tpx_yyy_xxzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yyy_xxzz_0[j];

            tlx_yyyy_xyyy_0[j] = pa_y[j] * tlx_yyy_xyyy_0[j] + 1.5 * fl1_fx * tlx_yy_xyyy_0[j] + 1.5 * fl1_fx * tlx_yyy_xyy_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xyyy_0[j];

            tly_yyyy_xyyy_0[j] = pa_y[j] * tly_yyy_xyyy_0[j] + 1.5 * fl1_fx * tly_yy_xyyy_0[j] + 1.5 * fl1_fx * tly_yyy_xyy_0[j];

            tlz_yyyy_xyyy_0[j] = pa_y[j] * tlz_yyy_xyyy_0[j] + 1.5 * fl1_fx * tlz_yy_xyyy_0[j] + 1.5 * fl1_fx * tlz_yyy_xyy_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xyyy_0[j];

            tlx_yyyy_xyyz_0[j] = pa_y[j] * tlx_yyy_xyyz_0[j] + 1.5 * fl1_fx * tlx_yy_xyyz_0[j] + fl1_fx * tlx_yyy_xyz_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xyyz_0[j];

            tly_yyyy_xyyz_0[j] = pa_y[j] * tly_yyy_xyyz_0[j] + 1.5 * fl1_fx * tly_yy_xyyz_0[j] + fl1_fx * tly_yyy_xyz_0[j];

            tlz_yyyy_xyyz_0[j] = pa_y[j] * tlz_yyy_xyyz_0[j] + 1.5 * fl1_fx * tlz_yy_xyyz_0[j] + fl1_fx * tlz_yyy_xyz_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xyyz_0[j];

            tlx_yyyy_xyzz_0[j] = pa_y[j] * tlx_yyy_xyzz_0[j] + 1.5 * fl1_fx * tlx_yy_xyzz_0[j] + 0.5 * fl1_fx * tlx_yyy_xzz_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_xyzz_0[j];

            tly_yyyy_xyzz_0[j] = pa_y[j] * tly_yyy_xyzz_0[j] + 1.5 * fl1_fx * tly_yy_xyzz_0[j] + 0.5 * fl1_fx * tly_yyy_xzz_0[j];

            tlz_yyyy_xyzz_0[j] = pa_y[j] * tlz_yyy_xyzz_0[j] + 1.5 * fl1_fx * tlz_yy_xyzz_0[j] + 0.5 * fl1_fx * tlz_yyy_xzz_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_xyzz_0[j];

            tlx_yyyy_xzzz_0[j] = pa_y[j] * tlx_yyy_xzzz_0[j] + 1.5 * fl1_fx * tlx_yy_xzzz_0[j] - 0.5 * fl1_fx * tpz_yyy_xzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yyy_xzzz_0[j];

            tly_yyyy_xzzz_0[j] = pa_y[j] * tly_yyy_xzzz_0[j] + 1.5 * fl1_fx * tly_yy_xzzz_0[j];

            tlz_yyyy_xzzz_0[j] = pa_y[j] * tlz_yyy_xzzz_0[j] + 1.5 * fl1_fx * tlz_yy_xzzz_0[j] + 0.5 * fl1_fx * tpx_yyy_xzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yyy_xzzz_0[j];

            tlx_yyyy_yyyy_0[j] = pa_y[j] * tlx_yyy_yyyy_0[j] + 1.5 * fl1_fx * tlx_yy_yyyy_0[j] + 2.0 * fl1_fx * tlx_yyy_yyy_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yyyy_0[j];

            tly_yyyy_yyyy_0[j] = pa_y[j] * tly_yyy_yyyy_0[j] + 1.5 * fl1_fx * tly_yy_yyyy_0[j] + 2.0 * fl1_fx * tly_yyy_yyy_0[j];

            tlz_yyyy_yyyy_0[j] = pa_y[j] * tlz_yyy_yyyy_0[j] + 1.5 * fl1_fx * tlz_yy_yyyy_0[j] + 2.0 * fl1_fx * tlz_yyy_yyy_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yyyy_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_483_531(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 101);

        auto tly_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 101);

        auto tlz_yyy_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tlx_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 102);

        auto tly_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 102);

        auto tlz_yyy_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tlx_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 103);

        auto tly_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 103);

        auto tlz_yyy_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tlx_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 104);

        auto tly_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 104);

        auto tlz_yyy_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tlx_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 105);

        auto tly_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 105);

        auto tlz_yyz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tlx_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 106);

        auto tly_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 106);

        auto tlz_yyz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tlx_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 107);

        auto tly_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 107);

        auto tlz_yyz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tlx_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 108);

        auto tly_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 108);

        auto tlz_yyz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tlx_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 109);

        auto tly_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 109);

        auto tlz_yyz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tlx_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 110);

        auto tly_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 110);

        auto tlz_yyz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tlx_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 111);

        auto tly_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 111);

        auto tlz_yyz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tlx_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 112);

        auto tly_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 112);

        auto tlz_yyz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tlx_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 113);

        auto tly_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 113);

        auto tlz_yyz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tlx_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 114);

        auto tly_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 114);

        auto tlz_yyz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tlx_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 115);

        auto tly_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 115);

        auto tlz_yyz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tlx_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 116);

        auto tly_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 116);

        auto tlz_yyz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tlx_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 56);

        auto tly_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 56);

        auto tlz_yy_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 56);

        auto tlx_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 57);

        auto tly_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 57);

        auto tlz_yy_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 57);

        auto tlx_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 58);

        auto tly_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 58);

        auto tlz_yy_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 58);

        auto tlx_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 59);

        auto tly_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 59);

        auto tlz_yy_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 59);

        auto tlx_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 60);

        auto tly_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 60);

        auto tlz_yz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 60);

        auto tlx_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 61);

        auto tly_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 61);

        auto tlz_yz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 61);

        auto tlx_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 62);

        auto tly_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 62);

        auto tlz_yz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 62);

        auto tlx_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 63);

        auto tly_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 63);

        auto tlz_yz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 63);

        auto tlx_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 64);

        auto tly_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 64);

        auto tlz_yz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 64);

        auto tlx_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 65);

        auto tly_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 65);

        auto tlz_yz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 65);

        auto tlx_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 66);

        auto tly_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 66);

        auto tlz_yz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 66);

        auto tlx_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 67);

        auto tly_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 67);

        auto tlz_yz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 67);

        auto tlx_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 68);

        auto tly_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 68);

        auto tlz_yz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 68);

        auto tlx_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 69);

        auto tly_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 69);

        auto tlz_yz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 69);

        auto tlx_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 70);

        auto tly_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 70);

        auto tlz_yz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 70);

        auto tlx_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 71);

        auto tly_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 71);

        auto tlz_yz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 71);

        auto tlx_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 67);

        auto tly_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 67);

        auto tlz_yyy_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 67);

        auto tlx_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 68);

        auto tly_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 68);

        auto tlz_yyy_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 68);

        auto tlx_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 69);

        auto tly_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 69);

        auto tlz_yyy_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 69);

        auto tlx_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 70);

        auto tly_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 70);

        auto tlz_yyz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 70);

        auto tlx_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 71);

        auto tly_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 71);

        auto tlz_yyz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 71);

        auto tlx_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 72);

        auto tly_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 72);

        auto tlz_yyz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 72);

        auto tlx_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 73);

        auto tly_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 73);

        auto tlz_yyz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 73);

        auto tlx_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 74);

        auto tly_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 74);

        auto tlz_yyz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 74);

        auto tlx_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 75);

        auto tly_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 75);

        auto tlz_yyz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 75);

        auto tlx_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 76);

        auto tly_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 76);

        auto tlz_yyz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 76);

        auto tlx_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 77);

        auto tly_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 77);

        auto tlz_yyz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 77);

        auto tpx_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 101);

        auto tpz_yyy_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tpx_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 102);

        auto tpz_yyy_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tpx_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 103);

        auto tpz_yyy_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tpx_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 104);

        auto tpz_yyy_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tpx_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 105);

        auto tpz_yyz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tpx_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 106);

        auto tpz_yyz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tpx_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 107);

        auto tpz_yyz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tpx_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 108);

        auto tpz_yyz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tpx_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 109);

        auto tpz_yyz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tpx_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 110);

        auto tpz_yyz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tpx_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 111);

        auto tpz_yyz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tpx_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 112);

        auto tpz_yyz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tpx_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 113);

        auto tpz_yyz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tpx_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 114);

        auto tpz_yyz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tpx_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 115);

        auto tpz_yyz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tpx_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 116);

        auto tpz_yyz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 116);

        auto tdx_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 101);

        auto tdz_yyy_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 101);

        auto tdx_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 102);

        auto tdz_yyy_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 102);

        auto tdx_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 103);

        auto tdz_yyy_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 103);

        auto tdx_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 104);

        auto tdz_yyy_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 104);

        auto tdx_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 105);

        auto tdz_yyz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 105);

        auto tdx_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 106);

        auto tdz_yyz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 106);

        auto tdx_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 107);

        auto tdz_yyz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 107);

        auto tdx_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 108);

        auto tdz_yyz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 108);

        auto tdx_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 109);

        auto tdz_yyz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 109);

        auto tdx_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 110);

        auto tdz_yyz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 110);

        auto tdx_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 111);

        auto tdz_yyz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 111);

        auto tdx_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 112);

        auto tdz_yyz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 112);

        auto tdx_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 113);

        auto tdz_yyz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 113);

        auto tdx_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 114);

        auto tdz_yyz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 114);

        auto tdx_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 115);

        auto tdz_yyz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 115);

        auto tdx_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 116);

        auto tdz_yyz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 116);

        // set up pointers to integrals

        auto tlx_yyyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 161);

        auto tly_yyyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 161);

        auto tlz_yyyy_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 161);

        auto tlx_yyyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 162);

        auto tly_yyyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 162);

        auto tlz_yyyy_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 162);

        auto tlx_yyyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 163);

        auto tly_yyyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 163);

        auto tlz_yyyy_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 163);

        auto tlx_yyyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 164);

        auto tly_yyyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 164);

        auto tlz_yyyy_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 164);

        auto tlx_yyyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 165);

        auto tly_yyyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 165);

        auto tlz_yyyz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 165);

        auto tlx_yyyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 166);

        auto tly_yyyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 166);

        auto tlz_yyyz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 166);

        auto tlx_yyyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 167);

        auto tly_yyyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 167);

        auto tlz_yyyz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 167);

        auto tlx_yyyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 168);

        auto tly_yyyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 168);

        auto tlz_yyyz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 168);

        auto tlx_yyyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 169);

        auto tly_yyyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 169);

        auto tlz_yyyz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 169);

        auto tlx_yyyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 170);

        auto tly_yyyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 170);

        auto tlz_yyyz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 170);

        auto tlx_yyyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 171);

        auto tly_yyyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 171);

        auto tlz_yyyz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 171);

        auto tlx_yyyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 172);

        auto tly_yyyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 172);

        auto tlz_yyyz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 172);

        auto tlx_yyyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 173);

        auto tly_yyyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 173);

        auto tlz_yyyz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 173);

        auto tlx_yyyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 174);

        auto tly_yyyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 174);

        auto tlz_yyyz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 174);

        auto tlx_yyyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 175);

        auto tly_yyyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 175);

        auto tlz_yyyz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 175);

        auto tlx_yyyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 176);

        auto tly_yyyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 176);

        auto tlz_yyyz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 176);

        // Batch of Integrals (483,531)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yyy_yyyz_0, tdx_yyy_yyzz_0, tdx_yyy_yzzz_0, \
                                     tdx_yyy_zzzz_0, tdx_yyz_xxxx_0, tdx_yyz_xxxy_0, tdx_yyz_xxxz_0, tdx_yyz_xxyy_0, \
                                     tdx_yyz_xxyz_0, tdx_yyz_xxzz_0, tdx_yyz_xyyy_0, tdx_yyz_xyyz_0, tdx_yyz_xyzz_0, \
                                     tdx_yyz_xzzz_0, tdx_yyz_yyyy_0, tdx_yyz_yyyz_0, tdz_yyy_yyyz_0, tdz_yyy_yyzz_0, \
                                     tdz_yyy_yzzz_0, tdz_yyy_zzzz_0, tdz_yyz_xxxx_0, tdz_yyz_xxxy_0, tdz_yyz_xxxz_0, \
                                     tdz_yyz_xxyy_0, tdz_yyz_xxyz_0, tdz_yyz_xxzz_0, tdz_yyz_xyyy_0, tdz_yyz_xyyz_0, \
                                     tdz_yyz_xyzz_0, tdz_yyz_xzzz_0, tdz_yyz_yyyy_0, tdz_yyz_yyyz_0, tlx_yy_yyyz_0, \
                                     tlx_yy_yyzz_0, tlx_yy_yzzz_0, tlx_yy_zzzz_0, tlx_yyy_yyyz_0, tlx_yyy_yyz_0, \
                                     tlx_yyy_yyzz_0, tlx_yyy_yzz_0, tlx_yyy_yzzz_0, tlx_yyy_zzz_0, tlx_yyy_zzzz_0, \
                                     tlx_yyyy_yyyz_0, tlx_yyyy_yyzz_0, tlx_yyyy_yzzz_0, tlx_yyyy_zzzz_0, tlx_yyyz_xxxx_0, \
                                     tlx_yyyz_xxxy_0, tlx_yyyz_xxxz_0, tlx_yyyz_xxyy_0, tlx_yyyz_xxyz_0, tlx_yyyz_xxzz_0, \
                                     tlx_yyyz_xyyy_0, tlx_yyyz_xyyz_0, tlx_yyyz_xyzz_0, tlx_yyyz_xzzz_0, tlx_yyyz_yyyy_0, \
                                     tlx_yyyz_yyyz_0, tlx_yyz_xxx_0, tlx_yyz_xxxx_0, tlx_yyz_xxxy_0, tlx_yyz_xxxz_0, \
                                     tlx_yyz_xxy_0, tlx_yyz_xxyy_0, tlx_yyz_xxyz_0, tlx_yyz_xxz_0, tlx_yyz_xxzz_0, \
                                     tlx_yyz_xyy_0, tlx_yyz_xyyy_0, tlx_yyz_xyyz_0, tlx_yyz_xyz_0, tlx_yyz_xyzz_0, \
                                     tlx_yyz_xzz_0, tlx_yyz_xzzz_0, tlx_yyz_yyy_0, tlx_yyz_yyyy_0, tlx_yyz_yyyz_0, \
                                     tlx_yyz_yyz_0, tlx_yz_xxxx_0, tlx_yz_xxxy_0, tlx_yz_xxxz_0, tlx_yz_xxyy_0, \
                                     tlx_yz_xxyz_0, tlx_yz_xxzz_0, tlx_yz_xyyy_0, tlx_yz_xyyz_0, tlx_yz_xyzz_0, \
                                     tlx_yz_xzzz_0, tlx_yz_yyyy_0, tlx_yz_yyyz_0, tly_yy_yyyz_0, tly_yy_yyzz_0, \
                                     tly_yy_yzzz_0, tly_yy_zzzz_0, tly_yyy_yyyz_0, tly_yyy_yyz_0, tly_yyy_yyzz_0, \
                                     tly_yyy_yzz_0, tly_yyy_yzzz_0, tly_yyy_zzz_0, tly_yyy_zzzz_0, tly_yyyy_yyyz_0, \
                                     tly_yyyy_yyzz_0, tly_yyyy_yzzz_0, tly_yyyy_zzzz_0, tly_yyyz_xxxx_0, tly_yyyz_xxxy_0, \
                                     tly_yyyz_xxxz_0, tly_yyyz_xxyy_0, tly_yyyz_xxyz_0, tly_yyyz_xxzz_0, tly_yyyz_xyyy_0, \
                                     tly_yyyz_xyyz_0, tly_yyyz_xyzz_0, tly_yyyz_xzzz_0, tly_yyyz_yyyy_0, tly_yyyz_yyyz_0, \
                                     tly_yyz_xxx_0, tly_yyz_xxxx_0, tly_yyz_xxxy_0, tly_yyz_xxxz_0, tly_yyz_xxy_0, \
                                     tly_yyz_xxyy_0, tly_yyz_xxyz_0, tly_yyz_xxz_0, tly_yyz_xxzz_0, tly_yyz_xyy_0, \
                                     tly_yyz_xyyy_0, tly_yyz_xyyz_0, tly_yyz_xyz_0, tly_yyz_xyzz_0, tly_yyz_xzz_0, \
                                     tly_yyz_xzzz_0, tly_yyz_yyy_0, tly_yyz_yyyy_0, tly_yyz_yyyz_0, tly_yyz_yyz_0, \
                                     tly_yz_xxxx_0, tly_yz_xxxy_0, tly_yz_xxxz_0, tly_yz_xxyy_0, tly_yz_xxyz_0, \
                                     tly_yz_xxzz_0, tly_yz_xyyy_0, tly_yz_xyyz_0, tly_yz_xyzz_0, tly_yz_xzzz_0, \
                                     tly_yz_yyyy_0, tly_yz_yyyz_0, tlz_yy_yyyz_0, tlz_yy_yyzz_0, tlz_yy_yzzz_0, \
                                     tlz_yy_zzzz_0, tlz_yyy_yyyz_0, tlz_yyy_yyz_0, tlz_yyy_yyzz_0, tlz_yyy_yzz_0, \
                                     tlz_yyy_yzzz_0, tlz_yyy_zzz_0, tlz_yyy_zzzz_0, tlz_yyyy_yyyz_0, tlz_yyyy_yyzz_0, \
                                     tlz_yyyy_yzzz_0, tlz_yyyy_zzzz_0, tlz_yyyz_xxxx_0, tlz_yyyz_xxxy_0, tlz_yyyz_xxxz_0, \
                                     tlz_yyyz_xxyy_0, tlz_yyyz_xxyz_0, tlz_yyyz_xxzz_0, tlz_yyyz_xyyy_0, tlz_yyyz_xyyz_0, \
                                     tlz_yyyz_xyzz_0, tlz_yyyz_xzzz_0, tlz_yyyz_yyyy_0, tlz_yyyz_yyyz_0, tlz_yyz_xxx_0, \
                                     tlz_yyz_xxxx_0, tlz_yyz_xxxy_0, tlz_yyz_xxxz_0, tlz_yyz_xxy_0, tlz_yyz_xxyy_0, \
                                     tlz_yyz_xxyz_0, tlz_yyz_xxz_0, tlz_yyz_xxzz_0, tlz_yyz_xyy_0, tlz_yyz_xyyy_0, \
                                     tlz_yyz_xyyz_0, tlz_yyz_xyz_0, tlz_yyz_xyzz_0, tlz_yyz_xzz_0, tlz_yyz_xzzz_0, \
                                     tlz_yyz_yyy_0, tlz_yyz_yyyy_0, tlz_yyz_yyyz_0, tlz_yyz_yyz_0, tlz_yz_xxxx_0, \
                                     tlz_yz_xxxy_0, tlz_yz_xxxz_0, tlz_yz_xxyy_0, tlz_yz_xxyz_0, tlz_yz_xxzz_0, \
                                     tlz_yz_xyyy_0, tlz_yz_xyyz_0, tlz_yz_xyzz_0, tlz_yz_xzzz_0, tlz_yz_yyyy_0, \
                                     tlz_yz_yyyz_0, tpx_yyy_yyyz_0, tpx_yyy_yyzz_0, tpx_yyy_yzzz_0, tpx_yyy_zzzz_0, \
                                     tpx_yyz_xxxx_0, tpx_yyz_xxxy_0, tpx_yyz_xxxz_0, tpx_yyz_xxyy_0, tpx_yyz_xxyz_0, \
                                     tpx_yyz_xxzz_0, tpx_yyz_xyyy_0, tpx_yyz_xyyz_0, tpx_yyz_xyzz_0, tpx_yyz_xzzz_0, \
                                     tpx_yyz_yyyy_0, tpx_yyz_yyyz_0, tpz_yyy_yyyz_0, tpz_yyy_yyzz_0, tpz_yyy_yzzz_0, \
                                     tpz_yyy_zzzz_0, tpz_yyz_xxxx_0, tpz_yyz_xxxy_0, tpz_yyz_xxxz_0, tpz_yyz_xxyy_0, \
                                     tpz_yyz_xxyz_0, tpz_yyz_xxzz_0, tpz_yyz_xyyy_0, tpz_yyz_xyyz_0, tpz_yyz_xyzz_0, \
                                     tpz_yyz_xzzz_0, tpz_yyz_yyyy_0, tpz_yyz_yyyz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yyyy_yyyz_0[j] = pa_y[j] * tlx_yyy_yyyz_0[j] + 1.5 * fl1_fx * tlx_yy_yyyz_0[j] + 1.5 * fl1_fx * tlx_yyy_yyz_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yyyz_0[j];

            tly_yyyy_yyyz_0[j] = pa_y[j] * tly_yyy_yyyz_0[j] + 1.5 * fl1_fx * tly_yy_yyyz_0[j] + 1.5 * fl1_fx * tly_yyy_yyz_0[j];

            tlz_yyyy_yyyz_0[j] = pa_y[j] * tlz_yyy_yyyz_0[j] + 1.5 * fl1_fx * tlz_yy_yyyz_0[j] + 1.5 * fl1_fx * tlz_yyy_yyz_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yyyz_0[j];

            tlx_yyyy_yyzz_0[j] = pa_y[j] * tlx_yyy_yyzz_0[j] + 1.5 * fl1_fx * tlx_yy_yyzz_0[j] + fl1_fx * tlx_yyy_yzz_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yyzz_0[j];

            tly_yyyy_yyzz_0[j] = pa_y[j] * tly_yyy_yyzz_0[j] + 1.5 * fl1_fx * tly_yy_yyzz_0[j] + fl1_fx * tly_yyy_yzz_0[j];

            tlz_yyyy_yyzz_0[j] = pa_y[j] * tlz_yyy_yyzz_0[j] + 1.5 * fl1_fx * tlz_yy_yyzz_0[j] + fl1_fx * tlz_yyy_yzz_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yyzz_0[j];

            tlx_yyyy_yzzz_0[j] = pa_y[j] * tlx_yyy_yzzz_0[j] + 1.5 * fl1_fx * tlx_yy_yzzz_0[j] + 0.5 * fl1_fx * tlx_yyy_zzz_0[j] -
                                 0.5 * fl1_fx * tpz_yyy_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_yyy_yzzz_0[j];

            tly_yyyy_yzzz_0[j] = pa_y[j] * tly_yyy_yzzz_0[j] + 1.5 * fl1_fx * tly_yy_yzzz_0[j] + 0.5 * fl1_fx * tly_yyy_zzz_0[j];

            tlz_yyyy_yzzz_0[j] = pa_y[j] * tlz_yyy_yzzz_0[j] + 1.5 * fl1_fx * tlz_yy_yzzz_0[j] + 0.5 * fl1_fx * tlz_yyy_zzz_0[j] +
                                 0.5 * fl1_fx * tpx_yyy_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_yyy_yzzz_0[j];

            tlx_yyyy_zzzz_0[j] = pa_y[j] * tlx_yyy_zzzz_0[j] + 1.5 * fl1_fx * tlx_yy_zzzz_0[j] - 0.5 * fl1_fx * tpz_yyy_zzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yyy_zzzz_0[j];

            tly_yyyy_zzzz_0[j] = pa_y[j] * tly_yyy_zzzz_0[j] + 1.5 * fl1_fx * tly_yy_zzzz_0[j];

            tlz_yyyy_zzzz_0[j] = pa_y[j] * tlz_yyy_zzzz_0[j] + 1.5 * fl1_fx * tlz_yy_zzzz_0[j] + 0.5 * fl1_fx * tpx_yyy_zzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yyy_zzzz_0[j];

            tlx_yyyz_xxxx_0[j] =
                pa_y[j] * tlx_yyz_xxxx_0[j] + fl1_fx * tlx_yz_xxxx_0[j] - 0.5 * fl1_fx * tpz_yyz_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxxx_0[j];

            tly_yyyz_xxxx_0[j] = pa_y[j] * tly_yyz_xxxx_0[j] + fl1_fx * tly_yz_xxxx_0[j];

            tlz_yyyz_xxxx_0[j] =
                pa_y[j] * tlz_yyz_xxxx_0[j] + fl1_fx * tlz_yz_xxxx_0[j] + 0.5 * fl1_fx * tpx_yyz_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxxx_0[j];

            tlx_yyyz_xxxy_0[j] = pa_y[j] * tlx_yyz_xxxy_0[j] + fl1_fx * tlx_yz_xxxy_0[j] + 0.5 * fl1_fx * tlx_yyz_xxx_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxxy_0[j];

            tly_yyyz_xxxy_0[j] = pa_y[j] * tly_yyz_xxxy_0[j] + fl1_fx * tly_yz_xxxy_0[j] + 0.5 * fl1_fx * tly_yyz_xxx_0[j];

            tlz_yyyz_xxxy_0[j] = pa_y[j] * tlz_yyz_xxxy_0[j] + fl1_fx * tlz_yz_xxxy_0[j] + 0.5 * fl1_fx * tlz_yyz_xxx_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxxy_0[j];

            tlx_yyyz_xxxz_0[j] =
                pa_y[j] * tlx_yyz_xxxz_0[j] + fl1_fx * tlx_yz_xxxz_0[j] - 0.5 * fl1_fx * tpz_yyz_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxxz_0[j];

            tly_yyyz_xxxz_0[j] = pa_y[j] * tly_yyz_xxxz_0[j] + fl1_fx * tly_yz_xxxz_0[j];

            tlz_yyyz_xxxz_0[j] =
                pa_y[j] * tlz_yyz_xxxz_0[j] + fl1_fx * tlz_yz_xxxz_0[j] + 0.5 * fl1_fx * tpx_yyz_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxxz_0[j];

            tlx_yyyz_xxyy_0[j] = pa_y[j] * tlx_yyz_xxyy_0[j] + fl1_fx * tlx_yz_xxyy_0[j] + fl1_fx * tlx_yyz_xxy_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxyy_0[j];

            tly_yyyz_xxyy_0[j] = pa_y[j] * tly_yyz_xxyy_0[j] + fl1_fx * tly_yz_xxyy_0[j] + fl1_fx * tly_yyz_xxy_0[j];

            tlz_yyyz_xxyy_0[j] = pa_y[j] * tlz_yyz_xxyy_0[j] + fl1_fx * tlz_yz_xxyy_0[j] + fl1_fx * tlz_yyz_xxy_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxyy_0[j];

            tlx_yyyz_xxyz_0[j] = pa_y[j] * tlx_yyz_xxyz_0[j] + fl1_fx * tlx_yz_xxyz_0[j] + 0.5 * fl1_fx * tlx_yyz_xxz_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxyz_0[j];

            tly_yyyz_xxyz_0[j] = pa_y[j] * tly_yyz_xxyz_0[j] + fl1_fx * tly_yz_xxyz_0[j] + 0.5 * fl1_fx * tly_yyz_xxz_0[j];

            tlz_yyyz_xxyz_0[j] = pa_y[j] * tlz_yyz_xxyz_0[j] + fl1_fx * tlz_yz_xxyz_0[j] + 0.5 * fl1_fx * tlz_yyz_xxz_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxyz_0[j];

            tlx_yyyz_xxzz_0[j] =
                pa_y[j] * tlx_yyz_xxzz_0[j] + fl1_fx * tlx_yz_xxzz_0[j] - 0.5 * fl1_fx * tpz_yyz_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xxzz_0[j];

            tly_yyyz_xxzz_0[j] = pa_y[j] * tly_yyz_xxzz_0[j] + fl1_fx * tly_yz_xxzz_0[j];

            tlz_yyyz_xxzz_0[j] =
                pa_y[j] * tlz_yyz_xxzz_0[j] + fl1_fx * tlz_yz_xxzz_0[j] + 0.5 * fl1_fx * tpx_yyz_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xxzz_0[j];

            tlx_yyyz_xyyy_0[j] = pa_y[j] * tlx_yyz_xyyy_0[j] + fl1_fx * tlx_yz_xyyy_0[j] + 1.5 * fl1_fx * tlx_yyz_xyy_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xyyy_0[j];

            tly_yyyz_xyyy_0[j] = pa_y[j] * tly_yyz_xyyy_0[j] + fl1_fx * tly_yz_xyyy_0[j] + 1.5 * fl1_fx * tly_yyz_xyy_0[j];

            tlz_yyyz_xyyy_0[j] = pa_y[j] * tlz_yyz_xyyy_0[j] + fl1_fx * tlz_yz_xyyy_0[j] + 1.5 * fl1_fx * tlz_yyz_xyy_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xyyy_0[j];

            tlx_yyyz_xyyz_0[j] = pa_y[j] * tlx_yyz_xyyz_0[j] + fl1_fx * tlx_yz_xyyz_0[j] + fl1_fx * tlx_yyz_xyz_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xyyz_0[j];

            tly_yyyz_xyyz_0[j] = pa_y[j] * tly_yyz_xyyz_0[j] + fl1_fx * tly_yz_xyyz_0[j] + fl1_fx * tly_yyz_xyz_0[j];

            tlz_yyyz_xyyz_0[j] = pa_y[j] * tlz_yyz_xyyz_0[j] + fl1_fx * tlz_yz_xyyz_0[j] + fl1_fx * tlz_yyz_xyz_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xyyz_0[j];

            tlx_yyyz_xyzz_0[j] = pa_y[j] * tlx_yyz_xyzz_0[j] + fl1_fx * tlx_yz_xyzz_0[j] + 0.5 * fl1_fx * tlx_yyz_xzz_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xyzz_0[j];

            tly_yyyz_xyzz_0[j] = pa_y[j] * tly_yyz_xyzz_0[j] + fl1_fx * tly_yz_xyzz_0[j] + 0.5 * fl1_fx * tly_yyz_xzz_0[j];

            tlz_yyyz_xyzz_0[j] = pa_y[j] * tlz_yyz_xyzz_0[j] + fl1_fx * tlz_yz_xyzz_0[j] + 0.5 * fl1_fx * tlz_yyz_xzz_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xyzz_0[j];

            tlx_yyyz_xzzz_0[j] =
                pa_y[j] * tlx_yyz_xzzz_0[j] + fl1_fx * tlx_yz_xzzz_0[j] - 0.5 * fl1_fx * tpz_yyz_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_xzzz_0[j];

            tly_yyyz_xzzz_0[j] = pa_y[j] * tly_yyz_xzzz_0[j] + fl1_fx * tly_yz_xzzz_0[j];

            tlz_yyyz_xzzz_0[j] =
                pa_y[j] * tlz_yyz_xzzz_0[j] + fl1_fx * tlz_yz_xzzz_0[j] + 0.5 * fl1_fx * tpx_yyz_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_xzzz_0[j];

            tlx_yyyz_yyyy_0[j] = pa_y[j] * tlx_yyz_yyyy_0[j] + fl1_fx * tlx_yz_yyyy_0[j] + 2.0 * fl1_fx * tlx_yyz_yyy_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yyyy_0[j];

            tly_yyyz_yyyy_0[j] = pa_y[j] * tly_yyz_yyyy_0[j] + fl1_fx * tly_yz_yyyy_0[j] + 2.0 * fl1_fx * tly_yyz_yyy_0[j];

            tlz_yyyz_yyyy_0[j] = pa_y[j] * tlz_yyz_yyyy_0[j] + fl1_fx * tlz_yz_yyyy_0[j] + 2.0 * fl1_fx * tlz_yyz_yyy_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yyyy_0[j];

            tlx_yyyz_yyyz_0[j] = pa_y[j] * tlx_yyz_yyyz_0[j] + fl1_fx * tlx_yz_yyyz_0[j] + 1.5 * fl1_fx * tlx_yyz_yyz_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yyyz_0[j];

            tly_yyyz_yyyz_0[j] = pa_y[j] * tly_yyz_yyyz_0[j] + fl1_fx * tly_yz_yyyz_0[j] + 1.5 * fl1_fx * tly_yyz_yyz_0[j];

            tlz_yyyz_yyyz_0[j] = pa_y[j] * tlz_yyz_yyyz_0[j] + fl1_fx * tlz_yz_yyyz_0[j] + 1.5 * fl1_fx * tlz_yyz_yyz_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yyyz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_531_579(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 117);

        auto tly_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 117);

        auto tlz_yyz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tlx_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 118);

        auto tly_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 118);

        auto tlz_yyz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tlx_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 119);

        auto tly_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 119);

        auto tlz_yyz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tlx_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 120);

        auto tly_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 120);

        auto tlz_yzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tlx_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 121);

        auto tly_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 121);

        auto tlz_yzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tlx_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 122);

        auto tly_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 122);

        auto tlz_yzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tlx_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 123);

        auto tly_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 123);

        auto tlz_yzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tlx_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 124);

        auto tly_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 124);

        auto tlz_yzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tlx_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 125);

        auto tly_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 125);

        auto tlz_yzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tlx_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 126);

        auto tly_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 126);

        auto tlz_yzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tlx_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 127);

        auto tly_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 127);

        auto tlz_yzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tlx_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 128);

        auto tly_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 128);

        auto tlz_yzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tlx_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 129);

        auto tly_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 129);

        auto tlz_yzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tlx_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 130);

        auto tly_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 130);

        auto tlz_yzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tlx_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 131);

        auto tly_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 131);

        auto tlz_yzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tlx_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 132);

        auto tly_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 132);

        auto tlz_yzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tlx_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 72);

        auto tly_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 72);

        auto tlz_yz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 72);

        auto tlx_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 73);

        auto tly_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 73);

        auto tlz_yz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 73);

        auto tlx_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 74);

        auto tly_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 74);

        auto tlz_yz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 74);

        auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75);

        auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76);

        auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77);

        auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78);

        auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79);

        auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80);

        auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81);

        auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82);

        auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83);

        auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84);

        auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85);

        auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86);

        auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87);

        auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tlx_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 78);

        auto tly_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 78);

        auto tlz_yyz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 78);

        auto tlx_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 79);

        auto tly_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 79);

        auto tlz_yyz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 79);

        auto tlx_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 80);

        auto tly_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 80);

        auto tlz_yzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 80);

        auto tlx_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 81);

        auto tly_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 81);

        auto tlz_yzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 81);

        auto tlx_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 82);

        auto tly_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 82);

        auto tlz_yzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 82);

        auto tlx_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 83);

        auto tly_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 83);

        auto tlz_yzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 83);

        auto tlx_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 84);

        auto tly_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 84);

        auto tlz_yzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 84);

        auto tlx_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 85);

        auto tly_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 85);

        auto tlz_yzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 85);

        auto tlx_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 86);

        auto tly_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 86);

        auto tlz_yzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 86);

        auto tlx_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 87);

        auto tly_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 87);

        auto tlz_yzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 87);

        auto tlx_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 88);

        auto tly_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 88);

        auto tlz_yzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 88);

        auto tpx_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 117);

        auto tpz_yyz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tpx_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 118);

        auto tpz_yyz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tpx_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 119);

        auto tpz_yyz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tpx_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 120);

        auto tpz_yzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tpx_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 121);

        auto tpz_yzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tpx_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 122);

        auto tpz_yzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tpx_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 123);

        auto tpz_yzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tpx_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 124);

        auto tpz_yzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tpx_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 125);

        auto tpz_yzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tpx_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 126);

        auto tpz_yzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tpx_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 127);

        auto tpz_yzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tpx_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 128);

        auto tpz_yzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tpx_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 129);

        auto tpz_yzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tpx_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 130);

        auto tpz_yzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tpx_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 131);

        auto tpz_yzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tpx_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 132);

        auto tpz_yzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 132);

        auto tdx_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 117);

        auto tdz_yyz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 117);

        auto tdx_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 118);

        auto tdz_yyz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 118);

        auto tdx_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 119);

        auto tdz_yyz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 119);

        auto tdx_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 120);

        auto tdz_yzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 120);

        auto tdx_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 121);

        auto tdz_yzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 121);

        auto tdx_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 122);

        auto tdz_yzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 122);

        auto tdx_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 123);

        auto tdz_yzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 123);

        auto tdx_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 124);

        auto tdz_yzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 124);

        auto tdx_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 125);

        auto tdz_yzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 125);

        auto tdx_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 126);

        auto tdz_yzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 126);

        auto tdx_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 127);

        auto tdz_yzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 127);

        auto tdx_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 128);

        auto tdz_yzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 128);

        auto tdx_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 129);

        auto tdz_yzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 129);

        auto tdx_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 130);

        auto tdz_yzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 130);

        auto tdx_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 131);

        auto tdz_yzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 131);

        auto tdx_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 132);

        auto tdz_yzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 132);

        // set up pointers to integrals

        auto tlx_yyyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 177);

        auto tly_yyyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 177);

        auto tlz_yyyz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 177);

        auto tlx_yyyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 178);

        auto tly_yyyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 178);

        auto tlz_yyyz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 178);

        auto tlx_yyyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 179);

        auto tly_yyyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 179);

        auto tlz_yyyz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 179);

        auto tlx_yyzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 180);

        auto tly_yyzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 180);

        auto tlz_yyzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 180);

        auto tlx_yyzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 181);

        auto tly_yyzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 181);

        auto tlz_yyzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 181);

        auto tlx_yyzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 182);

        auto tly_yyzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 182);

        auto tlz_yyzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 182);

        auto tlx_yyzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 183);

        auto tly_yyzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 183);

        auto tlz_yyzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 183);

        auto tlx_yyzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 184);

        auto tly_yyzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 184);

        auto tlz_yyzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 184);

        auto tlx_yyzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 185);

        auto tly_yyzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 185);

        auto tlz_yyzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 185);

        auto tlx_yyzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 186);

        auto tly_yyzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 186);

        auto tlz_yyzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 186);

        auto tlx_yyzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 187);

        auto tly_yyzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 187);

        auto tlz_yyzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 187);

        auto tlx_yyzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 188);

        auto tly_yyzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 188);

        auto tlz_yyzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 188);

        auto tlx_yyzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 189);

        auto tly_yyzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 189);

        auto tlz_yyzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 189);

        auto tlx_yyzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 190);

        auto tly_yyzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 190);

        auto tlz_yyzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 190);

        auto tlx_yyzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 191);

        auto tly_yyzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 191);

        auto tlz_yyzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 191);

        auto tlx_yyzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 192);

        auto tly_yyzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 192);

        auto tlz_yyzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 192);

        // Batch of Integrals (531,579)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yyz_yyzz_0, tdx_yyz_yzzz_0, tdx_yyz_zzzz_0, \
                                     tdx_yzz_xxxx_0, tdx_yzz_xxxy_0, tdx_yzz_xxxz_0, tdx_yzz_xxyy_0, tdx_yzz_xxyz_0, \
                                     tdx_yzz_xxzz_0, tdx_yzz_xyyy_0, tdx_yzz_xyyz_0, tdx_yzz_xyzz_0, tdx_yzz_xzzz_0, \
                                     tdx_yzz_yyyy_0, tdx_yzz_yyyz_0, tdx_yzz_yyzz_0, tdz_yyz_yyzz_0, tdz_yyz_yzzz_0, \
                                     tdz_yyz_zzzz_0, tdz_yzz_xxxx_0, tdz_yzz_xxxy_0, tdz_yzz_xxxz_0, tdz_yzz_xxyy_0, \
                                     tdz_yzz_xxyz_0, tdz_yzz_xxzz_0, tdz_yzz_xyyy_0, tdz_yzz_xyyz_0, tdz_yzz_xyzz_0, \
                                     tdz_yzz_xzzz_0, tdz_yzz_yyyy_0, tdz_yzz_yyyz_0, tdz_yzz_yyzz_0, tlx_yyyz_yyzz_0, \
                                     tlx_yyyz_yzzz_0, tlx_yyyz_zzzz_0, tlx_yyz_yyzz_0, tlx_yyz_yzz_0, tlx_yyz_yzzz_0, \
                                     tlx_yyz_zzz_0, tlx_yyz_zzzz_0, tlx_yyzz_xxxx_0, tlx_yyzz_xxxy_0, tlx_yyzz_xxxz_0, \
                                     tlx_yyzz_xxyy_0, tlx_yyzz_xxyz_0, tlx_yyzz_xxzz_0, tlx_yyzz_xyyy_0, tlx_yyzz_xyyz_0, \
                                     tlx_yyzz_xyzz_0, tlx_yyzz_xzzz_0, tlx_yyzz_yyyy_0, tlx_yyzz_yyyz_0, tlx_yyzz_yyzz_0, \
                                     tlx_yz_yyzz_0, tlx_yz_yzzz_0, tlx_yz_zzzz_0, tlx_yzz_xxx_0, tlx_yzz_xxxx_0, \
                                     tlx_yzz_xxxy_0, tlx_yzz_xxxz_0, tlx_yzz_xxy_0, tlx_yzz_xxyy_0, tlx_yzz_xxyz_0, \
                                     tlx_yzz_xxz_0, tlx_yzz_xxzz_0, tlx_yzz_xyy_0, tlx_yzz_xyyy_0, tlx_yzz_xyyz_0, \
                                     tlx_yzz_xyz_0, tlx_yzz_xyzz_0, tlx_yzz_xzz_0, tlx_yzz_xzzz_0, tlx_yzz_yyy_0, \
                                     tlx_yzz_yyyy_0, tlx_yzz_yyyz_0, tlx_yzz_yyz_0, tlx_yzz_yyzz_0, tlx_yzz_yzz_0, \
                                     tlx_zz_xxxx_0, tlx_zz_xxxy_0, tlx_zz_xxxz_0, tlx_zz_xxyy_0, tlx_zz_xxyz_0, \
                                     tlx_zz_xxzz_0, tlx_zz_xyyy_0, tlx_zz_xyyz_0, tlx_zz_xyzz_0, tlx_zz_xzzz_0, \
                                     tlx_zz_yyyy_0, tlx_zz_yyyz_0, tlx_zz_yyzz_0, tly_yyyz_yyzz_0, tly_yyyz_yzzz_0, \
                                     tly_yyyz_zzzz_0, tly_yyz_yyzz_0, tly_yyz_yzz_0, tly_yyz_yzzz_0, tly_yyz_zzz_0, \
                                     tly_yyz_zzzz_0, tly_yyzz_xxxx_0, tly_yyzz_xxxy_0, tly_yyzz_xxxz_0, tly_yyzz_xxyy_0, \
                                     tly_yyzz_xxyz_0, tly_yyzz_xxzz_0, tly_yyzz_xyyy_0, tly_yyzz_xyyz_0, tly_yyzz_xyzz_0, \
                                     tly_yyzz_xzzz_0, tly_yyzz_yyyy_0, tly_yyzz_yyyz_0, tly_yyzz_yyzz_0, tly_yz_yyzz_0, \
                                     tly_yz_yzzz_0, tly_yz_zzzz_0, tly_yzz_xxx_0, tly_yzz_xxxx_0, tly_yzz_xxxy_0, \
                                     tly_yzz_xxxz_0, tly_yzz_xxy_0, tly_yzz_xxyy_0, tly_yzz_xxyz_0, tly_yzz_xxz_0, \
                                     tly_yzz_xxzz_0, tly_yzz_xyy_0, tly_yzz_xyyy_0, tly_yzz_xyyz_0, tly_yzz_xyz_0, \
                                     tly_yzz_xyzz_0, tly_yzz_xzz_0, tly_yzz_xzzz_0, tly_yzz_yyy_0, tly_yzz_yyyy_0, \
                                     tly_yzz_yyyz_0, tly_yzz_yyz_0, tly_yzz_yyzz_0, tly_yzz_yzz_0, tly_zz_xxxx_0, \
                                     tly_zz_xxxy_0, tly_zz_xxxz_0, tly_zz_xxyy_0, tly_zz_xxyz_0, tly_zz_xxzz_0, \
                                     tly_zz_xyyy_0, tly_zz_xyyz_0, tly_zz_xyzz_0, tly_zz_xzzz_0, tly_zz_yyyy_0, \
                                     tly_zz_yyyz_0, tly_zz_yyzz_0, tlz_yyyz_yyzz_0, tlz_yyyz_yzzz_0, tlz_yyyz_zzzz_0, \
                                     tlz_yyz_yyzz_0, tlz_yyz_yzz_0, tlz_yyz_yzzz_0, tlz_yyz_zzz_0, tlz_yyz_zzzz_0, \
                                     tlz_yyzz_xxxx_0, tlz_yyzz_xxxy_0, tlz_yyzz_xxxz_0, tlz_yyzz_xxyy_0, tlz_yyzz_xxyz_0, \
                                     tlz_yyzz_xxzz_0, tlz_yyzz_xyyy_0, tlz_yyzz_xyyz_0, tlz_yyzz_xyzz_0, tlz_yyzz_xzzz_0, \
                                     tlz_yyzz_yyyy_0, tlz_yyzz_yyyz_0, tlz_yyzz_yyzz_0, tlz_yz_yyzz_0, tlz_yz_yzzz_0, \
                                     tlz_yz_zzzz_0, tlz_yzz_xxx_0, tlz_yzz_xxxx_0, tlz_yzz_xxxy_0, tlz_yzz_xxxz_0, \
                                     tlz_yzz_xxy_0, tlz_yzz_xxyy_0, tlz_yzz_xxyz_0, tlz_yzz_xxz_0, tlz_yzz_xxzz_0, \
                                     tlz_yzz_xyy_0, tlz_yzz_xyyy_0, tlz_yzz_xyyz_0, tlz_yzz_xyz_0, tlz_yzz_xyzz_0, \
                                     tlz_yzz_xzz_0, tlz_yzz_xzzz_0, tlz_yzz_yyy_0, tlz_yzz_yyyy_0, tlz_yzz_yyyz_0, \
                                     tlz_yzz_yyz_0, tlz_yzz_yyzz_0, tlz_yzz_yzz_0, tlz_zz_xxxx_0, tlz_zz_xxxy_0, \
                                     tlz_zz_xxxz_0, tlz_zz_xxyy_0, tlz_zz_xxyz_0, tlz_zz_xxzz_0, tlz_zz_xyyy_0, \
                                     tlz_zz_xyyz_0, tlz_zz_xyzz_0, tlz_zz_xzzz_0, tlz_zz_yyyy_0, tlz_zz_yyyz_0, \
                                     tlz_zz_yyzz_0, tpx_yyz_yyzz_0, tpx_yyz_yzzz_0, tpx_yyz_zzzz_0, tpx_yzz_xxxx_0, \
                                     tpx_yzz_xxxy_0, tpx_yzz_xxxz_0, tpx_yzz_xxyy_0, tpx_yzz_xxyz_0, tpx_yzz_xxzz_0, \
                                     tpx_yzz_xyyy_0, tpx_yzz_xyyz_0, tpx_yzz_xyzz_0, tpx_yzz_xzzz_0, tpx_yzz_yyyy_0, \
                                     tpx_yzz_yyyz_0, tpx_yzz_yyzz_0, tpz_yyz_yyzz_0, tpz_yyz_yzzz_0, tpz_yyz_zzzz_0, \
                                     tpz_yzz_xxxx_0, tpz_yzz_xxxy_0, tpz_yzz_xxxz_0, tpz_yzz_xxyy_0, tpz_yzz_xxyz_0, \
                                     tpz_yzz_xxzz_0, tpz_yzz_xyyy_0, tpz_yzz_xyyz_0, tpz_yzz_xyzz_0, tpz_yzz_xzzz_0, \
                                     tpz_yzz_yyyy_0, tpz_yzz_yyyz_0, tpz_yzz_yyzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yyyz_yyzz_0[j] = pa_y[j] * tlx_yyz_yyzz_0[j] + fl1_fx * tlx_yz_yyzz_0[j] + fl1_fx * tlx_yyz_yzz_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yyzz_0[j];

            tly_yyyz_yyzz_0[j] = pa_y[j] * tly_yyz_yyzz_0[j] + fl1_fx * tly_yz_yyzz_0[j] + fl1_fx * tly_yyz_yzz_0[j];

            tlz_yyyz_yyzz_0[j] = pa_y[j] * tlz_yyz_yyzz_0[j] + fl1_fx * tlz_yz_yyzz_0[j] + fl1_fx * tlz_yyz_yzz_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yyzz_0[j];

            tlx_yyyz_yzzz_0[j] = pa_y[j] * tlx_yyz_yzzz_0[j] + fl1_fx * tlx_yz_yzzz_0[j] + 0.5 * fl1_fx * tlx_yyz_zzz_0[j] -
                                 0.5 * fl1_fx * tpz_yyz_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_yzzz_0[j];

            tly_yyyz_yzzz_0[j] = pa_y[j] * tly_yyz_yzzz_0[j] + fl1_fx * tly_yz_yzzz_0[j] + 0.5 * fl1_fx * tly_yyz_zzz_0[j];

            tlz_yyyz_yzzz_0[j] = pa_y[j] * tlz_yyz_yzzz_0[j] + fl1_fx * tlz_yz_yzzz_0[j] + 0.5 * fl1_fx * tlz_yyz_zzz_0[j] +
                                 0.5 * fl1_fx * tpx_yyz_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_yzzz_0[j];

            tlx_yyyz_zzzz_0[j] =
                pa_y[j] * tlx_yyz_zzzz_0[j] + fl1_fx * tlx_yz_zzzz_0[j] - 0.5 * fl1_fx * tpz_yyz_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_yyz_zzzz_0[j];

            tly_yyyz_zzzz_0[j] = pa_y[j] * tly_yyz_zzzz_0[j] + fl1_fx * tly_yz_zzzz_0[j];

            tlz_yyyz_zzzz_0[j] =
                pa_y[j] * tlz_yyz_zzzz_0[j] + fl1_fx * tlz_yz_zzzz_0[j] + 0.5 * fl1_fx * tpx_yyz_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_yyz_zzzz_0[j];

            tlx_yyzz_xxxx_0[j] = pa_y[j] * tlx_yzz_xxxx_0[j] + 0.5 * fl1_fx * tlx_zz_xxxx_0[j] - 0.5 * fl1_fx * tpz_yzz_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yzz_xxxx_0[j];

            tly_yyzz_xxxx_0[j] = pa_y[j] * tly_yzz_xxxx_0[j] + 0.5 * fl1_fx * tly_zz_xxxx_0[j];

            tlz_yyzz_xxxx_0[j] = pa_y[j] * tlz_yzz_xxxx_0[j] + 0.5 * fl1_fx * tlz_zz_xxxx_0[j] + 0.5 * fl1_fx * tpx_yzz_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yzz_xxxx_0[j];

            tlx_yyzz_xxxy_0[j] = pa_y[j] * tlx_yzz_xxxy_0[j] + 0.5 * fl1_fx * tlx_zz_xxxy_0[j] + 0.5 * fl1_fx * tlx_yzz_xxx_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_xxxy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xxxy_0[j];

            tly_yyzz_xxxy_0[j] = pa_y[j] * tly_yzz_xxxy_0[j] + 0.5 * fl1_fx * tly_zz_xxxy_0[j] + 0.5 * fl1_fx * tly_yzz_xxx_0[j];

            tlz_yyzz_xxxy_0[j] = pa_y[j] * tlz_yzz_xxxy_0[j] + 0.5 * fl1_fx * tlz_zz_xxxy_0[j] + 0.5 * fl1_fx * tlz_yzz_xxx_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_xxxy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xxxy_0[j];

            tlx_yyzz_xxxz_0[j] = pa_y[j] * tlx_yzz_xxxz_0[j] + 0.5 * fl1_fx * tlx_zz_xxxz_0[j] - 0.5 * fl1_fx * tpz_yzz_xxxz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yzz_xxxz_0[j];

            tly_yyzz_xxxz_0[j] = pa_y[j] * tly_yzz_xxxz_0[j] + 0.5 * fl1_fx * tly_zz_xxxz_0[j];

            tlz_yyzz_xxxz_0[j] = pa_y[j] * tlz_yzz_xxxz_0[j] + 0.5 * fl1_fx * tlz_zz_xxxz_0[j] + 0.5 * fl1_fx * tpx_yzz_xxxz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yzz_xxxz_0[j];

            tlx_yyzz_xxyy_0[j] = pa_y[j] * tlx_yzz_xxyy_0[j] + 0.5 * fl1_fx * tlx_zz_xxyy_0[j] + fl1_fx * tlx_yzz_xxy_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xxyy_0[j];

            tly_yyzz_xxyy_0[j] = pa_y[j] * tly_yzz_xxyy_0[j] + 0.5 * fl1_fx * tly_zz_xxyy_0[j] + fl1_fx * tly_yzz_xxy_0[j];

            tlz_yyzz_xxyy_0[j] = pa_y[j] * tlz_yzz_xxyy_0[j] + 0.5 * fl1_fx * tlz_zz_xxyy_0[j] + fl1_fx * tlz_yzz_xxy_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xxyy_0[j];

            tlx_yyzz_xxyz_0[j] = pa_y[j] * tlx_yzz_xxyz_0[j] + 0.5 * fl1_fx * tlx_zz_xxyz_0[j] + 0.5 * fl1_fx * tlx_yzz_xxz_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_xxyz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xxyz_0[j];

            tly_yyzz_xxyz_0[j] = pa_y[j] * tly_yzz_xxyz_0[j] + 0.5 * fl1_fx * tly_zz_xxyz_0[j] + 0.5 * fl1_fx * tly_yzz_xxz_0[j];

            tlz_yyzz_xxyz_0[j] = pa_y[j] * tlz_yzz_xxyz_0[j] + 0.5 * fl1_fx * tlz_zz_xxyz_0[j] + 0.5 * fl1_fx * tlz_yzz_xxz_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_xxyz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xxyz_0[j];

            tlx_yyzz_xxzz_0[j] = pa_y[j] * tlx_yzz_xxzz_0[j] + 0.5 * fl1_fx * tlx_zz_xxzz_0[j] - 0.5 * fl1_fx * tpz_yzz_xxzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yzz_xxzz_0[j];

            tly_yyzz_xxzz_0[j] = pa_y[j] * tly_yzz_xxzz_0[j] + 0.5 * fl1_fx * tly_zz_xxzz_0[j];

            tlz_yyzz_xxzz_0[j] = pa_y[j] * tlz_yzz_xxzz_0[j] + 0.5 * fl1_fx * tlz_zz_xxzz_0[j] + 0.5 * fl1_fx * tpx_yzz_xxzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yzz_xxzz_0[j];

            tlx_yyzz_xyyy_0[j] = pa_y[j] * tlx_yzz_xyyy_0[j] + 0.5 * fl1_fx * tlx_zz_xyyy_0[j] + 1.5 * fl1_fx * tlx_yzz_xyy_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_xyyy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xyyy_0[j];

            tly_yyzz_xyyy_0[j] = pa_y[j] * tly_yzz_xyyy_0[j] + 0.5 * fl1_fx * tly_zz_xyyy_0[j] + 1.5 * fl1_fx * tly_yzz_xyy_0[j];

            tlz_yyzz_xyyy_0[j] = pa_y[j] * tlz_yzz_xyyy_0[j] + 0.5 * fl1_fx * tlz_zz_xyyy_0[j] + 1.5 * fl1_fx * tlz_yzz_xyy_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_xyyy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xyyy_0[j];

            tlx_yyzz_xyyz_0[j] = pa_y[j] * tlx_yzz_xyyz_0[j] + 0.5 * fl1_fx * tlx_zz_xyyz_0[j] + fl1_fx * tlx_yzz_xyz_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xyyz_0[j];

            tly_yyzz_xyyz_0[j] = pa_y[j] * tly_yzz_xyyz_0[j] + 0.5 * fl1_fx * tly_zz_xyyz_0[j] + fl1_fx * tly_yzz_xyz_0[j];

            tlz_yyzz_xyyz_0[j] = pa_y[j] * tlz_yzz_xyyz_0[j] + 0.5 * fl1_fx * tlz_zz_xyyz_0[j] + fl1_fx * tlz_yzz_xyz_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xyyz_0[j];

            tlx_yyzz_xyzz_0[j] = pa_y[j] * tlx_yzz_xyzz_0[j] + 0.5 * fl1_fx * tlx_zz_xyzz_0[j] + 0.5 * fl1_fx * tlx_yzz_xzz_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_xyzz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_xyzz_0[j];

            tly_yyzz_xyzz_0[j] = pa_y[j] * tly_yzz_xyzz_0[j] + 0.5 * fl1_fx * tly_zz_xyzz_0[j] + 0.5 * fl1_fx * tly_yzz_xzz_0[j];

            tlz_yyzz_xyzz_0[j] = pa_y[j] * tlz_yzz_xyzz_0[j] + 0.5 * fl1_fx * tlz_zz_xyzz_0[j] + 0.5 * fl1_fx * tlz_yzz_xzz_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_xyzz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_xyzz_0[j];

            tlx_yyzz_xzzz_0[j] = pa_y[j] * tlx_yzz_xzzz_0[j] + 0.5 * fl1_fx * tlx_zz_xzzz_0[j] - 0.5 * fl1_fx * tpz_yzz_xzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yzz_xzzz_0[j];

            tly_yyzz_xzzz_0[j] = pa_y[j] * tly_yzz_xzzz_0[j] + 0.5 * fl1_fx * tly_zz_xzzz_0[j];

            tlz_yyzz_xzzz_0[j] = pa_y[j] * tlz_yzz_xzzz_0[j] + 0.5 * fl1_fx * tlz_zz_xzzz_0[j] + 0.5 * fl1_fx * tpx_yzz_xzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yzz_xzzz_0[j];

            tlx_yyzz_yyyy_0[j] = pa_y[j] * tlx_yzz_yyyy_0[j] + 0.5 * fl1_fx * tlx_zz_yyyy_0[j] + 2.0 * fl1_fx * tlx_yzz_yyy_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_yyyy_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yyyy_0[j];

            tly_yyzz_yyyy_0[j] = pa_y[j] * tly_yzz_yyyy_0[j] + 0.5 * fl1_fx * tly_zz_yyyy_0[j] + 2.0 * fl1_fx * tly_yzz_yyy_0[j];

            tlz_yyzz_yyyy_0[j] = pa_y[j] * tlz_yzz_yyyy_0[j] + 0.5 * fl1_fx * tlz_zz_yyyy_0[j] + 2.0 * fl1_fx * tlz_yzz_yyy_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_yyyy_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yyyy_0[j];

            tlx_yyzz_yyyz_0[j] = pa_y[j] * tlx_yzz_yyyz_0[j] + 0.5 * fl1_fx * tlx_zz_yyyz_0[j] + 1.5 * fl1_fx * tlx_yzz_yyz_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_yyyz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yyyz_0[j];

            tly_yyzz_yyyz_0[j] = pa_y[j] * tly_yzz_yyyz_0[j] + 0.5 * fl1_fx * tly_zz_yyyz_0[j] + 1.5 * fl1_fx * tly_yzz_yyz_0[j];

            tlz_yyzz_yyyz_0[j] = pa_y[j] * tlz_yzz_yyyz_0[j] + 0.5 * fl1_fx * tlz_zz_yyyz_0[j] + 1.5 * fl1_fx * tlz_yzz_yyz_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_yyyz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yyyz_0[j];

            tlx_yyzz_yyzz_0[j] = pa_y[j] * tlx_yzz_yyzz_0[j] + 0.5 * fl1_fx * tlx_zz_yyzz_0[j] + fl1_fx * tlx_yzz_yzz_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yyzz_0[j];

            tly_yyzz_yyzz_0[j] = pa_y[j] * tly_yzz_yyzz_0[j] + 0.5 * fl1_fx * tly_zz_yyzz_0[j] + fl1_fx * tly_yzz_yzz_0[j];

            tlz_yyzz_yyzz_0[j] = pa_y[j] * tlz_yzz_yyzz_0[j] + 0.5 * fl1_fx * tlz_zz_yyzz_0[j] + fl1_fx * tlz_yzz_yzz_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yyzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_579_627(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 133);

        auto tly_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 133);

        auto tlz_yzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tlx_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 134);

        auto tly_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 134);

        auto tlz_yzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tlx_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 135);

        auto tly_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tlz_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tlx_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 136);

        auto tly_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tlz_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tlx_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 137);

        auto tly_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tlz_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tlx_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 138);

        auto tly_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tlz_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tlx_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 139);

        auto tly_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tlz_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tlx_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 140);

        auto tly_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tlz_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tlx_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 141);

        auto tly_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tlz_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tlx_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 142);

        auto tly_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tlz_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tlx_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 143);

        auto tly_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tlz_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tlx_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 144);

        auto tly_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tlz_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tlx_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 145);

        auto tly_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tlz_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tlx_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 146);

        auto tly_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tlz_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tlx_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 147);

        auto tly_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tlz_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tlx_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 148);

        auto tly_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tlz_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88);

        auto tly_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tlz_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tlx_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 89);

        auto tly_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tlz_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tlx_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 89);

        auto tly_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 89);

        auto tlz_yzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 89);

        auto tlx_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 90);

        auto tly_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tlz_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tlx_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 91);

        auto tly_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tlz_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tlx_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 92);

        auto tly_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tlz_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tlx_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 93);

        auto tly_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tlz_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tlx_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 94);

        auto tly_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tlz_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tlx_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 95);

        auto tly_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tlz_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tlx_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 96);

        auto tly_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tlz_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tlx_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 97);

        auto tly_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tlz_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tlx_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 98);

        auto tly_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tlz_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tlx_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 99);

        auto tly_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tlz_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto tpx_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 133);

        auto tpz_yzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tpx_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 134);

        auto tpz_yzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tpx_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 135);

        auto tpz_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tpx_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 136);

        auto tpz_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tpx_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 137);

        auto tpz_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tpx_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 138);

        auto tpz_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tpx_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 139);

        auto tpz_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tpx_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 140);

        auto tpz_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tpx_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 141);

        auto tpz_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tpx_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 142);

        auto tpz_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tpx_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 143);

        auto tpz_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tpx_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 144);

        auto tpz_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tpx_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 145);

        auto tpz_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tpx_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 146);

        auto tpz_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tpx_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 147);

        auto tpz_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tpx_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 148);

        auto tpz_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tdx_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 133);

        auto tdz_yzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 133);

        auto tdx_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 134);

        auto tdz_yzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 134);

        auto tdx_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 135);

        auto tdz_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tdx_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 136);

        auto tdz_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tdx_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 137);

        auto tdz_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tdx_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 138);

        auto tdz_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tdx_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 139);

        auto tdz_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tdx_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 140);

        auto tdz_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tdx_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 141);

        auto tdz_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tdx_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 142);

        auto tdz_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tdx_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 143);

        auto tdz_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tdx_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 144);

        auto tdz_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tdx_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 145);

        auto tdz_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tdx_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 146);

        auto tdz_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tdx_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 147);

        auto tdz_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tdx_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 148);

        auto tdz_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 148);

        // set up pointers to integrals

        auto tlx_yyzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 193);

        auto tly_yyzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 193);

        auto tlz_yyzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 193);

        auto tlx_yyzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 194);

        auto tly_yyzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 194);

        auto tlz_yyzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 194);

        auto tlx_yzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 195);

        auto tly_yzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 195);

        auto tlz_yzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 195);

        auto tlx_yzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 196);

        auto tly_yzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 196);

        auto tlz_yzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 196);

        auto tlx_yzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 197);

        auto tly_yzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 197);

        auto tlz_yzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 197);

        auto tlx_yzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 198);

        auto tly_yzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 198);

        auto tlz_yzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 198);

        auto tlx_yzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 199);

        auto tly_yzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 199);

        auto tlz_yzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 199);

        auto tlx_yzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 200);

        auto tly_yzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 200);

        auto tlz_yzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 200);

        auto tlx_yzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 201);

        auto tly_yzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 201);

        auto tlz_yzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 201);

        auto tlx_yzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 202);

        auto tly_yzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 202);

        auto tlz_yzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 202);

        auto tlx_yzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 203);

        auto tly_yzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 203);

        auto tlz_yzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 203);

        auto tlx_yzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 204);

        auto tly_yzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 204);

        auto tlz_yzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 204);

        auto tlx_yzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 205);

        auto tly_yzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 205);

        auto tlz_yzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 205);

        auto tlx_yzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 206);

        auto tly_yzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 206);

        auto tlz_yzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 206);

        auto tlx_yzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 207);

        auto tly_yzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 207);

        auto tlz_yzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 207);

        auto tlx_yzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 208);

        auto tly_yzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 208);

        auto tlz_yzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 208);

        // Batch of Integrals (579,627)

        #pragma omp simd aligned(fgb, fx, pa_y, tdx_yzz_yzzz_0, tdx_yzz_zzzz_0, tdx_zzz_xxxx_0, \
                                     tdx_zzz_xxxy_0, tdx_zzz_xxxz_0, tdx_zzz_xxyy_0, tdx_zzz_xxyz_0, tdx_zzz_xxzz_0, \
                                     tdx_zzz_xyyy_0, tdx_zzz_xyyz_0, tdx_zzz_xyzz_0, tdx_zzz_xzzz_0, tdx_zzz_yyyy_0, \
                                     tdx_zzz_yyyz_0, tdx_zzz_yyzz_0, tdx_zzz_yzzz_0, tdz_yzz_yzzz_0, tdz_yzz_zzzz_0, \
                                     tdz_zzz_xxxx_0, tdz_zzz_xxxy_0, tdz_zzz_xxxz_0, tdz_zzz_xxyy_0, tdz_zzz_xxyz_0, \
                                     tdz_zzz_xxzz_0, tdz_zzz_xyyy_0, tdz_zzz_xyyz_0, tdz_zzz_xyzz_0, tdz_zzz_xzzz_0, \
                                     tdz_zzz_yyyy_0, tdz_zzz_yyyz_0, tdz_zzz_yyzz_0, tdz_zzz_yzzz_0, tlx_yyzz_yzzz_0, \
                                     tlx_yyzz_zzzz_0, tlx_yzz_yzzz_0, tlx_yzz_zzz_0, tlx_yzz_zzzz_0, tlx_yzzz_xxxx_0, \
                                     tlx_yzzz_xxxy_0, tlx_yzzz_xxxz_0, tlx_yzzz_xxyy_0, tlx_yzzz_xxyz_0, tlx_yzzz_xxzz_0, \
                                     tlx_yzzz_xyyy_0, tlx_yzzz_xyyz_0, tlx_yzzz_xyzz_0, tlx_yzzz_xzzz_0, tlx_yzzz_yyyy_0, \
                                     tlx_yzzz_yyyz_0, tlx_yzzz_yyzz_0, tlx_yzzz_yzzz_0, tlx_zz_yzzz_0, tlx_zz_zzzz_0, \
                                     tlx_zzz_xxx_0, tlx_zzz_xxxx_0, tlx_zzz_xxxy_0, tlx_zzz_xxxz_0, tlx_zzz_xxy_0, \
                                     tlx_zzz_xxyy_0, tlx_zzz_xxyz_0, tlx_zzz_xxz_0, tlx_zzz_xxzz_0, tlx_zzz_xyy_0, \
                                     tlx_zzz_xyyy_0, tlx_zzz_xyyz_0, tlx_zzz_xyz_0, tlx_zzz_xyzz_0, tlx_zzz_xzz_0, \
                                     tlx_zzz_xzzz_0, tlx_zzz_yyy_0, tlx_zzz_yyyy_0, tlx_zzz_yyyz_0, tlx_zzz_yyz_0, \
                                     tlx_zzz_yyzz_0, tlx_zzz_yzz_0, tlx_zzz_yzzz_0, tlx_zzz_zzz_0, tly_yyzz_yzzz_0, \
                                     tly_yyzz_zzzz_0, tly_yzz_yzzz_0, tly_yzz_zzz_0, tly_yzz_zzzz_0, tly_yzzz_xxxx_0, \
                                     tly_yzzz_xxxy_0, tly_yzzz_xxxz_0, tly_yzzz_xxyy_0, tly_yzzz_xxyz_0, tly_yzzz_xxzz_0, \
                                     tly_yzzz_xyyy_0, tly_yzzz_xyyz_0, tly_yzzz_xyzz_0, tly_yzzz_xzzz_0, tly_yzzz_yyyy_0, \
                                     tly_yzzz_yyyz_0, tly_yzzz_yyzz_0, tly_yzzz_yzzz_0, tly_zz_yzzz_0, tly_zz_zzzz_0, \
                                     tly_zzz_xxx_0, tly_zzz_xxxx_0, tly_zzz_xxxy_0, tly_zzz_xxxz_0, tly_zzz_xxy_0, \
                                     tly_zzz_xxyy_0, tly_zzz_xxyz_0, tly_zzz_xxz_0, tly_zzz_xxzz_0, tly_zzz_xyy_0, \
                                     tly_zzz_xyyy_0, tly_zzz_xyyz_0, tly_zzz_xyz_0, tly_zzz_xyzz_0, tly_zzz_xzz_0, \
                                     tly_zzz_xzzz_0, tly_zzz_yyy_0, tly_zzz_yyyy_0, tly_zzz_yyyz_0, tly_zzz_yyz_0, \
                                     tly_zzz_yyzz_0, tly_zzz_yzz_0, tly_zzz_yzzz_0, tly_zzz_zzz_0, tlz_yyzz_yzzz_0, \
                                     tlz_yyzz_zzzz_0, tlz_yzz_yzzz_0, tlz_yzz_zzz_0, tlz_yzz_zzzz_0, tlz_yzzz_xxxx_0, \
                                     tlz_yzzz_xxxy_0, tlz_yzzz_xxxz_0, tlz_yzzz_xxyy_0, tlz_yzzz_xxyz_0, tlz_yzzz_xxzz_0, \
                                     tlz_yzzz_xyyy_0, tlz_yzzz_xyyz_0, tlz_yzzz_xyzz_0, tlz_yzzz_xzzz_0, tlz_yzzz_yyyy_0, \
                                     tlz_yzzz_yyyz_0, tlz_yzzz_yyzz_0, tlz_yzzz_yzzz_0, tlz_zz_yzzz_0, tlz_zz_zzzz_0, \
                                     tlz_zzz_xxx_0, tlz_zzz_xxxx_0, tlz_zzz_xxxy_0, tlz_zzz_xxxz_0, tlz_zzz_xxy_0, \
                                     tlz_zzz_xxyy_0, tlz_zzz_xxyz_0, tlz_zzz_xxz_0, tlz_zzz_xxzz_0, tlz_zzz_xyy_0, \
                                     tlz_zzz_xyyy_0, tlz_zzz_xyyz_0, tlz_zzz_xyz_0, tlz_zzz_xyzz_0, tlz_zzz_xzz_0, \
                                     tlz_zzz_xzzz_0, tlz_zzz_yyy_0, tlz_zzz_yyyy_0, tlz_zzz_yyyz_0, tlz_zzz_yyz_0, \
                                     tlz_zzz_yyzz_0, tlz_zzz_yzz_0, tlz_zzz_yzzz_0, tlz_zzz_zzz_0, tpx_yzz_yzzz_0, \
                                     tpx_yzz_zzzz_0, tpx_zzz_xxxx_0, tpx_zzz_xxxy_0, tpx_zzz_xxxz_0, tpx_zzz_xxyy_0, \
                                     tpx_zzz_xxyz_0, tpx_zzz_xxzz_0, tpx_zzz_xyyy_0, tpx_zzz_xyyz_0, tpx_zzz_xyzz_0, \
                                     tpx_zzz_xzzz_0, tpx_zzz_yyyy_0, tpx_zzz_yyyz_0, tpx_zzz_yyzz_0, tpx_zzz_yzzz_0, \
                                     tpz_yzz_yzzz_0, tpz_yzz_zzzz_0, tpz_zzz_xxxx_0, tpz_zzz_xxxy_0, tpz_zzz_xxxz_0, \
                                     tpz_zzz_xxyy_0, tpz_zzz_xxyz_0, tpz_zzz_xxzz_0, tpz_zzz_xyyy_0, tpz_zzz_xyyz_0, \
                                     tpz_zzz_xyzz_0, tpz_zzz_xzzz_0, tpz_zzz_yyyy_0, tpz_zzz_yyyz_0, tpz_zzz_yyzz_0, \
                                     tpz_zzz_yzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yyzz_yzzz_0[j] = pa_y[j] * tlx_yzz_yzzz_0[j] + 0.5 * fl1_fx * tlx_zz_yzzz_0[j] + 0.5 * fl1_fx * tlx_yzz_zzz_0[j] -
                                 0.5 * fl1_fx * tpz_yzz_yzzz_0[j] - fl1_fx * fl1_fgb * tdz_yzz_yzzz_0[j];

            tly_yyzz_yzzz_0[j] = pa_y[j] * tly_yzz_yzzz_0[j] + 0.5 * fl1_fx * tly_zz_yzzz_0[j] + 0.5 * fl1_fx * tly_yzz_zzz_0[j];

            tlz_yyzz_yzzz_0[j] = pa_y[j] * tlz_yzz_yzzz_0[j] + 0.5 * fl1_fx * tlz_zz_yzzz_0[j] + 0.5 * fl1_fx * tlz_yzz_zzz_0[j] +
                                 0.5 * fl1_fx * tpx_yzz_yzzz_0[j] + fl1_fx * fl1_fgb * tdx_yzz_yzzz_0[j];

            tlx_yyzz_zzzz_0[j] = pa_y[j] * tlx_yzz_zzzz_0[j] + 0.5 * fl1_fx * tlx_zz_zzzz_0[j] - 0.5 * fl1_fx * tpz_yzz_zzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_yzz_zzzz_0[j];

            tly_yyzz_zzzz_0[j] = pa_y[j] * tly_yzz_zzzz_0[j] + 0.5 * fl1_fx * tly_zz_zzzz_0[j];

            tlz_yyzz_zzzz_0[j] = pa_y[j] * tlz_yzz_zzzz_0[j] + 0.5 * fl1_fx * tlz_zz_zzzz_0[j] + 0.5 * fl1_fx * tpx_yzz_zzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_yzz_zzzz_0[j];

            tlx_yzzz_xxxx_0[j] = pa_y[j] * tlx_zzz_xxxx_0[j] - 0.5 * fl1_fx * tpz_zzz_xxxx_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxxx_0[j];

            tly_yzzz_xxxx_0[j] = pa_y[j] * tly_zzz_xxxx_0[j];

            tlz_yzzz_xxxx_0[j] = pa_y[j] * tlz_zzz_xxxx_0[j] + 0.5 * fl1_fx * tpx_zzz_xxxx_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxxx_0[j];

            tlx_yzzz_xxxy_0[j] = pa_y[j] * tlx_zzz_xxxy_0[j] + 0.5 * fl1_fx * tlx_zzz_xxx_0[j] - 0.5 * fl1_fx * tpz_zzz_xxxy_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_xxxy_0[j];

            tly_yzzz_xxxy_0[j] = pa_y[j] * tly_zzz_xxxy_0[j] + 0.5 * fl1_fx * tly_zzz_xxx_0[j];

            tlz_yzzz_xxxy_0[j] = pa_y[j] * tlz_zzz_xxxy_0[j] + 0.5 * fl1_fx * tlz_zzz_xxx_0[j] + 0.5 * fl1_fx * tpx_zzz_xxxy_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_xxxy_0[j];

            tlx_yzzz_xxxz_0[j] = pa_y[j] * tlx_zzz_xxxz_0[j] - 0.5 * fl1_fx * tpz_zzz_xxxz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxxz_0[j];

            tly_yzzz_xxxz_0[j] = pa_y[j] * tly_zzz_xxxz_0[j];

            tlz_yzzz_xxxz_0[j] = pa_y[j] * tlz_zzz_xxxz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxxz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxxz_0[j];

            tlx_yzzz_xxyy_0[j] =
                pa_y[j] * tlx_zzz_xxyy_0[j] + fl1_fx * tlx_zzz_xxy_0[j] - 0.5 * fl1_fx * tpz_zzz_xxyy_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxyy_0[j];

            tly_yzzz_xxyy_0[j] = pa_y[j] * tly_zzz_xxyy_0[j] + fl1_fx * tly_zzz_xxy_0[j];

            tlz_yzzz_xxyy_0[j] =
                pa_y[j] * tlz_zzz_xxyy_0[j] + fl1_fx * tlz_zzz_xxy_0[j] + 0.5 * fl1_fx * tpx_zzz_xxyy_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxyy_0[j];

            tlx_yzzz_xxyz_0[j] = pa_y[j] * tlx_zzz_xxyz_0[j] + 0.5 * fl1_fx * tlx_zzz_xxz_0[j] - 0.5 * fl1_fx * tpz_zzz_xxyz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_xxyz_0[j];

            tly_yzzz_xxyz_0[j] = pa_y[j] * tly_zzz_xxyz_0[j] + 0.5 * fl1_fx * tly_zzz_xxz_0[j];

            tlz_yzzz_xxyz_0[j] = pa_y[j] * tlz_zzz_xxyz_0[j] + 0.5 * fl1_fx * tlz_zzz_xxz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxyz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_xxyz_0[j];

            tlx_yzzz_xxzz_0[j] = pa_y[j] * tlx_zzz_xxzz_0[j] - 0.5 * fl1_fx * tpz_zzz_xxzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xxzz_0[j];

            tly_yzzz_xxzz_0[j] = pa_y[j] * tly_zzz_xxzz_0[j];

            tlz_yzzz_xxzz_0[j] = pa_y[j] * tlz_zzz_xxzz_0[j] + 0.5 * fl1_fx * tpx_zzz_xxzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xxzz_0[j];

            tlx_yzzz_xyyy_0[j] = pa_y[j] * tlx_zzz_xyyy_0[j] + 1.5 * fl1_fx * tlx_zzz_xyy_0[j] - 0.5 * fl1_fx * tpz_zzz_xyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_xyyy_0[j];

            tly_yzzz_xyyy_0[j] = pa_y[j] * tly_zzz_xyyy_0[j] + 1.5 * fl1_fx * tly_zzz_xyy_0[j];

            tlz_yzzz_xyyy_0[j] = pa_y[j] * tlz_zzz_xyyy_0[j] + 1.5 * fl1_fx * tlz_zzz_xyy_0[j] + 0.5 * fl1_fx * tpx_zzz_xyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_xyyy_0[j];

            tlx_yzzz_xyyz_0[j] =
                pa_y[j] * tlx_zzz_xyyz_0[j] + fl1_fx * tlx_zzz_xyz_0[j] - 0.5 * fl1_fx * tpz_zzz_xyyz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xyyz_0[j];

            tly_yzzz_xyyz_0[j] = pa_y[j] * tly_zzz_xyyz_0[j] + fl1_fx * tly_zzz_xyz_0[j];

            tlz_yzzz_xyyz_0[j] =
                pa_y[j] * tlz_zzz_xyyz_0[j] + fl1_fx * tlz_zzz_xyz_0[j] + 0.5 * fl1_fx * tpx_zzz_xyyz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xyyz_0[j];

            tlx_yzzz_xyzz_0[j] = pa_y[j] * tlx_zzz_xyzz_0[j] + 0.5 * fl1_fx * tlx_zzz_xzz_0[j] - 0.5 * fl1_fx * tpz_zzz_xyzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_xyzz_0[j];

            tly_yzzz_xyzz_0[j] = pa_y[j] * tly_zzz_xyzz_0[j] + 0.5 * fl1_fx * tly_zzz_xzz_0[j];

            tlz_yzzz_xyzz_0[j] = pa_y[j] * tlz_zzz_xyzz_0[j] + 0.5 * fl1_fx * tlz_zzz_xzz_0[j] + 0.5 * fl1_fx * tpx_zzz_xyzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_xyzz_0[j];

            tlx_yzzz_xzzz_0[j] = pa_y[j] * tlx_zzz_xzzz_0[j] - 0.5 * fl1_fx * tpz_zzz_xzzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_xzzz_0[j];

            tly_yzzz_xzzz_0[j] = pa_y[j] * tly_zzz_xzzz_0[j];

            tlz_yzzz_xzzz_0[j] = pa_y[j] * tlz_zzz_xzzz_0[j] + 0.5 * fl1_fx * tpx_zzz_xzzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_xzzz_0[j];

            tlx_yzzz_yyyy_0[j] = pa_y[j] * tlx_zzz_yyyy_0[j] + 2.0 * fl1_fx * tlx_zzz_yyy_0[j] - 0.5 * fl1_fx * tpz_zzz_yyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_yyyy_0[j];

            tly_yzzz_yyyy_0[j] = pa_y[j] * tly_zzz_yyyy_0[j] + 2.0 * fl1_fx * tly_zzz_yyy_0[j];

            tlz_yzzz_yyyy_0[j] = pa_y[j] * tlz_zzz_yyyy_0[j] + 2.0 * fl1_fx * tlz_zzz_yyy_0[j] + 0.5 * fl1_fx * tpx_zzz_yyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_yyyy_0[j];

            tlx_yzzz_yyyz_0[j] = pa_y[j] * tlx_zzz_yyyz_0[j] + 1.5 * fl1_fx * tlx_zzz_yyz_0[j] - 0.5 * fl1_fx * tpz_zzz_yyyz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_yyyz_0[j];

            tly_yzzz_yyyz_0[j] = pa_y[j] * tly_zzz_yyyz_0[j] + 1.5 * fl1_fx * tly_zzz_yyz_0[j];

            tlz_yzzz_yyyz_0[j] = pa_y[j] * tlz_zzz_yyyz_0[j] + 1.5 * fl1_fx * tlz_zzz_yyz_0[j] + 0.5 * fl1_fx * tpx_zzz_yyyz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_yyyz_0[j];

            tlx_yzzz_yyzz_0[j] =
                pa_y[j] * tlx_zzz_yyzz_0[j] + fl1_fx * tlx_zzz_yzz_0[j] - 0.5 * fl1_fx * tpz_zzz_yyzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_yyzz_0[j];

            tly_yzzz_yyzz_0[j] = pa_y[j] * tly_zzz_yyzz_0[j] + fl1_fx * tly_zzz_yzz_0[j];

            tlz_yzzz_yyzz_0[j] =
                pa_y[j] * tlz_zzz_yyzz_0[j] + fl1_fx * tlz_zzz_yzz_0[j] + 0.5 * fl1_fx * tpx_zzz_yyzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_yyzz_0[j];

            tlx_yzzz_yzzz_0[j] = pa_y[j] * tlx_zzz_yzzz_0[j] + 0.5 * fl1_fx * tlx_zzz_zzz_0[j] - 0.5 * fl1_fx * tpz_zzz_yzzz_0[j] -
                                 fl1_fx * fl1_fgb * tdz_zzz_yzzz_0[j];

            tly_yzzz_yzzz_0[j] = pa_y[j] * tly_zzz_yzzz_0[j] + 0.5 * fl1_fx * tly_zzz_zzz_0[j];

            tlz_yzzz_yzzz_0[j] = pa_y[j] * tlz_zzz_yzzz_0[j] + 0.5 * fl1_fx * tlz_zzz_zzz_0[j] + 0.5 * fl1_fx * tpx_zzz_yzzz_0[j] +
                                 fl1_fx * fl1_fgb * tdx_zzz_yzzz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGG_627_675(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_l_4_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tlx_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 135);

        auto tly_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tlz_zzz_xxxx_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 135);

        auto tlx_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 136);

        auto tly_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tlz_zzz_xxxy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 136);

        auto tlx_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 137);

        auto tly_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tlz_zzz_xxxz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 137);

        auto tlx_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 138);

        auto tly_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tlz_zzz_xxyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 138);

        auto tlx_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 139);

        auto tly_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tlz_zzz_xxyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 139);

        auto tlx_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 140);

        auto tly_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tlz_zzz_xxzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 140);

        auto tlx_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 141);

        auto tly_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tlz_zzz_xyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 141);

        auto tlx_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 142);

        auto tly_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tlz_zzz_xyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 142);

        auto tlx_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 143);

        auto tly_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tlz_zzz_xyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 143);

        auto tlx_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 144);

        auto tly_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tlz_zzz_xzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 144);

        auto tlx_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 145);

        auto tly_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tlz_zzz_yyyy_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 145);

        auto tlx_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 146);

        auto tly_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tlz_zzz_yyyz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 146);

        auto tlx_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 147);

        auto tly_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tlz_zzz_yyzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 147);

        auto tlx_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 148);

        auto tly_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tlz_zzz_yzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 148);

        auto tlx_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * idx + 149);

        auto tly_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tlz_zzz_zzzz_0 = primBuffer.data(pidx_l_3_4_m0 + 300 * bdim + 150 * idx + 149);

        auto tlx_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 75);

        auto tly_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 75);

        auto tlz_zz_xxxx_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 75);

        auto tlx_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 76);

        auto tly_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 76);

        auto tlz_zz_xxxy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 76);

        auto tlx_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 77);

        auto tly_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 77);

        auto tlz_zz_xxxz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 77);

        auto tlx_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 78);

        auto tly_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 78);

        auto tlz_zz_xxyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 78);

        auto tlx_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 79);

        auto tly_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 79);

        auto tlz_zz_xxyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 79);

        auto tlx_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 80);

        auto tly_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 80);

        auto tlz_zz_xxzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 80);

        auto tlx_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 81);

        auto tly_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 81);

        auto tlz_zz_xyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 81);

        auto tlx_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 82);

        auto tly_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 82);

        auto tlz_zz_xyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 82);

        auto tlx_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 83);

        auto tly_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 83);

        auto tlz_zz_xyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 83);

        auto tlx_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 84);

        auto tly_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 84);

        auto tlz_zz_xzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 84);

        auto tlx_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 85);

        auto tly_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 85);

        auto tlz_zz_yyyy_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 85);

        auto tlx_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 86);

        auto tly_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 86);

        auto tlz_zz_yyyz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 86);

        auto tlx_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 87);

        auto tly_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 87);

        auto tlz_zz_yyzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 87);

        auto tlx_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 88);

        auto tly_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 88);

        auto tlz_zz_yzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 88);

        auto tlx_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * idx + 89);

        auto tly_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 90 * bdim + 90 * idx + 89);

        auto tlz_zz_zzzz_0 = primBuffer.data(pidx_l_2_4_m0 + 180 * bdim + 90 * idx + 89);

        auto tlx_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 90);

        auto tly_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 90);

        auto tlz_zzz_xxx_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 90);

        auto tlx_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 91);

        auto tly_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 91);

        auto tlz_zzz_xxy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 91);

        auto tlx_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 92);

        auto tly_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 92);

        auto tlz_zzz_xxz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 92);

        auto tlx_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 93);

        auto tly_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 93);

        auto tlz_zzz_xyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 93);

        auto tlx_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 94);

        auto tly_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 94);

        auto tlz_zzz_xyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 94);

        auto tlx_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 95);

        auto tly_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 95);

        auto tlz_zzz_xzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 95);

        auto tlx_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 96);

        auto tly_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 96);

        auto tlz_zzz_yyy_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 96);

        auto tlx_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 97);

        auto tly_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 97);

        auto tlz_zzz_yyz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 97);

        auto tlx_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 98);

        auto tly_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 98);

        auto tlz_zzz_yzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 98);

        auto tlx_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * idx + 99);

        auto tly_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 100 * bdim + 100 * idx + 99);

        auto tlz_zzz_zzz_0 = primBuffer.data(pidx_l_3_3_m0 + 200 * bdim + 100 * idx + 99);

        auto tpx_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 135);

        auto tpy_zzz_xxxx_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tpx_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 136);

        auto tpy_zzz_xxxy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tpx_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 137);

        auto tpy_zzz_xxxz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tpx_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 138);

        auto tpy_zzz_xxyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tpx_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 139);

        auto tpy_zzz_xxyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tpx_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 140);

        auto tpy_zzz_xxzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tpx_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 141);

        auto tpy_zzz_xyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tpx_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 142);

        auto tpy_zzz_xyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tpx_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 143);

        auto tpy_zzz_xyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tpx_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 144);

        auto tpy_zzz_xzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tpx_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 145);

        auto tpy_zzz_yyyy_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tpx_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 146);

        auto tpy_zzz_yyyz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tpx_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 147);

        auto tpy_zzz_yyzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tpx_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 148);

        auto tpy_zzz_yzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tpx_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * idx + 149);

        auto tpy_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tpz_zzz_zzzz_0 = primBuffer.data(pidx_p_3_4_m0 + 300 * bdim + 150 * idx + 149);

        auto tdx_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 135);

        auto tdy_zzz_xxxx_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 135);

        auto tdx_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 136);

        auto tdy_zzz_xxxy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 136);

        auto tdx_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 137);

        auto tdy_zzz_xxxz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 137);

        auto tdx_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 138);

        auto tdy_zzz_xxyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 138);

        auto tdx_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 139);

        auto tdy_zzz_xxyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 139);

        auto tdx_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 140);

        auto tdy_zzz_xxzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 140);

        auto tdx_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 141);

        auto tdy_zzz_xyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 141);

        auto tdx_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 142);

        auto tdy_zzz_xyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 142);

        auto tdx_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 143);

        auto tdy_zzz_xyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 143);

        auto tdx_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 144);

        auto tdy_zzz_xzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 144);

        auto tdx_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 145);

        auto tdy_zzz_yyyy_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 145);

        auto tdx_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 146);

        auto tdy_zzz_yyyz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 146);

        auto tdx_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 147);

        auto tdy_zzz_yyzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 147);

        auto tdx_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 148);

        auto tdy_zzz_yzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 148);

        auto tdx_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * idx + 149);

        auto tdy_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 150 * bdim + 150 * idx + 149);

        auto tdz_zzz_zzzz_0 = primBuffer.data(pidx_d_3_4_m0 + 300 * bdim + 150 * idx + 149);

        // set up pointers to integrals

        auto tlx_yzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 209);

        auto tly_yzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 209);

        auto tlz_yzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 209);

        auto tlx_zzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 210);

        auto tly_zzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 210);

        auto tlz_zzzz_xxxx_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 210);

        auto tlx_zzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 211);

        auto tly_zzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 211);

        auto tlz_zzzz_xxxy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 211);

        auto tlx_zzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 212);

        auto tly_zzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 212);

        auto tlz_zzzz_xxxz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 212);

        auto tlx_zzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 213);

        auto tly_zzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 213);

        auto tlz_zzzz_xxyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 213);

        auto tlx_zzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 214);

        auto tly_zzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 214);

        auto tlz_zzzz_xxyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 214);

        auto tlx_zzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 215);

        auto tly_zzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 215);

        auto tlz_zzzz_xxzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 215);

        auto tlx_zzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 216);

        auto tly_zzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 216);

        auto tlz_zzzz_xyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 216);

        auto tlx_zzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 217);

        auto tly_zzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 217);

        auto tlz_zzzz_xyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 217);

        auto tlx_zzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 218);

        auto tly_zzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 218);

        auto tlz_zzzz_xyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 218);

        auto tlx_zzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 219);

        auto tly_zzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 219);

        auto tlz_zzzz_xzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 219);

        auto tlx_zzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 220);

        auto tly_zzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 220);

        auto tlz_zzzz_yyyy_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 220);

        auto tlx_zzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 221);

        auto tly_zzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 221);

        auto tlz_zzzz_yyyz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 221);

        auto tlx_zzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 222);

        auto tly_zzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 222);

        auto tlz_zzzz_yyzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 222);

        auto tlx_zzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 223);

        auto tly_zzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 223);

        auto tlz_zzzz_yzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 223);

        auto tlx_zzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * idx + 224);

        auto tly_zzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 225 * bdim + 225 * idx + 224);

        auto tlz_zzzz_zzzz_0 = primBuffer.data(pidx_l_4_4_m0 + 450 * bdim + 225 * idx + 224);

        // Batch of Integrals (627,675)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tdx_zzz_xxxx_0, tdx_zzz_xxxy_0, tdx_zzz_xxxz_0, \
                                     tdx_zzz_xxyy_0, tdx_zzz_xxyz_0, tdx_zzz_xxzz_0, tdx_zzz_xyyy_0, tdx_zzz_xyyz_0, \
                                     tdx_zzz_xyzz_0, tdx_zzz_xzzz_0, tdx_zzz_yyyy_0, tdx_zzz_yyyz_0, tdx_zzz_yyzz_0, \
                                     tdx_zzz_yzzz_0, tdx_zzz_zzzz_0, tdy_zzz_xxxx_0, tdy_zzz_xxxy_0, tdy_zzz_xxxz_0, \
                                     tdy_zzz_xxyy_0, tdy_zzz_xxyz_0, tdy_zzz_xxzz_0, tdy_zzz_xyyy_0, tdy_zzz_xyyz_0, \
                                     tdy_zzz_xyzz_0, tdy_zzz_xzzz_0, tdy_zzz_yyyy_0, tdy_zzz_yyyz_0, tdy_zzz_yyzz_0, \
                                     tdy_zzz_yzzz_0, tdy_zzz_zzzz_0, tdz_zzz_zzzz_0, tlx_yzzz_zzzz_0, tlx_zz_xxxx_0, \
                                     tlx_zz_xxxy_0, tlx_zz_xxxz_0, tlx_zz_xxyy_0, tlx_zz_xxyz_0, tlx_zz_xxzz_0, \
                                     tlx_zz_xyyy_0, tlx_zz_xyyz_0, tlx_zz_xyzz_0, tlx_zz_xzzz_0, tlx_zz_yyyy_0, \
                                     tlx_zz_yyyz_0, tlx_zz_yyzz_0, tlx_zz_yzzz_0, tlx_zz_zzzz_0, tlx_zzz_xxx_0, \
                                     tlx_zzz_xxxx_0, tlx_zzz_xxxy_0, tlx_zzz_xxxz_0, tlx_zzz_xxy_0, tlx_zzz_xxyy_0, \
                                     tlx_zzz_xxyz_0, tlx_zzz_xxz_0, tlx_zzz_xxzz_0, tlx_zzz_xyy_0, tlx_zzz_xyyy_0, \
                                     tlx_zzz_xyyz_0, tlx_zzz_xyz_0, tlx_zzz_xyzz_0, tlx_zzz_xzz_0, tlx_zzz_xzzz_0, \
                                     tlx_zzz_yyy_0, tlx_zzz_yyyy_0, tlx_zzz_yyyz_0, tlx_zzz_yyz_0, tlx_zzz_yyzz_0, \
                                     tlx_zzz_yzz_0, tlx_zzz_yzzz_0, tlx_zzz_zzz_0, tlx_zzz_zzzz_0, tlx_zzzz_xxxx_0, \
                                     tlx_zzzz_xxxy_0, tlx_zzzz_xxxz_0, tlx_zzzz_xxyy_0, tlx_zzzz_xxyz_0, tlx_zzzz_xxzz_0, \
                                     tlx_zzzz_xyyy_0, tlx_zzzz_xyyz_0, tlx_zzzz_xyzz_0, tlx_zzzz_xzzz_0, tlx_zzzz_yyyy_0, \
                                     tlx_zzzz_yyyz_0, tlx_zzzz_yyzz_0, tlx_zzzz_yzzz_0, tlx_zzzz_zzzz_0, tly_yzzz_zzzz_0, \
                                     tly_zz_xxxx_0, tly_zz_xxxy_0, tly_zz_xxxz_0, tly_zz_xxyy_0, tly_zz_xxyz_0, \
                                     tly_zz_xxzz_0, tly_zz_xyyy_0, tly_zz_xyyz_0, tly_zz_xyzz_0, tly_zz_xzzz_0, \
                                     tly_zz_yyyy_0, tly_zz_yyyz_0, tly_zz_yyzz_0, tly_zz_yzzz_0, tly_zz_zzzz_0, \
                                     tly_zzz_xxx_0, tly_zzz_xxxx_0, tly_zzz_xxxy_0, tly_zzz_xxxz_0, tly_zzz_xxy_0, \
                                     tly_zzz_xxyy_0, tly_zzz_xxyz_0, tly_zzz_xxz_0, tly_zzz_xxzz_0, tly_zzz_xyy_0, \
                                     tly_zzz_xyyy_0, tly_zzz_xyyz_0, tly_zzz_xyz_0, tly_zzz_xyzz_0, tly_zzz_xzz_0, \
                                     tly_zzz_xzzz_0, tly_zzz_yyy_0, tly_zzz_yyyy_0, tly_zzz_yyyz_0, tly_zzz_yyz_0, \
                                     tly_zzz_yyzz_0, tly_zzz_yzz_0, tly_zzz_yzzz_0, tly_zzz_zzz_0, tly_zzz_zzzz_0, \
                                     tly_zzzz_xxxx_0, tly_zzzz_xxxy_0, tly_zzzz_xxxz_0, tly_zzzz_xxyy_0, tly_zzzz_xxyz_0, \
                                     tly_zzzz_xxzz_0, tly_zzzz_xyyy_0, tly_zzzz_xyyz_0, tly_zzzz_xyzz_0, tly_zzzz_xzzz_0, \
                                     tly_zzzz_yyyy_0, tly_zzzz_yyyz_0, tly_zzzz_yyzz_0, tly_zzzz_yzzz_0, tly_zzzz_zzzz_0, \
                                     tlz_yzzz_zzzz_0, tlz_zz_xxxx_0, tlz_zz_xxxy_0, tlz_zz_xxxz_0, tlz_zz_xxyy_0, \
                                     tlz_zz_xxyz_0, tlz_zz_xxzz_0, tlz_zz_xyyy_0, tlz_zz_xyyz_0, tlz_zz_xyzz_0, \
                                     tlz_zz_xzzz_0, tlz_zz_yyyy_0, tlz_zz_yyyz_0, tlz_zz_yyzz_0, tlz_zz_yzzz_0, \
                                     tlz_zz_zzzz_0, tlz_zzz_xxx_0, tlz_zzz_xxxx_0, tlz_zzz_xxxy_0, tlz_zzz_xxxz_0, \
                                     tlz_zzz_xxy_0, tlz_zzz_xxyy_0, tlz_zzz_xxyz_0, tlz_zzz_xxz_0, tlz_zzz_xxzz_0, \
                                     tlz_zzz_xyy_0, tlz_zzz_xyyy_0, tlz_zzz_xyyz_0, tlz_zzz_xyz_0, tlz_zzz_xyzz_0, \
                                     tlz_zzz_xzz_0, tlz_zzz_xzzz_0, tlz_zzz_yyy_0, tlz_zzz_yyyy_0, tlz_zzz_yyyz_0, \
                                     tlz_zzz_yyz_0, tlz_zzz_yyzz_0, tlz_zzz_yzz_0, tlz_zzz_yzzz_0, tlz_zzz_zzz_0, \
                                     tlz_zzz_zzzz_0, tlz_zzzz_xxxx_0, tlz_zzzz_xxxy_0, tlz_zzzz_xxxz_0, tlz_zzzz_xxyy_0, \
                                     tlz_zzzz_xxyz_0, tlz_zzzz_xxzz_0, tlz_zzzz_xyyy_0, tlz_zzzz_xyyz_0, tlz_zzzz_xyzz_0, \
                                     tlz_zzzz_xzzz_0, tlz_zzzz_yyyy_0, tlz_zzzz_yyyz_0, tlz_zzzz_yyzz_0, tlz_zzzz_yzzz_0, \
                                     tlz_zzzz_zzzz_0, tpx_zzz_xxxx_0, tpx_zzz_xxxy_0, tpx_zzz_xxxz_0, tpx_zzz_xxyy_0, \
                                     tpx_zzz_xxyz_0, tpx_zzz_xxzz_0, tpx_zzz_xyyy_0, tpx_zzz_xyyz_0, tpx_zzz_xyzz_0, \
                                     tpx_zzz_xzzz_0, tpx_zzz_yyyy_0, tpx_zzz_yyyz_0, tpx_zzz_yyzz_0, tpx_zzz_yzzz_0, \
                                     tpx_zzz_zzzz_0, tpy_zzz_xxxx_0, tpy_zzz_xxxy_0, tpy_zzz_xxxz_0, tpy_zzz_xxyy_0, \
                                     tpy_zzz_xxyz_0, tpy_zzz_xxzz_0, tpy_zzz_xyyy_0, tpy_zzz_xyyz_0, tpy_zzz_xyzz_0, \
                                     tpy_zzz_xzzz_0, tpy_zzz_yyyy_0, tpy_zzz_yyyz_0, tpy_zzz_yyzz_0, tpy_zzz_yzzz_0, \
                                     tpy_zzz_zzzz_0, tpz_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_yzzz_zzzz_0[j] = pa_y[j] * tlx_zzz_zzzz_0[j] - 0.5 * fl1_fx * tpz_zzz_zzzz_0[j] - fl1_fx * fl1_fgb * tdz_zzz_zzzz_0[j];

            tly_yzzz_zzzz_0[j] = pa_y[j] * tly_zzz_zzzz_0[j];

            tlz_yzzz_zzzz_0[j] = pa_y[j] * tlz_zzz_zzzz_0[j] + 0.5 * fl1_fx * tpx_zzz_zzzz_0[j] + fl1_fx * fl1_fgb * tdx_zzz_zzzz_0[j];

            tlx_zzzz_xxxx_0[j] = pa_z[j] * tlx_zzz_xxxx_0[j] + 1.5 * fl1_fx * tlx_zz_xxxx_0[j] + 0.5 * fl1_fx * tpy_zzz_xxxx_0[j] +
                                 fl1_fx * fl1_fgb * tdy_zzz_xxxx_0[j];

            tly_zzzz_xxxx_0[j] = pa_z[j] * tly_zzz_xxxx_0[j] + 1.5 * fl1_fx * tly_zz_xxxx_0[j] - 0.5 * fl1_fx * tpx_zzz_xxxx_0[j] -
                                 fl1_fx * fl1_fgb * tdx_zzz_xxxx_0[j];

            tlz_zzzz_xxxx_0[j] = pa_z[j] * tlz_zzz_xxxx_0[j] + 1.5 * fl1_fx * tlz_zz_xxxx_0[j];

            tlx_zzzz_xxxy_0[j] = pa_z[j] * tlx_zzz_xxxy_0[j] + 1.5 * fl1_fx * tlx_zz_xxxy_0[j] + 0.5 * fl1_fx * tpy_zzz_xxxy_0[j] +
                                 fl1_fx * fl1_fgb * tdy_zzz_xxxy_0[j];

            tly_zzzz_xxxy_0[j] = pa_z[j] * tly_zzz_xxxy_0[j] + 1.5 * fl1_fx * tly_zz_xxxy_0[j] - 0.5 * fl1_fx * tpx_zzz_xxxy_0[j] -
                                 fl1_fx * fl1_fgb * tdx_zzz_xxxy_0[j];

            tlz_zzzz_xxxy_0[j] = pa_z[j] * tlz_zzz_xxxy_0[j] + 1.5 * fl1_fx * tlz_zz_xxxy_0[j];

            tlx_zzzz_xxxz_0[j] = pa_z[j] * tlx_zzz_xxxz_0[j] + 1.5 * fl1_fx * tlx_zz_xxxz_0[j] + 0.5 * fl1_fx * tlx_zzz_xxx_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_xxxz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xxxz_0[j];

            tly_zzzz_xxxz_0[j] = pa_z[j] * tly_zzz_xxxz_0[j] + 1.5 * fl1_fx * tly_zz_xxxz_0[j] + 0.5 * fl1_fx * tly_zzz_xxx_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_xxxz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xxxz_0[j];

            tlz_zzzz_xxxz_0[j] = pa_z[j] * tlz_zzz_xxxz_0[j] + 1.5 * fl1_fx * tlz_zz_xxxz_0[j] + 0.5 * fl1_fx * tlz_zzz_xxx_0[j];

            tlx_zzzz_xxyy_0[j] = pa_z[j] * tlx_zzz_xxyy_0[j] + 1.5 * fl1_fx * tlx_zz_xxyy_0[j] + 0.5 * fl1_fx * tpy_zzz_xxyy_0[j] +
                                 fl1_fx * fl1_fgb * tdy_zzz_xxyy_0[j];

            tly_zzzz_xxyy_0[j] = pa_z[j] * tly_zzz_xxyy_0[j] + 1.5 * fl1_fx * tly_zz_xxyy_0[j] - 0.5 * fl1_fx * tpx_zzz_xxyy_0[j] -
                                 fl1_fx * fl1_fgb * tdx_zzz_xxyy_0[j];

            tlz_zzzz_xxyy_0[j] = pa_z[j] * tlz_zzz_xxyy_0[j] + 1.5 * fl1_fx * tlz_zz_xxyy_0[j];

            tlx_zzzz_xxyz_0[j] = pa_z[j] * tlx_zzz_xxyz_0[j] + 1.5 * fl1_fx * tlx_zz_xxyz_0[j] + 0.5 * fl1_fx * tlx_zzz_xxy_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_xxyz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xxyz_0[j];

            tly_zzzz_xxyz_0[j] = pa_z[j] * tly_zzz_xxyz_0[j] + 1.5 * fl1_fx * tly_zz_xxyz_0[j] + 0.5 * fl1_fx * tly_zzz_xxy_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_xxyz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xxyz_0[j];

            tlz_zzzz_xxyz_0[j] = pa_z[j] * tlz_zzz_xxyz_0[j] + 1.5 * fl1_fx * tlz_zz_xxyz_0[j] + 0.5 * fl1_fx * tlz_zzz_xxy_0[j];

            tlx_zzzz_xxzz_0[j] = pa_z[j] * tlx_zzz_xxzz_0[j] + 1.5 * fl1_fx * tlx_zz_xxzz_0[j] + fl1_fx * tlx_zzz_xxz_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_xxzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xxzz_0[j];

            tly_zzzz_xxzz_0[j] = pa_z[j] * tly_zzz_xxzz_0[j] + 1.5 * fl1_fx * tly_zz_xxzz_0[j] + fl1_fx * tly_zzz_xxz_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_xxzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xxzz_0[j];

            tlz_zzzz_xxzz_0[j] = pa_z[j] * tlz_zzz_xxzz_0[j] + 1.5 * fl1_fx * tlz_zz_xxzz_0[j] + fl1_fx * tlz_zzz_xxz_0[j];

            tlx_zzzz_xyyy_0[j] = pa_z[j] * tlx_zzz_xyyy_0[j] + 1.5 * fl1_fx * tlx_zz_xyyy_0[j] + 0.5 * fl1_fx * tpy_zzz_xyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdy_zzz_xyyy_0[j];

            tly_zzzz_xyyy_0[j] = pa_z[j] * tly_zzz_xyyy_0[j] + 1.5 * fl1_fx * tly_zz_xyyy_0[j] - 0.5 * fl1_fx * tpx_zzz_xyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdx_zzz_xyyy_0[j];

            tlz_zzzz_xyyy_0[j] = pa_z[j] * tlz_zzz_xyyy_0[j] + 1.5 * fl1_fx * tlz_zz_xyyy_0[j];

            tlx_zzzz_xyyz_0[j] = pa_z[j] * tlx_zzz_xyyz_0[j] + 1.5 * fl1_fx * tlx_zz_xyyz_0[j] + 0.5 * fl1_fx * tlx_zzz_xyy_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_xyyz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xyyz_0[j];

            tly_zzzz_xyyz_0[j] = pa_z[j] * tly_zzz_xyyz_0[j] + 1.5 * fl1_fx * tly_zz_xyyz_0[j] + 0.5 * fl1_fx * tly_zzz_xyy_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_xyyz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xyyz_0[j];

            tlz_zzzz_xyyz_0[j] = pa_z[j] * tlz_zzz_xyyz_0[j] + 1.5 * fl1_fx * tlz_zz_xyyz_0[j] + 0.5 * fl1_fx * tlz_zzz_xyy_0[j];

            tlx_zzzz_xyzz_0[j] = pa_z[j] * tlx_zzz_xyzz_0[j] + 1.5 * fl1_fx * tlx_zz_xyzz_0[j] + fl1_fx * tlx_zzz_xyz_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_xyzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xyzz_0[j];

            tly_zzzz_xyzz_0[j] = pa_z[j] * tly_zzz_xyzz_0[j] + 1.5 * fl1_fx * tly_zz_xyzz_0[j] + fl1_fx * tly_zzz_xyz_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_xyzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xyzz_0[j];

            tlz_zzzz_xyzz_0[j] = pa_z[j] * tlz_zzz_xyzz_0[j] + 1.5 * fl1_fx * tlz_zz_xyzz_0[j] + fl1_fx * tlz_zzz_xyz_0[j];

            tlx_zzzz_xzzz_0[j] = pa_z[j] * tlx_zzz_xzzz_0[j] + 1.5 * fl1_fx * tlx_zz_xzzz_0[j] + 1.5 * fl1_fx * tlx_zzz_xzz_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_xzzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_xzzz_0[j];

            tly_zzzz_xzzz_0[j] = pa_z[j] * tly_zzz_xzzz_0[j] + 1.5 * fl1_fx * tly_zz_xzzz_0[j] + 1.5 * fl1_fx * tly_zzz_xzz_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_xzzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_xzzz_0[j];

            tlz_zzzz_xzzz_0[j] = pa_z[j] * tlz_zzz_xzzz_0[j] + 1.5 * fl1_fx * tlz_zz_xzzz_0[j] + 1.5 * fl1_fx * tlz_zzz_xzz_0[j];

            tlx_zzzz_yyyy_0[j] = pa_z[j] * tlx_zzz_yyyy_0[j] + 1.5 * fl1_fx * tlx_zz_yyyy_0[j] + 0.5 * fl1_fx * tpy_zzz_yyyy_0[j] +
                                 fl1_fx * fl1_fgb * tdy_zzz_yyyy_0[j];

            tly_zzzz_yyyy_0[j] = pa_z[j] * tly_zzz_yyyy_0[j] + 1.5 * fl1_fx * tly_zz_yyyy_0[j] - 0.5 * fl1_fx * tpx_zzz_yyyy_0[j] -
                                 fl1_fx * fl1_fgb * tdx_zzz_yyyy_0[j];

            tlz_zzzz_yyyy_0[j] = pa_z[j] * tlz_zzz_yyyy_0[j] + 1.5 * fl1_fx * tlz_zz_yyyy_0[j];

            tlx_zzzz_yyyz_0[j] = pa_z[j] * tlx_zzz_yyyz_0[j] + 1.5 * fl1_fx * tlx_zz_yyyz_0[j] + 0.5 * fl1_fx * tlx_zzz_yyy_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_yyyz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yyyz_0[j];

            tly_zzzz_yyyz_0[j] = pa_z[j] * tly_zzz_yyyz_0[j] + 1.5 * fl1_fx * tly_zz_yyyz_0[j] + 0.5 * fl1_fx * tly_zzz_yyy_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_yyyz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yyyz_0[j];

            tlz_zzzz_yyyz_0[j] = pa_z[j] * tlz_zzz_yyyz_0[j] + 1.5 * fl1_fx * tlz_zz_yyyz_0[j] + 0.5 * fl1_fx * tlz_zzz_yyy_0[j];

            tlx_zzzz_yyzz_0[j] = pa_z[j] * tlx_zzz_yyzz_0[j] + 1.5 * fl1_fx * tlx_zz_yyzz_0[j] + fl1_fx * tlx_zzz_yyz_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_yyzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yyzz_0[j];

            tly_zzzz_yyzz_0[j] = pa_z[j] * tly_zzz_yyzz_0[j] + 1.5 * fl1_fx * tly_zz_yyzz_0[j] + fl1_fx * tly_zzz_yyz_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_yyzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yyzz_0[j];

            tlz_zzzz_yyzz_0[j] = pa_z[j] * tlz_zzz_yyzz_0[j] + 1.5 * fl1_fx * tlz_zz_yyzz_0[j] + fl1_fx * tlz_zzz_yyz_0[j];

            tlx_zzzz_yzzz_0[j] = pa_z[j] * tlx_zzz_yzzz_0[j] + 1.5 * fl1_fx * tlx_zz_yzzz_0[j] + 1.5 * fl1_fx * tlx_zzz_yzz_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_yzzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_yzzz_0[j];

            tly_zzzz_yzzz_0[j] = pa_z[j] * tly_zzz_yzzz_0[j] + 1.5 * fl1_fx * tly_zz_yzzz_0[j] + 1.5 * fl1_fx * tly_zzz_yzz_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_yzzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_yzzz_0[j];

            tlz_zzzz_yzzz_0[j] = pa_z[j] * tlz_zzz_yzzz_0[j] + 1.5 * fl1_fx * tlz_zz_yzzz_0[j] + 1.5 * fl1_fx * tlz_zzz_yzz_0[j];

            tlx_zzzz_zzzz_0[j] = pa_z[j] * tlx_zzz_zzzz_0[j] + 1.5 * fl1_fx * tlx_zz_zzzz_0[j] + 2.0 * fl1_fx * tlx_zzz_zzz_0[j] +
                                 0.5 * fl1_fx * tpy_zzz_zzzz_0[j] + fl1_fx * fl1_fgb * tdy_zzz_zzzz_0[j];

            tly_zzzz_zzzz_0[j] = pa_z[j] * tly_zzz_zzzz_0[j] + 1.5 * fl1_fx * tly_zz_zzzz_0[j] + 2.0 * fl1_fx * tly_zzz_zzz_0[j] -
                                 0.5 * fl1_fx * tpx_zzz_zzzz_0[j] - fl1_fx * fl1_fgb * tdx_zzz_zzzz_0[j];

            tlz_zzzz_zzzz_0[j] = pa_z[j] * tlz_zzz_zzzz_0[j] + 1.5 * fl1_fx * tlz_zz_zzzz_0[j] + 2.0 * fl1_fx * tlz_zzz_zzz_0[j];
        }

        idx++;
    }
}

}  // namespace amomrecfunc
