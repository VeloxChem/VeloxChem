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

#include "OverlapRecFuncForFF.hpp"

namespace ovlrecfunc {  // ovlrecfunc namespace

void
compOverlapForFF(CMemBlock2D<double>&       primBuffer,
                 const CRecursionMap&       recursionMap,
                 const CMemBlock2D<double>& osFactors,
                 const int32_t              nOSFactors,
                 const CMemBlock2D<double>& paDistances,
                 const CGtoBlock&           braGtoBlock,
                 const CGtoBlock&           ketGtoBlock,
                 const int32_t              iContrGto)
{
    ovlrecfunc::compOverlapForFF_0_50(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFF_50_100(primBuffer, recursionMap, osFactors, nOSFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compOverlapForFF_0_50(CMemBlock2D<double>&       primBuffer,
                      const CRecursionMap&       recursionMap,
                      const CMemBlock2D<double>& osFactors,
                      const int32_t              nOSFactors,
                      const CMemBlock2D<double>& paDistances,
                      const CGtoBlock&           braGtoBlock,
                      const CGtoBlock&           ketGtoBlock,
                      const int32_t              iContrGto)
{
    // Batch of Integrals (0,50)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

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

        auto ts_x_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx);

        auto ts_x_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 1);

        auto ts_x_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 2);

        auto ts_x_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 3);

        auto ts_x_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 4);

        auto ts_x_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 5);

        auto ts_x_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 6);

        auto ts_x_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 7);

        auto ts_x_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 8);

        auto ts_x_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 9);

        auto ts_y_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 10);

        auto ts_y_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 11);

        auto ts_y_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 12);

        auto ts_y_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 13);

        auto ts_y_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 14);

        auto ts_y_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 15);

        auto ts_y_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 16);

        auto ts_y_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 17);

        auto ts_y_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 18);

        auto ts_y_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 19);

        auto ts_z_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 20);

        auto ts_z_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 21);

        auto ts_z_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 22);

        auto ts_z_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 23);

        auto ts_z_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 24);

        auto ts_z_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 25);

        auto ts_z_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 26);

        auto ts_z_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 27);

        auto ts_z_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 28);

        auto ts_z_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 29);

        auto ts_xx_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx);

        auto ts_xx_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 1);

        auto ts_xx_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 2);

        auto ts_xx_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 3);

        auto ts_xx_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 4);

        auto ts_xx_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 5);

        auto ts_xy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 6);

        auto ts_xy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 7);

        auto ts_xy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 8);

        auto ts_xy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 9);

        auto ts_xy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 10);

        auto ts_xy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 11);

        auto ts_xz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 12);

        auto ts_xz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 13);

        auto ts_xz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 14);

        auto ts_xz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 15);

        auto ts_xz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 16);

        auto ts_xz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 17);

        auto ts_yy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 18);

        auto ts_yy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 19);

        auto ts_yy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 20);

        auto ts_yy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 21);

        auto ts_yy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 22);

        auto ts_yy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 23);

        auto ts_yz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 24);

        auto ts_yz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 25);

        auto ts_yz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 26);

        auto ts_yz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 27);

        auto ts_yz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 28);

        auto ts_yz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 29);

        // set up pointers to integrals

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

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fx, pa_x, ts_x_xxx_0, ts_x_xxy_0, ts_x_xxz_0, ts_x_xyy_0, ts_x_xyz_0, \
                                     ts_x_xzz_0, ts_x_yyy_0, ts_x_yyz_0, ts_x_yzz_0, ts_x_zzz_0, ts_xx_xx_0, \
                                     ts_xx_xxx_0, ts_xx_xxy_0, ts_xx_xxz_0, ts_xx_xy_0, ts_xx_xyy_0, ts_xx_xyz_0, \
                                     ts_xx_xz_0, ts_xx_xzz_0, ts_xx_yy_0, ts_xx_yyy_0, ts_xx_yyz_0, ts_xx_yz_0, \
                                     ts_xx_yzz_0, ts_xx_zz_0, ts_xx_zzz_0, ts_xxx_xxx_0, ts_xxx_xxy_0, ts_xxx_xxz_0, \
                                     ts_xxx_xyy_0, ts_xxx_xyz_0, ts_xxx_xzz_0, ts_xxx_yyy_0, ts_xxx_yyz_0, ts_xxx_yzz_0, \
                                     ts_xxx_zzz_0, ts_xxy_xxx_0, ts_xxy_xxy_0, ts_xxy_xxz_0, ts_xxy_xyy_0, ts_xxy_xyz_0, \
                                     ts_xxy_xzz_0, ts_xxy_yyy_0, ts_xxy_yyz_0, ts_xxy_yzz_0, ts_xxy_zzz_0, ts_xxz_xxx_0, \
                                     ts_xxz_xxy_0, ts_xxz_xxz_0, ts_xxz_xyy_0, ts_xxz_xyz_0, ts_xxz_xzz_0, ts_xxz_yyy_0, \
                                     ts_xxz_yyz_0, ts_xxz_yzz_0, ts_xxz_zzz_0, ts_xy_xx_0, ts_xy_xxx_0, ts_xy_xxy_0, \
                                     ts_xy_xxz_0, ts_xy_xy_0, ts_xy_xyy_0, ts_xy_xyz_0, ts_xy_xz_0, ts_xy_xzz_0, \
                                     ts_xy_yy_0, ts_xy_yyy_0, ts_xy_yyz_0, ts_xy_yz_0, ts_xy_yzz_0, ts_xy_zz_0, \
                                     ts_xy_zzz_0, ts_xyy_xxx_0, ts_xyy_xxy_0, ts_xyy_xxz_0, ts_xyy_xyy_0, ts_xyy_xyz_0, \
                                     ts_xyy_xzz_0, ts_xyy_yyy_0, ts_xyy_yyz_0, ts_xyy_yzz_0, ts_xyy_zzz_0, ts_xyz_xxx_0, \
                                     ts_xyz_xxy_0, ts_xyz_xxz_0, ts_xyz_xyy_0, ts_xyz_xyz_0, ts_xyz_xzz_0, ts_xyz_yyy_0, \
                                     ts_xyz_yyz_0, ts_xyz_yzz_0, ts_xyz_zzz_0, ts_xz_xx_0, ts_xz_xxx_0, ts_xz_xxy_0, \
                                     ts_xz_xxz_0, ts_xz_xy_0, ts_xz_xyy_0, ts_xz_xyz_0, ts_xz_xz_0, ts_xz_xzz_0, \
                                     ts_xz_yy_0, ts_xz_yyy_0, ts_xz_yyz_0, ts_xz_yz_0, ts_xz_yzz_0, ts_xz_zz_0, \
                                     ts_xz_zzz_0, ts_y_xxx_0, ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0, ts_y_xzz_0, \
                                     ts_y_yyy_0, ts_y_yyz_0, ts_y_yzz_0, ts_y_zzz_0, ts_yy_xx_0, ts_yy_xxx_0, \
                                     ts_yy_xxy_0, ts_yy_xxz_0, ts_yy_xy_0, ts_yy_xyy_0, ts_yy_xyz_0, ts_yy_xz_0, \
                                     ts_yy_xzz_0, ts_yy_yy_0, ts_yy_yyy_0, ts_yy_yyz_0, ts_yy_yz_0, ts_yy_yzz_0, \
                                     ts_yy_zz_0, ts_yy_zzz_0, ts_yz_xx_0, ts_yz_xxx_0, ts_yz_xxy_0, ts_yz_xxz_0, \
                                     ts_yz_xy_0, ts_yz_xyy_0, ts_yz_xyz_0, ts_yz_xz_0, ts_yz_xzz_0, ts_yz_yy_0, \
                                     ts_yz_yyy_0, ts_yz_yyz_0, ts_yz_yz_0, ts_yz_yzz_0, ts_yz_zz_0, ts_yz_zzz_0, \
                                     ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, \
                                     ts_z_yyz_0, ts_z_yzz_0, ts_z_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_xxx_xxx_0[j] = pa_x[j] * ts_xx_xxx_0[j] + fl1_fx * ts_x_xxx_0[j] + 1.5 * fl1_fx * ts_xx_xx_0[j];

            ts_xxx_xxy_0[j] = pa_x[j] * ts_xx_xxy_0[j] + fl1_fx * ts_x_xxy_0[j] + fl1_fx * ts_xx_xy_0[j];

            ts_xxx_xxz_0[j] = pa_x[j] * ts_xx_xxz_0[j] + fl1_fx * ts_x_xxz_0[j] + fl1_fx * ts_xx_xz_0[j];

            ts_xxx_xyy_0[j] = pa_x[j] * ts_xx_xyy_0[j] + fl1_fx * ts_x_xyy_0[j] + 0.5 * fl1_fx * ts_xx_yy_0[j];

            ts_xxx_xyz_0[j] = pa_x[j] * ts_xx_xyz_0[j] + fl1_fx * ts_x_xyz_0[j] + 0.5 * fl1_fx * ts_xx_yz_0[j];

            ts_xxx_xzz_0[j] = pa_x[j] * ts_xx_xzz_0[j] + fl1_fx * ts_x_xzz_0[j] + 0.5 * fl1_fx * ts_xx_zz_0[j];

            ts_xxx_yyy_0[j] = pa_x[j] * ts_xx_yyy_0[j] + fl1_fx * ts_x_yyy_0[j];

            ts_xxx_yyz_0[j] = pa_x[j] * ts_xx_yyz_0[j] + fl1_fx * ts_x_yyz_0[j];

            ts_xxx_yzz_0[j] = pa_x[j] * ts_xx_yzz_0[j] + fl1_fx * ts_x_yzz_0[j];

            ts_xxx_zzz_0[j] = pa_x[j] * ts_xx_zzz_0[j] + fl1_fx * ts_x_zzz_0[j];

            ts_xxy_xxx_0[j] = pa_x[j] * ts_xy_xxx_0[j] + 0.5 * fl1_fx * ts_y_xxx_0[j] + 1.5 * fl1_fx * ts_xy_xx_0[j];

            ts_xxy_xxy_0[j] = pa_x[j] * ts_xy_xxy_0[j] + 0.5 * fl1_fx * ts_y_xxy_0[j] + fl1_fx * ts_xy_xy_0[j];

            ts_xxy_xxz_0[j] = pa_x[j] * ts_xy_xxz_0[j] + 0.5 * fl1_fx * ts_y_xxz_0[j] + fl1_fx * ts_xy_xz_0[j];

            ts_xxy_xyy_0[j] = pa_x[j] * ts_xy_xyy_0[j] + 0.5 * fl1_fx * ts_y_xyy_0[j] + 0.5 * fl1_fx * ts_xy_yy_0[j];

            ts_xxy_xyz_0[j] = pa_x[j] * ts_xy_xyz_0[j] + 0.5 * fl1_fx * ts_y_xyz_0[j] + 0.5 * fl1_fx * ts_xy_yz_0[j];

            ts_xxy_xzz_0[j] = pa_x[j] * ts_xy_xzz_0[j] + 0.5 * fl1_fx * ts_y_xzz_0[j] + 0.5 * fl1_fx * ts_xy_zz_0[j];

            ts_xxy_yyy_0[j] = pa_x[j] * ts_xy_yyy_0[j] + 0.5 * fl1_fx * ts_y_yyy_0[j];

            ts_xxy_yyz_0[j] = pa_x[j] * ts_xy_yyz_0[j] + 0.5 * fl1_fx * ts_y_yyz_0[j];

            ts_xxy_yzz_0[j] = pa_x[j] * ts_xy_yzz_0[j] + 0.5 * fl1_fx * ts_y_yzz_0[j];

            ts_xxy_zzz_0[j] = pa_x[j] * ts_xy_zzz_0[j] + 0.5 * fl1_fx * ts_y_zzz_0[j];

            ts_xxz_xxx_0[j] = pa_x[j] * ts_xz_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxx_0[j] + 1.5 * fl1_fx * ts_xz_xx_0[j];

            ts_xxz_xxy_0[j] = pa_x[j] * ts_xz_xxy_0[j] + 0.5 * fl1_fx * ts_z_xxy_0[j] + fl1_fx * ts_xz_xy_0[j];

            ts_xxz_xxz_0[j] = pa_x[j] * ts_xz_xxz_0[j] + 0.5 * fl1_fx * ts_z_xxz_0[j] + fl1_fx * ts_xz_xz_0[j];

            ts_xxz_xyy_0[j] = pa_x[j] * ts_xz_xyy_0[j] + 0.5 * fl1_fx * ts_z_xyy_0[j] + 0.5 * fl1_fx * ts_xz_yy_0[j];

            ts_xxz_xyz_0[j] = pa_x[j] * ts_xz_xyz_0[j] + 0.5 * fl1_fx * ts_z_xyz_0[j] + 0.5 * fl1_fx * ts_xz_yz_0[j];

            ts_xxz_xzz_0[j] = pa_x[j] * ts_xz_xzz_0[j] + 0.5 * fl1_fx * ts_z_xzz_0[j] + 0.5 * fl1_fx * ts_xz_zz_0[j];

            ts_xxz_yyy_0[j] = pa_x[j] * ts_xz_yyy_0[j] + 0.5 * fl1_fx * ts_z_yyy_0[j];

            ts_xxz_yyz_0[j] = pa_x[j] * ts_xz_yyz_0[j] + 0.5 * fl1_fx * ts_z_yyz_0[j];

            ts_xxz_yzz_0[j] = pa_x[j] * ts_xz_yzz_0[j] + 0.5 * fl1_fx * ts_z_yzz_0[j];

            ts_xxz_zzz_0[j] = pa_x[j] * ts_xz_zzz_0[j] + 0.5 * fl1_fx * ts_z_zzz_0[j];

            ts_xyy_xxx_0[j] = pa_x[j] * ts_yy_xxx_0[j] + 1.5 * fl1_fx * ts_yy_xx_0[j];

            ts_xyy_xxy_0[j] = pa_x[j] * ts_yy_xxy_0[j] + fl1_fx * ts_yy_xy_0[j];

            ts_xyy_xxz_0[j] = pa_x[j] * ts_yy_xxz_0[j] + fl1_fx * ts_yy_xz_0[j];

            ts_xyy_xyy_0[j] = pa_x[j] * ts_yy_xyy_0[j] + 0.5 * fl1_fx * ts_yy_yy_0[j];

            ts_xyy_xyz_0[j] = pa_x[j] * ts_yy_xyz_0[j] + 0.5 * fl1_fx * ts_yy_yz_0[j];

            ts_xyy_xzz_0[j] = pa_x[j] * ts_yy_xzz_0[j] + 0.5 * fl1_fx * ts_yy_zz_0[j];

            ts_xyy_yyy_0[j] = pa_x[j] * ts_yy_yyy_0[j];

            ts_xyy_yyz_0[j] = pa_x[j] * ts_yy_yyz_0[j];

            ts_xyy_yzz_0[j] = pa_x[j] * ts_yy_yzz_0[j];

            ts_xyy_zzz_0[j] = pa_x[j] * ts_yy_zzz_0[j];

            ts_xyz_xxx_0[j] = pa_x[j] * ts_yz_xxx_0[j] + 1.5 * fl1_fx * ts_yz_xx_0[j];

            ts_xyz_xxy_0[j] = pa_x[j] * ts_yz_xxy_0[j] + fl1_fx * ts_yz_xy_0[j];

            ts_xyz_xxz_0[j] = pa_x[j] * ts_yz_xxz_0[j] + fl1_fx * ts_yz_xz_0[j];

            ts_xyz_xyy_0[j] = pa_x[j] * ts_yz_xyy_0[j] + 0.5 * fl1_fx * ts_yz_yy_0[j];

            ts_xyz_xyz_0[j] = pa_x[j] * ts_yz_xyz_0[j] + 0.5 * fl1_fx * ts_yz_yz_0[j];

            ts_xyz_xzz_0[j] = pa_x[j] * ts_yz_xzz_0[j] + 0.5 * fl1_fx * ts_yz_zz_0[j];

            ts_xyz_yyy_0[j] = pa_x[j] * ts_yz_yyy_0[j];

            ts_xyz_yyz_0[j] = pa_x[j] * ts_yz_yyz_0[j];

            ts_xyz_yzz_0[j] = pa_x[j] * ts_yz_yzz_0[j];

            ts_xyz_zzz_0[j] = pa_x[j] * ts_yz_zzz_0[j];
        }

        idx++;
    }
}

void
compOverlapForFF_50_100(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const int32_t              nOSFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    // Batch of Integrals (50,100)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_s_3_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_s_3_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(nOSFactors * idx);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto ts_yy_xxx_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 30);

        auto ts_yy_xxy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 31);

        auto ts_yy_xxz_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 32);

        auto ts_yy_xyy_0 = primBuffer.data(pidx_s_2_3_m0 + 60 * idx + 33);

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

        auto ts_y_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 10);

        auto ts_y_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 11);

        auto ts_y_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 12);

        auto ts_y_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 13);

        auto ts_y_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 14);

        auto ts_y_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 15);

        auto ts_y_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 16);

        auto ts_y_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 17);

        auto ts_y_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 18);

        auto ts_y_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 19);

        auto ts_z_xxx_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 20);

        auto ts_z_xxy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 21);

        auto ts_z_xxz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 22);

        auto ts_z_xyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 23);

        auto ts_z_xyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 24);

        auto ts_z_xzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 25);

        auto ts_z_yyy_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 26);

        auto ts_z_yyz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 27);

        auto ts_z_yzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 28);

        auto ts_z_zzz_0 = primBuffer.data(pidx_s_1_3_m0 + 30 * idx + 29);

        auto ts_yy_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 18);

        auto ts_yy_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 19);

        auto ts_yy_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 20);

        auto ts_yy_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 21);

        auto ts_yy_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 22);

        auto ts_yy_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 23);

        auto ts_yz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 24);

        auto ts_yz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 25);

        auto ts_yz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 26);

        auto ts_yz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 27);

        auto ts_yz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 28);

        auto ts_yz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 29);

        auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30);

        auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31);

        auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32);

        auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33);

        auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34);

        auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35);

        // set up pointers to integrals

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

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, ts_xzz_xxx_0, ts_xzz_xxy_0, ts_xzz_xxz_0, \
                                     ts_xzz_xyy_0, ts_xzz_xyz_0, ts_xzz_xzz_0, ts_xzz_yyy_0, ts_xzz_yyz_0, ts_xzz_yzz_0, \
                                     ts_xzz_zzz_0, ts_y_xxx_0, ts_y_xxy_0, ts_y_xxz_0, ts_y_xyy_0, ts_y_xyz_0, ts_y_xzz_0, \
                                     ts_y_yyy_0, ts_y_yyz_0, ts_y_yzz_0, ts_y_zzz_0, ts_yy_xx_0, ts_yy_xxx_0, \
                                     ts_yy_xxy_0, ts_yy_xxz_0, ts_yy_xy_0, ts_yy_xyy_0, ts_yy_xyz_0, ts_yy_xz_0, \
                                     ts_yy_xzz_0, ts_yy_yy_0, ts_yy_yyy_0, ts_yy_yyz_0, ts_yy_yz_0, ts_yy_yzz_0, \
                                     ts_yy_zz_0, ts_yy_zzz_0, ts_yyy_xxx_0, ts_yyy_xxy_0, ts_yyy_xxz_0, ts_yyy_xyy_0, \
                                     ts_yyy_xyz_0, ts_yyy_xzz_0, ts_yyy_yyy_0, ts_yyy_yyz_0, ts_yyy_yzz_0, ts_yyy_zzz_0, \
                                     ts_yyz_xxx_0, ts_yyz_xxy_0, ts_yyz_xxz_0, ts_yyz_xyy_0, ts_yyz_xyz_0, ts_yyz_xzz_0, \
                                     ts_yyz_yyy_0, ts_yyz_yyz_0, ts_yyz_yzz_0, ts_yyz_zzz_0, ts_yz_xx_0, ts_yz_xxx_0, \
                                     ts_yz_xxy_0, ts_yz_xxz_0, ts_yz_xy_0, ts_yz_xyy_0, ts_yz_xyz_0, ts_yz_xz_0, \
                                     ts_yz_xzz_0, ts_yz_yy_0, ts_yz_yyy_0, ts_yz_yyz_0, ts_yz_yz_0, ts_yz_yzz_0, \
                                     ts_yz_zz_0, ts_yz_zzz_0, ts_yzz_xxx_0, ts_yzz_xxy_0, ts_yzz_xxz_0, ts_yzz_xyy_0, \
                                     ts_yzz_xyz_0, ts_yzz_xzz_0, ts_yzz_yyy_0, ts_yzz_yyz_0, ts_yzz_yzz_0, ts_yzz_zzz_0, \
                                     ts_z_xxx_0, ts_z_xxy_0, ts_z_xxz_0, ts_z_xyy_0, ts_z_xyz_0, ts_z_xzz_0, ts_z_yyy_0, \
                                     ts_z_yyz_0, ts_z_yzz_0, ts_z_zzz_0, ts_zz_xx_0, ts_zz_xxx_0, ts_zz_xxy_0, \
                                     ts_zz_xxz_0, ts_zz_xy_0, ts_zz_xyy_0, ts_zz_xyz_0, ts_zz_xz_0, ts_zz_xzz_0, \
                                     ts_zz_yy_0, ts_zz_yyy_0, ts_zz_yyz_0, ts_zz_yz_0, ts_zz_yzz_0, ts_zz_zz_0, \
                                     ts_zz_zzz_0, ts_zzz_xxx_0, ts_zzz_xxy_0, ts_zzz_xxz_0, ts_zzz_xyy_0, ts_zzz_xyz_0, \
                                     ts_zzz_xzz_0, ts_zzz_yyy_0, ts_zzz_yyz_0, ts_zzz_yzz_0, ts_zzz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fx = fx[j];

            ts_xzz_xxx_0[j] = pa_x[j] * ts_zz_xxx_0[j] + 1.5 * fl1_fx * ts_zz_xx_0[j];

            ts_xzz_xxy_0[j] = pa_x[j] * ts_zz_xxy_0[j] + fl1_fx * ts_zz_xy_0[j];

            ts_xzz_xxz_0[j] = pa_x[j] * ts_zz_xxz_0[j] + fl1_fx * ts_zz_xz_0[j];

            ts_xzz_xyy_0[j] = pa_x[j] * ts_zz_xyy_0[j] + 0.5 * fl1_fx * ts_zz_yy_0[j];

            ts_xzz_xyz_0[j] = pa_x[j] * ts_zz_xyz_0[j] + 0.5 * fl1_fx * ts_zz_yz_0[j];

            ts_xzz_xzz_0[j] = pa_x[j] * ts_zz_xzz_0[j] + 0.5 * fl1_fx * ts_zz_zz_0[j];

            ts_xzz_yyy_0[j] = pa_x[j] * ts_zz_yyy_0[j];

            ts_xzz_yyz_0[j] = pa_x[j] * ts_zz_yyz_0[j];

            ts_xzz_yzz_0[j] = pa_x[j] * ts_zz_yzz_0[j];

            ts_xzz_zzz_0[j] = pa_x[j] * ts_zz_zzz_0[j];

            ts_yyy_xxx_0[j] = pa_y[j] * ts_yy_xxx_0[j] + fl1_fx * ts_y_xxx_0[j];

            ts_yyy_xxy_0[j] = pa_y[j] * ts_yy_xxy_0[j] + fl1_fx * ts_y_xxy_0[j] + 0.5 * fl1_fx * ts_yy_xx_0[j];

            ts_yyy_xxz_0[j] = pa_y[j] * ts_yy_xxz_0[j] + fl1_fx * ts_y_xxz_0[j];

            ts_yyy_xyy_0[j] = pa_y[j] * ts_yy_xyy_0[j] + fl1_fx * ts_y_xyy_0[j] + fl1_fx * ts_yy_xy_0[j];

            ts_yyy_xyz_0[j] = pa_y[j] * ts_yy_xyz_0[j] + fl1_fx * ts_y_xyz_0[j] + 0.5 * fl1_fx * ts_yy_xz_0[j];

            ts_yyy_xzz_0[j] = pa_y[j] * ts_yy_xzz_0[j] + fl1_fx * ts_y_xzz_0[j];

            ts_yyy_yyy_0[j] = pa_y[j] * ts_yy_yyy_0[j] + fl1_fx * ts_y_yyy_0[j] + 1.5 * fl1_fx * ts_yy_yy_0[j];

            ts_yyy_yyz_0[j] = pa_y[j] * ts_yy_yyz_0[j] + fl1_fx * ts_y_yyz_0[j] + fl1_fx * ts_yy_yz_0[j];

            ts_yyy_yzz_0[j] = pa_y[j] * ts_yy_yzz_0[j] + fl1_fx * ts_y_yzz_0[j] + 0.5 * fl1_fx * ts_yy_zz_0[j];

            ts_yyy_zzz_0[j] = pa_y[j] * ts_yy_zzz_0[j] + fl1_fx * ts_y_zzz_0[j];

            ts_yyz_xxx_0[j] = pa_y[j] * ts_yz_xxx_0[j] + 0.5 * fl1_fx * ts_z_xxx_0[j];

            ts_yyz_xxy_0[j] = pa_y[j] * ts_yz_xxy_0[j] + 0.5 * fl1_fx * ts_z_xxy_0[j] + 0.5 * fl1_fx * ts_yz_xx_0[j];

            ts_yyz_xxz_0[j] = pa_y[j] * ts_yz_xxz_0[j] + 0.5 * fl1_fx * ts_z_xxz_0[j];

            ts_yyz_xyy_0[j] = pa_y[j] * ts_yz_xyy_0[j] + 0.5 * fl1_fx * ts_z_xyy_0[j] + fl1_fx * ts_yz_xy_0[j];

            ts_yyz_xyz_0[j] = pa_y[j] * ts_yz_xyz_0[j] + 0.5 * fl1_fx * ts_z_xyz_0[j] + 0.5 * fl1_fx * ts_yz_xz_0[j];

            ts_yyz_xzz_0[j] = pa_y[j] * ts_yz_xzz_0[j] + 0.5 * fl1_fx * ts_z_xzz_0[j];

            ts_yyz_yyy_0[j] = pa_y[j] * ts_yz_yyy_0[j] + 0.5 * fl1_fx * ts_z_yyy_0[j] + 1.5 * fl1_fx * ts_yz_yy_0[j];

            ts_yyz_yyz_0[j] = pa_y[j] * ts_yz_yyz_0[j] + 0.5 * fl1_fx * ts_z_yyz_0[j] + fl1_fx * ts_yz_yz_0[j];

            ts_yyz_yzz_0[j] = pa_y[j] * ts_yz_yzz_0[j] + 0.5 * fl1_fx * ts_z_yzz_0[j] + 0.5 * fl1_fx * ts_yz_zz_0[j];

            ts_yyz_zzz_0[j] = pa_y[j] * ts_yz_zzz_0[j] + 0.5 * fl1_fx * ts_z_zzz_0[j];

            ts_yzz_xxx_0[j] = pa_y[j] * ts_zz_xxx_0[j];

            ts_yzz_xxy_0[j] = pa_y[j] * ts_zz_xxy_0[j] + 0.5 * fl1_fx * ts_zz_xx_0[j];

            ts_yzz_xxz_0[j] = pa_y[j] * ts_zz_xxz_0[j];

            ts_yzz_xyy_0[j] = pa_y[j] * ts_zz_xyy_0[j] + fl1_fx * ts_zz_xy_0[j];

            ts_yzz_xyz_0[j] = pa_y[j] * ts_zz_xyz_0[j] + 0.5 * fl1_fx * ts_zz_xz_0[j];

            ts_yzz_xzz_0[j] = pa_y[j] * ts_zz_xzz_0[j];

            ts_yzz_yyy_0[j] = pa_y[j] * ts_zz_yyy_0[j] + 1.5 * fl1_fx * ts_zz_yy_0[j];

            ts_yzz_yyz_0[j] = pa_y[j] * ts_zz_yyz_0[j] + fl1_fx * ts_zz_yz_0[j];

            ts_yzz_yzz_0[j] = pa_y[j] * ts_zz_yzz_0[j] + 0.5 * fl1_fx * ts_zz_zz_0[j];

            ts_yzz_zzz_0[j] = pa_y[j] * ts_zz_zzz_0[j];

            ts_zzz_xxx_0[j] = pa_z[j] * ts_zz_xxx_0[j] + fl1_fx * ts_z_xxx_0[j];

            ts_zzz_xxy_0[j] = pa_z[j] * ts_zz_xxy_0[j] + fl1_fx * ts_z_xxy_0[j];

            ts_zzz_xxz_0[j] = pa_z[j] * ts_zz_xxz_0[j] + fl1_fx * ts_z_xxz_0[j] + 0.5 * fl1_fx * ts_zz_xx_0[j];

            ts_zzz_xyy_0[j] = pa_z[j] * ts_zz_xyy_0[j] + fl1_fx * ts_z_xyy_0[j];

            ts_zzz_xyz_0[j] = pa_z[j] * ts_zz_xyz_0[j] + fl1_fx * ts_z_xyz_0[j] + 0.5 * fl1_fx * ts_zz_xy_0[j];

            ts_zzz_xzz_0[j] = pa_z[j] * ts_zz_xzz_0[j] + fl1_fx * ts_z_xzz_0[j] + fl1_fx * ts_zz_xz_0[j];

            ts_zzz_yyy_0[j] = pa_z[j] * ts_zz_yyy_0[j] + fl1_fx * ts_z_yyy_0[j];

            ts_zzz_yyz_0[j] = pa_z[j] * ts_zz_yyz_0[j] + fl1_fx * ts_z_yyz_0[j] + 0.5 * fl1_fx * ts_zz_yy_0[j];

            ts_zzz_yzz_0[j] = pa_z[j] * ts_zz_yzz_0[j] + fl1_fx * ts_z_yzz_0[j] + fl1_fx * ts_zz_yz_0[j];

            ts_zzz_zzz_0[j] = pa_z[j] * ts_zz_zzz_0[j] + fl1_fx * ts_z_zzz_0[j] + 1.5 * fl1_fx * ts_zz_zz_0[j];
        }

        idx++;
    }
}

}  // namespace ovlrecfunc
