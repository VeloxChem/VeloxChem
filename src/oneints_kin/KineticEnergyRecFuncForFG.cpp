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

#include "KineticEnergyRecFuncForFG.hpp"

namespace kinrecfunc {  // kinrecfunc namespace

void
compKineticEnergyForFG(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    kinrecfunc::compKineticEnergyForFG_0_50(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFG_50_100(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFG_100_150(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compKineticEnergyForFG_0_50(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fz = osFactors.data(4 * idx + 1);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        // set up pointers to auxilary integrals

        auto tt_xx_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx);

        auto tt_xx_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 1);

        auto tt_xx_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 2);

        auto tt_xx_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 3);

        auto tt_xx_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 4);

        auto tt_xx_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 5);

        auto tt_xx_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 6);

        auto tt_xx_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 7);

        auto tt_xx_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 8);

        auto tt_xx_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 9);

        auto tt_xx_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 10);

        auto tt_xx_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 11);

        auto tt_xx_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 12);

        auto tt_xx_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 13);

        auto tt_xx_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 14);

        auto tt_xy_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 15);

        auto tt_xy_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 16);

        auto tt_xy_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 17);

        auto tt_xy_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 18);

        auto tt_xy_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 19);

        auto tt_xy_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 20);

        auto tt_xy_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 21);

        auto tt_xy_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 22);

        auto tt_xy_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 23);

        auto tt_xy_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 24);

        auto tt_xy_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 25);

        auto tt_xy_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 26);

        auto tt_xy_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 27);

        auto tt_xy_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 28);

        auto tt_xy_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 29);

        auto tt_xz_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 30);

        auto tt_xz_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 31);

        auto tt_xz_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 32);

        auto tt_xz_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 33);

        auto tt_xz_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 34);

        auto tt_xz_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 35);

        auto tt_xz_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 36);

        auto tt_xz_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 37);

        auto tt_xz_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 38);

        auto tt_xz_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 39);

        auto tt_xz_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 40);

        auto tt_xz_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 41);

        auto tt_xz_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 42);

        auto tt_xz_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 43);

        auto tt_xz_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 44);

        auto tt_yy_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 45);

        auto tt_yy_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 46);

        auto tt_yy_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 47);

        auto tt_yy_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 48);

        auto tt_yy_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 49);

        auto tt_x_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx);

        auto tt_x_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 1);

        auto tt_x_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 2);

        auto tt_x_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 3);

        auto tt_x_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 4);

        auto tt_x_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 5);

        auto tt_x_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 6);

        auto tt_x_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 7);

        auto tt_x_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 8);

        auto tt_x_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 9);

        auto tt_x_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 10);

        auto tt_x_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 11);

        auto tt_x_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 12);

        auto tt_x_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 13);

        auto tt_x_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 14);

        auto tt_y_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 15);

        auto tt_y_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 16);

        auto tt_y_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 17);

        auto tt_y_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 18);

        auto tt_y_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 19);

        auto tt_y_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 20);

        auto tt_y_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 21);

        auto tt_y_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 22);

        auto tt_y_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 23);

        auto tt_y_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 24);

        auto tt_y_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 25);

        auto tt_y_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 26);

        auto tt_y_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 27);

        auto tt_y_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 28);

        auto tt_y_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 29);

        auto tt_z_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 30);

        auto tt_z_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 31);

        auto tt_z_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 32);

        auto tt_z_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 33);

        auto tt_z_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 34);

        auto tt_z_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 35);

        auto tt_z_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 36);

        auto tt_z_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 37);

        auto tt_z_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 38);

        auto tt_z_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 39);

        auto tt_z_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 40);

        auto tt_z_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 41);

        auto tt_z_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 42);

        auto tt_z_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 43);

        auto tt_z_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 44);

        auto tt_xx_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx);

        auto tt_xx_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 1);

        auto tt_xx_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 2);

        auto tt_xx_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 3);

        auto tt_xx_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 4);

        auto tt_xx_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 5);

        auto tt_xx_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 6);

        auto tt_xx_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 7);

        auto tt_xx_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 8);

        auto tt_xx_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 9);

        auto tt_xy_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 10);

        auto tt_xy_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 11);

        auto tt_xy_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 12);

        auto tt_xy_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 13);

        auto tt_xy_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 14);

        auto tt_xy_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 15);

        auto tt_xy_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 16);

        auto tt_xy_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 17);

        auto tt_xy_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 18);

        auto tt_xy_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 19);

        auto tt_xz_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 20);

        auto tt_xz_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 21);

        auto tt_xz_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 22);

        auto tt_xz_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 23);

        auto tt_xz_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 24);

        auto tt_xz_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 25);

        auto tt_xz_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 26);

        auto tt_xz_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 27);

        auto tt_xz_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 28);

        auto tt_xz_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 29);

        auto tt_yy_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 30);

        auto tt_yy_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 31);

        auto tt_yy_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 32);

        auto tt_yy_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 33);

        auto tt_yy_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 34);

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

        auto ts_xyy_xxyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 49);

        auto ts_x_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx);

        auto ts_x_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 1);

        auto ts_x_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 2);

        auto ts_x_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 3);

        auto ts_x_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 4);

        auto ts_x_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 5);

        auto ts_x_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 6);

        auto ts_x_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 7);

        auto ts_x_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 8);

        auto ts_x_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 9);

        auto ts_x_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 10);

        auto ts_x_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 11);

        auto ts_x_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 12);

        auto ts_x_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 13);

        auto ts_x_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 14);

        auto ts_y_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 15);

        auto ts_y_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 16);

        auto ts_y_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 17);

        auto ts_y_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 18);

        auto ts_y_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 19);

        auto ts_y_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 20);

        auto ts_y_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 21);

        auto ts_y_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 22);

        auto ts_y_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 23);

        auto ts_y_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 24);

        auto ts_y_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 25);

        auto ts_y_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 26);

        auto ts_y_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 27);

        auto ts_y_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 28);

        auto ts_y_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 29);

        auto ts_z_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 30);

        auto ts_z_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 31);

        auto ts_z_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 32);

        auto ts_z_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 33);

        auto ts_z_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 34);

        auto ts_z_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 35);

        auto ts_z_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 36);

        auto ts_z_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 37);

        auto ts_z_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 38);

        auto ts_z_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 39);

        auto ts_z_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 40);

        auto ts_z_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 41);

        auto ts_z_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 42);

        auto ts_z_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 43);

        auto ts_z_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 44);

        // set up pointers to integrals

        auto tt_xxx_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx);

        auto tt_xxx_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 1);

        auto tt_xxx_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 2);

        auto tt_xxx_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 3);

        auto tt_xxx_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 4);

        auto tt_xxx_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 5);

        auto tt_xxx_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 6);

        auto tt_xxx_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 7);

        auto tt_xxx_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 8);

        auto tt_xxx_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 9);

        auto tt_xxx_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 10);

        auto tt_xxx_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 11);

        auto tt_xxx_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 12);

        auto tt_xxx_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 13);

        auto tt_xxx_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 14);

        auto tt_xxy_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 15);

        auto tt_xxy_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 16);

        auto tt_xxy_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 17);

        auto tt_xxy_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 18);

        auto tt_xxy_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 19);

        auto tt_xxy_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 20);

        auto tt_xxy_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 21);

        auto tt_xxy_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 22);

        auto tt_xxy_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 23);

        auto tt_xxy_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 24);

        auto tt_xxy_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 25);

        auto tt_xxy_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 26);

        auto tt_xxy_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 27);

        auto tt_xxy_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 28);

        auto tt_xxy_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 29);

        auto tt_xxz_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 30);

        auto tt_xxz_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 31);

        auto tt_xxz_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 32);

        auto tt_xxz_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 33);

        auto tt_xxz_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 34);

        auto tt_xxz_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 35);

        auto tt_xxz_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 36);

        auto tt_xxz_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 37);

        auto tt_xxz_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 38);

        auto tt_xxz_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 39);

        auto tt_xxz_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 40);

        auto tt_xxz_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 41);

        auto tt_xxz_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 42);

        auto tt_xxz_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 43);

        auto tt_xxz_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 44);

        auto tt_xyy_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 45);

        auto tt_xyy_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 46);

        auto tt_xyy_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 47);

        auto tt_xyy_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 48);

        auto tt_xyy_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 49);

        // Batch of Integrals (0,50)

        #pragma omp simd aligned(fga, fx, fz, pa_x, ts_x_xxxx_0, ts_x_xxxy_0, ts_x_xxxz_0, ts_x_xxyy_0, \
                                     ts_x_xxyz_0, ts_x_xxzz_0, ts_x_xyyy_0, ts_x_xyyz_0, ts_x_xyzz_0, ts_x_xzzz_0, \
                                     ts_x_yyyy_0, ts_x_yyyz_0, ts_x_yyzz_0, ts_x_yzzz_0, ts_x_zzzz_0, ts_xxx_xxxx_0, \
                                     ts_xxx_xxxy_0, ts_xxx_xxxz_0, ts_xxx_xxyy_0, ts_xxx_xxyz_0, ts_xxx_xxzz_0, \
                                     ts_xxx_xyyy_0, ts_xxx_xyyz_0, ts_xxx_xyzz_0, ts_xxx_xzzz_0, ts_xxx_yyyy_0, \
                                     ts_xxx_yyyz_0, ts_xxx_yyzz_0, ts_xxx_yzzz_0, ts_xxx_zzzz_0, ts_xxy_xxxx_0, \
                                     ts_xxy_xxxy_0, ts_xxy_xxxz_0, ts_xxy_xxyy_0, ts_xxy_xxyz_0, ts_xxy_xxzz_0, \
                                     ts_xxy_xyyy_0, ts_xxy_xyyz_0, ts_xxy_xyzz_0, ts_xxy_xzzz_0, ts_xxy_yyyy_0, \
                                     ts_xxy_yyyz_0, ts_xxy_yyzz_0, ts_xxy_yzzz_0, ts_xxy_zzzz_0, ts_xxz_xxxx_0, \
                                     ts_xxz_xxxy_0, ts_xxz_xxxz_0, ts_xxz_xxyy_0, ts_xxz_xxyz_0, ts_xxz_xxzz_0, \
                                     ts_xxz_xyyy_0, ts_xxz_xyyz_0, ts_xxz_xyzz_0, ts_xxz_xzzz_0, ts_xxz_yyyy_0, \
                                     ts_xxz_yyyz_0, ts_xxz_yyzz_0, ts_xxz_yzzz_0, ts_xxz_zzzz_0, ts_xyy_xxxx_0, \
                                     ts_xyy_xxxy_0, ts_xyy_xxxz_0, ts_xyy_xxyy_0, ts_xyy_xxyz_0, ts_y_xxxx_0, \
                                     ts_y_xxxy_0, ts_y_xxxz_0, ts_y_xxyy_0, ts_y_xxyz_0, ts_y_xxzz_0, ts_y_xyyy_0, \
                                     ts_y_xyyz_0, ts_y_xyzz_0, ts_y_xzzz_0, ts_y_yyyy_0, ts_y_yyyz_0, ts_y_yyzz_0, \
                                     ts_y_yzzz_0, ts_y_zzzz_0, ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, \
                                     ts_z_xxyz_0, ts_z_xxzz_0, ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, \
                                     ts_z_yyyy_0, ts_z_yyyz_0, ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0, tt_x_xxxx_0, \
                                     tt_x_xxxy_0, tt_x_xxxz_0, tt_x_xxyy_0, tt_x_xxyz_0, tt_x_xxzz_0, tt_x_xyyy_0, \
                                     tt_x_xyyz_0, tt_x_xyzz_0, tt_x_xzzz_0, tt_x_yyyy_0, tt_x_yyyz_0, tt_x_yyzz_0, \
                                     tt_x_yzzz_0, tt_x_zzzz_0, tt_xx_xxx_0, tt_xx_xxxx_0, tt_xx_xxxy_0, tt_xx_xxxz_0, \
                                     tt_xx_xxy_0, tt_xx_xxyy_0, tt_xx_xxyz_0, tt_xx_xxz_0, tt_xx_xxzz_0, tt_xx_xyy_0, \
                                     tt_xx_xyyy_0, tt_xx_xyyz_0, tt_xx_xyz_0, tt_xx_xyzz_0, tt_xx_xzz_0, tt_xx_xzzz_0, \
                                     tt_xx_yyy_0, tt_xx_yyyy_0, tt_xx_yyyz_0, tt_xx_yyz_0, tt_xx_yyzz_0, tt_xx_yzz_0, \
                                     tt_xx_yzzz_0, tt_xx_zzz_0, tt_xx_zzzz_0, tt_xxx_xxxx_0, tt_xxx_xxxy_0, \
                                     tt_xxx_xxxz_0, tt_xxx_xxyy_0, tt_xxx_xxyz_0, tt_xxx_xxzz_0, tt_xxx_xyyy_0, \
                                     tt_xxx_xyyz_0, tt_xxx_xyzz_0, tt_xxx_xzzz_0, tt_xxx_yyyy_0, tt_xxx_yyyz_0, \
                                     tt_xxx_yyzz_0, tt_xxx_yzzz_0, tt_xxx_zzzz_0, tt_xxy_xxxx_0, tt_xxy_xxxy_0, \
                                     tt_xxy_xxxz_0, tt_xxy_xxyy_0, tt_xxy_xxyz_0, tt_xxy_xxzz_0, tt_xxy_xyyy_0, \
                                     tt_xxy_xyyz_0, tt_xxy_xyzz_0, tt_xxy_xzzz_0, tt_xxy_yyyy_0, tt_xxy_yyyz_0, \
                                     tt_xxy_yyzz_0, tt_xxy_yzzz_0, tt_xxy_zzzz_0, tt_xxz_xxxx_0, tt_xxz_xxxy_0, \
                                     tt_xxz_xxxz_0, tt_xxz_xxyy_0, tt_xxz_xxyz_0, tt_xxz_xxzz_0, tt_xxz_xyyy_0, \
                                     tt_xxz_xyyz_0, tt_xxz_xyzz_0, tt_xxz_xzzz_0, tt_xxz_yyyy_0, tt_xxz_yyyz_0, \
                                     tt_xxz_yyzz_0, tt_xxz_yzzz_0, tt_xxz_zzzz_0, tt_xy_xxx_0, tt_xy_xxxx_0, \
                                     tt_xy_xxxy_0, tt_xy_xxxz_0, tt_xy_xxy_0, tt_xy_xxyy_0, tt_xy_xxyz_0, tt_xy_xxz_0, \
                                     tt_xy_xxzz_0, tt_xy_xyy_0, tt_xy_xyyy_0, tt_xy_xyyz_0, tt_xy_xyz_0, tt_xy_xyzz_0, \
                                     tt_xy_xzz_0, tt_xy_xzzz_0, tt_xy_yyy_0, tt_xy_yyyy_0, tt_xy_yyyz_0, tt_xy_yyz_0, \
                                     tt_xy_yyzz_0, tt_xy_yzz_0, tt_xy_yzzz_0, tt_xy_zzz_0, tt_xy_zzzz_0, tt_xyy_xxxx_0, \
                                     tt_xyy_xxxy_0, tt_xyy_xxxz_0, tt_xyy_xxyy_0, tt_xyy_xxyz_0, tt_xz_xxx_0, \
                                     tt_xz_xxxx_0, tt_xz_xxxy_0, tt_xz_xxxz_0, tt_xz_xxy_0, tt_xz_xxyy_0, tt_xz_xxyz_0, \
                                     tt_xz_xxz_0, tt_xz_xxzz_0, tt_xz_xyy_0, tt_xz_xyyy_0, tt_xz_xyyz_0, tt_xz_xyz_0, \
                                     tt_xz_xyzz_0, tt_xz_xzz_0, tt_xz_xzzz_0, tt_xz_yyy_0, tt_xz_yyyy_0, tt_xz_yyyz_0, \
                                     tt_xz_yyz_0, tt_xz_yyzz_0, tt_xz_yzz_0, tt_xz_yzzz_0, tt_xz_zzz_0, tt_xz_zzzz_0, \
                                     tt_y_xxxx_0, tt_y_xxxy_0, tt_y_xxxz_0, tt_y_xxyy_0, tt_y_xxyz_0, tt_y_xxzz_0, \
                                     tt_y_xyyy_0, tt_y_xyyz_0, tt_y_xyzz_0, tt_y_xzzz_0, tt_y_yyyy_0, tt_y_yyyz_0, \
                                     tt_y_yyzz_0, tt_y_yzzz_0, tt_y_zzzz_0, tt_yy_xxx_0, tt_yy_xxxx_0, tt_yy_xxxy_0, \
                                     tt_yy_xxxz_0, tt_yy_xxy_0, tt_yy_xxyy_0, tt_yy_xxyz_0, tt_yy_xxz_0, tt_yy_xyy_0, \
                                     tt_yy_xyz_0, tt_z_xxxx_0, tt_z_xxxy_0, tt_z_xxxz_0, tt_z_xxyy_0, tt_z_xxyz_0, \
                                     tt_z_xxzz_0, tt_z_xyyy_0, tt_z_xyyz_0, tt_z_xyzz_0, tt_z_xzzz_0, tt_z_yyyy_0, \
                                     tt_z_yyyz_0, tt_z_yyzz_0, tt_z_yzzz_0, tt_z_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xxx_xxxx_0[j] = pa_x[j] * tt_xx_xxxx_0[j] + fl1_fx * tt_x_xxxx_0[j] + 2.0 * fl1_fx * tt_xx_xxx_0[j] + 2.0 * fl1_fz * ts_xxx_xxxx_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xxxx_0[j];

            tt_xxx_xxxy_0[j] = pa_x[j] * tt_xx_xxxy_0[j] + fl1_fx * tt_x_xxxy_0[j] + 1.5 * fl1_fx * tt_xx_xxy_0[j] + 2.0 * fl1_fz * ts_xxx_xxxy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xxxy_0[j];

            tt_xxx_xxxz_0[j] = pa_x[j] * tt_xx_xxxz_0[j] + fl1_fx * tt_x_xxxz_0[j] + 1.5 * fl1_fx * tt_xx_xxz_0[j] + 2.0 * fl1_fz * ts_xxx_xxxz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xxxz_0[j];

            tt_xxx_xxyy_0[j] = pa_x[j] * tt_xx_xxyy_0[j] + fl1_fx * tt_x_xxyy_0[j] + fl1_fx * tt_xx_xyy_0[j] + 2.0 * fl1_fz * ts_xxx_xxyy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xxyy_0[j];

            tt_xxx_xxyz_0[j] = pa_x[j] * tt_xx_xxyz_0[j] + fl1_fx * tt_x_xxyz_0[j] + fl1_fx * tt_xx_xyz_0[j] + 2.0 * fl1_fz * ts_xxx_xxyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xxyz_0[j];

            tt_xxx_xxzz_0[j] = pa_x[j] * tt_xx_xxzz_0[j] + fl1_fx * tt_x_xxzz_0[j] + fl1_fx * tt_xx_xzz_0[j] + 2.0 * fl1_fz * ts_xxx_xxzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xxzz_0[j];

            tt_xxx_xyyy_0[j] = pa_x[j] * tt_xx_xyyy_0[j] + fl1_fx * tt_x_xyyy_0[j] + 0.5 * fl1_fx * tt_xx_yyy_0[j] + 2.0 * fl1_fz * ts_xxx_xyyy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xyyy_0[j];

            tt_xxx_xyyz_0[j] = pa_x[j] * tt_xx_xyyz_0[j] + fl1_fx * tt_x_xyyz_0[j] + 0.5 * fl1_fx * tt_xx_yyz_0[j] + 2.0 * fl1_fz * ts_xxx_xyyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xyyz_0[j];

            tt_xxx_xyzz_0[j] = pa_x[j] * tt_xx_xyzz_0[j] + fl1_fx * tt_x_xyzz_0[j] + 0.5 * fl1_fx * tt_xx_yzz_0[j] + 2.0 * fl1_fz * ts_xxx_xyzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xyzz_0[j];

            tt_xxx_xzzz_0[j] = pa_x[j] * tt_xx_xzzz_0[j] + fl1_fx * tt_x_xzzz_0[j] + 0.5 * fl1_fx * tt_xx_zzz_0[j] + 2.0 * fl1_fz * ts_xxx_xzzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_x_xzzz_0[j];

            tt_xxx_yyyy_0[j] =
                pa_x[j] * tt_xx_yyyy_0[j] + fl1_fx * tt_x_yyyy_0[j] + 2.0 * fl1_fz * ts_xxx_yyyy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_yyyy_0[j];

            tt_xxx_yyyz_0[j] =
                pa_x[j] * tt_xx_yyyz_0[j] + fl1_fx * tt_x_yyyz_0[j] + 2.0 * fl1_fz * ts_xxx_yyyz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_yyyz_0[j];

            tt_xxx_yyzz_0[j] =
                pa_x[j] * tt_xx_yyzz_0[j] + fl1_fx * tt_x_yyzz_0[j] + 2.0 * fl1_fz * ts_xxx_yyzz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_yyzz_0[j];

            tt_xxx_yzzz_0[j] =
                pa_x[j] * tt_xx_yzzz_0[j] + fl1_fx * tt_x_yzzz_0[j] + 2.0 * fl1_fz * ts_xxx_yzzz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_yzzz_0[j];

            tt_xxx_zzzz_0[j] =
                pa_x[j] * tt_xx_zzzz_0[j] + fl1_fx * tt_x_zzzz_0[j] + 2.0 * fl1_fz * ts_xxx_zzzz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_zzzz_0[j];

            tt_xxy_xxxx_0[j] = pa_x[j] * tt_xy_xxxx_0[j] + 0.5 * fl1_fx * tt_y_xxxx_0[j] + 2.0 * fl1_fx * tt_xy_xxx_0[j] +
                               2.0 * fl1_fz * ts_xxy_xxxx_0[j] - fl1_fz * fl1_fga * ts_y_xxxx_0[j];

            tt_xxy_xxxy_0[j] = pa_x[j] * tt_xy_xxxy_0[j] + 0.5 * fl1_fx * tt_y_xxxy_0[j] + 1.5 * fl1_fx * tt_xy_xxy_0[j] +
                               2.0 * fl1_fz * ts_xxy_xxxy_0[j] - fl1_fz * fl1_fga * ts_y_xxxy_0[j];

            tt_xxy_xxxz_0[j] = pa_x[j] * tt_xy_xxxz_0[j] + 0.5 * fl1_fx * tt_y_xxxz_0[j] + 1.5 * fl1_fx * tt_xy_xxz_0[j] +
                               2.0 * fl1_fz * ts_xxy_xxxz_0[j] - fl1_fz * fl1_fga * ts_y_xxxz_0[j];

            tt_xxy_xxyy_0[j] = pa_x[j] * tt_xy_xxyy_0[j] + 0.5 * fl1_fx * tt_y_xxyy_0[j] + fl1_fx * tt_xy_xyy_0[j] + 2.0 * fl1_fz * ts_xxy_xxyy_0[j] -
                               fl1_fz * fl1_fga * ts_y_xxyy_0[j];

            tt_xxy_xxyz_0[j] = pa_x[j] * tt_xy_xxyz_0[j] + 0.5 * fl1_fx * tt_y_xxyz_0[j] + fl1_fx * tt_xy_xyz_0[j] + 2.0 * fl1_fz * ts_xxy_xxyz_0[j] -
                               fl1_fz * fl1_fga * ts_y_xxyz_0[j];

            tt_xxy_xxzz_0[j] = pa_x[j] * tt_xy_xxzz_0[j] + 0.5 * fl1_fx * tt_y_xxzz_0[j] + fl1_fx * tt_xy_xzz_0[j] + 2.0 * fl1_fz * ts_xxy_xxzz_0[j] -
                               fl1_fz * fl1_fga * ts_y_xxzz_0[j];

            tt_xxy_xyyy_0[j] = pa_x[j] * tt_xy_xyyy_0[j] + 0.5 * fl1_fx * tt_y_xyyy_0[j] + 0.5 * fl1_fx * tt_xy_yyy_0[j] +
                               2.0 * fl1_fz * ts_xxy_xyyy_0[j] - fl1_fz * fl1_fga * ts_y_xyyy_0[j];

            tt_xxy_xyyz_0[j] = pa_x[j] * tt_xy_xyyz_0[j] + 0.5 * fl1_fx * tt_y_xyyz_0[j] + 0.5 * fl1_fx * tt_xy_yyz_0[j] +
                               2.0 * fl1_fz * ts_xxy_xyyz_0[j] - fl1_fz * fl1_fga * ts_y_xyyz_0[j];

            tt_xxy_xyzz_0[j] = pa_x[j] * tt_xy_xyzz_0[j] + 0.5 * fl1_fx * tt_y_xyzz_0[j] + 0.5 * fl1_fx * tt_xy_yzz_0[j] +
                               2.0 * fl1_fz * ts_xxy_xyzz_0[j] - fl1_fz * fl1_fga * ts_y_xyzz_0[j];

            tt_xxy_xzzz_0[j] = pa_x[j] * tt_xy_xzzz_0[j] + 0.5 * fl1_fx * tt_y_xzzz_0[j] + 0.5 * fl1_fx * tt_xy_zzz_0[j] +
                               2.0 * fl1_fz * ts_xxy_xzzz_0[j] - fl1_fz * fl1_fga * ts_y_xzzz_0[j];

            tt_xxy_yyyy_0[j] =
                pa_x[j] * tt_xy_yyyy_0[j] + 0.5 * fl1_fx * tt_y_yyyy_0[j] + 2.0 * fl1_fz * ts_xxy_yyyy_0[j] - fl1_fz * fl1_fga * ts_y_yyyy_0[j];

            tt_xxy_yyyz_0[j] =
                pa_x[j] * tt_xy_yyyz_0[j] + 0.5 * fl1_fx * tt_y_yyyz_0[j] + 2.0 * fl1_fz * ts_xxy_yyyz_0[j] - fl1_fz * fl1_fga * ts_y_yyyz_0[j];

            tt_xxy_yyzz_0[j] =
                pa_x[j] * tt_xy_yyzz_0[j] + 0.5 * fl1_fx * tt_y_yyzz_0[j] + 2.0 * fl1_fz * ts_xxy_yyzz_0[j] - fl1_fz * fl1_fga * ts_y_yyzz_0[j];

            tt_xxy_yzzz_0[j] =
                pa_x[j] * tt_xy_yzzz_0[j] + 0.5 * fl1_fx * tt_y_yzzz_0[j] + 2.0 * fl1_fz * ts_xxy_yzzz_0[j] - fl1_fz * fl1_fga * ts_y_yzzz_0[j];

            tt_xxy_zzzz_0[j] =
                pa_x[j] * tt_xy_zzzz_0[j] + 0.5 * fl1_fx * tt_y_zzzz_0[j] + 2.0 * fl1_fz * ts_xxy_zzzz_0[j] - fl1_fz * fl1_fga * ts_y_zzzz_0[j];

            tt_xxz_xxxx_0[j] = pa_x[j] * tt_xz_xxxx_0[j] + 0.5 * fl1_fx * tt_z_xxxx_0[j] + 2.0 * fl1_fx * tt_xz_xxx_0[j] +
                               2.0 * fl1_fz * ts_xxz_xxxx_0[j] - fl1_fz * fl1_fga * ts_z_xxxx_0[j];

            tt_xxz_xxxy_0[j] = pa_x[j] * tt_xz_xxxy_0[j] + 0.5 * fl1_fx * tt_z_xxxy_0[j] + 1.5 * fl1_fx * tt_xz_xxy_0[j] +
                               2.0 * fl1_fz * ts_xxz_xxxy_0[j] - fl1_fz * fl1_fga * ts_z_xxxy_0[j];

            tt_xxz_xxxz_0[j] = pa_x[j] * tt_xz_xxxz_0[j] + 0.5 * fl1_fx * tt_z_xxxz_0[j] + 1.5 * fl1_fx * tt_xz_xxz_0[j] +
                               2.0 * fl1_fz * ts_xxz_xxxz_0[j] - fl1_fz * fl1_fga * ts_z_xxxz_0[j];

            tt_xxz_xxyy_0[j] = pa_x[j] * tt_xz_xxyy_0[j] + 0.5 * fl1_fx * tt_z_xxyy_0[j] + fl1_fx * tt_xz_xyy_0[j] + 2.0 * fl1_fz * ts_xxz_xxyy_0[j] -
                               fl1_fz * fl1_fga * ts_z_xxyy_0[j];

            tt_xxz_xxyz_0[j] = pa_x[j] * tt_xz_xxyz_0[j] + 0.5 * fl1_fx * tt_z_xxyz_0[j] + fl1_fx * tt_xz_xyz_0[j] + 2.0 * fl1_fz * ts_xxz_xxyz_0[j] -
                               fl1_fz * fl1_fga * ts_z_xxyz_0[j];

            tt_xxz_xxzz_0[j] = pa_x[j] * tt_xz_xxzz_0[j] + 0.5 * fl1_fx * tt_z_xxzz_0[j] + fl1_fx * tt_xz_xzz_0[j] + 2.0 * fl1_fz * ts_xxz_xxzz_0[j] -
                               fl1_fz * fl1_fga * ts_z_xxzz_0[j];

            tt_xxz_xyyy_0[j] = pa_x[j] * tt_xz_xyyy_0[j] + 0.5 * fl1_fx * tt_z_xyyy_0[j] + 0.5 * fl1_fx * tt_xz_yyy_0[j] +
                               2.0 * fl1_fz * ts_xxz_xyyy_0[j] - fl1_fz * fl1_fga * ts_z_xyyy_0[j];

            tt_xxz_xyyz_0[j] = pa_x[j] * tt_xz_xyyz_0[j] + 0.5 * fl1_fx * tt_z_xyyz_0[j] + 0.5 * fl1_fx * tt_xz_yyz_0[j] +
                               2.0 * fl1_fz * ts_xxz_xyyz_0[j] - fl1_fz * fl1_fga * ts_z_xyyz_0[j];

            tt_xxz_xyzz_0[j] = pa_x[j] * tt_xz_xyzz_0[j] + 0.5 * fl1_fx * tt_z_xyzz_0[j] + 0.5 * fl1_fx * tt_xz_yzz_0[j] +
                               2.0 * fl1_fz * ts_xxz_xyzz_0[j] - fl1_fz * fl1_fga * ts_z_xyzz_0[j];

            tt_xxz_xzzz_0[j] = pa_x[j] * tt_xz_xzzz_0[j] + 0.5 * fl1_fx * tt_z_xzzz_0[j] + 0.5 * fl1_fx * tt_xz_zzz_0[j] +
                               2.0 * fl1_fz * ts_xxz_xzzz_0[j] - fl1_fz * fl1_fga * ts_z_xzzz_0[j];

            tt_xxz_yyyy_0[j] =
                pa_x[j] * tt_xz_yyyy_0[j] + 0.5 * fl1_fx * tt_z_yyyy_0[j] + 2.0 * fl1_fz * ts_xxz_yyyy_0[j] - fl1_fz * fl1_fga * ts_z_yyyy_0[j];

            tt_xxz_yyyz_0[j] =
                pa_x[j] * tt_xz_yyyz_0[j] + 0.5 * fl1_fx * tt_z_yyyz_0[j] + 2.0 * fl1_fz * ts_xxz_yyyz_0[j] - fl1_fz * fl1_fga * ts_z_yyyz_0[j];

            tt_xxz_yyzz_0[j] =
                pa_x[j] * tt_xz_yyzz_0[j] + 0.5 * fl1_fx * tt_z_yyzz_0[j] + 2.0 * fl1_fz * ts_xxz_yyzz_0[j] - fl1_fz * fl1_fga * ts_z_yyzz_0[j];

            tt_xxz_yzzz_0[j] =
                pa_x[j] * tt_xz_yzzz_0[j] + 0.5 * fl1_fx * tt_z_yzzz_0[j] + 2.0 * fl1_fz * ts_xxz_yzzz_0[j] - fl1_fz * fl1_fga * ts_z_yzzz_0[j];

            tt_xxz_zzzz_0[j] =
                pa_x[j] * tt_xz_zzzz_0[j] + 0.5 * fl1_fx * tt_z_zzzz_0[j] + 2.0 * fl1_fz * ts_xxz_zzzz_0[j] - fl1_fz * fl1_fga * ts_z_zzzz_0[j];

            tt_xyy_xxxx_0[j] = pa_x[j] * tt_yy_xxxx_0[j] + 2.0 * fl1_fx * tt_yy_xxx_0[j] + 2.0 * fl1_fz * ts_xyy_xxxx_0[j];

            tt_xyy_xxxy_0[j] = pa_x[j] * tt_yy_xxxy_0[j] + 1.5 * fl1_fx * tt_yy_xxy_0[j] + 2.0 * fl1_fz * ts_xyy_xxxy_0[j];

            tt_xyy_xxxz_0[j] = pa_x[j] * tt_yy_xxxz_0[j] + 1.5 * fl1_fx * tt_yy_xxz_0[j] + 2.0 * fl1_fz * ts_xyy_xxxz_0[j];

            tt_xyy_xxyy_0[j] = pa_x[j] * tt_yy_xxyy_0[j] + fl1_fx * tt_yy_xyy_0[j] + 2.0 * fl1_fz * ts_xyy_xxyy_0[j];

            tt_xyy_xxyz_0[j] = pa_x[j] * tt_yy_xxyz_0[j] + fl1_fx * tt_yy_xyz_0[j] + 2.0 * fl1_fz * ts_xyy_xxyz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForFG_50_100(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fz = osFactors.data(4 * idx + 1);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_x = paDistances.data(3 * idx);

        auto pa_y = paDistances.data(3 * idx + 1);

        // set up pointers to auxilary integrals

        auto tt_yy_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 45);

        auto tt_yy_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 46);

        auto tt_yy_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 47);

        auto tt_yy_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 48);

        auto tt_yy_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 49);

        auto tt_yy_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 50);

        auto tt_yy_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 51);

        auto tt_yy_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 52);

        auto tt_yy_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 53);

        auto tt_yy_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 54);

        auto tt_yy_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 55);

        auto tt_yy_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 56);

        auto tt_yy_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 57);

        auto tt_yy_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 58);

        auto tt_yy_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 59);

        auto tt_yz_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 60);

        auto tt_yz_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 61);

        auto tt_yz_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 62);

        auto tt_yz_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 63);

        auto tt_yz_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 64);

        auto tt_yz_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 65);

        auto tt_yz_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 66);

        auto tt_yz_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 67);

        auto tt_yz_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 68);

        auto tt_yz_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 69);

        auto tt_yz_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 70);

        auto tt_yz_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 71);

        auto tt_yz_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 72);

        auto tt_yz_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 73);

        auto tt_yz_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 74);

        auto tt_zz_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 75);

        auto tt_zz_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 76);

        auto tt_zz_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 77);

        auto tt_zz_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 78);

        auto tt_zz_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 79);

        auto tt_zz_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 80);

        auto tt_zz_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 81);

        auto tt_zz_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 82);

        auto tt_zz_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 83);

        auto tt_zz_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 84);

        auto tt_zz_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 85);

        auto tt_zz_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 86);

        auto tt_zz_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 87);

        auto tt_zz_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 88);

        auto tt_zz_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 89);

        auto tt_y_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 15);

        auto tt_y_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 16);

        auto tt_y_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 17);

        auto tt_y_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 18);

        auto tt_y_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 19);

        auto tt_y_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 20);

        auto tt_y_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 21);

        auto tt_y_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 22);

        auto tt_y_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 23);

        auto tt_y_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 24);

        auto tt_yy_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 30);

        auto tt_yy_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 31);

        auto tt_yy_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 32);

        auto tt_yy_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 33);

        auto tt_yy_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 34);

        auto tt_yy_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 35);

        auto tt_yy_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 36);

        auto tt_yy_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 37);

        auto tt_yy_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 38);

        auto tt_yy_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 39);

        auto tt_yz_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 40);

        auto tt_yz_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 41);

        auto tt_yz_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 42);

        auto tt_yz_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 43);

        auto tt_yz_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 44);

        auto tt_yz_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 45);

        auto tt_yz_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 46);

        auto tt_yz_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 47);

        auto tt_yz_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 48);

        auto tt_yz_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 49);

        auto tt_zz_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 50);

        auto tt_zz_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 51);

        auto tt_zz_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 52);

        auto tt_zz_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 53);

        auto tt_zz_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 54);

        auto tt_zz_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 55);

        auto tt_zz_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 56);

        auto tt_zz_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 57);

        auto tt_zz_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 58);

        auto tt_zz_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 59);

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

        auto ts_yyy_xyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 97);

        auto ts_yyy_xyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 98);

        auto ts_yyy_xzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 99);

        auto ts_y_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 15);

        auto ts_y_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 16);

        auto ts_y_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 17);

        auto ts_y_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 18);

        auto ts_y_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 19);

        auto ts_y_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 20);

        auto ts_y_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 21);

        auto ts_y_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 22);

        auto ts_y_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 23);

        auto ts_y_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 24);

        // set up pointers to integrals

        auto tt_xyy_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 50);

        auto tt_xyy_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 51);

        auto tt_xyy_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 52);

        auto tt_xyy_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 53);

        auto tt_xyy_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 54);

        auto tt_xyy_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 55);

        auto tt_xyy_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 56);

        auto tt_xyy_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 57);

        auto tt_xyy_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 58);

        auto tt_xyy_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 59);

        auto tt_xyz_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 60);

        auto tt_xyz_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 61);

        auto tt_xyz_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 62);

        auto tt_xyz_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 63);

        auto tt_xyz_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 64);

        auto tt_xyz_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 65);

        auto tt_xyz_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 66);

        auto tt_xyz_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 67);

        auto tt_xyz_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 68);

        auto tt_xyz_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 69);

        auto tt_xyz_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 70);

        auto tt_xyz_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 71);

        auto tt_xyz_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 72);

        auto tt_xyz_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 73);

        auto tt_xyz_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 74);

        auto tt_xzz_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 75);

        auto tt_xzz_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 76);

        auto tt_xzz_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 77);

        auto tt_xzz_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 78);

        auto tt_xzz_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 79);

        auto tt_xzz_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 80);

        auto tt_xzz_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 81);

        auto tt_xzz_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 82);

        auto tt_xzz_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 83);

        auto tt_xzz_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 84);

        auto tt_xzz_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 85);

        auto tt_xzz_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 86);

        auto tt_xzz_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 87);

        auto tt_xzz_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 88);

        auto tt_xzz_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 89);

        auto tt_yyy_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 90);

        auto tt_yyy_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 91);

        auto tt_yyy_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 92);

        auto tt_yyy_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 93);

        auto tt_yyy_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 94);

        auto tt_yyy_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 95);

        auto tt_yyy_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 96);

        auto tt_yyy_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 97);

        auto tt_yyy_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 98);

        auto tt_yyy_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 99);

        // Batch of Integrals (50,100)

        #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, ts_xyy_xxzz_0, ts_xyy_xyyy_0, ts_xyy_xyyz_0, \
                                     ts_xyy_xyzz_0, ts_xyy_xzzz_0, ts_xyy_yyyy_0, ts_xyy_yyyz_0, ts_xyy_yyzz_0, \
                                     ts_xyy_yzzz_0, ts_xyy_zzzz_0, ts_xyz_xxxx_0, ts_xyz_xxxy_0, ts_xyz_xxxz_0, \
                                     ts_xyz_xxyy_0, ts_xyz_xxyz_0, ts_xyz_xxzz_0, ts_xyz_xyyy_0, ts_xyz_xyyz_0, \
                                     ts_xyz_xyzz_0, ts_xyz_xzzz_0, ts_xyz_yyyy_0, ts_xyz_yyyz_0, ts_xyz_yyzz_0, \
                                     ts_xyz_yzzz_0, ts_xyz_zzzz_0, ts_xzz_xxxx_0, ts_xzz_xxxy_0, ts_xzz_xxxz_0, \
                                     ts_xzz_xxyy_0, ts_xzz_xxyz_0, ts_xzz_xxzz_0, ts_xzz_xyyy_0, ts_xzz_xyyz_0, \
                                     ts_xzz_xyzz_0, ts_xzz_xzzz_0, ts_xzz_yyyy_0, ts_xzz_yyyz_0, ts_xzz_yyzz_0, \
                                     ts_xzz_yzzz_0, ts_xzz_zzzz_0, ts_y_xxxx_0, ts_y_xxxy_0, ts_y_xxxz_0, ts_y_xxyy_0, \
                                     ts_y_xxyz_0, ts_y_xxzz_0, ts_y_xyyy_0, ts_y_xyyz_0, ts_y_xyzz_0, ts_y_xzzz_0, \
                                     ts_yyy_xxxx_0, ts_yyy_xxxy_0, ts_yyy_xxxz_0, ts_yyy_xxyy_0, ts_yyy_xxyz_0, \
                                     ts_yyy_xxzz_0, ts_yyy_xyyy_0, ts_yyy_xyyz_0, ts_yyy_xyzz_0, ts_yyy_xzzz_0, \
                                     tt_xyy_xxzz_0, tt_xyy_xyyy_0, tt_xyy_xyyz_0, tt_xyy_xyzz_0, tt_xyy_xzzz_0, \
                                     tt_xyy_yyyy_0, tt_xyy_yyyz_0, tt_xyy_yyzz_0, tt_xyy_yzzz_0, tt_xyy_zzzz_0, \
                                     tt_xyz_xxxx_0, tt_xyz_xxxy_0, tt_xyz_xxxz_0, tt_xyz_xxyy_0, tt_xyz_xxyz_0, \
                                     tt_xyz_xxzz_0, tt_xyz_xyyy_0, tt_xyz_xyyz_0, tt_xyz_xyzz_0, tt_xyz_xzzz_0, \
                                     tt_xyz_yyyy_0, tt_xyz_yyyz_0, tt_xyz_yyzz_0, tt_xyz_yzzz_0, tt_xyz_zzzz_0, \
                                     tt_xzz_xxxx_0, tt_xzz_xxxy_0, tt_xzz_xxxz_0, tt_xzz_xxyy_0, tt_xzz_xxyz_0, \
                                     tt_xzz_xxzz_0, tt_xzz_xyyy_0, tt_xzz_xyyz_0, tt_xzz_xyzz_0, tt_xzz_xzzz_0, \
                                     tt_xzz_yyyy_0, tt_xzz_yyyz_0, tt_xzz_yyzz_0, tt_xzz_yzzz_0, tt_xzz_zzzz_0, \
                                     tt_y_xxxx_0, tt_y_xxxy_0, tt_y_xxxz_0, tt_y_xxyy_0, tt_y_xxyz_0, tt_y_xxzz_0, \
                                     tt_y_xyyy_0, tt_y_xyyz_0, tt_y_xyzz_0, tt_y_xzzz_0, tt_yy_xxx_0, tt_yy_xxxx_0, \
                                     tt_yy_xxxy_0, tt_yy_xxxz_0, tt_yy_xxy_0, tt_yy_xxyy_0, tt_yy_xxyz_0, tt_yy_xxz_0, \
                                     tt_yy_xxzz_0, tt_yy_xyy_0, tt_yy_xyyy_0, tt_yy_xyyz_0, tt_yy_xyz_0, tt_yy_xyzz_0, \
                                     tt_yy_xzz_0, tt_yy_xzzz_0, tt_yy_yyy_0, tt_yy_yyyy_0, tt_yy_yyyz_0, tt_yy_yyz_0, \
                                     tt_yy_yyzz_0, tt_yy_yzz_0, tt_yy_yzzz_0, tt_yy_zzz_0, tt_yy_zzzz_0, tt_yyy_xxxx_0, \
                                     tt_yyy_xxxy_0, tt_yyy_xxxz_0, tt_yyy_xxyy_0, tt_yyy_xxyz_0, tt_yyy_xxzz_0, \
                                     tt_yyy_xyyy_0, tt_yyy_xyyz_0, tt_yyy_xyzz_0, tt_yyy_xzzz_0, tt_yz_xxx_0, \
                                     tt_yz_xxxx_0, tt_yz_xxxy_0, tt_yz_xxxz_0, tt_yz_xxy_0, tt_yz_xxyy_0, tt_yz_xxyz_0, \
                                     tt_yz_xxz_0, tt_yz_xxzz_0, tt_yz_xyy_0, tt_yz_xyyy_0, tt_yz_xyyz_0, tt_yz_xyz_0, \
                                     tt_yz_xyzz_0, tt_yz_xzz_0, tt_yz_xzzz_0, tt_yz_yyy_0, tt_yz_yyyy_0, tt_yz_yyyz_0, \
                                     tt_yz_yyz_0, tt_yz_yyzz_0, tt_yz_yzz_0, tt_yz_yzzz_0, tt_yz_zzz_0, tt_yz_zzzz_0, \
                                     tt_zz_xxx_0, tt_zz_xxxx_0, tt_zz_xxxy_0, tt_zz_xxxz_0, tt_zz_xxy_0, tt_zz_xxyy_0, \
                                     tt_zz_xxyz_0, tt_zz_xxz_0, tt_zz_xxzz_0, tt_zz_xyy_0, tt_zz_xyyy_0, tt_zz_xyyz_0, \
                                     tt_zz_xyz_0, tt_zz_xyzz_0, tt_zz_xzz_0, tt_zz_xzzz_0, tt_zz_yyy_0, tt_zz_yyyy_0, \
                                     tt_zz_yyyz_0, tt_zz_yyz_0, tt_zz_yyzz_0, tt_zz_yzz_0, tt_zz_yzzz_0, tt_zz_zzz_0, \
                                     tt_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xyy_xxzz_0[j] = pa_x[j] * tt_yy_xxzz_0[j] + fl1_fx * tt_yy_xzz_0[j] + 2.0 * fl1_fz * ts_xyy_xxzz_0[j];

            tt_xyy_xyyy_0[j] = pa_x[j] * tt_yy_xyyy_0[j] + 0.5 * fl1_fx * tt_yy_yyy_0[j] + 2.0 * fl1_fz * ts_xyy_xyyy_0[j];

            tt_xyy_xyyz_0[j] = pa_x[j] * tt_yy_xyyz_0[j] + 0.5 * fl1_fx * tt_yy_yyz_0[j] + 2.0 * fl1_fz * ts_xyy_xyyz_0[j];

            tt_xyy_xyzz_0[j] = pa_x[j] * tt_yy_xyzz_0[j] + 0.5 * fl1_fx * tt_yy_yzz_0[j] + 2.0 * fl1_fz * ts_xyy_xyzz_0[j];

            tt_xyy_xzzz_0[j] = pa_x[j] * tt_yy_xzzz_0[j] + 0.5 * fl1_fx * tt_yy_zzz_0[j] + 2.0 * fl1_fz * ts_xyy_xzzz_0[j];

            tt_xyy_yyyy_0[j] = pa_x[j] * tt_yy_yyyy_0[j] + 2.0 * fl1_fz * ts_xyy_yyyy_0[j];

            tt_xyy_yyyz_0[j] = pa_x[j] * tt_yy_yyyz_0[j] + 2.0 * fl1_fz * ts_xyy_yyyz_0[j];

            tt_xyy_yyzz_0[j] = pa_x[j] * tt_yy_yyzz_0[j] + 2.0 * fl1_fz * ts_xyy_yyzz_0[j];

            tt_xyy_yzzz_0[j] = pa_x[j] * tt_yy_yzzz_0[j] + 2.0 * fl1_fz * ts_xyy_yzzz_0[j];

            tt_xyy_zzzz_0[j] = pa_x[j] * tt_yy_zzzz_0[j] + 2.0 * fl1_fz * ts_xyy_zzzz_0[j];

            tt_xyz_xxxx_0[j] = pa_x[j] * tt_yz_xxxx_0[j] + 2.0 * fl1_fx * tt_yz_xxx_0[j] + 2.0 * fl1_fz * ts_xyz_xxxx_0[j];

            tt_xyz_xxxy_0[j] = pa_x[j] * tt_yz_xxxy_0[j] + 1.5 * fl1_fx * tt_yz_xxy_0[j] + 2.0 * fl1_fz * ts_xyz_xxxy_0[j];

            tt_xyz_xxxz_0[j] = pa_x[j] * tt_yz_xxxz_0[j] + 1.5 * fl1_fx * tt_yz_xxz_0[j] + 2.0 * fl1_fz * ts_xyz_xxxz_0[j];

            tt_xyz_xxyy_0[j] = pa_x[j] * tt_yz_xxyy_0[j] + fl1_fx * tt_yz_xyy_0[j] + 2.0 * fl1_fz * ts_xyz_xxyy_0[j];

            tt_xyz_xxyz_0[j] = pa_x[j] * tt_yz_xxyz_0[j] + fl1_fx * tt_yz_xyz_0[j] + 2.0 * fl1_fz * ts_xyz_xxyz_0[j];

            tt_xyz_xxzz_0[j] = pa_x[j] * tt_yz_xxzz_0[j] + fl1_fx * tt_yz_xzz_0[j] + 2.0 * fl1_fz * ts_xyz_xxzz_0[j];

            tt_xyz_xyyy_0[j] = pa_x[j] * tt_yz_xyyy_0[j] + 0.5 * fl1_fx * tt_yz_yyy_0[j] + 2.0 * fl1_fz * ts_xyz_xyyy_0[j];

            tt_xyz_xyyz_0[j] = pa_x[j] * tt_yz_xyyz_0[j] + 0.5 * fl1_fx * tt_yz_yyz_0[j] + 2.0 * fl1_fz * ts_xyz_xyyz_0[j];

            tt_xyz_xyzz_0[j] = pa_x[j] * tt_yz_xyzz_0[j] + 0.5 * fl1_fx * tt_yz_yzz_0[j] + 2.0 * fl1_fz * ts_xyz_xyzz_0[j];

            tt_xyz_xzzz_0[j] = pa_x[j] * tt_yz_xzzz_0[j] + 0.5 * fl1_fx * tt_yz_zzz_0[j] + 2.0 * fl1_fz * ts_xyz_xzzz_0[j];

            tt_xyz_yyyy_0[j] = pa_x[j] * tt_yz_yyyy_0[j] + 2.0 * fl1_fz * ts_xyz_yyyy_0[j];

            tt_xyz_yyyz_0[j] = pa_x[j] * tt_yz_yyyz_0[j] + 2.0 * fl1_fz * ts_xyz_yyyz_0[j];

            tt_xyz_yyzz_0[j] = pa_x[j] * tt_yz_yyzz_0[j] + 2.0 * fl1_fz * ts_xyz_yyzz_0[j];

            tt_xyz_yzzz_0[j] = pa_x[j] * tt_yz_yzzz_0[j] + 2.0 * fl1_fz * ts_xyz_yzzz_0[j];

            tt_xyz_zzzz_0[j] = pa_x[j] * tt_yz_zzzz_0[j] + 2.0 * fl1_fz * ts_xyz_zzzz_0[j];

            tt_xzz_xxxx_0[j] = pa_x[j] * tt_zz_xxxx_0[j] + 2.0 * fl1_fx * tt_zz_xxx_0[j] + 2.0 * fl1_fz * ts_xzz_xxxx_0[j];

            tt_xzz_xxxy_0[j] = pa_x[j] * tt_zz_xxxy_0[j] + 1.5 * fl1_fx * tt_zz_xxy_0[j] + 2.0 * fl1_fz * ts_xzz_xxxy_0[j];

            tt_xzz_xxxz_0[j] = pa_x[j] * tt_zz_xxxz_0[j] + 1.5 * fl1_fx * tt_zz_xxz_0[j] + 2.0 * fl1_fz * ts_xzz_xxxz_0[j];

            tt_xzz_xxyy_0[j] = pa_x[j] * tt_zz_xxyy_0[j] + fl1_fx * tt_zz_xyy_0[j] + 2.0 * fl1_fz * ts_xzz_xxyy_0[j];

            tt_xzz_xxyz_0[j] = pa_x[j] * tt_zz_xxyz_0[j] + fl1_fx * tt_zz_xyz_0[j] + 2.0 * fl1_fz * ts_xzz_xxyz_0[j];

            tt_xzz_xxzz_0[j] = pa_x[j] * tt_zz_xxzz_0[j] + fl1_fx * tt_zz_xzz_0[j] + 2.0 * fl1_fz * ts_xzz_xxzz_0[j];

            tt_xzz_xyyy_0[j] = pa_x[j] * tt_zz_xyyy_0[j] + 0.5 * fl1_fx * tt_zz_yyy_0[j] + 2.0 * fl1_fz * ts_xzz_xyyy_0[j];

            tt_xzz_xyyz_0[j] = pa_x[j] * tt_zz_xyyz_0[j] + 0.5 * fl1_fx * tt_zz_yyz_0[j] + 2.0 * fl1_fz * ts_xzz_xyyz_0[j];

            tt_xzz_xyzz_0[j] = pa_x[j] * tt_zz_xyzz_0[j] + 0.5 * fl1_fx * tt_zz_yzz_0[j] + 2.0 * fl1_fz * ts_xzz_xyzz_0[j];

            tt_xzz_xzzz_0[j] = pa_x[j] * tt_zz_xzzz_0[j] + 0.5 * fl1_fx * tt_zz_zzz_0[j] + 2.0 * fl1_fz * ts_xzz_xzzz_0[j];

            tt_xzz_yyyy_0[j] = pa_x[j] * tt_zz_yyyy_0[j] + 2.0 * fl1_fz * ts_xzz_yyyy_0[j];

            tt_xzz_yyyz_0[j] = pa_x[j] * tt_zz_yyyz_0[j] + 2.0 * fl1_fz * ts_xzz_yyyz_0[j];

            tt_xzz_yyzz_0[j] = pa_x[j] * tt_zz_yyzz_0[j] + 2.0 * fl1_fz * ts_xzz_yyzz_0[j];

            tt_xzz_yzzz_0[j] = pa_x[j] * tt_zz_yzzz_0[j] + 2.0 * fl1_fz * ts_xzz_yzzz_0[j];

            tt_xzz_zzzz_0[j] = pa_x[j] * tt_zz_zzzz_0[j] + 2.0 * fl1_fz * ts_xzz_zzzz_0[j];

            tt_yyy_xxxx_0[j] =
                pa_y[j] * tt_yy_xxxx_0[j] + fl1_fx * tt_y_xxxx_0[j] + 2.0 * fl1_fz * ts_yyy_xxxx_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_xxxx_0[j];

            tt_yyy_xxxy_0[j] = pa_y[j] * tt_yy_xxxy_0[j] + fl1_fx * tt_y_xxxy_0[j] + 0.5 * fl1_fx * tt_yy_xxx_0[j] + 2.0 * fl1_fz * ts_yyy_xxxy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_xxxy_0[j];

            tt_yyy_xxxz_0[j] =
                pa_y[j] * tt_yy_xxxz_0[j] + fl1_fx * tt_y_xxxz_0[j] + 2.0 * fl1_fz * ts_yyy_xxxz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_xxxz_0[j];

            tt_yyy_xxyy_0[j] = pa_y[j] * tt_yy_xxyy_0[j] + fl1_fx * tt_y_xxyy_0[j] + fl1_fx * tt_yy_xxy_0[j] + 2.0 * fl1_fz * ts_yyy_xxyy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_xxyy_0[j];

            tt_yyy_xxyz_0[j] = pa_y[j] * tt_yy_xxyz_0[j] + fl1_fx * tt_y_xxyz_0[j] + 0.5 * fl1_fx * tt_yy_xxz_0[j] + 2.0 * fl1_fz * ts_yyy_xxyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_xxyz_0[j];

            tt_yyy_xxzz_0[j] =
                pa_y[j] * tt_yy_xxzz_0[j] + fl1_fx * tt_y_xxzz_0[j] + 2.0 * fl1_fz * ts_yyy_xxzz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_xxzz_0[j];

            tt_yyy_xyyy_0[j] = pa_y[j] * tt_yy_xyyy_0[j] + fl1_fx * tt_y_xyyy_0[j] + 1.5 * fl1_fx * tt_yy_xyy_0[j] + 2.0 * fl1_fz * ts_yyy_xyyy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_xyyy_0[j];

            tt_yyy_xyyz_0[j] = pa_y[j] * tt_yy_xyyz_0[j] + fl1_fx * tt_y_xyyz_0[j] + fl1_fx * tt_yy_xyz_0[j] + 2.0 * fl1_fz * ts_yyy_xyyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_xyyz_0[j];

            tt_yyy_xyzz_0[j] = pa_y[j] * tt_yy_xyzz_0[j] + fl1_fx * tt_y_xyzz_0[j] + 0.5 * fl1_fx * tt_yy_xzz_0[j] + 2.0 * fl1_fz * ts_yyy_xyzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_xyzz_0[j];

            tt_yyy_xzzz_0[j] =
                pa_y[j] * tt_yy_xzzz_0[j] + fl1_fx * tt_y_xzzz_0[j] + 2.0 * fl1_fz * ts_yyy_xzzz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_xzzz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForFG_100_150(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_3_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fz = osFactors.data(4 * idx + 1);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_y = paDistances.data(3 * idx + 1);

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tt_yy_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 55);

        auto tt_yy_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 56);

        auto tt_yy_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 57);

        auto tt_yy_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 58);

        auto tt_yy_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 59);

        auto tt_yz_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 60);

        auto tt_yz_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 61);

        auto tt_yz_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 62);

        auto tt_yz_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 63);

        auto tt_yz_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 64);

        auto tt_yz_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 65);

        auto tt_yz_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 66);

        auto tt_yz_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 67);

        auto tt_yz_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 68);

        auto tt_yz_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 69);

        auto tt_yz_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 70);

        auto tt_yz_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 71);

        auto tt_yz_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 72);

        auto tt_yz_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 73);

        auto tt_yz_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 74);

        auto tt_zz_xxxx_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 75);

        auto tt_zz_xxxy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 76);

        auto tt_zz_xxxz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 77);

        auto tt_zz_xxyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 78);

        auto tt_zz_xxyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 79);

        auto tt_zz_xxzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 80);

        auto tt_zz_xyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 81);

        auto tt_zz_xyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 82);

        auto tt_zz_xyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 83);

        auto tt_zz_xzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 84);

        auto tt_zz_yyyy_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 85);

        auto tt_zz_yyyz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 86);

        auto tt_zz_yyzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 87);

        auto tt_zz_yzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 88);

        auto tt_zz_zzzz_0 = primBuffer.data(pidx_t_2_4_m0 + 90 * idx + 89);

        auto tt_y_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 25);

        auto tt_y_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 26);

        auto tt_y_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 27);

        auto tt_y_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 28);

        auto tt_y_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 29);

        auto tt_z_xxxx_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 30);

        auto tt_z_xxxy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 31);

        auto tt_z_xxxz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 32);

        auto tt_z_xxyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 33);

        auto tt_z_xxyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 34);

        auto tt_z_xxzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 35);

        auto tt_z_xyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 36);

        auto tt_z_xyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 37);

        auto tt_z_xyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 38);

        auto tt_z_xzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 39);

        auto tt_z_yyyy_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 40);

        auto tt_z_yyyz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 41);

        auto tt_z_yyzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 42);

        auto tt_z_yzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 43);

        auto tt_z_zzzz_0 = primBuffer.data(pidx_t_1_4_m0 + 45 * idx + 44);

        auto tt_yy_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 36);

        auto tt_yy_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 37);

        auto tt_yy_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 38);

        auto tt_yy_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 39);

        auto tt_yz_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 40);

        auto tt_yz_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 41);

        auto tt_yz_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 42);

        auto tt_yz_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 43);

        auto tt_yz_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 44);

        auto tt_yz_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 45);

        auto tt_yz_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 46);

        auto tt_yz_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 47);

        auto tt_yz_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 48);

        auto tt_yz_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 49);

        auto tt_zz_xxx_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 50);

        auto tt_zz_xxy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 51);

        auto tt_zz_xxz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 52);

        auto tt_zz_xyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 53);

        auto tt_zz_xyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 54);

        auto tt_zz_xzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 55);

        auto tt_zz_yyy_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 56);

        auto tt_zz_yyz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 57);

        auto tt_zz_yzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 58);

        auto tt_zz_zzz_0 = primBuffer.data(pidx_t_2_3_m0 + 60 * idx + 59);

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

        auto ts_zzz_yyyy_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 145);

        auto ts_zzz_yyyz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 146);

        auto ts_zzz_yyzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 147);

        auto ts_zzz_yzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 148);

        auto ts_zzz_zzzz_0 = primBuffer.data(pidx_s_3_4_m0 + 150 * idx + 149);

        auto ts_y_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 25);

        auto ts_y_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 26);

        auto ts_y_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 27);

        auto ts_y_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 28);

        auto ts_y_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 29);

        auto ts_z_xxxx_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 30);

        auto ts_z_xxxy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 31);

        auto ts_z_xxxz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 32);

        auto ts_z_xxyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 33);

        auto ts_z_xxyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 34);

        auto ts_z_xxzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 35);

        auto ts_z_xyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 36);

        auto ts_z_xyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 37);

        auto ts_z_xyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 38);

        auto ts_z_xzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 39);

        auto ts_z_yyyy_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 40);

        auto ts_z_yyyz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 41);

        auto ts_z_yyzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 42);

        auto ts_z_yzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 43);

        auto ts_z_zzzz_0 = primBuffer.data(pidx_s_1_4_m0 + 45 * idx + 44);

        // set up pointers to integrals

        auto tt_yyy_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 100);

        auto tt_yyy_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 101);

        auto tt_yyy_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 102);

        auto tt_yyy_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 103);

        auto tt_yyy_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 104);

        auto tt_yyz_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 105);

        auto tt_yyz_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 106);

        auto tt_yyz_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 107);

        auto tt_yyz_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 108);

        auto tt_yyz_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 109);

        auto tt_yyz_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 110);

        auto tt_yyz_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 111);

        auto tt_yyz_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 112);

        auto tt_yyz_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 113);

        auto tt_yyz_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 114);

        auto tt_yyz_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 115);

        auto tt_yyz_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 116);

        auto tt_yyz_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 117);

        auto tt_yyz_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 118);

        auto tt_yyz_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 119);

        auto tt_yzz_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 120);

        auto tt_yzz_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 121);

        auto tt_yzz_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 122);

        auto tt_yzz_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 123);

        auto tt_yzz_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 124);

        auto tt_yzz_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 125);

        auto tt_yzz_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 126);

        auto tt_yzz_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 127);

        auto tt_yzz_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 128);

        auto tt_yzz_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 129);

        auto tt_yzz_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 130);

        auto tt_yzz_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 131);

        auto tt_yzz_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 132);

        auto tt_yzz_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 133);

        auto tt_yzz_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 134);

        auto tt_zzz_xxxx_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 135);

        auto tt_zzz_xxxy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 136);

        auto tt_zzz_xxxz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 137);

        auto tt_zzz_xxyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 138);

        auto tt_zzz_xxyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 139);

        auto tt_zzz_xxzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 140);

        auto tt_zzz_xyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 141);

        auto tt_zzz_xyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 142);

        auto tt_zzz_xyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 143);

        auto tt_zzz_xzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 144);

        auto tt_zzz_yyyy_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 145);

        auto tt_zzz_yyyz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 146);

        auto tt_zzz_yyzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 147);

        auto tt_zzz_yzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 148);

        auto tt_zzz_zzzz_0 = primBuffer.data(pidx_t_3_4_m0 + 150 * idx + 149);

        // Batch of Integrals (100,150)

        #pragma omp simd aligned(fga, fx, fz, pa_y, pa_z, ts_y_yyyy_0, ts_y_yyyz_0, ts_y_yyzz_0, ts_y_yzzz_0, \
                                     ts_y_zzzz_0, ts_yyy_yyyy_0, ts_yyy_yyyz_0, ts_yyy_yyzz_0, ts_yyy_yzzz_0, \
                                     ts_yyy_zzzz_0, ts_yyz_xxxx_0, ts_yyz_xxxy_0, ts_yyz_xxxz_0, ts_yyz_xxyy_0, \
                                     ts_yyz_xxyz_0, ts_yyz_xxzz_0, ts_yyz_xyyy_0, ts_yyz_xyyz_0, ts_yyz_xyzz_0, \
                                     ts_yyz_xzzz_0, ts_yyz_yyyy_0, ts_yyz_yyyz_0, ts_yyz_yyzz_0, ts_yyz_yzzz_0, \
                                     ts_yyz_zzzz_0, ts_yzz_xxxx_0, ts_yzz_xxxy_0, ts_yzz_xxxz_0, ts_yzz_xxyy_0, \
                                     ts_yzz_xxyz_0, ts_yzz_xxzz_0, ts_yzz_xyyy_0, ts_yzz_xyyz_0, ts_yzz_xyzz_0, \
                                     ts_yzz_xzzz_0, ts_yzz_yyyy_0, ts_yzz_yyyz_0, ts_yzz_yyzz_0, ts_yzz_yzzz_0, \
                                     ts_yzz_zzzz_0, ts_z_xxxx_0, ts_z_xxxy_0, ts_z_xxxz_0, ts_z_xxyy_0, ts_z_xxyz_0, \
                                     ts_z_xxzz_0, ts_z_xyyy_0, ts_z_xyyz_0, ts_z_xyzz_0, ts_z_xzzz_0, ts_z_yyyy_0, \
                                     ts_z_yyyz_0, ts_z_yyzz_0, ts_z_yzzz_0, ts_z_zzzz_0, ts_zzz_xxxx_0, ts_zzz_xxxy_0, \
                                     ts_zzz_xxxz_0, ts_zzz_xxyy_0, ts_zzz_xxyz_0, ts_zzz_xxzz_0, ts_zzz_xyyy_0, \
                                     ts_zzz_xyyz_0, ts_zzz_xyzz_0, ts_zzz_xzzz_0, ts_zzz_yyyy_0, ts_zzz_yyyz_0, \
                                     ts_zzz_yyzz_0, ts_zzz_yzzz_0, ts_zzz_zzzz_0, tt_y_yyyy_0, tt_y_yyyz_0, tt_y_yyzz_0, \
                                     tt_y_yzzz_0, tt_y_zzzz_0, tt_yy_yyy_0, tt_yy_yyyy_0, tt_yy_yyyz_0, tt_yy_yyz_0, \
                                     tt_yy_yyzz_0, tt_yy_yzz_0, tt_yy_yzzz_0, tt_yy_zzz_0, tt_yy_zzzz_0, tt_yyy_yyyy_0, \
                                     tt_yyy_yyyz_0, tt_yyy_yyzz_0, tt_yyy_yzzz_0, tt_yyy_zzzz_0, tt_yyz_xxxx_0, \
                                     tt_yyz_xxxy_0, tt_yyz_xxxz_0, tt_yyz_xxyy_0, tt_yyz_xxyz_0, tt_yyz_xxzz_0, \
                                     tt_yyz_xyyy_0, tt_yyz_xyyz_0, tt_yyz_xyzz_0, tt_yyz_xzzz_0, tt_yyz_yyyy_0, \
                                     tt_yyz_yyyz_0, tt_yyz_yyzz_0, tt_yyz_yzzz_0, tt_yyz_zzzz_0, tt_yz_xxx_0, \
                                     tt_yz_xxxx_0, tt_yz_xxxy_0, tt_yz_xxxz_0, tt_yz_xxy_0, tt_yz_xxyy_0, tt_yz_xxyz_0, \
                                     tt_yz_xxz_0, tt_yz_xxzz_0, tt_yz_xyy_0, tt_yz_xyyy_0, tt_yz_xyyz_0, tt_yz_xyz_0, \
                                     tt_yz_xyzz_0, tt_yz_xzz_0, tt_yz_xzzz_0, tt_yz_yyy_0, tt_yz_yyyy_0, tt_yz_yyyz_0, \
                                     tt_yz_yyz_0, tt_yz_yyzz_0, tt_yz_yzz_0, tt_yz_yzzz_0, tt_yz_zzz_0, tt_yz_zzzz_0, \
                                     tt_yzz_xxxx_0, tt_yzz_xxxy_0, tt_yzz_xxxz_0, tt_yzz_xxyy_0, tt_yzz_xxyz_0, \
                                     tt_yzz_xxzz_0, tt_yzz_xyyy_0, tt_yzz_xyyz_0, tt_yzz_xyzz_0, tt_yzz_xzzz_0, \
                                     tt_yzz_yyyy_0, tt_yzz_yyyz_0, tt_yzz_yyzz_0, tt_yzz_yzzz_0, tt_yzz_zzzz_0, \
                                     tt_z_xxxx_0, tt_z_xxxy_0, tt_z_xxxz_0, tt_z_xxyy_0, tt_z_xxyz_0, tt_z_xxzz_0, \
                                     tt_z_xyyy_0, tt_z_xyyz_0, tt_z_xyzz_0, tt_z_xzzz_0, tt_z_yyyy_0, tt_z_yyyz_0, \
                                     tt_z_yyzz_0, tt_z_yzzz_0, tt_z_zzzz_0, tt_zz_xxx_0, tt_zz_xxxx_0, tt_zz_xxxy_0, \
                                     tt_zz_xxxz_0, tt_zz_xxy_0, tt_zz_xxyy_0, tt_zz_xxyz_0, tt_zz_xxz_0, tt_zz_xxzz_0, \
                                     tt_zz_xyy_0, tt_zz_xyyy_0, tt_zz_xyyz_0, tt_zz_xyz_0, tt_zz_xyzz_0, tt_zz_xzz_0, \
                                     tt_zz_xzzz_0, tt_zz_yyy_0, tt_zz_yyyy_0, tt_zz_yyyz_0, tt_zz_yyz_0, tt_zz_yyzz_0, \
                                     tt_zz_yzz_0, tt_zz_yzzz_0, tt_zz_zzz_0, tt_zz_zzzz_0, tt_zzz_xxxx_0, \
                                     tt_zzz_xxxy_0, tt_zzz_xxxz_0, tt_zzz_xxyy_0, tt_zzz_xxyz_0, tt_zzz_xxzz_0, \
                                     tt_zzz_xyyy_0, tt_zzz_xyyz_0, tt_zzz_xyzz_0, tt_zzz_xzzz_0, tt_zzz_yyyy_0, \
                                     tt_zzz_yyyz_0, tt_zzz_yyzz_0, tt_zzz_yzzz_0, tt_zzz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_yyy_yyyy_0[j] = pa_y[j] * tt_yy_yyyy_0[j] + fl1_fx * tt_y_yyyy_0[j] + 2.0 * fl1_fx * tt_yy_yyy_0[j] + 2.0 * fl1_fz * ts_yyy_yyyy_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_yyyy_0[j];

            tt_yyy_yyyz_0[j] = pa_y[j] * tt_yy_yyyz_0[j] + fl1_fx * tt_y_yyyz_0[j] + 1.5 * fl1_fx * tt_yy_yyz_0[j] + 2.0 * fl1_fz * ts_yyy_yyyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_yyyz_0[j];

            tt_yyy_yyzz_0[j] = pa_y[j] * tt_yy_yyzz_0[j] + fl1_fx * tt_y_yyzz_0[j] + fl1_fx * tt_yy_yzz_0[j] + 2.0 * fl1_fz * ts_yyy_yyzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_yyzz_0[j];

            tt_yyy_yzzz_0[j] = pa_y[j] * tt_yy_yzzz_0[j] + fl1_fx * tt_y_yzzz_0[j] + 0.5 * fl1_fx * tt_yy_zzz_0[j] + 2.0 * fl1_fz * ts_yyy_yzzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_y_yzzz_0[j];

            tt_yyy_zzzz_0[j] =
                pa_y[j] * tt_yy_zzzz_0[j] + fl1_fx * tt_y_zzzz_0[j] + 2.0 * fl1_fz * ts_yyy_zzzz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_zzzz_0[j];

            tt_yyz_xxxx_0[j] =
                pa_y[j] * tt_yz_xxxx_0[j] + 0.5 * fl1_fx * tt_z_xxxx_0[j] + 2.0 * fl1_fz * ts_yyz_xxxx_0[j] - fl1_fz * fl1_fga * ts_z_xxxx_0[j];

            tt_yyz_xxxy_0[j] = pa_y[j] * tt_yz_xxxy_0[j] + 0.5 * fl1_fx * tt_z_xxxy_0[j] + 0.5 * fl1_fx * tt_yz_xxx_0[j] +
                               2.0 * fl1_fz * ts_yyz_xxxy_0[j] - fl1_fz * fl1_fga * ts_z_xxxy_0[j];

            tt_yyz_xxxz_0[j] =
                pa_y[j] * tt_yz_xxxz_0[j] + 0.5 * fl1_fx * tt_z_xxxz_0[j] + 2.0 * fl1_fz * ts_yyz_xxxz_0[j] - fl1_fz * fl1_fga * ts_z_xxxz_0[j];

            tt_yyz_xxyy_0[j] = pa_y[j] * tt_yz_xxyy_0[j] + 0.5 * fl1_fx * tt_z_xxyy_0[j] + fl1_fx * tt_yz_xxy_0[j] + 2.0 * fl1_fz * ts_yyz_xxyy_0[j] -
                               fl1_fz * fl1_fga * ts_z_xxyy_0[j];

            tt_yyz_xxyz_0[j] = pa_y[j] * tt_yz_xxyz_0[j] + 0.5 * fl1_fx * tt_z_xxyz_0[j] + 0.5 * fl1_fx * tt_yz_xxz_0[j] +
                               2.0 * fl1_fz * ts_yyz_xxyz_0[j] - fl1_fz * fl1_fga * ts_z_xxyz_0[j];

            tt_yyz_xxzz_0[j] =
                pa_y[j] * tt_yz_xxzz_0[j] + 0.5 * fl1_fx * tt_z_xxzz_0[j] + 2.0 * fl1_fz * ts_yyz_xxzz_0[j] - fl1_fz * fl1_fga * ts_z_xxzz_0[j];

            tt_yyz_xyyy_0[j] = pa_y[j] * tt_yz_xyyy_0[j] + 0.5 * fl1_fx * tt_z_xyyy_0[j] + 1.5 * fl1_fx * tt_yz_xyy_0[j] +
                               2.0 * fl1_fz * ts_yyz_xyyy_0[j] - fl1_fz * fl1_fga * ts_z_xyyy_0[j];

            tt_yyz_xyyz_0[j] = pa_y[j] * tt_yz_xyyz_0[j] + 0.5 * fl1_fx * tt_z_xyyz_0[j] + fl1_fx * tt_yz_xyz_0[j] + 2.0 * fl1_fz * ts_yyz_xyyz_0[j] -
                               fl1_fz * fl1_fga * ts_z_xyyz_0[j];

            tt_yyz_xyzz_0[j] = pa_y[j] * tt_yz_xyzz_0[j] + 0.5 * fl1_fx * tt_z_xyzz_0[j] + 0.5 * fl1_fx * tt_yz_xzz_0[j] +
                               2.0 * fl1_fz * ts_yyz_xyzz_0[j] - fl1_fz * fl1_fga * ts_z_xyzz_0[j];

            tt_yyz_xzzz_0[j] =
                pa_y[j] * tt_yz_xzzz_0[j] + 0.5 * fl1_fx * tt_z_xzzz_0[j] + 2.0 * fl1_fz * ts_yyz_xzzz_0[j] - fl1_fz * fl1_fga * ts_z_xzzz_0[j];

            tt_yyz_yyyy_0[j] = pa_y[j] * tt_yz_yyyy_0[j] + 0.5 * fl1_fx * tt_z_yyyy_0[j] + 2.0 * fl1_fx * tt_yz_yyy_0[j] +
                               2.0 * fl1_fz * ts_yyz_yyyy_0[j] - fl1_fz * fl1_fga * ts_z_yyyy_0[j];

            tt_yyz_yyyz_0[j] = pa_y[j] * tt_yz_yyyz_0[j] + 0.5 * fl1_fx * tt_z_yyyz_0[j] + 1.5 * fl1_fx * tt_yz_yyz_0[j] +
                               2.0 * fl1_fz * ts_yyz_yyyz_0[j] - fl1_fz * fl1_fga * ts_z_yyyz_0[j];

            tt_yyz_yyzz_0[j] = pa_y[j] * tt_yz_yyzz_0[j] + 0.5 * fl1_fx * tt_z_yyzz_0[j] + fl1_fx * tt_yz_yzz_0[j] + 2.0 * fl1_fz * ts_yyz_yyzz_0[j] -
                               fl1_fz * fl1_fga * ts_z_yyzz_0[j];

            tt_yyz_yzzz_0[j] = pa_y[j] * tt_yz_yzzz_0[j] + 0.5 * fl1_fx * tt_z_yzzz_0[j] + 0.5 * fl1_fx * tt_yz_zzz_0[j] +
                               2.0 * fl1_fz * ts_yyz_yzzz_0[j] - fl1_fz * fl1_fga * ts_z_yzzz_0[j];

            tt_yyz_zzzz_0[j] =
                pa_y[j] * tt_yz_zzzz_0[j] + 0.5 * fl1_fx * tt_z_zzzz_0[j] + 2.0 * fl1_fz * ts_yyz_zzzz_0[j] - fl1_fz * fl1_fga * ts_z_zzzz_0[j];

            tt_yzz_xxxx_0[j] = pa_y[j] * tt_zz_xxxx_0[j] + 2.0 * fl1_fz * ts_yzz_xxxx_0[j];

            tt_yzz_xxxy_0[j] = pa_y[j] * tt_zz_xxxy_0[j] + 0.5 * fl1_fx * tt_zz_xxx_0[j] + 2.0 * fl1_fz * ts_yzz_xxxy_0[j];

            tt_yzz_xxxz_0[j] = pa_y[j] * tt_zz_xxxz_0[j] + 2.0 * fl1_fz * ts_yzz_xxxz_0[j];

            tt_yzz_xxyy_0[j] = pa_y[j] * tt_zz_xxyy_0[j] + fl1_fx * tt_zz_xxy_0[j] + 2.0 * fl1_fz * ts_yzz_xxyy_0[j];

            tt_yzz_xxyz_0[j] = pa_y[j] * tt_zz_xxyz_0[j] + 0.5 * fl1_fx * tt_zz_xxz_0[j] + 2.0 * fl1_fz * ts_yzz_xxyz_0[j];

            tt_yzz_xxzz_0[j] = pa_y[j] * tt_zz_xxzz_0[j] + 2.0 * fl1_fz * ts_yzz_xxzz_0[j];

            tt_yzz_xyyy_0[j] = pa_y[j] * tt_zz_xyyy_0[j] + 1.5 * fl1_fx * tt_zz_xyy_0[j] + 2.0 * fl1_fz * ts_yzz_xyyy_0[j];

            tt_yzz_xyyz_0[j] = pa_y[j] * tt_zz_xyyz_0[j] + fl1_fx * tt_zz_xyz_0[j] + 2.0 * fl1_fz * ts_yzz_xyyz_0[j];

            tt_yzz_xyzz_0[j] = pa_y[j] * tt_zz_xyzz_0[j] + 0.5 * fl1_fx * tt_zz_xzz_0[j] + 2.0 * fl1_fz * ts_yzz_xyzz_0[j];

            tt_yzz_xzzz_0[j] = pa_y[j] * tt_zz_xzzz_0[j] + 2.0 * fl1_fz * ts_yzz_xzzz_0[j];

            tt_yzz_yyyy_0[j] = pa_y[j] * tt_zz_yyyy_0[j] + 2.0 * fl1_fx * tt_zz_yyy_0[j] + 2.0 * fl1_fz * ts_yzz_yyyy_0[j];

            tt_yzz_yyyz_0[j] = pa_y[j] * tt_zz_yyyz_0[j] + 1.5 * fl1_fx * tt_zz_yyz_0[j] + 2.0 * fl1_fz * ts_yzz_yyyz_0[j];

            tt_yzz_yyzz_0[j] = pa_y[j] * tt_zz_yyzz_0[j] + fl1_fx * tt_zz_yzz_0[j] + 2.0 * fl1_fz * ts_yzz_yyzz_0[j];

            tt_yzz_yzzz_0[j] = pa_y[j] * tt_zz_yzzz_0[j] + 0.5 * fl1_fx * tt_zz_zzz_0[j] + 2.0 * fl1_fz * ts_yzz_yzzz_0[j];

            tt_yzz_zzzz_0[j] = pa_y[j] * tt_zz_zzzz_0[j] + 2.0 * fl1_fz * ts_yzz_zzzz_0[j];

            tt_zzz_xxxx_0[j] =
                pa_z[j] * tt_zz_xxxx_0[j] + fl1_fx * tt_z_xxxx_0[j] + 2.0 * fl1_fz * ts_zzz_xxxx_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_xxxx_0[j];

            tt_zzz_xxxy_0[j] =
                pa_z[j] * tt_zz_xxxy_0[j] + fl1_fx * tt_z_xxxy_0[j] + 2.0 * fl1_fz * ts_zzz_xxxy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_xxxy_0[j];

            tt_zzz_xxxz_0[j] = pa_z[j] * tt_zz_xxxz_0[j] + fl1_fx * tt_z_xxxz_0[j] + 0.5 * fl1_fx * tt_zz_xxx_0[j] + 2.0 * fl1_fz * ts_zzz_xxxz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_xxxz_0[j];

            tt_zzz_xxyy_0[j] =
                pa_z[j] * tt_zz_xxyy_0[j] + fl1_fx * tt_z_xxyy_0[j] + 2.0 * fl1_fz * ts_zzz_xxyy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_xxyy_0[j];

            tt_zzz_xxyz_0[j] = pa_z[j] * tt_zz_xxyz_0[j] + fl1_fx * tt_z_xxyz_0[j] + 0.5 * fl1_fx * tt_zz_xxy_0[j] + 2.0 * fl1_fz * ts_zzz_xxyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_xxyz_0[j];

            tt_zzz_xxzz_0[j] = pa_z[j] * tt_zz_xxzz_0[j] + fl1_fx * tt_z_xxzz_0[j] + fl1_fx * tt_zz_xxz_0[j] + 2.0 * fl1_fz * ts_zzz_xxzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_xxzz_0[j];

            tt_zzz_xyyy_0[j] =
                pa_z[j] * tt_zz_xyyy_0[j] + fl1_fx * tt_z_xyyy_0[j] + 2.0 * fl1_fz * ts_zzz_xyyy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_xyyy_0[j];

            tt_zzz_xyyz_0[j] = pa_z[j] * tt_zz_xyyz_0[j] + fl1_fx * tt_z_xyyz_0[j] + 0.5 * fl1_fx * tt_zz_xyy_0[j] + 2.0 * fl1_fz * ts_zzz_xyyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_xyyz_0[j];

            tt_zzz_xyzz_0[j] = pa_z[j] * tt_zz_xyzz_0[j] + fl1_fx * tt_z_xyzz_0[j] + fl1_fx * tt_zz_xyz_0[j] + 2.0 * fl1_fz * ts_zzz_xyzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_xyzz_0[j];

            tt_zzz_xzzz_0[j] = pa_z[j] * tt_zz_xzzz_0[j] + fl1_fx * tt_z_xzzz_0[j] + 1.5 * fl1_fx * tt_zz_xzz_0[j] + 2.0 * fl1_fz * ts_zzz_xzzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_xzzz_0[j];

            tt_zzz_yyyy_0[j] =
                pa_z[j] * tt_zz_yyyy_0[j] + fl1_fx * tt_z_yyyy_0[j] + 2.0 * fl1_fz * ts_zzz_yyyy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_yyyy_0[j];

            tt_zzz_yyyz_0[j] = pa_z[j] * tt_zz_yyyz_0[j] + fl1_fx * tt_z_yyyz_0[j] + 0.5 * fl1_fx * tt_zz_yyy_0[j] + 2.0 * fl1_fz * ts_zzz_yyyz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_yyyz_0[j];

            tt_zzz_yyzz_0[j] = pa_z[j] * tt_zz_yyzz_0[j] + fl1_fx * tt_z_yyzz_0[j] + fl1_fx * tt_zz_yyz_0[j] + 2.0 * fl1_fz * ts_zzz_yyzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_yyzz_0[j];

            tt_zzz_yzzz_0[j] = pa_z[j] * tt_zz_yzzz_0[j] + fl1_fx * tt_z_yzzz_0[j] + 1.5 * fl1_fx * tt_zz_yzz_0[j] + 2.0 * fl1_fz * ts_zzz_yzzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_yzzz_0[j];

            tt_zzz_zzzz_0[j] = pa_z[j] * tt_zz_zzzz_0[j] + fl1_fx * tt_z_zzzz_0[j] + 2.0 * fl1_fx * tt_zz_zzz_0[j] + 2.0 * fl1_fz * ts_zzz_zzzz_0[j] -
                               2.0 * fl1_fz * fl1_fga * ts_z_zzzz_0[j];
        }

        idx++;
    }
}

}  // namespace kinrecfunc
