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

#include "NuclearPotentialRecFuncForFG.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

void
compNuclearPotentialForFG(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForFG_0_50(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFG_50_100(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFG_100_150(
        primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForFG_0_50(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto ta_xx_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx);

            auto ta_xx_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 1);

            auto ta_xx_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 2);

            auto ta_xx_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 3);

            auto ta_xx_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 4);

            auto ta_xx_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 5);

            auto ta_xx_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 6);

            auto ta_xx_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 7);

            auto ta_xx_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 8);

            auto ta_xx_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 9);

            auto ta_xx_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 10);

            auto ta_xx_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 11);

            auto ta_xx_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 12);

            auto ta_xx_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 13);

            auto ta_xx_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 14);

            auto ta_xy_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 15);

            auto ta_xy_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 16);

            auto ta_xy_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 17);

            auto ta_xy_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 18);

            auto ta_xy_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 19);

            auto ta_xy_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 20);

            auto ta_xy_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 21);

            auto ta_xy_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 22);

            auto ta_xy_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 23);

            auto ta_xy_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 24);

            auto ta_xy_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 25);

            auto ta_xy_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 26);

            auto ta_xy_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 27);

            auto ta_xy_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 28);

            auto ta_xy_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 29);

            auto ta_xz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 30);

            auto ta_xz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 31);

            auto ta_xz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 32);

            auto ta_xz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 33);

            auto ta_xz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 34);

            auto ta_xz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 35);

            auto ta_xz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 36);

            auto ta_xz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 37);

            auto ta_xz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 38);

            auto ta_xz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 39);

            auto ta_xz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 40);

            auto ta_xz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 41);

            auto ta_xz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 42);

            auto ta_xz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 43);

            auto ta_xz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 44);

            auto ta_yy_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 45);

            auto ta_yy_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 46);

            auto ta_yy_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 47);

            auto ta_yy_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 48);

            auto ta_yy_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 49);

            auto ta_xx_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx);

            auto ta_xx_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 1);

            auto ta_xx_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 2);

            auto ta_xx_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 3);

            auto ta_xx_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 4);

            auto ta_xx_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 5);

            auto ta_xx_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 6);

            auto ta_xx_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 7);

            auto ta_xx_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 8);

            auto ta_xx_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 9);

            auto ta_xx_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 10);

            auto ta_xx_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 11);

            auto ta_xx_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 12);

            auto ta_xx_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 13);

            auto ta_xx_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 14);

            auto ta_xy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 15);

            auto ta_xy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 16);

            auto ta_xy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 17);

            auto ta_xy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 18);

            auto ta_xy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 19);

            auto ta_xy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 20);

            auto ta_xy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 21);

            auto ta_xy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 22);

            auto ta_xy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 23);

            auto ta_xy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 24);

            auto ta_xy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 25);

            auto ta_xy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 26);

            auto ta_xy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 27);

            auto ta_xy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 28);

            auto ta_xy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 29);

            auto ta_xz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 30);

            auto ta_xz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 31);

            auto ta_xz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 32);

            auto ta_xz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 33);

            auto ta_xz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 34);

            auto ta_xz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 35);

            auto ta_xz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 36);

            auto ta_xz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 37);

            auto ta_xz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 38);

            auto ta_xz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 39);

            auto ta_xz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 40);

            auto ta_xz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 41);

            auto ta_xz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 42);

            auto ta_xz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 43);

            auto ta_xz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 44);

            auto ta_yy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 45);

            auto ta_yy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 46);

            auto ta_yy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 47);

            auto ta_yy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 48);

            auto ta_yy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 49);

            auto ta_x_xxxx_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx);

            auto ta_x_xxxy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 1);

            auto ta_x_xxxz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 2);

            auto ta_x_xxyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 3);

            auto ta_x_xxyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 4);

            auto ta_x_xxzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 5);

            auto ta_x_xyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 6);

            auto ta_x_xyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 7);

            auto ta_x_xyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 8);

            auto ta_x_xzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 9);

            auto ta_x_yyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 10);

            auto ta_x_yyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 11);

            auto ta_x_yyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 12);

            auto ta_x_yzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 13);

            auto ta_x_zzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 14);

            auto ta_y_xxxx_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 15);

            auto ta_y_xxxy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 16);

            auto ta_y_xxxz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 17);

            auto ta_y_xxyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 18);

            auto ta_y_xxyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 19);

            auto ta_y_xxzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 20);

            auto ta_y_xyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 21);

            auto ta_y_xyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 22);

            auto ta_y_xyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 23);

            auto ta_y_xzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 24);

            auto ta_y_yyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 25);

            auto ta_y_yyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 26);

            auto ta_y_yyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 27);

            auto ta_y_yzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 28);

            auto ta_y_zzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 29);

            auto ta_z_xxxx_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 30);

            auto ta_z_xxxy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 31);

            auto ta_z_xxxz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 32);

            auto ta_z_xxyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 33);

            auto ta_z_xxyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 34);

            auto ta_z_xxzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 35);

            auto ta_z_xyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 36);

            auto ta_z_xyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 37);

            auto ta_z_xyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 38);

            auto ta_z_xzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 39);

            auto ta_z_yyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 40);

            auto ta_z_yyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 41);

            auto ta_z_yyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 42);

            auto ta_z_yzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 43);

            auto ta_z_zzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 44);

            auto ta_x_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx);

            auto ta_x_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 1);

            auto ta_x_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 2);

            auto ta_x_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 3);

            auto ta_x_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 4);

            auto ta_x_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 5);

            auto ta_x_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 6);

            auto ta_x_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 7);

            auto ta_x_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 8);

            auto ta_x_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 9);

            auto ta_x_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 10);

            auto ta_x_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 11);

            auto ta_x_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 12);

            auto ta_x_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 13);

            auto ta_x_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 14);

            auto ta_y_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 15);

            auto ta_y_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 16);

            auto ta_y_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 17);

            auto ta_y_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 18);

            auto ta_y_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 19);

            auto ta_y_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 20);

            auto ta_y_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 21);

            auto ta_y_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 22);

            auto ta_y_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 23);

            auto ta_y_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 24);

            auto ta_y_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 25);

            auto ta_y_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 26);

            auto ta_y_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 27);

            auto ta_y_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 28);

            auto ta_y_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 29);

            auto ta_z_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 30);

            auto ta_z_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 31);

            auto ta_z_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 32);

            auto ta_z_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 33);

            auto ta_z_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 34);

            auto ta_z_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 35);

            auto ta_z_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 36);

            auto ta_z_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 37);

            auto ta_z_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 38);

            auto ta_z_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 39);

            auto ta_z_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 40);

            auto ta_z_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 41);

            auto ta_z_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 42);

            auto ta_z_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 43);

            auto ta_z_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 44);

            auto ta_xx_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx);

            auto ta_xx_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 1);

            auto ta_xx_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 2);

            auto ta_xx_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 3);

            auto ta_xx_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 4);

            auto ta_xx_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 5);

            auto ta_xx_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 6);

            auto ta_xx_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 7);

            auto ta_xx_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 8);

            auto ta_xx_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 9);

            auto ta_xy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 10);

            auto ta_xy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 11);

            auto ta_xy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 12);

            auto ta_xy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 13);

            auto ta_xy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 14);

            auto ta_xy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 15);

            auto ta_xy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 16);

            auto ta_xy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 17);

            auto ta_xy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 18);

            auto ta_xy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 19);

            auto ta_xz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 20);

            auto ta_xz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 21);

            auto ta_xz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 22);

            auto ta_xz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 23);

            auto ta_xz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 24);

            auto ta_xz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 25);

            auto ta_xz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 26);

            auto ta_xz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 27);

            auto ta_xz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 28);

            auto ta_xz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 29);

            auto ta_yy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 30);

            auto ta_yy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 31);

            auto ta_yy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 32);

            auto ta_yy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 33);

            auto ta_yy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 34);

            auto ta_xx_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx);

            auto ta_xx_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 1);

            auto ta_xx_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 2);

            auto ta_xx_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 3);

            auto ta_xx_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 4);

            auto ta_xx_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 5);

            auto ta_xx_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 6);

            auto ta_xx_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 7);

            auto ta_xx_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 8);

            auto ta_xx_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 9);

            auto ta_xy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 10);

            auto ta_xy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 11);

            auto ta_xy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 12);

            auto ta_xy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 13);

            auto ta_xy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 14);

            auto ta_xy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 15);

            auto ta_xy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 16);

            auto ta_xy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 17);

            auto ta_xy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 18);

            auto ta_xy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 19);

            auto ta_xz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 20);

            auto ta_xz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 21);

            auto ta_xz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 22);

            auto ta_xz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 23);

            auto ta_xz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 24);

            auto ta_xz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 25);

            auto ta_xz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 26);

            auto ta_xz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 27);

            auto ta_xz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 28);

            auto ta_xz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 29);

            auto ta_yy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 30);

            auto ta_yy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 31);

            auto ta_yy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 32);

            auto ta_yy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 33);

            auto ta_yy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 34);

            // set up pointers to integrals

            auto ta_xxx_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx);

            auto ta_xxx_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 1);

            auto ta_xxx_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 2);

            auto ta_xxx_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 3);

            auto ta_xxx_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 4);

            auto ta_xxx_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 5);

            auto ta_xxx_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 6);

            auto ta_xxx_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 7);

            auto ta_xxx_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 8);

            auto ta_xxx_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 9);

            auto ta_xxx_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 10);

            auto ta_xxx_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 11);

            auto ta_xxx_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 12);

            auto ta_xxx_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 13);

            auto ta_xxx_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 14);

            auto ta_xxy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 15);

            auto ta_xxy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 16);

            auto ta_xxy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 17);

            auto ta_xxy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 18);

            auto ta_xxy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 19);

            auto ta_xxy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 20);

            auto ta_xxy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 21);

            auto ta_xxy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 22);

            auto ta_xxy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 23);

            auto ta_xxy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 24);

            auto ta_xxy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 25);

            auto ta_xxy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 26);

            auto ta_xxy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 27);

            auto ta_xxy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 28);

            auto ta_xxy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 29);

            auto ta_xxz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 30);

            auto ta_xxz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 31);

            auto ta_xxz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 32);

            auto ta_xxz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 33);

            auto ta_xxz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 34);

            auto ta_xxz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 35);

            auto ta_xxz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 36);

            auto ta_xxz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 37);

            auto ta_xxz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 38);

            auto ta_xxz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 39);

            auto ta_xxz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 40);

            auto ta_xxz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 41);

            auto ta_xxz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 42);

            auto ta_xxz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 43);

            auto ta_xxz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 44);

            auto ta_xyy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 45);

            auto ta_xyy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 46);

            auto ta_xyy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 47);

            auto ta_xyy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 48);

            auto ta_xyy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 49);

            // Batch of Integrals (0,50)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_xxxx_0, ta_x_xxxx_1, ta_x_xxxy_0, ta_x_xxxy_1, \
                                         ta_x_xxxz_0, ta_x_xxxz_1, ta_x_xxyy_0, ta_x_xxyy_1, ta_x_xxyz_0, ta_x_xxyz_1, \
                                         ta_x_xxzz_0, ta_x_xxzz_1, ta_x_xyyy_0, ta_x_xyyy_1, ta_x_xyyz_0, ta_x_xyyz_1, \
                                         ta_x_xyzz_0, ta_x_xyzz_1, ta_x_xzzz_0, ta_x_xzzz_1, ta_x_yyyy_0, ta_x_yyyy_1, \
                                         ta_x_yyyz_0, ta_x_yyyz_1, ta_x_yyzz_0, ta_x_yyzz_1, ta_x_yzzz_0, ta_x_yzzz_1, \
                                         ta_x_zzzz_0, ta_x_zzzz_1, ta_xx_xxx_0, ta_xx_xxx_1, ta_xx_xxxx_0, ta_xx_xxxx_1, \
                                         ta_xx_xxxy_0, ta_xx_xxxy_1, ta_xx_xxxz_0, ta_xx_xxxz_1, ta_xx_xxy_0, ta_xx_xxy_1, \
                                         ta_xx_xxyy_0, ta_xx_xxyy_1, ta_xx_xxyz_0, ta_xx_xxyz_1, ta_xx_xxz_0, ta_xx_xxz_1, \
                                         ta_xx_xxzz_0, ta_xx_xxzz_1, ta_xx_xyy_0, ta_xx_xyy_1, ta_xx_xyyy_0, ta_xx_xyyy_1, \
                                         ta_xx_xyyz_0, ta_xx_xyyz_1, ta_xx_xyz_0, ta_xx_xyz_1, ta_xx_xyzz_0, ta_xx_xyzz_1, \
                                         ta_xx_xzz_0, ta_xx_xzz_1, ta_xx_xzzz_0, ta_xx_xzzz_1, ta_xx_yyy_0, ta_xx_yyy_1, \
                                         ta_xx_yyyy_0, ta_xx_yyyy_1, ta_xx_yyyz_0, ta_xx_yyyz_1, ta_xx_yyz_0, ta_xx_yyz_1, \
                                         ta_xx_yyzz_0, ta_xx_yyzz_1, ta_xx_yzz_0, ta_xx_yzz_1, ta_xx_yzzz_0, ta_xx_yzzz_1, \
                                         ta_xx_zzz_0, ta_xx_zzz_1, ta_xx_zzzz_0, ta_xx_zzzz_1, ta_xxx_xxxx_0, \
                                         ta_xxx_xxxy_0, ta_xxx_xxxz_0, ta_xxx_xxyy_0, ta_xxx_xxyz_0, ta_xxx_xxzz_0, \
                                         ta_xxx_xyyy_0, ta_xxx_xyyz_0, ta_xxx_xyzz_0, ta_xxx_xzzz_0, ta_xxx_yyyy_0, \
                                         ta_xxx_yyyz_0, ta_xxx_yyzz_0, ta_xxx_yzzz_0, ta_xxx_zzzz_0, ta_xxy_xxxx_0, \
                                         ta_xxy_xxxy_0, ta_xxy_xxxz_0, ta_xxy_xxyy_0, ta_xxy_xxyz_0, ta_xxy_xxzz_0, \
                                         ta_xxy_xyyy_0, ta_xxy_xyyz_0, ta_xxy_xyzz_0, ta_xxy_xzzz_0, ta_xxy_yyyy_0, \
                                         ta_xxy_yyyz_0, ta_xxy_yyzz_0, ta_xxy_yzzz_0, ta_xxy_zzzz_0, ta_xxz_xxxx_0, \
                                         ta_xxz_xxxy_0, ta_xxz_xxxz_0, ta_xxz_xxyy_0, ta_xxz_xxyz_0, ta_xxz_xxzz_0, \
                                         ta_xxz_xyyy_0, ta_xxz_xyyz_0, ta_xxz_xyzz_0, ta_xxz_xzzz_0, ta_xxz_yyyy_0, \
                                         ta_xxz_yyyz_0, ta_xxz_yyzz_0, ta_xxz_yzzz_0, ta_xxz_zzzz_0, ta_xy_xxx_0, \
                                         ta_xy_xxx_1, ta_xy_xxxx_0, ta_xy_xxxx_1, ta_xy_xxxy_0, ta_xy_xxxy_1, ta_xy_xxxz_0, \
                                         ta_xy_xxxz_1, ta_xy_xxy_0, ta_xy_xxy_1, ta_xy_xxyy_0, ta_xy_xxyy_1, ta_xy_xxyz_0, \
                                         ta_xy_xxyz_1, ta_xy_xxz_0, ta_xy_xxz_1, ta_xy_xxzz_0, ta_xy_xxzz_1, ta_xy_xyy_0, \
                                         ta_xy_xyy_1, ta_xy_xyyy_0, ta_xy_xyyy_1, ta_xy_xyyz_0, ta_xy_xyyz_1, ta_xy_xyz_0, \
                                         ta_xy_xyz_1, ta_xy_xyzz_0, ta_xy_xyzz_1, ta_xy_xzz_0, ta_xy_xzz_1, ta_xy_xzzz_0, \
                                         ta_xy_xzzz_1, ta_xy_yyy_0, ta_xy_yyy_1, ta_xy_yyyy_0, ta_xy_yyyy_1, ta_xy_yyyz_0, \
                                         ta_xy_yyyz_1, ta_xy_yyz_0, ta_xy_yyz_1, ta_xy_yyzz_0, ta_xy_yyzz_1, ta_xy_yzz_0, \
                                         ta_xy_yzz_1, ta_xy_yzzz_0, ta_xy_yzzz_1, ta_xy_zzz_0, ta_xy_zzz_1, ta_xy_zzzz_0, \
                                         ta_xy_zzzz_1, ta_xyy_xxxx_0, ta_xyy_xxxy_0, ta_xyy_xxxz_0, ta_xyy_xxyy_0, \
                                         ta_xyy_xxyz_0, ta_xz_xxx_0, ta_xz_xxx_1, ta_xz_xxxx_0, ta_xz_xxxx_1, ta_xz_xxxy_0, \
                                         ta_xz_xxxy_1, ta_xz_xxxz_0, ta_xz_xxxz_1, ta_xz_xxy_0, ta_xz_xxy_1, ta_xz_xxyy_0, \
                                         ta_xz_xxyy_1, ta_xz_xxyz_0, ta_xz_xxyz_1, ta_xz_xxz_0, ta_xz_xxz_1, ta_xz_xxzz_0, \
                                         ta_xz_xxzz_1, ta_xz_xyy_0, ta_xz_xyy_1, ta_xz_xyyy_0, ta_xz_xyyy_1, ta_xz_xyyz_0, \
                                         ta_xz_xyyz_1, ta_xz_xyz_0, ta_xz_xyz_1, ta_xz_xyzz_0, ta_xz_xyzz_1, ta_xz_xzz_0, \
                                         ta_xz_xzz_1, ta_xz_xzzz_0, ta_xz_xzzz_1, ta_xz_yyy_0, ta_xz_yyy_1, ta_xz_yyyy_0, \
                                         ta_xz_yyyy_1, ta_xz_yyyz_0, ta_xz_yyyz_1, ta_xz_yyz_0, ta_xz_yyz_1, ta_xz_yyzz_0, \
                                         ta_xz_yyzz_1, ta_xz_yzz_0, ta_xz_yzz_1, ta_xz_yzzz_0, ta_xz_yzzz_1, ta_xz_zzz_0, \
                                         ta_xz_zzz_1, ta_xz_zzzz_0, ta_xz_zzzz_1, ta_y_xxxx_0, ta_y_xxxx_1, ta_y_xxxy_0, \
                                         ta_y_xxxy_1, ta_y_xxxz_0, ta_y_xxxz_1, ta_y_xxyy_0, ta_y_xxyy_1, ta_y_xxyz_0, \
                                         ta_y_xxyz_1, ta_y_xxzz_0, ta_y_xxzz_1, ta_y_xyyy_0, ta_y_xyyy_1, ta_y_xyyz_0, \
                                         ta_y_xyyz_1, ta_y_xyzz_0, ta_y_xyzz_1, ta_y_xzzz_0, ta_y_xzzz_1, ta_y_yyyy_0, \
                                         ta_y_yyyy_1, ta_y_yyyz_0, ta_y_yyyz_1, ta_y_yyzz_0, ta_y_yyzz_1, ta_y_yzzz_0, \
                                         ta_y_yzzz_1, ta_y_zzzz_0, ta_y_zzzz_1, ta_yy_xxx_0, ta_yy_xxx_1, ta_yy_xxxx_0, \
                                         ta_yy_xxxx_1, ta_yy_xxxy_0, ta_yy_xxxy_1, ta_yy_xxxz_0, ta_yy_xxxz_1, ta_yy_xxy_0, \
                                         ta_yy_xxy_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyz_0, ta_yy_xxyz_1, ta_yy_xxz_0, \
                                         ta_yy_xxz_1, ta_yy_xyy_0, ta_yy_xyy_1, ta_yy_xyz_0, ta_yy_xyz_1, ta_z_xxxx_0, \
                                         ta_z_xxxx_1, ta_z_xxxy_0, ta_z_xxxy_1, ta_z_xxxz_0, ta_z_xxxz_1, ta_z_xxyy_0, \
                                         ta_z_xxyy_1, ta_z_xxyz_0, ta_z_xxyz_1, ta_z_xxzz_0, ta_z_xxzz_1, ta_z_xyyy_0, \
                                         ta_z_xyyy_1, ta_z_xyyz_0, ta_z_xyyz_1, ta_z_xyzz_0, ta_z_xyzz_1, ta_z_xzzz_0, \
                                         ta_z_xzzz_1, ta_z_yyyy_0, ta_z_yyyy_1, ta_z_yyyz_0, ta_z_yyyz_1, ta_z_yyzz_0, \
                                         ta_z_yyzz_1, ta_z_yzzz_0, ta_z_yzzz_1, ta_z_zzzz_0, ta_z_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxx_xxxx_0[j] = pa_x[j] * ta_xx_xxxx_0[j] - pc_x[j] * ta_xx_xxxx_1[j] + fl1_fx * ta_x_xxxx_0[j] - fl1_fx * ta_x_xxxx_1[j] +
                                   2.0 * fl1_fx * ta_xx_xxx_0[j] - 2.0 * fl1_fx * ta_xx_xxx_1[j];

                ta_xxx_xxxy_0[j] = pa_x[j] * ta_xx_xxxy_0[j] - pc_x[j] * ta_xx_xxxy_1[j] + fl1_fx * ta_x_xxxy_0[j] - fl1_fx * ta_x_xxxy_1[j] +
                                   1.5 * fl1_fx * ta_xx_xxy_0[j] - 1.5 * fl1_fx * ta_xx_xxy_1[j];

                ta_xxx_xxxz_0[j] = pa_x[j] * ta_xx_xxxz_0[j] - pc_x[j] * ta_xx_xxxz_1[j] + fl1_fx * ta_x_xxxz_0[j] - fl1_fx * ta_x_xxxz_1[j] +
                                   1.5 * fl1_fx * ta_xx_xxz_0[j] - 1.5 * fl1_fx * ta_xx_xxz_1[j];

                ta_xxx_xxyy_0[j] = pa_x[j] * ta_xx_xxyy_0[j] - pc_x[j] * ta_xx_xxyy_1[j] + fl1_fx * ta_x_xxyy_0[j] - fl1_fx * ta_x_xxyy_1[j] +
                                   fl1_fx * ta_xx_xyy_0[j] - fl1_fx * ta_xx_xyy_1[j];

                ta_xxx_xxyz_0[j] = pa_x[j] * ta_xx_xxyz_0[j] - pc_x[j] * ta_xx_xxyz_1[j] + fl1_fx * ta_x_xxyz_0[j] - fl1_fx * ta_x_xxyz_1[j] +
                                   fl1_fx * ta_xx_xyz_0[j] - fl1_fx * ta_xx_xyz_1[j];

                ta_xxx_xxzz_0[j] = pa_x[j] * ta_xx_xxzz_0[j] - pc_x[j] * ta_xx_xxzz_1[j] + fl1_fx * ta_x_xxzz_0[j] - fl1_fx * ta_x_xxzz_1[j] +
                                   fl1_fx * ta_xx_xzz_0[j] - fl1_fx * ta_xx_xzz_1[j];

                ta_xxx_xyyy_0[j] = pa_x[j] * ta_xx_xyyy_0[j] - pc_x[j] * ta_xx_xyyy_1[j] + fl1_fx * ta_x_xyyy_0[j] - fl1_fx * ta_x_xyyy_1[j] +
                                   0.5 * fl1_fx * ta_xx_yyy_0[j] - 0.5 * fl1_fx * ta_xx_yyy_1[j];

                ta_xxx_xyyz_0[j] = pa_x[j] * ta_xx_xyyz_0[j] - pc_x[j] * ta_xx_xyyz_1[j] + fl1_fx * ta_x_xyyz_0[j] - fl1_fx * ta_x_xyyz_1[j] +
                                   0.5 * fl1_fx * ta_xx_yyz_0[j] - 0.5 * fl1_fx * ta_xx_yyz_1[j];

                ta_xxx_xyzz_0[j] = pa_x[j] * ta_xx_xyzz_0[j] - pc_x[j] * ta_xx_xyzz_1[j] + fl1_fx * ta_x_xyzz_0[j] - fl1_fx * ta_x_xyzz_1[j] +
                                   0.5 * fl1_fx * ta_xx_yzz_0[j] - 0.5 * fl1_fx * ta_xx_yzz_1[j];

                ta_xxx_xzzz_0[j] = pa_x[j] * ta_xx_xzzz_0[j] - pc_x[j] * ta_xx_xzzz_1[j] + fl1_fx * ta_x_xzzz_0[j] - fl1_fx * ta_x_xzzz_1[j] +
                                   0.5 * fl1_fx * ta_xx_zzz_0[j] - 0.5 * fl1_fx * ta_xx_zzz_1[j];

                ta_xxx_yyyy_0[j] = pa_x[j] * ta_xx_yyyy_0[j] - pc_x[j] * ta_xx_yyyy_1[j] + fl1_fx * ta_x_yyyy_0[j] - fl1_fx * ta_x_yyyy_1[j];

                ta_xxx_yyyz_0[j] = pa_x[j] * ta_xx_yyyz_0[j] - pc_x[j] * ta_xx_yyyz_1[j] + fl1_fx * ta_x_yyyz_0[j] - fl1_fx * ta_x_yyyz_1[j];

                ta_xxx_yyzz_0[j] = pa_x[j] * ta_xx_yyzz_0[j] - pc_x[j] * ta_xx_yyzz_1[j] + fl1_fx * ta_x_yyzz_0[j] - fl1_fx * ta_x_yyzz_1[j];

                ta_xxx_yzzz_0[j] = pa_x[j] * ta_xx_yzzz_0[j] - pc_x[j] * ta_xx_yzzz_1[j] + fl1_fx * ta_x_yzzz_0[j] - fl1_fx * ta_x_yzzz_1[j];

                ta_xxx_zzzz_0[j] = pa_x[j] * ta_xx_zzzz_0[j] - pc_x[j] * ta_xx_zzzz_1[j] + fl1_fx * ta_x_zzzz_0[j] - fl1_fx * ta_x_zzzz_1[j];

                ta_xxy_xxxx_0[j] = pa_x[j] * ta_xy_xxxx_0[j] - pc_x[j] * ta_xy_xxxx_1[j] + 0.5 * fl1_fx * ta_y_xxxx_0[j] -
                                   0.5 * fl1_fx * ta_y_xxxx_1[j] + 2.0 * fl1_fx * ta_xy_xxx_0[j] - 2.0 * fl1_fx * ta_xy_xxx_1[j];

                ta_xxy_xxxy_0[j] = pa_x[j] * ta_xy_xxxy_0[j] - pc_x[j] * ta_xy_xxxy_1[j] + 0.5 * fl1_fx * ta_y_xxxy_0[j] -
                                   0.5 * fl1_fx * ta_y_xxxy_1[j] + 1.5 * fl1_fx * ta_xy_xxy_0[j] - 1.5 * fl1_fx * ta_xy_xxy_1[j];

                ta_xxy_xxxz_0[j] = pa_x[j] * ta_xy_xxxz_0[j] - pc_x[j] * ta_xy_xxxz_1[j] + 0.5 * fl1_fx * ta_y_xxxz_0[j] -
                                   0.5 * fl1_fx * ta_y_xxxz_1[j] + 1.5 * fl1_fx * ta_xy_xxz_0[j] - 1.5 * fl1_fx * ta_xy_xxz_1[j];

                ta_xxy_xxyy_0[j] = pa_x[j] * ta_xy_xxyy_0[j] - pc_x[j] * ta_xy_xxyy_1[j] + 0.5 * fl1_fx * ta_y_xxyy_0[j] -
                                   0.5 * fl1_fx * ta_y_xxyy_1[j] + fl1_fx * ta_xy_xyy_0[j] - fl1_fx * ta_xy_xyy_1[j];

                ta_xxy_xxyz_0[j] = pa_x[j] * ta_xy_xxyz_0[j] - pc_x[j] * ta_xy_xxyz_1[j] + 0.5 * fl1_fx * ta_y_xxyz_0[j] -
                                   0.5 * fl1_fx * ta_y_xxyz_1[j] + fl1_fx * ta_xy_xyz_0[j] - fl1_fx * ta_xy_xyz_1[j];

                ta_xxy_xxzz_0[j] = pa_x[j] * ta_xy_xxzz_0[j] - pc_x[j] * ta_xy_xxzz_1[j] + 0.5 * fl1_fx * ta_y_xxzz_0[j] -
                                   0.5 * fl1_fx * ta_y_xxzz_1[j] + fl1_fx * ta_xy_xzz_0[j] - fl1_fx * ta_xy_xzz_1[j];

                ta_xxy_xyyy_0[j] = pa_x[j] * ta_xy_xyyy_0[j] - pc_x[j] * ta_xy_xyyy_1[j] + 0.5 * fl1_fx * ta_y_xyyy_0[j] -
                                   0.5 * fl1_fx * ta_y_xyyy_1[j] + 0.5 * fl1_fx * ta_xy_yyy_0[j] - 0.5 * fl1_fx * ta_xy_yyy_1[j];

                ta_xxy_xyyz_0[j] = pa_x[j] * ta_xy_xyyz_0[j] - pc_x[j] * ta_xy_xyyz_1[j] + 0.5 * fl1_fx * ta_y_xyyz_0[j] -
                                   0.5 * fl1_fx * ta_y_xyyz_1[j] + 0.5 * fl1_fx * ta_xy_yyz_0[j] - 0.5 * fl1_fx * ta_xy_yyz_1[j];

                ta_xxy_xyzz_0[j] = pa_x[j] * ta_xy_xyzz_0[j] - pc_x[j] * ta_xy_xyzz_1[j] + 0.5 * fl1_fx * ta_y_xyzz_0[j] -
                                   0.5 * fl1_fx * ta_y_xyzz_1[j] + 0.5 * fl1_fx * ta_xy_yzz_0[j] - 0.5 * fl1_fx * ta_xy_yzz_1[j];

                ta_xxy_xzzz_0[j] = pa_x[j] * ta_xy_xzzz_0[j] - pc_x[j] * ta_xy_xzzz_1[j] + 0.5 * fl1_fx * ta_y_xzzz_0[j] -
                                   0.5 * fl1_fx * ta_y_xzzz_1[j] + 0.5 * fl1_fx * ta_xy_zzz_0[j] - 0.5 * fl1_fx * ta_xy_zzz_1[j];

                ta_xxy_yyyy_0[j] =
                    pa_x[j] * ta_xy_yyyy_0[j] - pc_x[j] * ta_xy_yyyy_1[j] + 0.5 * fl1_fx * ta_y_yyyy_0[j] - 0.5 * fl1_fx * ta_y_yyyy_1[j];

                ta_xxy_yyyz_0[j] =
                    pa_x[j] * ta_xy_yyyz_0[j] - pc_x[j] * ta_xy_yyyz_1[j] + 0.5 * fl1_fx * ta_y_yyyz_0[j] - 0.5 * fl1_fx * ta_y_yyyz_1[j];

                ta_xxy_yyzz_0[j] =
                    pa_x[j] * ta_xy_yyzz_0[j] - pc_x[j] * ta_xy_yyzz_1[j] + 0.5 * fl1_fx * ta_y_yyzz_0[j] - 0.5 * fl1_fx * ta_y_yyzz_1[j];

                ta_xxy_yzzz_0[j] =
                    pa_x[j] * ta_xy_yzzz_0[j] - pc_x[j] * ta_xy_yzzz_1[j] + 0.5 * fl1_fx * ta_y_yzzz_0[j] - 0.5 * fl1_fx * ta_y_yzzz_1[j];

                ta_xxy_zzzz_0[j] =
                    pa_x[j] * ta_xy_zzzz_0[j] - pc_x[j] * ta_xy_zzzz_1[j] + 0.5 * fl1_fx * ta_y_zzzz_0[j] - 0.5 * fl1_fx * ta_y_zzzz_1[j];

                ta_xxz_xxxx_0[j] = pa_x[j] * ta_xz_xxxx_0[j] - pc_x[j] * ta_xz_xxxx_1[j] + 0.5 * fl1_fx * ta_z_xxxx_0[j] -
                                   0.5 * fl1_fx * ta_z_xxxx_1[j] + 2.0 * fl1_fx * ta_xz_xxx_0[j] - 2.0 * fl1_fx * ta_xz_xxx_1[j];

                ta_xxz_xxxy_0[j] = pa_x[j] * ta_xz_xxxy_0[j] - pc_x[j] * ta_xz_xxxy_1[j] + 0.5 * fl1_fx * ta_z_xxxy_0[j] -
                                   0.5 * fl1_fx * ta_z_xxxy_1[j] + 1.5 * fl1_fx * ta_xz_xxy_0[j] - 1.5 * fl1_fx * ta_xz_xxy_1[j];

                ta_xxz_xxxz_0[j] = pa_x[j] * ta_xz_xxxz_0[j] - pc_x[j] * ta_xz_xxxz_1[j] + 0.5 * fl1_fx * ta_z_xxxz_0[j] -
                                   0.5 * fl1_fx * ta_z_xxxz_1[j] + 1.5 * fl1_fx * ta_xz_xxz_0[j] - 1.5 * fl1_fx * ta_xz_xxz_1[j];

                ta_xxz_xxyy_0[j] = pa_x[j] * ta_xz_xxyy_0[j] - pc_x[j] * ta_xz_xxyy_1[j] + 0.5 * fl1_fx * ta_z_xxyy_0[j] -
                                   0.5 * fl1_fx * ta_z_xxyy_1[j] + fl1_fx * ta_xz_xyy_0[j] - fl1_fx * ta_xz_xyy_1[j];

                ta_xxz_xxyz_0[j] = pa_x[j] * ta_xz_xxyz_0[j] - pc_x[j] * ta_xz_xxyz_1[j] + 0.5 * fl1_fx * ta_z_xxyz_0[j] -
                                   0.5 * fl1_fx * ta_z_xxyz_1[j] + fl1_fx * ta_xz_xyz_0[j] - fl1_fx * ta_xz_xyz_1[j];

                ta_xxz_xxzz_0[j] = pa_x[j] * ta_xz_xxzz_0[j] - pc_x[j] * ta_xz_xxzz_1[j] + 0.5 * fl1_fx * ta_z_xxzz_0[j] -
                                   0.5 * fl1_fx * ta_z_xxzz_1[j] + fl1_fx * ta_xz_xzz_0[j] - fl1_fx * ta_xz_xzz_1[j];

                ta_xxz_xyyy_0[j] = pa_x[j] * ta_xz_xyyy_0[j] - pc_x[j] * ta_xz_xyyy_1[j] + 0.5 * fl1_fx * ta_z_xyyy_0[j] -
                                   0.5 * fl1_fx * ta_z_xyyy_1[j] + 0.5 * fl1_fx * ta_xz_yyy_0[j] - 0.5 * fl1_fx * ta_xz_yyy_1[j];

                ta_xxz_xyyz_0[j] = pa_x[j] * ta_xz_xyyz_0[j] - pc_x[j] * ta_xz_xyyz_1[j] + 0.5 * fl1_fx * ta_z_xyyz_0[j] -
                                   0.5 * fl1_fx * ta_z_xyyz_1[j] + 0.5 * fl1_fx * ta_xz_yyz_0[j] - 0.5 * fl1_fx * ta_xz_yyz_1[j];

                ta_xxz_xyzz_0[j] = pa_x[j] * ta_xz_xyzz_0[j] - pc_x[j] * ta_xz_xyzz_1[j] + 0.5 * fl1_fx * ta_z_xyzz_0[j] -
                                   0.5 * fl1_fx * ta_z_xyzz_1[j] + 0.5 * fl1_fx * ta_xz_yzz_0[j] - 0.5 * fl1_fx * ta_xz_yzz_1[j];

                ta_xxz_xzzz_0[j] = pa_x[j] * ta_xz_xzzz_0[j] - pc_x[j] * ta_xz_xzzz_1[j] + 0.5 * fl1_fx * ta_z_xzzz_0[j] -
                                   0.5 * fl1_fx * ta_z_xzzz_1[j] + 0.5 * fl1_fx * ta_xz_zzz_0[j] - 0.5 * fl1_fx * ta_xz_zzz_1[j];

                ta_xxz_yyyy_0[j] =
                    pa_x[j] * ta_xz_yyyy_0[j] - pc_x[j] * ta_xz_yyyy_1[j] + 0.5 * fl1_fx * ta_z_yyyy_0[j] - 0.5 * fl1_fx * ta_z_yyyy_1[j];

                ta_xxz_yyyz_0[j] =
                    pa_x[j] * ta_xz_yyyz_0[j] - pc_x[j] * ta_xz_yyyz_1[j] + 0.5 * fl1_fx * ta_z_yyyz_0[j] - 0.5 * fl1_fx * ta_z_yyyz_1[j];

                ta_xxz_yyzz_0[j] =
                    pa_x[j] * ta_xz_yyzz_0[j] - pc_x[j] * ta_xz_yyzz_1[j] + 0.5 * fl1_fx * ta_z_yyzz_0[j] - 0.5 * fl1_fx * ta_z_yyzz_1[j];

                ta_xxz_yzzz_0[j] =
                    pa_x[j] * ta_xz_yzzz_0[j] - pc_x[j] * ta_xz_yzzz_1[j] + 0.5 * fl1_fx * ta_z_yzzz_0[j] - 0.5 * fl1_fx * ta_z_yzzz_1[j];

                ta_xxz_zzzz_0[j] =
                    pa_x[j] * ta_xz_zzzz_0[j] - pc_x[j] * ta_xz_zzzz_1[j] + 0.5 * fl1_fx * ta_z_zzzz_0[j] - 0.5 * fl1_fx * ta_z_zzzz_1[j];

                ta_xyy_xxxx_0[j] =
                    pa_x[j] * ta_yy_xxxx_0[j] - pc_x[j] * ta_yy_xxxx_1[j] + 2.0 * fl1_fx * ta_yy_xxx_0[j] - 2.0 * fl1_fx * ta_yy_xxx_1[j];

                ta_xyy_xxxy_0[j] =
                    pa_x[j] * ta_yy_xxxy_0[j] - pc_x[j] * ta_yy_xxxy_1[j] + 1.5 * fl1_fx * ta_yy_xxy_0[j] - 1.5 * fl1_fx * ta_yy_xxy_1[j];

                ta_xyy_xxxz_0[j] =
                    pa_x[j] * ta_yy_xxxz_0[j] - pc_x[j] * ta_yy_xxxz_1[j] + 1.5 * fl1_fx * ta_yy_xxz_0[j] - 1.5 * fl1_fx * ta_yy_xxz_1[j];

                ta_xyy_xxyy_0[j] = pa_x[j] * ta_yy_xxyy_0[j] - pc_x[j] * ta_yy_xxyy_1[j] + fl1_fx * ta_yy_xyy_0[j] - fl1_fx * ta_yy_xyy_1[j];

                ta_xyy_xxyz_0[j] = pa_x[j] * ta_yy_xxyz_0[j] - pc_x[j] * ta_yy_xxyz_1[j] + fl1_fx * ta_yy_xyz_0[j] - fl1_fx * ta_yy_xyz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFG_50_100(CMemBlock2D<double>&       primBuffer,
                                 const CRecursionMap&       recursionMap,
                                 const CMemBlock2D<double>& osFactors,
                                 const CMemBlock2D<double>& paDistances,
                                 const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto ta_yy_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 45);

            auto ta_yy_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 46);

            auto ta_yy_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 47);

            auto ta_yy_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 48);

            auto ta_yy_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 49);

            auto ta_yy_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 50);

            auto ta_yy_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 51);

            auto ta_yy_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 52);

            auto ta_yy_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 53);

            auto ta_yy_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 54);

            auto ta_yy_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 55);

            auto ta_yy_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 56);

            auto ta_yy_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 57);

            auto ta_yy_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 58);

            auto ta_yy_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 59);

            auto ta_yz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 60);

            auto ta_yz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 61);

            auto ta_yz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 62);

            auto ta_yz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 63);

            auto ta_yz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 64);

            auto ta_yz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 65);

            auto ta_yz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 66);

            auto ta_yz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 67);

            auto ta_yz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 68);

            auto ta_yz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 69);

            auto ta_yz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 70);

            auto ta_yz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 71);

            auto ta_yz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 72);

            auto ta_yz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 73);

            auto ta_yz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 74);

            auto ta_zz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 75);

            auto ta_zz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 76);

            auto ta_zz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 77);

            auto ta_zz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 78);

            auto ta_zz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 79);

            auto ta_zz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 80);

            auto ta_zz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 81);

            auto ta_zz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 82);

            auto ta_zz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 83);

            auto ta_zz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 84);

            auto ta_zz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 85);

            auto ta_zz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 86);

            auto ta_zz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 87);

            auto ta_zz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 88);

            auto ta_zz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 89);

            auto ta_yy_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 45);

            auto ta_yy_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 46);

            auto ta_yy_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 47);

            auto ta_yy_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 48);

            auto ta_yy_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 49);

            auto ta_yy_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 50);

            auto ta_yy_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 51);

            auto ta_yy_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 52);

            auto ta_yy_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 53);

            auto ta_yy_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 54);

            auto ta_yy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 55);

            auto ta_yy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 56);

            auto ta_yy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 57);

            auto ta_yy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 58);

            auto ta_yy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 59);

            auto ta_yz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 60);

            auto ta_yz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 61);

            auto ta_yz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 62);

            auto ta_yz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 63);

            auto ta_yz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 64);

            auto ta_yz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 65);

            auto ta_yz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 66);

            auto ta_yz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 67);

            auto ta_yz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 68);

            auto ta_yz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 69);

            auto ta_yz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 70);

            auto ta_yz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 71);

            auto ta_yz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 72);

            auto ta_yz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 73);

            auto ta_yz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 74);

            auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75);

            auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76);

            auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77);

            auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78);

            auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79);

            auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80);

            auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81);

            auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82);

            auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83);

            auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84);

            auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85);

            auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86);

            auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87);

            auto ta_zz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 88);

            auto ta_zz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 89);

            auto ta_y_xxxx_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 15);

            auto ta_y_xxxy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 16);

            auto ta_y_xxxz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 17);

            auto ta_y_xxyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 18);

            auto ta_y_xxyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 19);

            auto ta_y_xxzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 20);

            auto ta_y_xyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 21);

            auto ta_y_xyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 22);

            auto ta_y_xyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 23);

            auto ta_y_xzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 24);

            auto ta_y_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 15);

            auto ta_y_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 16);

            auto ta_y_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 17);

            auto ta_y_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 18);

            auto ta_y_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 19);

            auto ta_y_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 20);

            auto ta_y_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 21);

            auto ta_y_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 22);

            auto ta_y_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 23);

            auto ta_y_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 24);

            auto ta_yy_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 30);

            auto ta_yy_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 31);

            auto ta_yy_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 32);

            auto ta_yy_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 33);

            auto ta_yy_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 34);

            auto ta_yy_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 35);

            auto ta_yy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 36);

            auto ta_yy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 37);

            auto ta_yy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 38);

            auto ta_yy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 39);

            auto ta_yz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 40);

            auto ta_yz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 41);

            auto ta_yz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 42);

            auto ta_yz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 43);

            auto ta_yz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 44);

            auto ta_yz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 45);

            auto ta_yz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 46);

            auto ta_yz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 47);

            auto ta_yz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 48);

            auto ta_yz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 49);

            auto ta_zz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 50);

            auto ta_zz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 51);

            auto ta_zz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 52);

            auto ta_zz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 53);

            auto ta_zz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 54);

            auto ta_zz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 55);

            auto ta_zz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 56);

            auto ta_zz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 57);

            auto ta_zz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 58);

            auto ta_zz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 59);

            auto ta_yy_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 30);

            auto ta_yy_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 31);

            auto ta_yy_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 32);

            auto ta_yy_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 33);

            auto ta_yy_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 34);

            auto ta_yy_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 35);

            auto ta_yy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 36);

            auto ta_yy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 37);

            auto ta_yy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 38);

            auto ta_yy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 39);

            auto ta_yz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 40);

            auto ta_yz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 41);

            auto ta_yz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 42);

            auto ta_yz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 43);

            auto ta_yz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 44);

            auto ta_yz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 45);

            auto ta_yz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 46);

            auto ta_yz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 47);

            auto ta_yz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 48);

            auto ta_yz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 49);

            auto ta_zz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 50);

            auto ta_zz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 51);

            auto ta_zz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 52);

            auto ta_zz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 53);

            auto ta_zz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 54);

            auto ta_zz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 55);

            auto ta_zz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 56);

            auto ta_zz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 57);

            auto ta_zz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 58);

            auto ta_zz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 59);

            // set up pointers to integrals

            auto ta_xyy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 50);

            auto ta_xyy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 51);

            auto ta_xyy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 52);

            auto ta_xyy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 53);

            auto ta_xyy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 54);

            auto ta_xyy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 55);

            auto ta_xyy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 56);

            auto ta_xyy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 57);

            auto ta_xyy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 58);

            auto ta_xyy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 59);

            auto ta_xyz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 60);

            auto ta_xyz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 61);

            auto ta_xyz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 62);

            auto ta_xyz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 63);

            auto ta_xyz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 64);

            auto ta_xyz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 65);

            auto ta_xyz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 66);

            auto ta_xyz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 67);

            auto ta_xyz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 68);

            auto ta_xyz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 69);

            auto ta_xyz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 70);

            auto ta_xyz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 71);

            auto ta_xyz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 72);

            auto ta_xyz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 73);

            auto ta_xyz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 74);

            auto ta_xzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 75);

            auto ta_xzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 76);

            auto ta_xzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 77);

            auto ta_xzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 78);

            auto ta_xzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 79);

            auto ta_xzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 80);

            auto ta_xzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 81);

            auto ta_xzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 82);

            auto ta_xzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 83);

            auto ta_xzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 84);

            auto ta_xzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 85);

            auto ta_xzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 86);

            auto ta_xzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 87);

            auto ta_xzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 88);

            auto ta_xzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 89);

            auto ta_yyy_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 90);

            auto ta_yyy_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 91);

            auto ta_yyy_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 92);

            auto ta_yyy_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 93);

            auto ta_yyy_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 94);

            auto ta_yyy_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 95);

            auto ta_yyy_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 96);

            auto ta_yyy_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 97);

            auto ta_yyy_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 98);

            auto ta_yyy_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 99);

            // Batch of Integrals (50,100)

            #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_xyy_xxzz_0, ta_xyy_xyyy_0, ta_xyy_xyyz_0, \
                                         ta_xyy_xyzz_0, ta_xyy_xzzz_0, ta_xyy_yyyy_0, ta_xyy_yyyz_0, ta_xyy_yyzz_0, \
                                         ta_xyy_yzzz_0, ta_xyy_zzzz_0, ta_xyz_xxxx_0, ta_xyz_xxxy_0, ta_xyz_xxxz_0, \
                                         ta_xyz_xxyy_0, ta_xyz_xxyz_0, ta_xyz_xxzz_0, ta_xyz_xyyy_0, ta_xyz_xyyz_0, \
                                         ta_xyz_xyzz_0, ta_xyz_xzzz_0, ta_xyz_yyyy_0, ta_xyz_yyyz_0, ta_xyz_yyzz_0, \
                                         ta_xyz_yzzz_0, ta_xyz_zzzz_0, ta_xzz_xxxx_0, ta_xzz_xxxy_0, ta_xzz_xxxz_0, \
                                         ta_xzz_xxyy_0, ta_xzz_xxyz_0, ta_xzz_xxzz_0, ta_xzz_xyyy_0, ta_xzz_xyyz_0, \
                                         ta_xzz_xyzz_0, ta_xzz_xzzz_0, ta_xzz_yyyy_0, ta_xzz_yyyz_0, ta_xzz_yyzz_0, \
                                         ta_xzz_yzzz_0, ta_xzz_zzzz_0, ta_y_xxxx_0, ta_y_xxxx_1, ta_y_xxxy_0, ta_y_xxxy_1, \
                                         ta_y_xxxz_0, ta_y_xxxz_1, ta_y_xxyy_0, ta_y_xxyy_1, ta_y_xxyz_0, ta_y_xxyz_1, \
                                         ta_y_xxzz_0, ta_y_xxzz_1, ta_y_xyyy_0, ta_y_xyyy_1, ta_y_xyyz_0, ta_y_xyyz_1, \
                                         ta_y_xyzz_0, ta_y_xyzz_1, ta_y_xzzz_0, ta_y_xzzz_1, ta_yy_xxx_0, ta_yy_xxx_1, \
                                         ta_yy_xxxx_0, ta_yy_xxxx_1, ta_yy_xxxy_0, ta_yy_xxxy_1, ta_yy_xxxz_0, ta_yy_xxxz_1, \
                                         ta_yy_xxy_0, ta_yy_xxy_1, ta_yy_xxyy_0, ta_yy_xxyy_1, ta_yy_xxyz_0, ta_yy_xxyz_1, \
                                         ta_yy_xxz_0, ta_yy_xxz_1, ta_yy_xxzz_0, ta_yy_xxzz_1, ta_yy_xyy_0, ta_yy_xyy_1, \
                                         ta_yy_xyyy_0, ta_yy_xyyy_1, ta_yy_xyyz_0, ta_yy_xyyz_1, ta_yy_xyz_0, ta_yy_xyz_1, \
                                         ta_yy_xyzz_0, ta_yy_xyzz_1, ta_yy_xzz_0, ta_yy_xzz_1, ta_yy_xzzz_0, ta_yy_xzzz_1, \
                                         ta_yy_yyy_0, ta_yy_yyy_1, ta_yy_yyyy_0, ta_yy_yyyy_1, ta_yy_yyyz_0, ta_yy_yyyz_1, \
                                         ta_yy_yyz_0, ta_yy_yyz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yzz_0, ta_yy_yzz_1, \
                                         ta_yy_yzzz_0, ta_yy_yzzz_1, ta_yy_zzz_0, ta_yy_zzz_1, ta_yy_zzzz_0, ta_yy_zzzz_1, \
                                         ta_yyy_xxxx_0, ta_yyy_xxxy_0, ta_yyy_xxxz_0, ta_yyy_xxyy_0, ta_yyy_xxyz_0, \
                                         ta_yyy_xxzz_0, ta_yyy_xyyy_0, ta_yyy_xyyz_0, ta_yyy_xyzz_0, ta_yyy_xzzz_0, \
                                         ta_yz_xxx_0, ta_yz_xxx_1, ta_yz_xxxx_0, ta_yz_xxxx_1, ta_yz_xxxy_0, ta_yz_xxxy_1, \
                                         ta_yz_xxxz_0, ta_yz_xxxz_1, ta_yz_xxy_0, ta_yz_xxy_1, ta_yz_xxyy_0, ta_yz_xxyy_1, \
                                         ta_yz_xxyz_0, ta_yz_xxyz_1, ta_yz_xxz_0, ta_yz_xxz_1, ta_yz_xxzz_0, ta_yz_xxzz_1, \
                                         ta_yz_xyy_0, ta_yz_xyy_1, ta_yz_xyyy_0, ta_yz_xyyy_1, ta_yz_xyyz_0, ta_yz_xyyz_1, \
                                         ta_yz_xyz_0, ta_yz_xyz_1, ta_yz_xyzz_0, ta_yz_xyzz_1, ta_yz_xzz_0, ta_yz_xzz_1, \
                                         ta_yz_xzzz_0, ta_yz_xzzz_1, ta_yz_yyy_0, ta_yz_yyy_1, ta_yz_yyyy_0, ta_yz_yyyy_1, \
                                         ta_yz_yyyz_0, ta_yz_yyyz_1, ta_yz_yyz_0, ta_yz_yyz_1, ta_yz_yyzz_0, ta_yz_yyzz_1, \
                                         ta_yz_yzz_0, ta_yz_yzz_1, ta_yz_yzzz_0, ta_yz_yzzz_1, ta_yz_zzz_0, ta_yz_zzz_1, \
                                         ta_yz_zzzz_0, ta_yz_zzzz_1, ta_zz_xxx_0, ta_zz_xxx_1, ta_zz_xxxx_0, ta_zz_xxxx_1, \
                                         ta_zz_xxxy_0, ta_zz_xxxy_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxy_0, ta_zz_xxy_1, \
                                         ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyz_0, ta_zz_xxyz_1, ta_zz_xxz_0, ta_zz_xxz_1, \
                                         ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xyy_0, ta_zz_xyy_1, ta_zz_xyyy_0, ta_zz_xyyy_1, \
                                         ta_zz_xyyz_0, ta_zz_xyyz_1, ta_zz_xyz_0, ta_zz_xyz_1, ta_zz_xyzz_0, ta_zz_xyzz_1, \
                                         ta_zz_xzz_0, ta_zz_xzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_yyy_0, ta_zz_yyy_1, \
                                         ta_zz_yyyy_0, ta_zz_yyyy_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyz_0, ta_zz_yyz_1, \
                                         ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yzz_0, ta_zz_yzz_1, ta_zz_yzzz_0, ta_zz_yzzz_1, \
                                         ta_zz_zzz_0, ta_zz_zzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xyy_xxzz_0[j] = pa_x[j] * ta_yy_xxzz_0[j] - pc_x[j] * ta_yy_xxzz_1[j] + fl1_fx * ta_yy_xzz_0[j] - fl1_fx * ta_yy_xzz_1[j];

                ta_xyy_xyyy_0[j] =
                    pa_x[j] * ta_yy_xyyy_0[j] - pc_x[j] * ta_yy_xyyy_1[j] + 0.5 * fl1_fx * ta_yy_yyy_0[j] - 0.5 * fl1_fx * ta_yy_yyy_1[j];

                ta_xyy_xyyz_0[j] =
                    pa_x[j] * ta_yy_xyyz_0[j] - pc_x[j] * ta_yy_xyyz_1[j] + 0.5 * fl1_fx * ta_yy_yyz_0[j] - 0.5 * fl1_fx * ta_yy_yyz_1[j];

                ta_xyy_xyzz_0[j] =
                    pa_x[j] * ta_yy_xyzz_0[j] - pc_x[j] * ta_yy_xyzz_1[j] + 0.5 * fl1_fx * ta_yy_yzz_0[j] - 0.5 * fl1_fx * ta_yy_yzz_1[j];

                ta_xyy_xzzz_0[j] =
                    pa_x[j] * ta_yy_xzzz_0[j] - pc_x[j] * ta_yy_xzzz_1[j] + 0.5 * fl1_fx * ta_yy_zzz_0[j] - 0.5 * fl1_fx * ta_yy_zzz_1[j];

                ta_xyy_yyyy_0[j] = pa_x[j] * ta_yy_yyyy_0[j] - pc_x[j] * ta_yy_yyyy_1[j];

                ta_xyy_yyyz_0[j] = pa_x[j] * ta_yy_yyyz_0[j] - pc_x[j] * ta_yy_yyyz_1[j];

                ta_xyy_yyzz_0[j] = pa_x[j] * ta_yy_yyzz_0[j] - pc_x[j] * ta_yy_yyzz_1[j];

                ta_xyy_yzzz_0[j] = pa_x[j] * ta_yy_yzzz_0[j] - pc_x[j] * ta_yy_yzzz_1[j];

                ta_xyy_zzzz_0[j] = pa_x[j] * ta_yy_zzzz_0[j] - pc_x[j] * ta_yy_zzzz_1[j];

                ta_xyz_xxxx_0[j] =
                    pa_x[j] * ta_yz_xxxx_0[j] - pc_x[j] * ta_yz_xxxx_1[j] + 2.0 * fl1_fx * ta_yz_xxx_0[j] - 2.0 * fl1_fx * ta_yz_xxx_1[j];

                ta_xyz_xxxy_0[j] =
                    pa_x[j] * ta_yz_xxxy_0[j] - pc_x[j] * ta_yz_xxxy_1[j] + 1.5 * fl1_fx * ta_yz_xxy_0[j] - 1.5 * fl1_fx * ta_yz_xxy_1[j];

                ta_xyz_xxxz_0[j] =
                    pa_x[j] * ta_yz_xxxz_0[j] - pc_x[j] * ta_yz_xxxz_1[j] + 1.5 * fl1_fx * ta_yz_xxz_0[j] - 1.5 * fl1_fx * ta_yz_xxz_1[j];

                ta_xyz_xxyy_0[j] = pa_x[j] * ta_yz_xxyy_0[j] - pc_x[j] * ta_yz_xxyy_1[j] + fl1_fx * ta_yz_xyy_0[j] - fl1_fx * ta_yz_xyy_1[j];

                ta_xyz_xxyz_0[j] = pa_x[j] * ta_yz_xxyz_0[j] - pc_x[j] * ta_yz_xxyz_1[j] + fl1_fx * ta_yz_xyz_0[j] - fl1_fx * ta_yz_xyz_1[j];

                ta_xyz_xxzz_0[j] = pa_x[j] * ta_yz_xxzz_0[j] - pc_x[j] * ta_yz_xxzz_1[j] + fl1_fx * ta_yz_xzz_0[j] - fl1_fx * ta_yz_xzz_1[j];

                ta_xyz_xyyy_0[j] =
                    pa_x[j] * ta_yz_xyyy_0[j] - pc_x[j] * ta_yz_xyyy_1[j] + 0.5 * fl1_fx * ta_yz_yyy_0[j] - 0.5 * fl1_fx * ta_yz_yyy_1[j];

                ta_xyz_xyyz_0[j] =
                    pa_x[j] * ta_yz_xyyz_0[j] - pc_x[j] * ta_yz_xyyz_1[j] + 0.5 * fl1_fx * ta_yz_yyz_0[j] - 0.5 * fl1_fx * ta_yz_yyz_1[j];

                ta_xyz_xyzz_0[j] =
                    pa_x[j] * ta_yz_xyzz_0[j] - pc_x[j] * ta_yz_xyzz_1[j] + 0.5 * fl1_fx * ta_yz_yzz_0[j] - 0.5 * fl1_fx * ta_yz_yzz_1[j];

                ta_xyz_xzzz_0[j] =
                    pa_x[j] * ta_yz_xzzz_0[j] - pc_x[j] * ta_yz_xzzz_1[j] + 0.5 * fl1_fx * ta_yz_zzz_0[j] - 0.5 * fl1_fx * ta_yz_zzz_1[j];

                ta_xyz_yyyy_0[j] = pa_x[j] * ta_yz_yyyy_0[j] - pc_x[j] * ta_yz_yyyy_1[j];

                ta_xyz_yyyz_0[j] = pa_x[j] * ta_yz_yyyz_0[j] - pc_x[j] * ta_yz_yyyz_1[j];

                ta_xyz_yyzz_0[j] = pa_x[j] * ta_yz_yyzz_0[j] - pc_x[j] * ta_yz_yyzz_1[j];

                ta_xyz_yzzz_0[j] = pa_x[j] * ta_yz_yzzz_0[j] - pc_x[j] * ta_yz_yzzz_1[j];

                ta_xyz_zzzz_0[j] = pa_x[j] * ta_yz_zzzz_0[j] - pc_x[j] * ta_yz_zzzz_1[j];

                ta_xzz_xxxx_0[j] =
                    pa_x[j] * ta_zz_xxxx_0[j] - pc_x[j] * ta_zz_xxxx_1[j] + 2.0 * fl1_fx * ta_zz_xxx_0[j] - 2.0 * fl1_fx * ta_zz_xxx_1[j];

                ta_xzz_xxxy_0[j] =
                    pa_x[j] * ta_zz_xxxy_0[j] - pc_x[j] * ta_zz_xxxy_1[j] + 1.5 * fl1_fx * ta_zz_xxy_0[j] - 1.5 * fl1_fx * ta_zz_xxy_1[j];

                ta_xzz_xxxz_0[j] =
                    pa_x[j] * ta_zz_xxxz_0[j] - pc_x[j] * ta_zz_xxxz_1[j] + 1.5 * fl1_fx * ta_zz_xxz_0[j] - 1.5 * fl1_fx * ta_zz_xxz_1[j];

                ta_xzz_xxyy_0[j] = pa_x[j] * ta_zz_xxyy_0[j] - pc_x[j] * ta_zz_xxyy_1[j] + fl1_fx * ta_zz_xyy_0[j] - fl1_fx * ta_zz_xyy_1[j];

                ta_xzz_xxyz_0[j] = pa_x[j] * ta_zz_xxyz_0[j] - pc_x[j] * ta_zz_xxyz_1[j] + fl1_fx * ta_zz_xyz_0[j] - fl1_fx * ta_zz_xyz_1[j];

                ta_xzz_xxzz_0[j] = pa_x[j] * ta_zz_xxzz_0[j] - pc_x[j] * ta_zz_xxzz_1[j] + fl1_fx * ta_zz_xzz_0[j] - fl1_fx * ta_zz_xzz_1[j];

                ta_xzz_xyyy_0[j] =
                    pa_x[j] * ta_zz_xyyy_0[j] - pc_x[j] * ta_zz_xyyy_1[j] + 0.5 * fl1_fx * ta_zz_yyy_0[j] - 0.5 * fl1_fx * ta_zz_yyy_1[j];

                ta_xzz_xyyz_0[j] =
                    pa_x[j] * ta_zz_xyyz_0[j] - pc_x[j] * ta_zz_xyyz_1[j] + 0.5 * fl1_fx * ta_zz_yyz_0[j] - 0.5 * fl1_fx * ta_zz_yyz_1[j];

                ta_xzz_xyzz_0[j] =
                    pa_x[j] * ta_zz_xyzz_0[j] - pc_x[j] * ta_zz_xyzz_1[j] + 0.5 * fl1_fx * ta_zz_yzz_0[j] - 0.5 * fl1_fx * ta_zz_yzz_1[j];

                ta_xzz_xzzz_0[j] =
                    pa_x[j] * ta_zz_xzzz_0[j] - pc_x[j] * ta_zz_xzzz_1[j] + 0.5 * fl1_fx * ta_zz_zzz_0[j] - 0.5 * fl1_fx * ta_zz_zzz_1[j];

                ta_xzz_yyyy_0[j] = pa_x[j] * ta_zz_yyyy_0[j] - pc_x[j] * ta_zz_yyyy_1[j];

                ta_xzz_yyyz_0[j] = pa_x[j] * ta_zz_yyyz_0[j] - pc_x[j] * ta_zz_yyyz_1[j];

                ta_xzz_yyzz_0[j] = pa_x[j] * ta_zz_yyzz_0[j] - pc_x[j] * ta_zz_yyzz_1[j];

                ta_xzz_yzzz_0[j] = pa_x[j] * ta_zz_yzzz_0[j] - pc_x[j] * ta_zz_yzzz_1[j];

                ta_xzz_zzzz_0[j] = pa_x[j] * ta_zz_zzzz_0[j] - pc_x[j] * ta_zz_zzzz_1[j];

                ta_yyy_xxxx_0[j] = pa_y[j] * ta_yy_xxxx_0[j] - pc_y[j] * ta_yy_xxxx_1[j] + fl1_fx * ta_y_xxxx_0[j] - fl1_fx * ta_y_xxxx_1[j];

                ta_yyy_xxxy_0[j] = pa_y[j] * ta_yy_xxxy_0[j] - pc_y[j] * ta_yy_xxxy_1[j] + fl1_fx * ta_y_xxxy_0[j] - fl1_fx * ta_y_xxxy_1[j] +
                                   0.5 * fl1_fx * ta_yy_xxx_0[j] - 0.5 * fl1_fx * ta_yy_xxx_1[j];

                ta_yyy_xxxz_0[j] = pa_y[j] * ta_yy_xxxz_0[j] - pc_y[j] * ta_yy_xxxz_1[j] + fl1_fx * ta_y_xxxz_0[j] - fl1_fx * ta_y_xxxz_1[j];

                ta_yyy_xxyy_0[j] = pa_y[j] * ta_yy_xxyy_0[j] - pc_y[j] * ta_yy_xxyy_1[j] + fl1_fx * ta_y_xxyy_0[j] - fl1_fx * ta_y_xxyy_1[j] +
                                   fl1_fx * ta_yy_xxy_0[j] - fl1_fx * ta_yy_xxy_1[j];

                ta_yyy_xxyz_0[j] = pa_y[j] * ta_yy_xxyz_0[j] - pc_y[j] * ta_yy_xxyz_1[j] + fl1_fx * ta_y_xxyz_0[j] - fl1_fx * ta_y_xxyz_1[j] +
                                   0.5 * fl1_fx * ta_yy_xxz_0[j] - 0.5 * fl1_fx * ta_yy_xxz_1[j];

                ta_yyy_xxzz_0[j] = pa_y[j] * ta_yy_xxzz_0[j] - pc_y[j] * ta_yy_xxzz_1[j] + fl1_fx * ta_y_xxzz_0[j] - fl1_fx * ta_y_xxzz_1[j];

                ta_yyy_xyyy_0[j] = pa_y[j] * ta_yy_xyyy_0[j] - pc_y[j] * ta_yy_xyyy_1[j] + fl1_fx * ta_y_xyyy_0[j] - fl1_fx * ta_y_xyyy_1[j] +
                                   1.5 * fl1_fx * ta_yy_xyy_0[j] - 1.5 * fl1_fx * ta_yy_xyy_1[j];

                ta_yyy_xyyz_0[j] = pa_y[j] * ta_yy_xyyz_0[j] - pc_y[j] * ta_yy_xyyz_1[j] + fl1_fx * ta_y_xyyz_0[j] - fl1_fx * ta_y_xyyz_1[j] +
                                   fl1_fx * ta_yy_xyz_0[j] - fl1_fx * ta_yy_xyz_1[j];

                ta_yyy_xyzz_0[j] = pa_y[j] * ta_yy_xyzz_0[j] - pc_y[j] * ta_yy_xyzz_1[j] + fl1_fx * ta_y_xyzz_0[j] - fl1_fx * ta_y_xyzz_1[j] +
                                   0.5 * fl1_fx * ta_yy_xzz_0[j] - 0.5 * fl1_fx * ta_yy_xzz_1[j];

                ta_yyy_xzzz_0[j] = pa_y[j] * ta_yy_xzzz_0[j] - pc_y[j] * ta_yy_xzzz_1[j] + fl1_fx * ta_y_xzzz_0[j] - fl1_fx * ta_y_xzzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFG_100_150(CMemBlock2D<double>&       primBuffer,
                                  const CRecursionMap&       recursionMap,
                                  const CMemBlock2D<double>& osFactors,
                                  const CMemBlock2D<double>& paDistances,
                                  const CMemBlock2D<double>& pcDistances,
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

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ta_yy_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 55);

            auto ta_yy_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 56);

            auto ta_yy_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 57);

            auto ta_yy_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 58);

            auto ta_yy_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 59);

            auto ta_yz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 60);

            auto ta_yz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 61);

            auto ta_yz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 62);

            auto ta_yz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 63);

            auto ta_yz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 64);

            auto ta_yz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 65);

            auto ta_yz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 66);

            auto ta_yz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 67);

            auto ta_yz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 68);

            auto ta_yz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 69);

            auto ta_yz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 70);

            auto ta_yz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 71);

            auto ta_yz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 72);

            auto ta_yz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 73);

            auto ta_yz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 74);

            auto ta_zz_xxxx_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 75);

            auto ta_zz_xxxy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 76);

            auto ta_zz_xxxz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 77);

            auto ta_zz_xxyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 78);

            auto ta_zz_xxyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 79);

            auto ta_zz_xxzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 80);

            auto ta_zz_xyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 81);

            auto ta_zz_xyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 82);

            auto ta_zz_xyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 83);

            auto ta_zz_xzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 84);

            auto ta_zz_yyyy_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 85);

            auto ta_zz_yyyz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 86);

            auto ta_zz_yyzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 87);

            auto ta_zz_yzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 88);

            auto ta_zz_zzzz_0 = primBuffer.data(pidx_a_2_4_m0 + 90 * idx + 89);

            auto ta_yy_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 55);

            auto ta_yy_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 56);

            auto ta_yy_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 57);

            auto ta_yy_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 58);

            auto ta_yy_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 59);

            auto ta_yz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 60);

            auto ta_yz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 61);

            auto ta_yz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 62);

            auto ta_yz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 63);

            auto ta_yz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 64);

            auto ta_yz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 65);

            auto ta_yz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 66);

            auto ta_yz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 67);

            auto ta_yz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 68);

            auto ta_yz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 69);

            auto ta_yz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 70);

            auto ta_yz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 71);

            auto ta_yz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 72);

            auto ta_yz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 73);

            auto ta_yz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 74);

            auto ta_zz_xxxx_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 75);

            auto ta_zz_xxxy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 76);

            auto ta_zz_xxxz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 77);

            auto ta_zz_xxyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 78);

            auto ta_zz_xxyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 79);

            auto ta_zz_xxzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 80);

            auto ta_zz_xyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 81);

            auto ta_zz_xyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 82);

            auto ta_zz_xyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 83);

            auto ta_zz_xzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 84);

            auto ta_zz_yyyy_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 85);

            auto ta_zz_yyyz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 86);

            auto ta_zz_yyzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 87);

            auto ta_zz_yzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 88);

            auto ta_zz_zzzz_1 = primBuffer.data(pidx_a_2_4_m1 + 90 * idx + 89);

            auto ta_y_yyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 25);

            auto ta_y_yyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 26);

            auto ta_y_yyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 27);

            auto ta_y_yzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 28);

            auto ta_y_zzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 29);

            auto ta_z_xxxx_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 30);

            auto ta_z_xxxy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 31);

            auto ta_z_xxxz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 32);

            auto ta_z_xxyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 33);

            auto ta_z_xxyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 34);

            auto ta_z_xxzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 35);

            auto ta_z_xyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 36);

            auto ta_z_xyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 37);

            auto ta_z_xyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 38);

            auto ta_z_xzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 39);

            auto ta_z_yyyy_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 40);

            auto ta_z_yyyz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 41);

            auto ta_z_yyzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 42);

            auto ta_z_yzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 43);

            auto ta_z_zzzz_0 = primBuffer.data(pidx_a_1_4_m0 + 45 * idx + 44);

            auto ta_y_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 25);

            auto ta_y_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 26);

            auto ta_y_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 27);

            auto ta_y_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 28);

            auto ta_y_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 29);

            auto ta_z_xxxx_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 30);

            auto ta_z_xxxy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 31);

            auto ta_z_xxxz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 32);

            auto ta_z_xxyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 33);

            auto ta_z_xxyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 34);

            auto ta_z_xxzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 35);

            auto ta_z_xyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 36);

            auto ta_z_xyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 37);

            auto ta_z_xyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 38);

            auto ta_z_xzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 39);

            auto ta_z_yyyy_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 40);

            auto ta_z_yyyz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 41);

            auto ta_z_yyzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 42);

            auto ta_z_yzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 43);

            auto ta_z_zzzz_1 = primBuffer.data(pidx_a_1_4_m1 + 45 * idx + 44);

            auto ta_yy_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 36);

            auto ta_yy_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 37);

            auto ta_yy_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 38);

            auto ta_yy_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 39);

            auto ta_yz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 40);

            auto ta_yz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 41);

            auto ta_yz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 42);

            auto ta_yz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 43);

            auto ta_yz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 44);

            auto ta_yz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 45);

            auto ta_yz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 46);

            auto ta_yz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 47);

            auto ta_yz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 48);

            auto ta_yz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 49);

            auto ta_zz_xxx_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 50);

            auto ta_zz_xxy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 51);

            auto ta_zz_xxz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 52);

            auto ta_zz_xyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 53);

            auto ta_zz_xyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 54);

            auto ta_zz_xzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 55);

            auto ta_zz_yyy_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 56);

            auto ta_zz_yyz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 57);

            auto ta_zz_yzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 58);

            auto ta_zz_zzz_0 = primBuffer.data(pidx_a_2_3_m0 + 60 * idx + 59);

            auto ta_yy_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 36);

            auto ta_yy_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 37);

            auto ta_yy_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 38);

            auto ta_yy_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 39);

            auto ta_yz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 40);

            auto ta_yz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 41);

            auto ta_yz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 42);

            auto ta_yz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 43);

            auto ta_yz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 44);

            auto ta_yz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 45);

            auto ta_yz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 46);

            auto ta_yz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 47);

            auto ta_yz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 48);

            auto ta_yz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 49);

            auto ta_zz_xxx_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 50);

            auto ta_zz_xxy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 51);

            auto ta_zz_xxz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 52);

            auto ta_zz_xyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 53);

            auto ta_zz_xyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 54);

            auto ta_zz_xzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 55);

            auto ta_zz_yyy_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 56);

            auto ta_zz_yyz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 57);

            auto ta_zz_yzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 58);

            auto ta_zz_zzz_1 = primBuffer.data(pidx_a_2_3_m1 + 60 * idx + 59);

            // set up pointers to integrals

            auto ta_yyy_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 100);

            auto ta_yyy_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 101);

            auto ta_yyy_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 102);

            auto ta_yyy_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 103);

            auto ta_yyy_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 104);

            auto ta_yyz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 105);

            auto ta_yyz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 106);

            auto ta_yyz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 107);

            auto ta_yyz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 108);

            auto ta_yyz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 109);

            auto ta_yyz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 110);

            auto ta_yyz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 111);

            auto ta_yyz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 112);

            auto ta_yyz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 113);

            auto ta_yyz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 114);

            auto ta_yyz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 115);

            auto ta_yyz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 116);

            auto ta_yyz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 117);

            auto ta_yyz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 118);

            auto ta_yyz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 119);

            auto ta_yzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 120);

            auto ta_yzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 121);

            auto ta_yzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 122);

            auto ta_yzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 123);

            auto ta_yzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 124);

            auto ta_yzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 125);

            auto ta_yzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 126);

            auto ta_yzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 127);

            auto ta_yzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 128);

            auto ta_yzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 129);

            auto ta_yzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 130);

            auto ta_yzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 131);

            auto ta_yzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 132);

            auto ta_yzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 133);

            auto ta_yzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 134);

            auto ta_zzz_xxxx_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 135);

            auto ta_zzz_xxxy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 136);

            auto ta_zzz_xxxz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 137);

            auto ta_zzz_xxyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 138);

            auto ta_zzz_xxyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 139);

            auto ta_zzz_xxzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 140);

            auto ta_zzz_xyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 141);

            auto ta_zzz_xyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 142);

            auto ta_zzz_xyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 143);

            auto ta_zzz_xzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 144);

            auto ta_zzz_yyyy_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 145);

            auto ta_zzz_yyyz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 146);

            auto ta_zzz_yyzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 147);

            auto ta_zzz_yzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 148);

            auto ta_zzz_zzzz_0 = primBuffer.data(pidx_a_3_4_m0 + 150 * idx + 149);

            // Batch of Integrals (100,150)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_y_yyyy_0, ta_y_yyyy_1, ta_y_yyyz_0, \
                                         ta_y_yyyz_1, ta_y_yyzz_0, ta_y_yyzz_1, ta_y_yzzz_0, ta_y_yzzz_1, ta_y_zzzz_0, \
                                         ta_y_zzzz_1, ta_yy_yyy_0, ta_yy_yyy_1, ta_yy_yyyy_0, ta_yy_yyyy_1, ta_yy_yyyz_0, \
                                         ta_yy_yyyz_1, ta_yy_yyz_0, ta_yy_yyz_1, ta_yy_yyzz_0, ta_yy_yyzz_1, ta_yy_yzz_0, \
                                         ta_yy_yzz_1, ta_yy_yzzz_0, ta_yy_yzzz_1, ta_yy_zzz_0, ta_yy_zzz_1, ta_yy_zzzz_0, \
                                         ta_yy_zzzz_1, ta_yyy_yyyy_0, ta_yyy_yyyz_0, ta_yyy_yyzz_0, ta_yyy_yzzz_0, \
                                         ta_yyy_zzzz_0, ta_yyz_xxxx_0, ta_yyz_xxxy_0, ta_yyz_xxxz_0, ta_yyz_xxyy_0, \
                                         ta_yyz_xxyz_0, ta_yyz_xxzz_0, ta_yyz_xyyy_0, ta_yyz_xyyz_0, ta_yyz_xyzz_0, \
                                         ta_yyz_xzzz_0, ta_yyz_yyyy_0, ta_yyz_yyyz_0, ta_yyz_yyzz_0, ta_yyz_yzzz_0, \
                                         ta_yyz_zzzz_0, ta_yz_xxx_0, ta_yz_xxx_1, ta_yz_xxxx_0, ta_yz_xxxx_1, ta_yz_xxxy_0, \
                                         ta_yz_xxxy_1, ta_yz_xxxz_0, ta_yz_xxxz_1, ta_yz_xxy_0, ta_yz_xxy_1, ta_yz_xxyy_0, \
                                         ta_yz_xxyy_1, ta_yz_xxyz_0, ta_yz_xxyz_1, ta_yz_xxz_0, ta_yz_xxz_1, ta_yz_xxzz_0, \
                                         ta_yz_xxzz_1, ta_yz_xyy_0, ta_yz_xyy_1, ta_yz_xyyy_0, ta_yz_xyyy_1, ta_yz_xyyz_0, \
                                         ta_yz_xyyz_1, ta_yz_xyz_0, ta_yz_xyz_1, ta_yz_xyzz_0, ta_yz_xyzz_1, ta_yz_xzz_0, \
                                         ta_yz_xzz_1, ta_yz_xzzz_0, ta_yz_xzzz_1, ta_yz_yyy_0, ta_yz_yyy_1, ta_yz_yyyy_0, \
                                         ta_yz_yyyy_1, ta_yz_yyyz_0, ta_yz_yyyz_1, ta_yz_yyz_0, ta_yz_yyz_1, ta_yz_yyzz_0, \
                                         ta_yz_yyzz_1, ta_yz_yzz_0, ta_yz_yzz_1, ta_yz_yzzz_0, ta_yz_yzzz_1, ta_yz_zzz_0, \
                                         ta_yz_zzz_1, ta_yz_zzzz_0, ta_yz_zzzz_1, ta_yzz_xxxx_0, ta_yzz_xxxy_0, \
                                         ta_yzz_xxxz_0, ta_yzz_xxyy_0, ta_yzz_xxyz_0, ta_yzz_xxzz_0, ta_yzz_xyyy_0, \
                                         ta_yzz_xyyz_0, ta_yzz_xyzz_0, ta_yzz_xzzz_0, ta_yzz_yyyy_0, ta_yzz_yyyz_0, \
                                         ta_yzz_yyzz_0, ta_yzz_yzzz_0, ta_yzz_zzzz_0, ta_z_xxxx_0, ta_z_xxxx_1, ta_z_xxxy_0, \
                                         ta_z_xxxy_1, ta_z_xxxz_0, ta_z_xxxz_1, ta_z_xxyy_0, ta_z_xxyy_1, ta_z_xxyz_0, \
                                         ta_z_xxyz_1, ta_z_xxzz_0, ta_z_xxzz_1, ta_z_xyyy_0, ta_z_xyyy_1, ta_z_xyyz_0, \
                                         ta_z_xyyz_1, ta_z_xyzz_0, ta_z_xyzz_1, ta_z_xzzz_0, ta_z_xzzz_1, ta_z_yyyy_0, \
                                         ta_z_yyyy_1, ta_z_yyyz_0, ta_z_yyyz_1, ta_z_yyzz_0, ta_z_yyzz_1, ta_z_yzzz_0, \
                                         ta_z_yzzz_1, ta_z_zzzz_0, ta_z_zzzz_1, ta_zz_xxx_0, ta_zz_xxx_1, ta_zz_xxxx_0, \
                                         ta_zz_xxxx_1, ta_zz_xxxy_0, ta_zz_xxxy_1, ta_zz_xxxz_0, ta_zz_xxxz_1, ta_zz_xxy_0, \
                                         ta_zz_xxy_1, ta_zz_xxyy_0, ta_zz_xxyy_1, ta_zz_xxyz_0, ta_zz_xxyz_1, ta_zz_xxz_0, \
                                         ta_zz_xxz_1, ta_zz_xxzz_0, ta_zz_xxzz_1, ta_zz_xyy_0, ta_zz_xyy_1, ta_zz_xyyy_0, \
                                         ta_zz_xyyy_1, ta_zz_xyyz_0, ta_zz_xyyz_1, ta_zz_xyz_0, ta_zz_xyz_1, ta_zz_xyzz_0, \
                                         ta_zz_xyzz_1, ta_zz_xzz_0, ta_zz_xzz_1, ta_zz_xzzz_0, ta_zz_xzzz_1, ta_zz_yyy_0, \
                                         ta_zz_yyy_1, ta_zz_yyyy_0, ta_zz_yyyy_1, ta_zz_yyyz_0, ta_zz_yyyz_1, ta_zz_yyz_0, \
                                         ta_zz_yyz_1, ta_zz_yyzz_0, ta_zz_yyzz_1, ta_zz_yzz_0, ta_zz_yzz_1, ta_zz_yzzz_0, \
                                         ta_zz_yzzz_1, ta_zz_zzz_0, ta_zz_zzz_1, ta_zz_zzzz_0, ta_zz_zzzz_1, ta_zzz_xxxx_0, \
                                         ta_zzz_xxxy_0, ta_zzz_xxxz_0, ta_zzz_xxyy_0, ta_zzz_xxyz_0, ta_zzz_xxzz_0, \
                                         ta_zzz_xyyy_0, ta_zzz_xyyz_0, ta_zzz_xyzz_0, ta_zzz_xzzz_0, ta_zzz_yyyy_0, \
                                         ta_zzz_yyyz_0, ta_zzz_yyzz_0, ta_zzz_yzzz_0, ta_zzz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_yyy_yyyy_0[j] = pa_y[j] * ta_yy_yyyy_0[j] - pc_y[j] * ta_yy_yyyy_1[j] + fl1_fx * ta_y_yyyy_0[j] - fl1_fx * ta_y_yyyy_1[j] +
                                   2.0 * fl1_fx * ta_yy_yyy_0[j] - 2.0 * fl1_fx * ta_yy_yyy_1[j];

                ta_yyy_yyyz_0[j] = pa_y[j] * ta_yy_yyyz_0[j] - pc_y[j] * ta_yy_yyyz_1[j] + fl1_fx * ta_y_yyyz_0[j] - fl1_fx * ta_y_yyyz_1[j] +
                                   1.5 * fl1_fx * ta_yy_yyz_0[j] - 1.5 * fl1_fx * ta_yy_yyz_1[j];

                ta_yyy_yyzz_0[j] = pa_y[j] * ta_yy_yyzz_0[j] - pc_y[j] * ta_yy_yyzz_1[j] + fl1_fx * ta_y_yyzz_0[j] - fl1_fx * ta_y_yyzz_1[j] +
                                   fl1_fx * ta_yy_yzz_0[j] - fl1_fx * ta_yy_yzz_1[j];

                ta_yyy_yzzz_0[j] = pa_y[j] * ta_yy_yzzz_0[j] - pc_y[j] * ta_yy_yzzz_1[j] + fl1_fx * ta_y_yzzz_0[j] - fl1_fx * ta_y_yzzz_1[j] +
                                   0.5 * fl1_fx * ta_yy_zzz_0[j] - 0.5 * fl1_fx * ta_yy_zzz_1[j];

                ta_yyy_zzzz_0[j] = pa_y[j] * ta_yy_zzzz_0[j] - pc_y[j] * ta_yy_zzzz_1[j] + fl1_fx * ta_y_zzzz_0[j] - fl1_fx * ta_y_zzzz_1[j];

                ta_yyz_xxxx_0[j] =
                    pa_y[j] * ta_yz_xxxx_0[j] - pc_y[j] * ta_yz_xxxx_1[j] + 0.5 * fl1_fx * ta_z_xxxx_0[j] - 0.5 * fl1_fx * ta_z_xxxx_1[j];

                ta_yyz_xxxy_0[j] = pa_y[j] * ta_yz_xxxy_0[j] - pc_y[j] * ta_yz_xxxy_1[j] + 0.5 * fl1_fx * ta_z_xxxy_0[j] -
                                   0.5 * fl1_fx * ta_z_xxxy_1[j] + 0.5 * fl1_fx * ta_yz_xxx_0[j] - 0.5 * fl1_fx * ta_yz_xxx_1[j];

                ta_yyz_xxxz_0[j] =
                    pa_y[j] * ta_yz_xxxz_0[j] - pc_y[j] * ta_yz_xxxz_1[j] + 0.5 * fl1_fx * ta_z_xxxz_0[j] - 0.5 * fl1_fx * ta_z_xxxz_1[j];

                ta_yyz_xxyy_0[j] = pa_y[j] * ta_yz_xxyy_0[j] - pc_y[j] * ta_yz_xxyy_1[j] + 0.5 * fl1_fx * ta_z_xxyy_0[j] -
                                   0.5 * fl1_fx * ta_z_xxyy_1[j] + fl1_fx * ta_yz_xxy_0[j] - fl1_fx * ta_yz_xxy_1[j];

                ta_yyz_xxyz_0[j] = pa_y[j] * ta_yz_xxyz_0[j] - pc_y[j] * ta_yz_xxyz_1[j] + 0.5 * fl1_fx * ta_z_xxyz_0[j] -
                                   0.5 * fl1_fx * ta_z_xxyz_1[j] + 0.5 * fl1_fx * ta_yz_xxz_0[j] - 0.5 * fl1_fx * ta_yz_xxz_1[j];

                ta_yyz_xxzz_0[j] =
                    pa_y[j] * ta_yz_xxzz_0[j] - pc_y[j] * ta_yz_xxzz_1[j] + 0.5 * fl1_fx * ta_z_xxzz_0[j] - 0.5 * fl1_fx * ta_z_xxzz_1[j];

                ta_yyz_xyyy_0[j] = pa_y[j] * ta_yz_xyyy_0[j] - pc_y[j] * ta_yz_xyyy_1[j] + 0.5 * fl1_fx * ta_z_xyyy_0[j] -
                                   0.5 * fl1_fx * ta_z_xyyy_1[j] + 1.5 * fl1_fx * ta_yz_xyy_0[j] - 1.5 * fl1_fx * ta_yz_xyy_1[j];

                ta_yyz_xyyz_0[j] = pa_y[j] * ta_yz_xyyz_0[j] - pc_y[j] * ta_yz_xyyz_1[j] + 0.5 * fl1_fx * ta_z_xyyz_0[j] -
                                   0.5 * fl1_fx * ta_z_xyyz_1[j] + fl1_fx * ta_yz_xyz_0[j] - fl1_fx * ta_yz_xyz_1[j];

                ta_yyz_xyzz_0[j] = pa_y[j] * ta_yz_xyzz_0[j] - pc_y[j] * ta_yz_xyzz_1[j] + 0.5 * fl1_fx * ta_z_xyzz_0[j] -
                                   0.5 * fl1_fx * ta_z_xyzz_1[j] + 0.5 * fl1_fx * ta_yz_xzz_0[j] - 0.5 * fl1_fx * ta_yz_xzz_1[j];

                ta_yyz_xzzz_0[j] =
                    pa_y[j] * ta_yz_xzzz_0[j] - pc_y[j] * ta_yz_xzzz_1[j] + 0.5 * fl1_fx * ta_z_xzzz_0[j] - 0.5 * fl1_fx * ta_z_xzzz_1[j];

                ta_yyz_yyyy_0[j] = pa_y[j] * ta_yz_yyyy_0[j] - pc_y[j] * ta_yz_yyyy_1[j] + 0.5 * fl1_fx * ta_z_yyyy_0[j] -
                                   0.5 * fl1_fx * ta_z_yyyy_1[j] + 2.0 * fl1_fx * ta_yz_yyy_0[j] - 2.0 * fl1_fx * ta_yz_yyy_1[j];

                ta_yyz_yyyz_0[j] = pa_y[j] * ta_yz_yyyz_0[j] - pc_y[j] * ta_yz_yyyz_1[j] + 0.5 * fl1_fx * ta_z_yyyz_0[j] -
                                   0.5 * fl1_fx * ta_z_yyyz_1[j] + 1.5 * fl1_fx * ta_yz_yyz_0[j] - 1.5 * fl1_fx * ta_yz_yyz_1[j];

                ta_yyz_yyzz_0[j] = pa_y[j] * ta_yz_yyzz_0[j] - pc_y[j] * ta_yz_yyzz_1[j] + 0.5 * fl1_fx * ta_z_yyzz_0[j] -
                                   0.5 * fl1_fx * ta_z_yyzz_1[j] + fl1_fx * ta_yz_yzz_0[j] - fl1_fx * ta_yz_yzz_1[j];

                ta_yyz_yzzz_0[j] = pa_y[j] * ta_yz_yzzz_0[j] - pc_y[j] * ta_yz_yzzz_1[j] + 0.5 * fl1_fx * ta_z_yzzz_0[j] -
                                   0.5 * fl1_fx * ta_z_yzzz_1[j] + 0.5 * fl1_fx * ta_yz_zzz_0[j] - 0.5 * fl1_fx * ta_yz_zzz_1[j];

                ta_yyz_zzzz_0[j] =
                    pa_y[j] * ta_yz_zzzz_0[j] - pc_y[j] * ta_yz_zzzz_1[j] + 0.5 * fl1_fx * ta_z_zzzz_0[j] - 0.5 * fl1_fx * ta_z_zzzz_1[j];

                ta_yzz_xxxx_0[j] = pa_y[j] * ta_zz_xxxx_0[j] - pc_y[j] * ta_zz_xxxx_1[j];

                ta_yzz_xxxy_0[j] =
                    pa_y[j] * ta_zz_xxxy_0[j] - pc_y[j] * ta_zz_xxxy_1[j] + 0.5 * fl1_fx * ta_zz_xxx_0[j] - 0.5 * fl1_fx * ta_zz_xxx_1[j];

                ta_yzz_xxxz_0[j] = pa_y[j] * ta_zz_xxxz_0[j] - pc_y[j] * ta_zz_xxxz_1[j];

                ta_yzz_xxyy_0[j] = pa_y[j] * ta_zz_xxyy_0[j] - pc_y[j] * ta_zz_xxyy_1[j] + fl1_fx * ta_zz_xxy_0[j] - fl1_fx * ta_zz_xxy_1[j];

                ta_yzz_xxyz_0[j] =
                    pa_y[j] * ta_zz_xxyz_0[j] - pc_y[j] * ta_zz_xxyz_1[j] + 0.5 * fl1_fx * ta_zz_xxz_0[j] - 0.5 * fl1_fx * ta_zz_xxz_1[j];

                ta_yzz_xxzz_0[j] = pa_y[j] * ta_zz_xxzz_0[j] - pc_y[j] * ta_zz_xxzz_1[j];

                ta_yzz_xyyy_0[j] =
                    pa_y[j] * ta_zz_xyyy_0[j] - pc_y[j] * ta_zz_xyyy_1[j] + 1.5 * fl1_fx * ta_zz_xyy_0[j] - 1.5 * fl1_fx * ta_zz_xyy_1[j];

                ta_yzz_xyyz_0[j] = pa_y[j] * ta_zz_xyyz_0[j] - pc_y[j] * ta_zz_xyyz_1[j] + fl1_fx * ta_zz_xyz_0[j] - fl1_fx * ta_zz_xyz_1[j];

                ta_yzz_xyzz_0[j] =
                    pa_y[j] * ta_zz_xyzz_0[j] - pc_y[j] * ta_zz_xyzz_1[j] + 0.5 * fl1_fx * ta_zz_xzz_0[j] - 0.5 * fl1_fx * ta_zz_xzz_1[j];

                ta_yzz_xzzz_0[j] = pa_y[j] * ta_zz_xzzz_0[j] - pc_y[j] * ta_zz_xzzz_1[j];

                ta_yzz_yyyy_0[j] =
                    pa_y[j] * ta_zz_yyyy_0[j] - pc_y[j] * ta_zz_yyyy_1[j] + 2.0 * fl1_fx * ta_zz_yyy_0[j] - 2.0 * fl1_fx * ta_zz_yyy_1[j];

                ta_yzz_yyyz_0[j] =
                    pa_y[j] * ta_zz_yyyz_0[j] - pc_y[j] * ta_zz_yyyz_1[j] + 1.5 * fl1_fx * ta_zz_yyz_0[j] - 1.5 * fl1_fx * ta_zz_yyz_1[j];

                ta_yzz_yyzz_0[j] = pa_y[j] * ta_zz_yyzz_0[j] - pc_y[j] * ta_zz_yyzz_1[j] + fl1_fx * ta_zz_yzz_0[j] - fl1_fx * ta_zz_yzz_1[j];

                ta_yzz_yzzz_0[j] =
                    pa_y[j] * ta_zz_yzzz_0[j] - pc_y[j] * ta_zz_yzzz_1[j] + 0.5 * fl1_fx * ta_zz_zzz_0[j] - 0.5 * fl1_fx * ta_zz_zzz_1[j];

                ta_yzz_zzzz_0[j] = pa_y[j] * ta_zz_zzzz_0[j] - pc_y[j] * ta_zz_zzzz_1[j];

                ta_zzz_xxxx_0[j] = pa_z[j] * ta_zz_xxxx_0[j] - pc_z[j] * ta_zz_xxxx_1[j] + fl1_fx * ta_z_xxxx_0[j] - fl1_fx * ta_z_xxxx_1[j];

                ta_zzz_xxxy_0[j] = pa_z[j] * ta_zz_xxxy_0[j] - pc_z[j] * ta_zz_xxxy_1[j] + fl1_fx * ta_z_xxxy_0[j] - fl1_fx * ta_z_xxxy_1[j];

                ta_zzz_xxxz_0[j] = pa_z[j] * ta_zz_xxxz_0[j] - pc_z[j] * ta_zz_xxxz_1[j] + fl1_fx * ta_z_xxxz_0[j] - fl1_fx * ta_z_xxxz_1[j] +
                                   0.5 * fl1_fx * ta_zz_xxx_0[j] - 0.5 * fl1_fx * ta_zz_xxx_1[j];

                ta_zzz_xxyy_0[j] = pa_z[j] * ta_zz_xxyy_0[j] - pc_z[j] * ta_zz_xxyy_1[j] + fl1_fx * ta_z_xxyy_0[j] - fl1_fx * ta_z_xxyy_1[j];

                ta_zzz_xxyz_0[j] = pa_z[j] * ta_zz_xxyz_0[j] - pc_z[j] * ta_zz_xxyz_1[j] + fl1_fx * ta_z_xxyz_0[j] - fl1_fx * ta_z_xxyz_1[j] +
                                   0.5 * fl1_fx * ta_zz_xxy_0[j] - 0.5 * fl1_fx * ta_zz_xxy_1[j];

                ta_zzz_xxzz_0[j] = pa_z[j] * ta_zz_xxzz_0[j] - pc_z[j] * ta_zz_xxzz_1[j] + fl1_fx * ta_z_xxzz_0[j] - fl1_fx * ta_z_xxzz_1[j] +
                                   fl1_fx * ta_zz_xxz_0[j] - fl1_fx * ta_zz_xxz_1[j];

                ta_zzz_xyyy_0[j] = pa_z[j] * ta_zz_xyyy_0[j] - pc_z[j] * ta_zz_xyyy_1[j] + fl1_fx * ta_z_xyyy_0[j] - fl1_fx * ta_z_xyyy_1[j];

                ta_zzz_xyyz_0[j] = pa_z[j] * ta_zz_xyyz_0[j] - pc_z[j] * ta_zz_xyyz_1[j] + fl1_fx * ta_z_xyyz_0[j] - fl1_fx * ta_z_xyyz_1[j] +
                                   0.5 * fl1_fx * ta_zz_xyy_0[j] - 0.5 * fl1_fx * ta_zz_xyy_1[j];

                ta_zzz_xyzz_0[j] = pa_z[j] * ta_zz_xyzz_0[j] - pc_z[j] * ta_zz_xyzz_1[j] + fl1_fx * ta_z_xyzz_0[j] - fl1_fx * ta_z_xyzz_1[j] +
                                   fl1_fx * ta_zz_xyz_0[j] - fl1_fx * ta_zz_xyz_1[j];

                ta_zzz_xzzz_0[j] = pa_z[j] * ta_zz_xzzz_0[j] - pc_z[j] * ta_zz_xzzz_1[j] + fl1_fx * ta_z_xzzz_0[j] - fl1_fx * ta_z_xzzz_1[j] +
                                   1.5 * fl1_fx * ta_zz_xzz_0[j] - 1.5 * fl1_fx * ta_zz_xzz_1[j];

                ta_zzz_yyyy_0[j] = pa_z[j] * ta_zz_yyyy_0[j] - pc_z[j] * ta_zz_yyyy_1[j] + fl1_fx * ta_z_yyyy_0[j] - fl1_fx * ta_z_yyyy_1[j];

                ta_zzz_yyyz_0[j] = pa_z[j] * ta_zz_yyyz_0[j] - pc_z[j] * ta_zz_yyyz_1[j] + fl1_fx * ta_z_yyyz_0[j] - fl1_fx * ta_z_yyyz_1[j] +
                                   0.5 * fl1_fx * ta_zz_yyy_0[j] - 0.5 * fl1_fx * ta_zz_yyy_1[j];

                ta_zzz_yyzz_0[j] = pa_z[j] * ta_zz_yyzz_0[j] - pc_z[j] * ta_zz_yyzz_1[j] + fl1_fx * ta_z_yyzz_0[j] - fl1_fx * ta_z_yyzz_1[j] +
                                   fl1_fx * ta_zz_yyz_0[j] - fl1_fx * ta_zz_yyz_1[j];

                ta_zzz_yzzz_0[j] = pa_z[j] * ta_zz_yzzz_0[j] - pc_z[j] * ta_zz_yzzz_1[j] + fl1_fx * ta_z_yzzz_0[j] - fl1_fx * ta_z_yzzz_1[j] +
                                   1.5 * fl1_fx * ta_zz_yzz_0[j] - 1.5 * fl1_fx * ta_zz_yzz_1[j];

                ta_zzz_zzzz_0[j] = pa_z[j] * ta_zz_zzzz_0[j] - pc_z[j] * ta_zz_zzzz_1[j] + fl1_fx * ta_z_zzzz_0[j] - fl1_fx * ta_z_zzzz_1[j] +
                                   2.0 * fl1_fx * ta_zz_zzz_0[j] - 2.0 * fl1_fx * ta_zz_zzz_1[j];
            }

            idx++;
        }
    }
}

}  // namespace npotrecfunc
