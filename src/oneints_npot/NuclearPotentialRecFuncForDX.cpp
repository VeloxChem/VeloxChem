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

#include "NuclearPotentialRecFuncForDX.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

void
compNuclearPotentialForDD(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ta_x_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx);

            auto ta_x_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 1);

            auto ta_x_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 2);

            auto ta_x_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 3);

            auto ta_x_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 4);

            auto ta_x_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 5);

            auto ta_y_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 6);

            auto ta_y_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 7);

            auto ta_y_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 8);

            auto ta_y_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 9);

            auto ta_y_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 10);

            auto ta_y_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 11);

            auto ta_z_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 12);

            auto ta_z_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 13);

            auto ta_z_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 14);

            auto ta_z_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 15);

            auto ta_z_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 16);

            auto ta_z_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 17);

            auto ta_x_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx);

            auto ta_x_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 1);

            auto ta_x_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 2);

            auto ta_x_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 3);

            auto ta_x_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 4);

            auto ta_x_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 5);

            auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6);

            auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7);

            auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8);

            auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9);

            auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10);

            auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11);

            auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12);

            auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13);

            auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14);

            auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15);

            auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16);

            auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17);

            auto ta_0_xx_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx);

            auto ta_0_xy_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 1);

            auto ta_0_xz_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 2);

            auto ta_0_yy_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 3);

            auto ta_0_yz_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 4);

            auto ta_0_zz_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 5);

            auto ta_0_xx_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx);

            auto ta_0_xy_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 1);

            auto ta_0_xz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 2);

            auto ta_0_yy_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 3);

            auto ta_0_yz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 4);

            auto ta_0_zz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 5);

            auto ta_x_x_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx);

            auto ta_x_y_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 1);

            auto ta_x_z_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 2);

            auto ta_y_x_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 3);

            auto ta_y_y_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 4);

            auto ta_y_z_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 5);

            auto ta_z_x_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 6);

            auto ta_z_y_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 7);

            auto ta_z_z_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 8);

            auto ta_x_x_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx);

            auto ta_x_y_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 1);

            auto ta_x_z_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 2);

            auto ta_y_x_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 3);

            auto ta_y_y_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 4);

            auto ta_y_z_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 5);

            auto ta_z_x_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 6);

            auto ta_z_y_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 7);

            auto ta_z_z_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 8);

            // set up pointers to integrals

            auto ta_xx_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx);

            auto ta_xx_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 1);

            auto ta_xx_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 2);

            auto ta_xx_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 3);

            auto ta_xx_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 4);

            auto ta_xx_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 5);

            auto ta_xy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 6);

            auto ta_xy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 7);

            auto ta_xy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 8);

            auto ta_xy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 9);

            auto ta_xy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 10);

            auto ta_xy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 11);

            auto ta_xz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 12);

            auto ta_xz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 13);

            auto ta_xz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 14);

            auto ta_xz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 15);

            auto ta_xz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 16);

            auto ta_xz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 17);

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_zz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 30);

            auto ta_zz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 31);

            auto ta_zz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 32);

            auto ta_zz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 33);

            auto ta_zz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 34);

            auto ta_zz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 35);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_xx_0, ta_0_xx_1, ta_0_xy_0, \
                                         ta_0_xy_1, ta_0_xz_0, ta_0_xz_1, ta_0_yy_0, ta_0_yy_1, ta_0_yz_0, ta_0_yz_1, \
                                         ta_0_zz_0, ta_0_zz_1, ta_x_x_0, ta_x_x_1, ta_x_xx_0, ta_x_xx_1, ta_x_xy_0, \
                                         ta_x_xy_1, ta_x_xz_0, ta_x_xz_1, ta_x_y_0, ta_x_y_1, ta_x_yy_0, ta_x_yy_1, \
                                         ta_x_yz_0, ta_x_yz_1, ta_x_z_0, ta_x_z_1, ta_x_zz_0, ta_x_zz_1, ta_xx_xx_0, \
                                         ta_xx_xy_0, ta_xx_xz_0, ta_xx_yy_0, ta_xx_yz_0, ta_xx_zz_0, ta_xy_xx_0, ta_xy_xy_0, \
                                         ta_xy_xz_0, ta_xy_yy_0, ta_xy_yz_0, ta_xy_zz_0, ta_xz_xx_0, ta_xz_xy_0, ta_xz_xz_0, \
                                         ta_xz_yy_0, ta_xz_yz_0, ta_xz_zz_0, ta_y_x_0, ta_y_x_1, ta_y_xx_0, ta_y_xx_1, \
                                         ta_y_xy_0, ta_y_xy_1, ta_y_xz_0, ta_y_xz_1, ta_y_y_0, ta_y_y_1, ta_y_yy_0, \
                                         ta_y_yy_1, ta_y_yz_0, ta_y_yz_1, ta_y_z_0, ta_y_z_1, ta_y_zz_0, ta_y_zz_1, \
                                         ta_yy_xx_0, ta_yy_xy_0, ta_yy_xz_0, ta_yy_yy_0, ta_yy_yz_0, ta_yy_zz_0, ta_yz_xx_0, \
                                         ta_yz_xy_0, ta_yz_xz_0, ta_yz_yy_0, ta_yz_yz_0, ta_yz_zz_0, ta_z_x_0, ta_z_x_1, \
                                         ta_z_xx_0, ta_z_xx_1, ta_z_xy_0, ta_z_xy_1, ta_z_xz_0, ta_z_xz_1, ta_z_y_0, \
                                         ta_z_y_1, ta_z_yy_0, ta_z_yy_1, ta_z_yz_0, ta_z_yz_1, ta_z_z_0, ta_z_z_1, \
                                         ta_z_zz_0, ta_z_zz_1, ta_zz_xx_0, ta_zz_xy_0, ta_zz_xz_0, ta_zz_yy_0, ta_zz_yz_0, \
                                         ta_zz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xx_xx_0[j] = pa_x[j] * ta_x_xx_0[j] - pc_x[j] * ta_x_xx_1[j] + 0.5 * fl1_fx * ta_0_xx_0[j] - 0.5 * fl1_fx * ta_0_xx_1[j] +
                                fl1_fx * ta_x_x_0[j] - fl1_fx * ta_x_x_1[j];

                ta_xx_xy_0[j] = pa_x[j] * ta_x_xy_0[j] - pc_x[j] * ta_x_xy_1[j] + 0.5 * fl1_fx * ta_0_xy_0[j] - 0.5 * fl1_fx * ta_0_xy_1[j] +
                                0.5 * fl1_fx * ta_x_y_0[j] - 0.5 * fl1_fx * ta_x_y_1[j];

                ta_xx_xz_0[j] = pa_x[j] * ta_x_xz_0[j] - pc_x[j] * ta_x_xz_1[j] + 0.5 * fl1_fx * ta_0_xz_0[j] - 0.5 * fl1_fx * ta_0_xz_1[j] +
                                0.5 * fl1_fx * ta_x_z_0[j] - 0.5 * fl1_fx * ta_x_z_1[j];

                ta_xx_yy_0[j] = pa_x[j] * ta_x_yy_0[j] - pc_x[j] * ta_x_yy_1[j] + 0.5 * fl1_fx * ta_0_yy_0[j] - 0.5 * fl1_fx * ta_0_yy_1[j];

                ta_xx_yz_0[j] = pa_x[j] * ta_x_yz_0[j] - pc_x[j] * ta_x_yz_1[j] + 0.5 * fl1_fx * ta_0_yz_0[j] - 0.5 * fl1_fx * ta_0_yz_1[j];

                ta_xx_zz_0[j] = pa_x[j] * ta_x_zz_0[j] - pc_x[j] * ta_x_zz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j];

                ta_xy_xx_0[j] = pa_x[j] * ta_y_xx_0[j] - pc_x[j] * ta_y_xx_1[j] + fl1_fx * ta_y_x_0[j] - fl1_fx * ta_y_x_1[j];

                ta_xy_xy_0[j] = pa_x[j] * ta_y_xy_0[j] - pc_x[j] * ta_y_xy_1[j] + 0.5 * fl1_fx * ta_y_y_0[j] - 0.5 * fl1_fx * ta_y_y_1[j];

                ta_xy_xz_0[j] = pa_x[j] * ta_y_xz_0[j] - pc_x[j] * ta_y_xz_1[j] + 0.5 * fl1_fx * ta_y_z_0[j] - 0.5 * fl1_fx * ta_y_z_1[j];

                ta_xy_yy_0[j] = pa_x[j] * ta_y_yy_0[j] - pc_x[j] * ta_y_yy_1[j];

                ta_xy_yz_0[j] = pa_x[j] * ta_y_yz_0[j] - pc_x[j] * ta_y_yz_1[j];

                ta_xy_zz_0[j] = pa_x[j] * ta_y_zz_0[j] - pc_x[j] * ta_y_zz_1[j];

                ta_xz_xx_0[j] = pa_x[j] * ta_z_xx_0[j] - pc_x[j] * ta_z_xx_1[j] + fl1_fx * ta_z_x_0[j] - fl1_fx * ta_z_x_1[j];

                ta_xz_xy_0[j] = pa_x[j] * ta_z_xy_0[j] - pc_x[j] * ta_z_xy_1[j] + 0.5 * fl1_fx * ta_z_y_0[j] - 0.5 * fl1_fx * ta_z_y_1[j];

                ta_xz_xz_0[j] = pa_x[j] * ta_z_xz_0[j] - pc_x[j] * ta_z_xz_1[j] + 0.5 * fl1_fx * ta_z_z_0[j] - 0.5 * fl1_fx * ta_z_z_1[j];

                ta_xz_yy_0[j] = pa_x[j] * ta_z_yy_0[j] - pc_x[j] * ta_z_yy_1[j];

                ta_xz_yz_0[j] = pa_x[j] * ta_z_yz_0[j] - pc_x[j] * ta_z_yz_1[j];

                ta_xz_zz_0[j] = pa_x[j] * ta_z_zz_0[j] - pc_x[j] * ta_z_zz_1[j];

                ta_yy_xx_0[j] = pa_y[j] * ta_y_xx_0[j] - pc_y[j] * ta_y_xx_1[j] + 0.5 * fl1_fx * ta_0_xx_0[j] - 0.5 * fl1_fx * ta_0_xx_1[j];

                ta_yy_xy_0[j] = pa_y[j] * ta_y_xy_0[j] - pc_y[j] * ta_y_xy_1[j] + 0.5 * fl1_fx * ta_0_xy_0[j] - 0.5 * fl1_fx * ta_0_xy_1[j] +
                                0.5 * fl1_fx * ta_y_x_0[j] - 0.5 * fl1_fx * ta_y_x_1[j];

                ta_yy_xz_0[j] = pa_y[j] * ta_y_xz_0[j] - pc_y[j] * ta_y_xz_1[j] + 0.5 * fl1_fx * ta_0_xz_0[j] - 0.5 * fl1_fx * ta_0_xz_1[j];

                ta_yy_yy_0[j] = pa_y[j] * ta_y_yy_0[j] - pc_y[j] * ta_y_yy_1[j] + 0.5 * fl1_fx * ta_0_yy_0[j] - 0.5 * fl1_fx * ta_0_yy_1[j] +
                                fl1_fx * ta_y_y_0[j] - fl1_fx * ta_y_y_1[j];

                ta_yy_yz_0[j] = pa_y[j] * ta_y_yz_0[j] - pc_y[j] * ta_y_yz_1[j] + 0.5 * fl1_fx * ta_0_yz_0[j] - 0.5 * fl1_fx * ta_0_yz_1[j] +
                                0.5 * fl1_fx * ta_y_z_0[j] - 0.5 * fl1_fx * ta_y_z_1[j];

                ta_yy_zz_0[j] = pa_y[j] * ta_y_zz_0[j] - pc_y[j] * ta_y_zz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j];

                ta_yz_xx_0[j] = pa_y[j] * ta_z_xx_0[j] - pc_y[j] * ta_z_xx_1[j];

                ta_yz_xy_0[j] = pa_y[j] * ta_z_xy_0[j] - pc_y[j] * ta_z_xy_1[j] + 0.5 * fl1_fx * ta_z_x_0[j] - 0.5 * fl1_fx * ta_z_x_1[j];

                ta_yz_xz_0[j] = pa_y[j] * ta_z_xz_0[j] - pc_y[j] * ta_z_xz_1[j];

                ta_yz_yy_0[j] = pa_y[j] * ta_z_yy_0[j] - pc_y[j] * ta_z_yy_1[j] + fl1_fx * ta_z_y_0[j] - fl1_fx * ta_z_y_1[j];

                ta_yz_yz_0[j] = pa_y[j] * ta_z_yz_0[j] - pc_y[j] * ta_z_yz_1[j] + 0.5 * fl1_fx * ta_z_z_0[j] - 0.5 * fl1_fx * ta_z_z_1[j];

                ta_yz_zz_0[j] = pa_y[j] * ta_z_zz_0[j] - pc_y[j] * ta_z_zz_1[j];

                ta_zz_xx_0[j] = pa_z[j] * ta_z_xx_0[j] - pc_z[j] * ta_z_xx_1[j] + 0.5 * fl1_fx * ta_0_xx_0[j] - 0.5 * fl1_fx * ta_0_xx_1[j];

                ta_zz_xy_0[j] = pa_z[j] * ta_z_xy_0[j] - pc_z[j] * ta_z_xy_1[j] + 0.5 * fl1_fx * ta_0_xy_0[j] - 0.5 * fl1_fx * ta_0_xy_1[j];

                ta_zz_xz_0[j] = pa_z[j] * ta_z_xz_0[j] - pc_z[j] * ta_z_xz_1[j] + 0.5 * fl1_fx * ta_0_xz_0[j] - 0.5 * fl1_fx * ta_0_xz_1[j] +
                                0.5 * fl1_fx * ta_z_x_0[j] - 0.5 * fl1_fx * ta_z_x_1[j];

                ta_zz_yy_0[j] = pa_z[j] * ta_z_yy_0[j] - pc_z[j] * ta_z_yy_1[j] + 0.5 * fl1_fx * ta_0_yy_0[j] - 0.5 * fl1_fx * ta_0_yy_1[j];

                ta_zz_yz_0[j] = pa_z[j] * ta_z_yz_0[j] - pc_z[j] * ta_z_yz_1[j] + 0.5 * fl1_fx * ta_0_yz_0[j] - 0.5 * fl1_fx * ta_0_yz_1[j] +
                                0.5 * fl1_fx * ta_z_y_0[j] - 0.5 * fl1_fx * ta_z_y_1[j];

                ta_zz_zz_0[j] = pa_z[j] * ta_z_zz_0[j] - pc_z[j] * ta_z_zz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j] +
                                fl1_fx * ta_z_z_0[j] - fl1_fx * ta_z_z_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForDF(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForDF_0_30(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDF_30_60(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForDF_0_30(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (0,30)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_x_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx);

            auto ta_x_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 1);

            auto ta_x_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 2);

            auto ta_x_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 3);

            auto ta_x_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 4);

            auto ta_x_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 5);

            auto ta_x_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 6);

            auto ta_x_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 7);

            auto ta_x_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 8);

            auto ta_x_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 9);

            auto ta_y_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 10);

            auto ta_y_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 11);

            auto ta_y_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 12);

            auto ta_y_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 13);

            auto ta_y_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 14);

            auto ta_y_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 15);

            auto ta_y_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 16);

            auto ta_y_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 17);

            auto ta_y_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 18);

            auto ta_y_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 19);

            auto ta_z_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 20);

            auto ta_z_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 21);

            auto ta_z_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 22);

            auto ta_z_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 23);

            auto ta_z_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 24);

            auto ta_z_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 25);

            auto ta_z_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 26);

            auto ta_z_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 27);

            auto ta_z_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 28);

            auto ta_z_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 29);

            auto ta_x_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx);

            auto ta_x_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 1);

            auto ta_x_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 2);

            auto ta_x_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 3);

            auto ta_x_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 4);

            auto ta_x_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 5);

            auto ta_x_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 6);

            auto ta_x_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 7);

            auto ta_x_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 8);

            auto ta_x_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 9);

            auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10);

            auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11);

            auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12);

            auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13);

            auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14);

            auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15);

            auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16);

            auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17);

            auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18);

            auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19);

            auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20);

            auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21);

            auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22);

            auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23);

            auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24);

            auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25);

            auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26);

            auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27);

            auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28);

            auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29);

            auto ta_0_xxx_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx);

            auto ta_0_xxy_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 1);

            auto ta_0_xxz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 2);

            auto ta_0_xyy_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 3);

            auto ta_0_xyz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 4);

            auto ta_0_xzz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 5);

            auto ta_0_yyy_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 6);

            auto ta_0_yyz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 7);

            auto ta_0_yzz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 8);

            auto ta_0_zzz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 9);

            auto ta_0_xxx_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx);

            auto ta_0_xxy_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 1);

            auto ta_0_xxz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 2);

            auto ta_0_xyy_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 3);

            auto ta_0_xyz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 4);

            auto ta_0_xzz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 5);

            auto ta_0_yyy_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 6);

            auto ta_0_yyz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 7);

            auto ta_0_yzz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 8);

            auto ta_0_zzz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 9);

            auto ta_x_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx);

            auto ta_x_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 1);

            auto ta_x_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 2);

            auto ta_x_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 3);

            auto ta_x_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 4);

            auto ta_x_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 5);

            auto ta_y_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 6);

            auto ta_y_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 7);

            auto ta_y_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 8);

            auto ta_y_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 9);

            auto ta_y_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 10);

            auto ta_y_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 11);

            auto ta_z_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 12);

            auto ta_z_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 13);

            auto ta_z_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 14);

            auto ta_z_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 15);

            auto ta_z_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 16);

            auto ta_z_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 17);

            auto ta_x_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx);

            auto ta_x_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 1);

            auto ta_x_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 2);

            auto ta_x_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 3);

            auto ta_x_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 4);

            auto ta_x_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 5);

            auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6);

            auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7);

            auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8);

            auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9);

            auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10);

            auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11);

            auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12);

            auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13);

            auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14);

            auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15);

            auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16);

            auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17);

            // set up pointers to integrals

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

            // Batch of Integrals (0,30)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_0_xxx_0, ta_0_xxx_1, ta_0_xxy_0, ta_0_xxy_1, ta_0_xxz_0, \
                                         ta_0_xxz_1, ta_0_xyy_0, ta_0_xyy_1, ta_0_xyz_0, ta_0_xyz_1, ta_0_xzz_0, ta_0_xzz_1, \
                                         ta_0_yyy_0, ta_0_yyy_1, ta_0_yyz_0, ta_0_yyz_1, ta_0_yzz_0, ta_0_yzz_1, ta_0_zzz_0, \
                                         ta_0_zzz_1, ta_x_xx_0, ta_x_xx_1, ta_x_xxx_0, ta_x_xxx_1, ta_x_xxy_0, ta_x_xxy_1, \
                                         ta_x_xxz_0, ta_x_xxz_1, ta_x_xy_0, ta_x_xy_1, ta_x_xyy_0, ta_x_xyy_1, ta_x_xyz_0, \
                                         ta_x_xyz_1, ta_x_xz_0, ta_x_xz_1, ta_x_xzz_0, ta_x_xzz_1, ta_x_yy_0, ta_x_yy_1, \
                                         ta_x_yyy_0, ta_x_yyy_1, ta_x_yyz_0, ta_x_yyz_1, ta_x_yz_0, ta_x_yz_1, ta_x_yzz_0, \
                                         ta_x_yzz_1, ta_x_zz_0, ta_x_zz_1, ta_x_zzz_0, ta_x_zzz_1, ta_xx_xxx_0, ta_xx_xxy_0, \
                                         ta_xx_xxz_0, ta_xx_xyy_0, ta_xx_xyz_0, ta_xx_xzz_0, ta_xx_yyy_0, ta_xx_yyz_0, \
                                         ta_xx_yzz_0, ta_xx_zzz_0, ta_xy_xxx_0, ta_xy_xxy_0, ta_xy_xxz_0, ta_xy_xyy_0, \
                                         ta_xy_xyz_0, ta_xy_xzz_0, ta_xy_yyy_0, ta_xy_yyz_0, ta_xy_yzz_0, ta_xy_zzz_0, \
                                         ta_xz_xxx_0, ta_xz_xxy_0, ta_xz_xxz_0, ta_xz_xyy_0, ta_xz_xyz_0, ta_xz_xzz_0, \
                                         ta_xz_yyy_0, ta_xz_yyz_0, ta_xz_yzz_0, ta_xz_zzz_0, ta_y_xx_0, ta_y_xx_1, \
                                         ta_y_xxx_0, ta_y_xxx_1, ta_y_xxy_0, ta_y_xxy_1, ta_y_xxz_0, ta_y_xxz_1, ta_y_xy_0, \
                                         ta_y_xy_1, ta_y_xyy_0, ta_y_xyy_1, ta_y_xyz_0, ta_y_xyz_1, ta_y_xz_0, ta_y_xz_1, \
                                         ta_y_xzz_0, ta_y_xzz_1, ta_y_yy_0, ta_y_yy_1, ta_y_yyy_0, ta_y_yyy_1, ta_y_yyz_0, \
                                         ta_y_yyz_1, ta_y_yz_0, ta_y_yz_1, ta_y_yzz_0, ta_y_yzz_1, ta_y_zz_0, ta_y_zz_1, \
                                         ta_y_zzz_0, ta_y_zzz_1, ta_z_xx_0, ta_z_xx_1, ta_z_xxx_0, ta_z_xxx_1, ta_z_xxy_0, \
                                         ta_z_xxy_1, ta_z_xxz_0, ta_z_xxz_1, ta_z_xy_0, ta_z_xy_1, ta_z_xyy_0, ta_z_xyy_1, \
                                         ta_z_xyz_0, ta_z_xyz_1, ta_z_xz_0, ta_z_xz_1, ta_z_xzz_0, ta_z_xzz_1, ta_z_yy_0, \
                                         ta_z_yy_1, ta_z_yyy_0, ta_z_yyy_1, ta_z_yyz_0, ta_z_yyz_1, ta_z_yz_0, ta_z_yz_1, \
                                         ta_z_yzz_0, ta_z_yzz_1, ta_z_zz_0, ta_z_zz_1, ta_z_zzz_0, ta_z_zzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xx_xxx_0[j] = pa_x[j] * ta_x_xxx_0[j] - pc_x[j] * ta_x_xxx_1[j] + 0.5 * fl1_fx * ta_0_xxx_0[j] - 0.5 * fl1_fx * ta_0_xxx_1[j] +
                                 1.5 * fl1_fx * ta_x_xx_0[j] - 1.5 * fl1_fx * ta_x_xx_1[j];

                ta_xx_xxy_0[j] = pa_x[j] * ta_x_xxy_0[j] - pc_x[j] * ta_x_xxy_1[j] + 0.5 * fl1_fx * ta_0_xxy_0[j] - 0.5 * fl1_fx * ta_0_xxy_1[j] +
                                 fl1_fx * ta_x_xy_0[j] - fl1_fx * ta_x_xy_1[j];

                ta_xx_xxz_0[j] = pa_x[j] * ta_x_xxz_0[j] - pc_x[j] * ta_x_xxz_1[j] + 0.5 * fl1_fx * ta_0_xxz_0[j] - 0.5 * fl1_fx * ta_0_xxz_1[j] +
                                 fl1_fx * ta_x_xz_0[j] - fl1_fx * ta_x_xz_1[j];

                ta_xx_xyy_0[j] = pa_x[j] * ta_x_xyy_0[j] - pc_x[j] * ta_x_xyy_1[j] + 0.5 * fl1_fx * ta_0_xyy_0[j] - 0.5 * fl1_fx * ta_0_xyy_1[j] +
                                 0.5 * fl1_fx * ta_x_yy_0[j] - 0.5 * fl1_fx * ta_x_yy_1[j];

                ta_xx_xyz_0[j] = pa_x[j] * ta_x_xyz_0[j] - pc_x[j] * ta_x_xyz_1[j] + 0.5 * fl1_fx * ta_0_xyz_0[j] - 0.5 * fl1_fx * ta_0_xyz_1[j] +
                                 0.5 * fl1_fx * ta_x_yz_0[j] - 0.5 * fl1_fx * ta_x_yz_1[j];

                ta_xx_xzz_0[j] = pa_x[j] * ta_x_xzz_0[j] - pc_x[j] * ta_x_xzz_1[j] + 0.5 * fl1_fx * ta_0_xzz_0[j] - 0.5 * fl1_fx * ta_0_xzz_1[j] +
                                 0.5 * fl1_fx * ta_x_zz_0[j] - 0.5 * fl1_fx * ta_x_zz_1[j];

                ta_xx_yyy_0[j] = pa_x[j] * ta_x_yyy_0[j] - pc_x[j] * ta_x_yyy_1[j] + 0.5 * fl1_fx * ta_0_yyy_0[j] - 0.5 * fl1_fx * ta_0_yyy_1[j];

                ta_xx_yyz_0[j] = pa_x[j] * ta_x_yyz_0[j] - pc_x[j] * ta_x_yyz_1[j] + 0.5 * fl1_fx * ta_0_yyz_0[j] - 0.5 * fl1_fx * ta_0_yyz_1[j];

                ta_xx_yzz_0[j] = pa_x[j] * ta_x_yzz_0[j] - pc_x[j] * ta_x_yzz_1[j] + 0.5 * fl1_fx * ta_0_yzz_0[j] - 0.5 * fl1_fx * ta_0_yzz_1[j];

                ta_xx_zzz_0[j] = pa_x[j] * ta_x_zzz_0[j] - pc_x[j] * ta_x_zzz_1[j] + 0.5 * fl1_fx * ta_0_zzz_0[j] - 0.5 * fl1_fx * ta_0_zzz_1[j];

                ta_xy_xxx_0[j] = pa_x[j] * ta_y_xxx_0[j] - pc_x[j] * ta_y_xxx_1[j] + 1.5 * fl1_fx * ta_y_xx_0[j] - 1.5 * fl1_fx * ta_y_xx_1[j];

                ta_xy_xxy_0[j] = pa_x[j] * ta_y_xxy_0[j] - pc_x[j] * ta_y_xxy_1[j] + fl1_fx * ta_y_xy_0[j] - fl1_fx * ta_y_xy_1[j];

                ta_xy_xxz_0[j] = pa_x[j] * ta_y_xxz_0[j] - pc_x[j] * ta_y_xxz_1[j] + fl1_fx * ta_y_xz_0[j] - fl1_fx * ta_y_xz_1[j];

                ta_xy_xyy_0[j] = pa_x[j] * ta_y_xyy_0[j] - pc_x[j] * ta_y_xyy_1[j] + 0.5 * fl1_fx * ta_y_yy_0[j] - 0.5 * fl1_fx * ta_y_yy_1[j];

                ta_xy_xyz_0[j] = pa_x[j] * ta_y_xyz_0[j] - pc_x[j] * ta_y_xyz_1[j] + 0.5 * fl1_fx * ta_y_yz_0[j] - 0.5 * fl1_fx * ta_y_yz_1[j];

                ta_xy_xzz_0[j] = pa_x[j] * ta_y_xzz_0[j] - pc_x[j] * ta_y_xzz_1[j] + 0.5 * fl1_fx * ta_y_zz_0[j] - 0.5 * fl1_fx * ta_y_zz_1[j];

                ta_xy_yyy_0[j] = pa_x[j] * ta_y_yyy_0[j] - pc_x[j] * ta_y_yyy_1[j];

                ta_xy_yyz_0[j] = pa_x[j] * ta_y_yyz_0[j] - pc_x[j] * ta_y_yyz_1[j];

                ta_xy_yzz_0[j] = pa_x[j] * ta_y_yzz_0[j] - pc_x[j] * ta_y_yzz_1[j];

                ta_xy_zzz_0[j] = pa_x[j] * ta_y_zzz_0[j] - pc_x[j] * ta_y_zzz_1[j];

                ta_xz_xxx_0[j] = pa_x[j] * ta_z_xxx_0[j] - pc_x[j] * ta_z_xxx_1[j] + 1.5 * fl1_fx * ta_z_xx_0[j] - 1.5 * fl1_fx * ta_z_xx_1[j];

                ta_xz_xxy_0[j] = pa_x[j] * ta_z_xxy_0[j] - pc_x[j] * ta_z_xxy_1[j] + fl1_fx * ta_z_xy_0[j] - fl1_fx * ta_z_xy_1[j];

                ta_xz_xxz_0[j] = pa_x[j] * ta_z_xxz_0[j] - pc_x[j] * ta_z_xxz_1[j] + fl1_fx * ta_z_xz_0[j] - fl1_fx * ta_z_xz_1[j];

                ta_xz_xyy_0[j] = pa_x[j] * ta_z_xyy_0[j] - pc_x[j] * ta_z_xyy_1[j] + 0.5 * fl1_fx * ta_z_yy_0[j] - 0.5 * fl1_fx * ta_z_yy_1[j];

                ta_xz_xyz_0[j] = pa_x[j] * ta_z_xyz_0[j] - pc_x[j] * ta_z_xyz_1[j] + 0.5 * fl1_fx * ta_z_yz_0[j] - 0.5 * fl1_fx * ta_z_yz_1[j];

                ta_xz_xzz_0[j] = pa_x[j] * ta_z_xzz_0[j] - pc_x[j] * ta_z_xzz_1[j] + 0.5 * fl1_fx * ta_z_zz_0[j] - 0.5 * fl1_fx * ta_z_zz_1[j];

                ta_xz_yyy_0[j] = pa_x[j] * ta_z_yyy_0[j] - pc_x[j] * ta_z_yyy_1[j];

                ta_xz_yyz_0[j] = pa_x[j] * ta_z_yyz_0[j] - pc_x[j] * ta_z_yyz_1[j];

                ta_xz_yzz_0[j] = pa_x[j] * ta_z_yzz_0[j] - pc_x[j] * ta_z_yzz_1[j];

                ta_xz_zzz_0[j] = pa_x[j] * ta_z_zzz_0[j] - pc_x[j] * ta_z_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForDF_30_60(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (30,60)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_y_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 10);

            auto ta_y_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 11);

            auto ta_y_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 12);

            auto ta_y_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 13);

            auto ta_y_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 14);

            auto ta_y_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 15);

            auto ta_y_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 16);

            auto ta_y_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 17);

            auto ta_y_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 18);

            auto ta_y_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 19);

            auto ta_z_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 20);

            auto ta_z_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 21);

            auto ta_z_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 22);

            auto ta_z_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 23);

            auto ta_z_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 24);

            auto ta_z_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 25);

            auto ta_z_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 26);

            auto ta_z_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 27);

            auto ta_z_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 28);

            auto ta_z_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 29);

            auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10);

            auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11);

            auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12);

            auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13);

            auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14);

            auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15);

            auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16);

            auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17);

            auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18);

            auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19);

            auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20);

            auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21);

            auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22);

            auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23);

            auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24);

            auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25);

            auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26);

            auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27);

            auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28);

            auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29);

            auto ta_0_xxx_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx);

            auto ta_0_xxy_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 1);

            auto ta_0_xxz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 2);

            auto ta_0_xyy_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 3);

            auto ta_0_xyz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 4);

            auto ta_0_xzz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 5);

            auto ta_0_yyy_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 6);

            auto ta_0_yyz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 7);

            auto ta_0_yzz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 8);

            auto ta_0_zzz_0 = primBuffer.data(pidx_a_0_3_m0 + 10 * idx + 9);

            auto ta_0_xxx_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx);

            auto ta_0_xxy_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 1);

            auto ta_0_xxz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 2);

            auto ta_0_xyy_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 3);

            auto ta_0_xyz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 4);

            auto ta_0_xzz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 5);

            auto ta_0_yyy_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 6);

            auto ta_0_yyz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 7);

            auto ta_0_yzz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 8);

            auto ta_0_zzz_1 = primBuffer.data(pidx_a_0_3_m1 + 10 * idx + 9);

            auto ta_y_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 6);

            auto ta_y_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 7);

            auto ta_y_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 8);

            auto ta_y_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 9);

            auto ta_y_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 10);

            auto ta_y_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 11);

            auto ta_z_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 12);

            auto ta_z_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 13);

            auto ta_z_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 14);

            auto ta_z_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 15);

            auto ta_z_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 16);

            auto ta_z_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 17);

            auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6);

            auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7);

            auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8);

            auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9);

            auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10);

            auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11);

            auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12);

            auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13);

            auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14);

            auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15);

            auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16);

            auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17);

            // set up pointers to integrals

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

            // Batch of Integrals (30,60)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_0_xxx_0, ta_0_xxx_1, ta_0_xxy_0, ta_0_xxy_1, \
                                         ta_0_xxz_0, ta_0_xxz_1, ta_0_xyy_0, ta_0_xyy_1, ta_0_xyz_0, ta_0_xyz_1, ta_0_xzz_0, \
                                         ta_0_xzz_1, ta_0_yyy_0, ta_0_yyy_1, ta_0_yyz_0, ta_0_yyz_1, ta_0_yzz_0, ta_0_yzz_1, \
                                         ta_0_zzz_0, ta_0_zzz_1, ta_y_xx_0, ta_y_xx_1, ta_y_xxx_0, ta_y_xxx_1, ta_y_xxy_0, \
                                         ta_y_xxy_1, ta_y_xxz_0, ta_y_xxz_1, ta_y_xy_0, ta_y_xy_1, ta_y_xyy_0, ta_y_xyy_1, \
                                         ta_y_xyz_0, ta_y_xyz_1, ta_y_xz_0, ta_y_xz_1, ta_y_xzz_0, ta_y_xzz_1, ta_y_yy_0, \
                                         ta_y_yy_1, ta_y_yyy_0, ta_y_yyy_1, ta_y_yyz_0, ta_y_yyz_1, ta_y_yz_0, ta_y_yz_1, \
                                         ta_y_yzz_0, ta_y_yzz_1, ta_y_zz_0, ta_y_zz_1, ta_y_zzz_0, ta_y_zzz_1, ta_yy_xxx_0, \
                                         ta_yy_xxy_0, ta_yy_xxz_0, ta_yy_xyy_0, ta_yy_xyz_0, ta_yy_xzz_0, ta_yy_yyy_0, \
                                         ta_yy_yyz_0, ta_yy_yzz_0, ta_yy_zzz_0, ta_yz_xxx_0, ta_yz_xxy_0, ta_yz_xxz_0, \
                                         ta_yz_xyy_0, ta_yz_xyz_0, ta_yz_xzz_0, ta_yz_yyy_0, ta_yz_yyz_0, ta_yz_yzz_0, \
                                         ta_yz_zzz_0, ta_z_xx_0, ta_z_xx_1, ta_z_xxx_0, ta_z_xxx_1, ta_z_xxy_0, ta_z_xxy_1, \
                                         ta_z_xxz_0, ta_z_xxz_1, ta_z_xy_0, ta_z_xy_1, ta_z_xyy_0, ta_z_xyy_1, ta_z_xyz_0, \
                                         ta_z_xyz_1, ta_z_xz_0, ta_z_xz_1, ta_z_xzz_0, ta_z_xzz_1, ta_z_yy_0, ta_z_yy_1, \
                                         ta_z_yyy_0, ta_z_yyy_1, ta_z_yyz_0, ta_z_yyz_1, ta_z_yz_0, ta_z_yz_1, ta_z_yzz_0, \
                                         ta_z_yzz_1, ta_z_zz_0, ta_z_zz_1, ta_z_zzz_0, ta_z_zzz_1, ta_zz_xxx_0, ta_zz_xxy_0, \
                                         ta_zz_xxz_0, ta_zz_xyy_0, ta_zz_xyz_0, ta_zz_xzz_0, ta_zz_yyy_0, ta_zz_yyz_0, \
                                         ta_zz_yzz_0, ta_zz_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_yy_xxx_0[j] = pa_y[j] * ta_y_xxx_0[j] - pc_y[j] * ta_y_xxx_1[j] + 0.5 * fl1_fx * ta_0_xxx_0[j] - 0.5 * fl1_fx * ta_0_xxx_1[j];

                ta_yy_xxy_0[j] = pa_y[j] * ta_y_xxy_0[j] - pc_y[j] * ta_y_xxy_1[j] + 0.5 * fl1_fx * ta_0_xxy_0[j] - 0.5 * fl1_fx * ta_0_xxy_1[j] +
                                 0.5 * fl1_fx * ta_y_xx_0[j] - 0.5 * fl1_fx * ta_y_xx_1[j];

                ta_yy_xxz_0[j] = pa_y[j] * ta_y_xxz_0[j] - pc_y[j] * ta_y_xxz_1[j] + 0.5 * fl1_fx * ta_0_xxz_0[j] - 0.5 * fl1_fx * ta_0_xxz_1[j];

                ta_yy_xyy_0[j] = pa_y[j] * ta_y_xyy_0[j] - pc_y[j] * ta_y_xyy_1[j] + 0.5 * fl1_fx * ta_0_xyy_0[j] - 0.5 * fl1_fx * ta_0_xyy_1[j] +
                                 fl1_fx * ta_y_xy_0[j] - fl1_fx * ta_y_xy_1[j];

                ta_yy_xyz_0[j] = pa_y[j] * ta_y_xyz_0[j] - pc_y[j] * ta_y_xyz_1[j] + 0.5 * fl1_fx * ta_0_xyz_0[j] - 0.5 * fl1_fx * ta_0_xyz_1[j] +
                                 0.5 * fl1_fx * ta_y_xz_0[j] - 0.5 * fl1_fx * ta_y_xz_1[j];

                ta_yy_xzz_0[j] = pa_y[j] * ta_y_xzz_0[j] - pc_y[j] * ta_y_xzz_1[j] + 0.5 * fl1_fx * ta_0_xzz_0[j] - 0.5 * fl1_fx * ta_0_xzz_1[j];

                ta_yy_yyy_0[j] = pa_y[j] * ta_y_yyy_0[j] - pc_y[j] * ta_y_yyy_1[j] + 0.5 * fl1_fx * ta_0_yyy_0[j] - 0.5 * fl1_fx * ta_0_yyy_1[j] +
                                 1.5 * fl1_fx * ta_y_yy_0[j] - 1.5 * fl1_fx * ta_y_yy_1[j];

                ta_yy_yyz_0[j] = pa_y[j] * ta_y_yyz_0[j] - pc_y[j] * ta_y_yyz_1[j] + 0.5 * fl1_fx * ta_0_yyz_0[j] - 0.5 * fl1_fx * ta_0_yyz_1[j] +
                                 fl1_fx * ta_y_yz_0[j] - fl1_fx * ta_y_yz_1[j];

                ta_yy_yzz_0[j] = pa_y[j] * ta_y_yzz_0[j] - pc_y[j] * ta_y_yzz_1[j] + 0.5 * fl1_fx * ta_0_yzz_0[j] - 0.5 * fl1_fx * ta_0_yzz_1[j] +
                                 0.5 * fl1_fx * ta_y_zz_0[j] - 0.5 * fl1_fx * ta_y_zz_1[j];

                ta_yy_zzz_0[j] = pa_y[j] * ta_y_zzz_0[j] - pc_y[j] * ta_y_zzz_1[j] + 0.5 * fl1_fx * ta_0_zzz_0[j] - 0.5 * fl1_fx * ta_0_zzz_1[j];

                ta_yz_xxx_0[j] = pa_y[j] * ta_z_xxx_0[j] - pc_y[j] * ta_z_xxx_1[j];

                ta_yz_xxy_0[j] = pa_y[j] * ta_z_xxy_0[j] - pc_y[j] * ta_z_xxy_1[j] + 0.5 * fl1_fx * ta_z_xx_0[j] - 0.5 * fl1_fx * ta_z_xx_1[j];

                ta_yz_xxz_0[j] = pa_y[j] * ta_z_xxz_0[j] - pc_y[j] * ta_z_xxz_1[j];

                ta_yz_xyy_0[j] = pa_y[j] * ta_z_xyy_0[j] - pc_y[j] * ta_z_xyy_1[j] + fl1_fx * ta_z_xy_0[j] - fl1_fx * ta_z_xy_1[j];

                ta_yz_xyz_0[j] = pa_y[j] * ta_z_xyz_0[j] - pc_y[j] * ta_z_xyz_1[j] + 0.5 * fl1_fx * ta_z_xz_0[j] - 0.5 * fl1_fx * ta_z_xz_1[j];

                ta_yz_xzz_0[j] = pa_y[j] * ta_z_xzz_0[j] - pc_y[j] * ta_z_xzz_1[j];

                ta_yz_yyy_0[j] = pa_y[j] * ta_z_yyy_0[j] - pc_y[j] * ta_z_yyy_1[j] + 1.5 * fl1_fx * ta_z_yy_0[j] - 1.5 * fl1_fx * ta_z_yy_1[j];

                ta_yz_yyz_0[j] = pa_y[j] * ta_z_yyz_0[j] - pc_y[j] * ta_z_yyz_1[j] + fl1_fx * ta_z_yz_0[j] - fl1_fx * ta_z_yz_1[j];

                ta_yz_yzz_0[j] = pa_y[j] * ta_z_yzz_0[j] - pc_y[j] * ta_z_yzz_1[j] + 0.5 * fl1_fx * ta_z_zz_0[j] - 0.5 * fl1_fx * ta_z_zz_1[j];

                ta_yz_zzz_0[j] = pa_y[j] * ta_z_zzz_0[j] - pc_y[j] * ta_z_zzz_1[j];

                ta_zz_xxx_0[j] = pa_z[j] * ta_z_xxx_0[j] - pc_z[j] * ta_z_xxx_1[j] + 0.5 * fl1_fx * ta_0_xxx_0[j] - 0.5 * fl1_fx * ta_0_xxx_1[j];

                ta_zz_xxy_0[j] = pa_z[j] * ta_z_xxy_0[j] - pc_z[j] * ta_z_xxy_1[j] + 0.5 * fl1_fx * ta_0_xxy_0[j] - 0.5 * fl1_fx * ta_0_xxy_1[j];

                ta_zz_xxz_0[j] = pa_z[j] * ta_z_xxz_0[j] - pc_z[j] * ta_z_xxz_1[j] + 0.5 * fl1_fx * ta_0_xxz_0[j] - 0.5 * fl1_fx * ta_0_xxz_1[j] +
                                 0.5 * fl1_fx * ta_z_xx_0[j] - 0.5 * fl1_fx * ta_z_xx_1[j];

                ta_zz_xyy_0[j] = pa_z[j] * ta_z_xyy_0[j] - pc_z[j] * ta_z_xyy_1[j] + 0.5 * fl1_fx * ta_0_xyy_0[j] - 0.5 * fl1_fx * ta_0_xyy_1[j];

                ta_zz_xyz_0[j] = pa_z[j] * ta_z_xyz_0[j] - pc_z[j] * ta_z_xyz_1[j] + 0.5 * fl1_fx * ta_0_xyz_0[j] - 0.5 * fl1_fx * ta_0_xyz_1[j] +
                                 0.5 * fl1_fx * ta_z_xy_0[j] - 0.5 * fl1_fx * ta_z_xy_1[j];

                ta_zz_xzz_0[j] = pa_z[j] * ta_z_xzz_0[j] - pc_z[j] * ta_z_xzz_1[j] + 0.5 * fl1_fx * ta_0_xzz_0[j] - 0.5 * fl1_fx * ta_0_xzz_1[j] +
                                 fl1_fx * ta_z_xz_0[j] - fl1_fx * ta_z_xz_1[j];

                ta_zz_yyy_0[j] = pa_z[j] * ta_z_yyy_0[j] - pc_z[j] * ta_z_yyy_1[j] + 0.5 * fl1_fx * ta_0_yyy_0[j] - 0.5 * fl1_fx * ta_0_yyy_1[j];

                ta_zz_yyz_0[j] = pa_z[j] * ta_z_yyz_0[j] - pc_z[j] * ta_z_yyz_1[j] + 0.5 * fl1_fx * ta_0_yyz_0[j] - 0.5 * fl1_fx * ta_0_yyz_1[j] +
                                 0.5 * fl1_fx * ta_z_yy_0[j] - 0.5 * fl1_fx * ta_z_yy_1[j];

                ta_zz_yzz_0[j] = pa_z[j] * ta_z_yzz_0[j] - pc_z[j] * ta_z_yzz_1[j] + 0.5 * fl1_fx * ta_0_yzz_0[j] - 0.5 * fl1_fx * ta_0_yzz_1[j] +
                                 fl1_fx * ta_z_yz_0[j] - fl1_fx * ta_z_yz_1[j];

                ta_zz_zzz_0[j] = pa_z[j] * ta_z_zzz_0[j] - pc_z[j] * ta_z_zzz_1[j] + 0.5 * fl1_fx * ta_0_zzz_0[j] - 0.5 * fl1_fx * ta_0_zzz_1[j] +
                                 1.5 * fl1_fx * ta_z_zz_0[j] - 1.5 * fl1_fx * ta_z_zz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFD(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForFD_0_30(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFD_30_60(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForFD_0_30(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (0,30)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_xx_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx);

            auto ta_xx_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 1);

            auto ta_xx_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 2);

            auto ta_xx_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 3);

            auto ta_xx_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 4);

            auto ta_xx_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 5);

            auto ta_xy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 6);

            auto ta_xy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 7);

            auto ta_xy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 8);

            auto ta_xy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 9);

            auto ta_xy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 10);

            auto ta_xy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 11);

            auto ta_xz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 12);

            auto ta_xz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 13);

            auto ta_xz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 14);

            auto ta_xz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 15);

            auto ta_xz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 16);

            auto ta_xz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 17);

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_xx_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx);

            auto ta_xx_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 1);

            auto ta_xx_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 2);

            auto ta_xx_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 3);

            auto ta_xx_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 4);

            auto ta_xx_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 5);

            auto ta_xy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 6);

            auto ta_xy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 7);

            auto ta_xy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 8);

            auto ta_xy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 9);

            auto ta_xy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 10);

            auto ta_xy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 11);

            auto ta_xz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 12);

            auto ta_xz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 13);

            auto ta_xz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 14);

            auto ta_xz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 15);

            auto ta_xz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 16);

            auto ta_xz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 17);

            auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18);

            auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19);

            auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20);

            auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21);

            auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22);

            auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23);

            auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24);

            auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25);

            auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26);

            auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27);

            auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28);

            auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29);

            auto ta_x_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx);

            auto ta_x_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 1);

            auto ta_x_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 2);

            auto ta_x_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 3);

            auto ta_x_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 4);

            auto ta_x_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 5);

            auto ta_y_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 6);

            auto ta_y_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 7);

            auto ta_y_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 8);

            auto ta_y_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 9);

            auto ta_y_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 10);

            auto ta_y_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 11);

            auto ta_z_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 12);

            auto ta_z_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 13);

            auto ta_z_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 14);

            auto ta_z_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 15);

            auto ta_z_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 16);

            auto ta_z_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 17);

            auto ta_x_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx);

            auto ta_x_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 1);

            auto ta_x_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 2);

            auto ta_x_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 3);

            auto ta_x_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 4);

            auto ta_x_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 5);

            auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6);

            auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7);

            auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8);

            auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9);

            auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10);

            auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11);

            auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12);

            auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13);

            auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14);

            auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15);

            auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16);

            auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17);

            auto ta_xx_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx);

            auto ta_xx_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 1);

            auto ta_xx_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 2);

            auto ta_xy_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 3);

            auto ta_xy_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 4);

            auto ta_xy_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 5);

            auto ta_xz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 6);

            auto ta_xz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 7);

            auto ta_xz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 8);

            auto ta_yy_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 9);

            auto ta_yy_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 10);

            auto ta_yy_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 11);

            auto ta_yz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 12);

            auto ta_yz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 13);

            auto ta_yz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 14);

            auto ta_xx_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx);

            auto ta_xx_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 1);

            auto ta_xx_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 2);

            auto ta_xy_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 3);

            auto ta_xy_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 4);

            auto ta_xy_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 5);

            auto ta_xz_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 6);

            auto ta_xz_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 7);

            auto ta_xz_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 8);

            auto ta_yy_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 9);

            auto ta_yy_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 10);

            auto ta_yy_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 11);

            auto ta_yz_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 12);

            auto ta_yz_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 13);

            auto ta_yz_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 14);

            // set up pointers to integrals

            auto ta_xxx_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx);

            auto ta_xxx_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 1);

            auto ta_xxx_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 2);

            auto ta_xxx_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 3);

            auto ta_xxx_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 4);

            auto ta_xxx_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 5);

            auto ta_xxy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 6);

            auto ta_xxy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 7);

            auto ta_xxy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 8);

            auto ta_xxy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 9);

            auto ta_xxy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 10);

            auto ta_xxy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 11);

            auto ta_xxz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 12);

            auto ta_xxz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 13);

            auto ta_xxz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 14);

            auto ta_xxz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 15);

            auto ta_xxz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 16);

            auto ta_xxz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 17);

            auto ta_xyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 18);

            auto ta_xyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 19);

            auto ta_xyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 20);

            auto ta_xyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 21);

            auto ta_xyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 22);

            auto ta_xyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 23);

            auto ta_xyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 24);

            auto ta_xyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 25);

            auto ta_xyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 26);

            auto ta_xyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 27);

            auto ta_xyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 28);

            auto ta_xyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 29);

            // Batch of Integrals (0,30)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_xx_0, ta_x_xx_1, ta_x_xy_0, ta_x_xy_1, ta_x_xz_0, \
                                         ta_x_xz_1, ta_x_yy_0, ta_x_yy_1, ta_x_yz_0, ta_x_yz_1, ta_x_zz_0, ta_x_zz_1, \
                                         ta_xx_x_0, ta_xx_x_1, ta_xx_xx_0, ta_xx_xx_1, ta_xx_xy_0, ta_xx_xy_1, ta_xx_xz_0, \
                                         ta_xx_xz_1, ta_xx_y_0, ta_xx_y_1, ta_xx_yy_0, ta_xx_yy_1, ta_xx_yz_0, ta_xx_yz_1, \
                                         ta_xx_z_0, ta_xx_z_1, ta_xx_zz_0, ta_xx_zz_1, ta_xxx_xx_0, ta_xxx_xy_0, \
                                         ta_xxx_xz_0, ta_xxx_yy_0, ta_xxx_yz_0, ta_xxx_zz_0, ta_xxy_xx_0, ta_xxy_xy_0, \
                                         ta_xxy_xz_0, ta_xxy_yy_0, ta_xxy_yz_0, ta_xxy_zz_0, ta_xxz_xx_0, ta_xxz_xy_0, \
                                         ta_xxz_xz_0, ta_xxz_yy_0, ta_xxz_yz_0, ta_xxz_zz_0, ta_xy_x_0, ta_xy_x_1, \
                                         ta_xy_xx_0, ta_xy_xx_1, ta_xy_xy_0, ta_xy_xy_1, ta_xy_xz_0, ta_xy_xz_1, ta_xy_y_0, \
                                         ta_xy_y_1, ta_xy_yy_0, ta_xy_yy_1, ta_xy_yz_0, ta_xy_yz_1, ta_xy_z_0, ta_xy_z_1, \
                                         ta_xy_zz_0, ta_xy_zz_1, ta_xyy_xx_0, ta_xyy_xy_0, ta_xyy_xz_0, ta_xyy_yy_0, \
                                         ta_xyy_yz_0, ta_xyy_zz_0, ta_xyz_xx_0, ta_xyz_xy_0, ta_xyz_xz_0, ta_xyz_yy_0, \
                                         ta_xyz_yz_0, ta_xyz_zz_0, ta_xz_x_0, ta_xz_x_1, ta_xz_xx_0, ta_xz_xx_1, ta_xz_xy_0, \
                                         ta_xz_xy_1, ta_xz_xz_0, ta_xz_xz_1, ta_xz_y_0, ta_xz_y_1, ta_xz_yy_0, ta_xz_yy_1, \
                                         ta_xz_yz_0, ta_xz_yz_1, ta_xz_z_0, ta_xz_z_1, ta_xz_zz_0, ta_xz_zz_1, ta_y_xx_0, \
                                         ta_y_xx_1, ta_y_xy_0, ta_y_xy_1, ta_y_xz_0, ta_y_xz_1, ta_y_yy_0, ta_y_yy_1, \
                                         ta_y_yz_0, ta_y_yz_1, ta_y_zz_0, ta_y_zz_1, ta_yy_x_0, ta_yy_x_1, ta_yy_xx_0, \
                                         ta_yy_xx_1, ta_yy_xy_0, ta_yy_xy_1, ta_yy_xz_0, ta_yy_xz_1, ta_yy_y_0, ta_yy_y_1, \
                                         ta_yy_yy_0, ta_yy_yy_1, ta_yy_yz_0, ta_yy_yz_1, ta_yy_z_0, ta_yy_z_1, ta_yy_zz_0, \
                                         ta_yy_zz_1, ta_yz_x_0, ta_yz_x_1, ta_yz_xx_0, ta_yz_xx_1, ta_yz_xy_0, ta_yz_xy_1, \
                                         ta_yz_xz_0, ta_yz_xz_1, ta_yz_y_0, ta_yz_y_1, ta_yz_yy_0, ta_yz_yy_1, ta_yz_yz_0, \
                                         ta_yz_yz_1, ta_yz_z_0, ta_yz_z_1, ta_yz_zz_0, ta_yz_zz_1, ta_z_xx_0, ta_z_xx_1, \
                                         ta_z_xy_0, ta_z_xy_1, ta_z_xz_0, ta_z_xz_1, ta_z_yy_0, ta_z_yy_1, ta_z_yz_0, \
                                         ta_z_yz_1, ta_z_zz_0, ta_z_zz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxx_xx_0[j] = pa_x[j] * ta_xx_xx_0[j] - pc_x[j] * ta_xx_xx_1[j] + fl1_fx * ta_x_xx_0[j] - fl1_fx * ta_x_xx_1[j] +
                                 fl1_fx * ta_xx_x_0[j] - fl1_fx * ta_xx_x_1[j];

                ta_xxx_xy_0[j] = pa_x[j] * ta_xx_xy_0[j] - pc_x[j] * ta_xx_xy_1[j] + fl1_fx * ta_x_xy_0[j] - fl1_fx * ta_x_xy_1[j] +
                                 0.5 * fl1_fx * ta_xx_y_0[j] - 0.5 * fl1_fx * ta_xx_y_1[j];

                ta_xxx_xz_0[j] = pa_x[j] * ta_xx_xz_0[j] - pc_x[j] * ta_xx_xz_1[j] + fl1_fx * ta_x_xz_0[j] - fl1_fx * ta_x_xz_1[j] +
                                 0.5 * fl1_fx * ta_xx_z_0[j] - 0.5 * fl1_fx * ta_xx_z_1[j];

                ta_xxx_yy_0[j] = pa_x[j] * ta_xx_yy_0[j] - pc_x[j] * ta_xx_yy_1[j] + fl1_fx * ta_x_yy_0[j] - fl1_fx * ta_x_yy_1[j];

                ta_xxx_yz_0[j] = pa_x[j] * ta_xx_yz_0[j] - pc_x[j] * ta_xx_yz_1[j] + fl1_fx * ta_x_yz_0[j] - fl1_fx * ta_x_yz_1[j];

                ta_xxx_zz_0[j] = pa_x[j] * ta_xx_zz_0[j] - pc_x[j] * ta_xx_zz_1[j] + fl1_fx * ta_x_zz_0[j] - fl1_fx * ta_x_zz_1[j];

                ta_xxy_xx_0[j] = pa_x[j] * ta_xy_xx_0[j] - pc_x[j] * ta_xy_xx_1[j] + 0.5 * fl1_fx * ta_y_xx_0[j] - 0.5 * fl1_fx * ta_y_xx_1[j] +
                                 fl1_fx * ta_xy_x_0[j] - fl1_fx * ta_xy_x_1[j];

                ta_xxy_xy_0[j] = pa_x[j] * ta_xy_xy_0[j] - pc_x[j] * ta_xy_xy_1[j] + 0.5 * fl1_fx * ta_y_xy_0[j] - 0.5 * fl1_fx * ta_y_xy_1[j] +
                                 0.5 * fl1_fx * ta_xy_y_0[j] - 0.5 * fl1_fx * ta_xy_y_1[j];

                ta_xxy_xz_0[j] = pa_x[j] * ta_xy_xz_0[j] - pc_x[j] * ta_xy_xz_1[j] + 0.5 * fl1_fx * ta_y_xz_0[j] - 0.5 * fl1_fx * ta_y_xz_1[j] +
                                 0.5 * fl1_fx * ta_xy_z_0[j] - 0.5 * fl1_fx * ta_xy_z_1[j];

                ta_xxy_yy_0[j] = pa_x[j] * ta_xy_yy_0[j] - pc_x[j] * ta_xy_yy_1[j] + 0.5 * fl1_fx * ta_y_yy_0[j] - 0.5 * fl1_fx * ta_y_yy_1[j];

                ta_xxy_yz_0[j] = pa_x[j] * ta_xy_yz_0[j] - pc_x[j] * ta_xy_yz_1[j] + 0.5 * fl1_fx * ta_y_yz_0[j] - 0.5 * fl1_fx * ta_y_yz_1[j];

                ta_xxy_zz_0[j] = pa_x[j] * ta_xy_zz_0[j] - pc_x[j] * ta_xy_zz_1[j] + 0.5 * fl1_fx * ta_y_zz_0[j] - 0.5 * fl1_fx * ta_y_zz_1[j];

                ta_xxz_xx_0[j] = pa_x[j] * ta_xz_xx_0[j] - pc_x[j] * ta_xz_xx_1[j] + 0.5 * fl1_fx * ta_z_xx_0[j] - 0.5 * fl1_fx * ta_z_xx_1[j] +
                                 fl1_fx * ta_xz_x_0[j] - fl1_fx * ta_xz_x_1[j];

                ta_xxz_xy_0[j] = pa_x[j] * ta_xz_xy_0[j] - pc_x[j] * ta_xz_xy_1[j] + 0.5 * fl1_fx * ta_z_xy_0[j] - 0.5 * fl1_fx * ta_z_xy_1[j] +
                                 0.5 * fl1_fx * ta_xz_y_0[j] - 0.5 * fl1_fx * ta_xz_y_1[j];

                ta_xxz_xz_0[j] = pa_x[j] * ta_xz_xz_0[j] - pc_x[j] * ta_xz_xz_1[j] + 0.5 * fl1_fx * ta_z_xz_0[j] - 0.5 * fl1_fx * ta_z_xz_1[j] +
                                 0.5 * fl1_fx * ta_xz_z_0[j] - 0.5 * fl1_fx * ta_xz_z_1[j];

                ta_xxz_yy_0[j] = pa_x[j] * ta_xz_yy_0[j] - pc_x[j] * ta_xz_yy_1[j] + 0.5 * fl1_fx * ta_z_yy_0[j] - 0.5 * fl1_fx * ta_z_yy_1[j];

                ta_xxz_yz_0[j] = pa_x[j] * ta_xz_yz_0[j] - pc_x[j] * ta_xz_yz_1[j] + 0.5 * fl1_fx * ta_z_yz_0[j] - 0.5 * fl1_fx * ta_z_yz_1[j];

                ta_xxz_zz_0[j] = pa_x[j] * ta_xz_zz_0[j] - pc_x[j] * ta_xz_zz_1[j] + 0.5 * fl1_fx * ta_z_zz_0[j] - 0.5 * fl1_fx * ta_z_zz_1[j];

                ta_xyy_xx_0[j] = pa_x[j] * ta_yy_xx_0[j] - pc_x[j] * ta_yy_xx_1[j] + fl1_fx * ta_yy_x_0[j] - fl1_fx * ta_yy_x_1[j];

                ta_xyy_xy_0[j] = pa_x[j] * ta_yy_xy_0[j] - pc_x[j] * ta_yy_xy_1[j] + 0.5 * fl1_fx * ta_yy_y_0[j] - 0.5 * fl1_fx * ta_yy_y_1[j];

                ta_xyy_xz_0[j] = pa_x[j] * ta_yy_xz_0[j] - pc_x[j] * ta_yy_xz_1[j] + 0.5 * fl1_fx * ta_yy_z_0[j] - 0.5 * fl1_fx * ta_yy_z_1[j];

                ta_xyy_yy_0[j] = pa_x[j] * ta_yy_yy_0[j] - pc_x[j] * ta_yy_yy_1[j];

                ta_xyy_yz_0[j] = pa_x[j] * ta_yy_yz_0[j] - pc_x[j] * ta_yy_yz_1[j];

                ta_xyy_zz_0[j] = pa_x[j] * ta_yy_zz_0[j] - pc_x[j] * ta_yy_zz_1[j];

                ta_xyz_xx_0[j] = pa_x[j] * ta_yz_xx_0[j] - pc_x[j] * ta_yz_xx_1[j] + fl1_fx * ta_yz_x_0[j] - fl1_fx * ta_yz_x_1[j];

                ta_xyz_xy_0[j] = pa_x[j] * ta_yz_xy_0[j] - pc_x[j] * ta_yz_xy_1[j] + 0.5 * fl1_fx * ta_yz_y_0[j] - 0.5 * fl1_fx * ta_yz_y_1[j];

                ta_xyz_xz_0[j] = pa_x[j] * ta_yz_xz_0[j] - pc_x[j] * ta_yz_xz_1[j] + 0.5 * fl1_fx * ta_yz_z_0[j] - 0.5 * fl1_fx * ta_yz_z_1[j];

                ta_xyz_yy_0[j] = pa_x[j] * ta_yz_yy_0[j] - pc_x[j] * ta_yz_yy_1[j];

                ta_xyz_yz_0[j] = pa_x[j] * ta_yz_yz_0[j] - pc_x[j] * ta_yz_yz_1[j];

                ta_xyz_zz_0[j] = pa_x[j] * ta_yz_zz_0[j] - pc_x[j] * ta_yz_zz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFD_30_60(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (30,60)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_zz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 30);

            auto ta_zz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 31);

            auto ta_zz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 32);

            auto ta_zz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 33);

            auto ta_zz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 34);

            auto ta_zz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 35);

            auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18);

            auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19);

            auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20);

            auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21);

            auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22);

            auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23);

            auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24);

            auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25);

            auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26);

            auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27);

            auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28);

            auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29);

            auto ta_zz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 30);

            auto ta_zz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 31);

            auto ta_zz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 32);

            auto ta_zz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 33);

            auto ta_zz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 34);

            auto ta_zz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 35);

            auto ta_y_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 6);

            auto ta_y_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 7);

            auto ta_y_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 8);

            auto ta_y_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 9);

            auto ta_y_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 10);

            auto ta_y_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 11);

            auto ta_z_xx_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 12);

            auto ta_z_xy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 13);

            auto ta_z_xz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 14);

            auto ta_z_yy_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 15);

            auto ta_z_yz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 16);

            auto ta_z_zz_0 = primBuffer.data(pidx_a_1_2_m0 + 18 * idx + 17);

            auto ta_y_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 6);

            auto ta_y_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 7);

            auto ta_y_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 8);

            auto ta_y_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 9);

            auto ta_y_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 10);

            auto ta_y_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 11);

            auto ta_z_xx_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 12);

            auto ta_z_xy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 13);

            auto ta_z_xz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 14);

            auto ta_z_yy_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 15);

            auto ta_z_yz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 16);

            auto ta_z_zz_1 = primBuffer.data(pidx_a_1_2_m1 + 18 * idx + 17);

            auto ta_yy_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 9);

            auto ta_yy_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 10);

            auto ta_yy_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 11);

            auto ta_yz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 12);

            auto ta_yz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 13);

            auto ta_yz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 14);

            auto ta_zz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 15);

            auto ta_zz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 16);

            auto ta_zz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 17);

            auto ta_yy_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 9);

            auto ta_yy_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 10);

            auto ta_yy_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 11);

            auto ta_yz_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 12);

            auto ta_yz_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 13);

            auto ta_yz_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 14);

            auto ta_zz_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 15);

            auto ta_zz_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 16);

            auto ta_zz_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 17);

            // set up pointers to integrals

            auto ta_xzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 30);

            auto ta_xzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 31);

            auto ta_xzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 32);

            auto ta_xzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 33);

            auto ta_xzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 34);

            auto ta_xzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 35);

            auto ta_yyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 36);

            auto ta_yyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 37);

            auto ta_yyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 38);

            auto ta_yyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 39);

            auto ta_yyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 40);

            auto ta_yyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 41);

            auto ta_yyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 42);

            auto ta_yyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 43);

            auto ta_yyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 44);

            auto ta_yyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 45);

            auto ta_yyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 46);

            auto ta_yyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 47);

            auto ta_yzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 48);

            auto ta_yzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 49);

            auto ta_yzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 50);

            auto ta_yzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 51);

            auto ta_yzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 52);

            auto ta_yzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 53);

            auto ta_zzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 54);

            auto ta_zzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 55);

            auto ta_zzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 56);

            auto ta_zzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 57);

            auto ta_zzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 58);

            auto ta_zzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 59);

            // Batch of Integrals (30,60)

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xzz_xx_0, ta_xzz_xy_0, ta_xzz_xz_0, \
                                         ta_xzz_yy_0, ta_xzz_yz_0, ta_xzz_zz_0, ta_y_xx_0, ta_y_xx_1, ta_y_xy_0, ta_y_xy_1, \
                                         ta_y_xz_0, ta_y_xz_1, ta_y_yy_0, ta_y_yy_1, ta_y_yz_0, ta_y_yz_1, ta_y_zz_0, \
                                         ta_y_zz_1, ta_yy_x_0, ta_yy_x_1, ta_yy_xx_0, ta_yy_xx_1, ta_yy_xy_0, ta_yy_xy_1, \
                                         ta_yy_xz_0, ta_yy_xz_1, ta_yy_y_0, ta_yy_y_1, ta_yy_yy_0, ta_yy_yy_1, ta_yy_yz_0, \
                                         ta_yy_yz_1, ta_yy_z_0, ta_yy_z_1, ta_yy_zz_0, ta_yy_zz_1, ta_yyy_xx_0, ta_yyy_xy_0, \
                                         ta_yyy_xz_0, ta_yyy_yy_0, ta_yyy_yz_0, ta_yyy_zz_0, ta_yyz_xx_0, ta_yyz_xy_0, \
                                         ta_yyz_xz_0, ta_yyz_yy_0, ta_yyz_yz_0, ta_yyz_zz_0, ta_yz_x_0, ta_yz_x_1, \
                                         ta_yz_xx_0, ta_yz_xx_1, ta_yz_xy_0, ta_yz_xy_1, ta_yz_xz_0, ta_yz_xz_1, ta_yz_y_0, \
                                         ta_yz_y_1, ta_yz_yy_0, ta_yz_yy_1, ta_yz_yz_0, ta_yz_yz_1, ta_yz_z_0, ta_yz_z_1, \
                                         ta_yz_zz_0, ta_yz_zz_1, ta_yzz_xx_0, ta_yzz_xy_0, ta_yzz_xz_0, ta_yzz_yy_0, \
                                         ta_yzz_yz_0, ta_yzz_zz_0, ta_z_xx_0, ta_z_xx_1, ta_z_xy_0, ta_z_xy_1, ta_z_xz_0, \
                                         ta_z_xz_1, ta_z_yy_0, ta_z_yy_1, ta_z_yz_0, ta_z_yz_1, ta_z_zz_0, ta_z_zz_1, \
                                         ta_zz_x_0, ta_zz_x_1, ta_zz_xx_0, ta_zz_xx_1, ta_zz_xy_0, ta_zz_xy_1, ta_zz_xz_0, \
                                         ta_zz_xz_1, ta_zz_y_0, ta_zz_y_1, ta_zz_yy_0, ta_zz_yy_1, ta_zz_yz_0, ta_zz_yz_1, \
                                         ta_zz_z_0, ta_zz_z_1, ta_zz_zz_0, ta_zz_zz_1, ta_zzz_xx_0, ta_zzz_xy_0, \
                                         ta_zzz_xz_0, ta_zzz_yy_0, ta_zzz_yz_0, ta_zzz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xzz_xx_0[j] = pa_x[j] * ta_zz_xx_0[j] - pc_x[j] * ta_zz_xx_1[j] + fl1_fx * ta_zz_x_0[j] - fl1_fx * ta_zz_x_1[j];

                ta_xzz_xy_0[j] = pa_x[j] * ta_zz_xy_0[j] - pc_x[j] * ta_zz_xy_1[j] + 0.5 * fl1_fx * ta_zz_y_0[j] - 0.5 * fl1_fx * ta_zz_y_1[j];

                ta_xzz_xz_0[j] = pa_x[j] * ta_zz_xz_0[j] - pc_x[j] * ta_zz_xz_1[j] + 0.5 * fl1_fx * ta_zz_z_0[j] - 0.5 * fl1_fx * ta_zz_z_1[j];

                ta_xzz_yy_0[j] = pa_x[j] * ta_zz_yy_0[j] - pc_x[j] * ta_zz_yy_1[j];

                ta_xzz_yz_0[j] = pa_x[j] * ta_zz_yz_0[j] - pc_x[j] * ta_zz_yz_1[j];

                ta_xzz_zz_0[j] = pa_x[j] * ta_zz_zz_0[j] - pc_x[j] * ta_zz_zz_1[j];

                ta_yyy_xx_0[j] = pa_y[j] * ta_yy_xx_0[j] - pc_y[j] * ta_yy_xx_1[j] + fl1_fx * ta_y_xx_0[j] - fl1_fx * ta_y_xx_1[j];

                ta_yyy_xy_0[j] = pa_y[j] * ta_yy_xy_0[j] - pc_y[j] * ta_yy_xy_1[j] + fl1_fx * ta_y_xy_0[j] - fl1_fx * ta_y_xy_1[j] +
                                 0.5 * fl1_fx * ta_yy_x_0[j] - 0.5 * fl1_fx * ta_yy_x_1[j];

                ta_yyy_xz_0[j] = pa_y[j] * ta_yy_xz_0[j] - pc_y[j] * ta_yy_xz_1[j] + fl1_fx * ta_y_xz_0[j] - fl1_fx * ta_y_xz_1[j];

                ta_yyy_yy_0[j] = pa_y[j] * ta_yy_yy_0[j] - pc_y[j] * ta_yy_yy_1[j] + fl1_fx * ta_y_yy_0[j] - fl1_fx * ta_y_yy_1[j] +
                                 fl1_fx * ta_yy_y_0[j] - fl1_fx * ta_yy_y_1[j];

                ta_yyy_yz_0[j] = pa_y[j] * ta_yy_yz_0[j] - pc_y[j] * ta_yy_yz_1[j] + fl1_fx * ta_y_yz_0[j] - fl1_fx * ta_y_yz_1[j] +
                                 0.5 * fl1_fx * ta_yy_z_0[j] - 0.5 * fl1_fx * ta_yy_z_1[j];

                ta_yyy_zz_0[j] = pa_y[j] * ta_yy_zz_0[j] - pc_y[j] * ta_yy_zz_1[j] + fl1_fx * ta_y_zz_0[j] - fl1_fx * ta_y_zz_1[j];

                ta_yyz_xx_0[j] = pa_y[j] * ta_yz_xx_0[j] - pc_y[j] * ta_yz_xx_1[j] + 0.5 * fl1_fx * ta_z_xx_0[j] - 0.5 * fl1_fx * ta_z_xx_1[j];

                ta_yyz_xy_0[j] = pa_y[j] * ta_yz_xy_0[j] - pc_y[j] * ta_yz_xy_1[j] + 0.5 * fl1_fx * ta_z_xy_0[j] - 0.5 * fl1_fx * ta_z_xy_1[j] +
                                 0.5 * fl1_fx * ta_yz_x_0[j] - 0.5 * fl1_fx * ta_yz_x_1[j];

                ta_yyz_xz_0[j] = pa_y[j] * ta_yz_xz_0[j] - pc_y[j] * ta_yz_xz_1[j] + 0.5 * fl1_fx * ta_z_xz_0[j] - 0.5 * fl1_fx * ta_z_xz_1[j];

                ta_yyz_yy_0[j] = pa_y[j] * ta_yz_yy_0[j] - pc_y[j] * ta_yz_yy_1[j] + 0.5 * fl1_fx * ta_z_yy_0[j] - 0.5 * fl1_fx * ta_z_yy_1[j] +
                                 fl1_fx * ta_yz_y_0[j] - fl1_fx * ta_yz_y_1[j];

                ta_yyz_yz_0[j] = pa_y[j] * ta_yz_yz_0[j] - pc_y[j] * ta_yz_yz_1[j] + 0.5 * fl1_fx * ta_z_yz_0[j] - 0.5 * fl1_fx * ta_z_yz_1[j] +
                                 0.5 * fl1_fx * ta_yz_z_0[j] - 0.5 * fl1_fx * ta_yz_z_1[j];

                ta_yyz_zz_0[j] = pa_y[j] * ta_yz_zz_0[j] - pc_y[j] * ta_yz_zz_1[j] + 0.5 * fl1_fx * ta_z_zz_0[j] - 0.5 * fl1_fx * ta_z_zz_1[j];

                ta_yzz_xx_0[j] = pa_y[j] * ta_zz_xx_0[j] - pc_y[j] * ta_zz_xx_1[j];

                ta_yzz_xy_0[j] = pa_y[j] * ta_zz_xy_0[j] - pc_y[j] * ta_zz_xy_1[j] + 0.5 * fl1_fx * ta_zz_x_0[j] - 0.5 * fl1_fx * ta_zz_x_1[j];

                ta_yzz_xz_0[j] = pa_y[j] * ta_zz_xz_0[j] - pc_y[j] * ta_zz_xz_1[j];

                ta_yzz_yy_0[j] = pa_y[j] * ta_zz_yy_0[j] - pc_y[j] * ta_zz_yy_1[j] + fl1_fx * ta_zz_y_0[j] - fl1_fx * ta_zz_y_1[j];

                ta_yzz_yz_0[j] = pa_y[j] * ta_zz_yz_0[j] - pc_y[j] * ta_zz_yz_1[j] + 0.5 * fl1_fx * ta_zz_z_0[j] - 0.5 * fl1_fx * ta_zz_z_1[j];

                ta_yzz_zz_0[j] = pa_y[j] * ta_zz_zz_0[j] - pc_y[j] * ta_zz_zz_1[j];

                ta_zzz_xx_0[j] = pa_z[j] * ta_zz_xx_0[j] - pc_z[j] * ta_zz_xx_1[j] + fl1_fx * ta_z_xx_0[j] - fl1_fx * ta_z_xx_1[j];

                ta_zzz_xy_0[j] = pa_z[j] * ta_zz_xy_0[j] - pc_z[j] * ta_zz_xy_1[j] + fl1_fx * ta_z_xy_0[j] - fl1_fx * ta_z_xy_1[j];

                ta_zzz_xz_0[j] = pa_z[j] * ta_zz_xz_0[j] - pc_z[j] * ta_zz_xz_1[j] + fl1_fx * ta_z_xz_0[j] - fl1_fx * ta_z_xz_1[j] +
                                 0.5 * fl1_fx * ta_zz_x_0[j] - 0.5 * fl1_fx * ta_zz_x_1[j];

                ta_zzz_yy_0[j] = pa_z[j] * ta_zz_yy_0[j] - pc_z[j] * ta_zz_yy_1[j] + fl1_fx * ta_z_yy_0[j] - fl1_fx * ta_z_yy_1[j];

                ta_zzz_yz_0[j] = pa_z[j] * ta_zz_yz_0[j] - pc_z[j] * ta_zz_yz_1[j] + fl1_fx * ta_z_yz_0[j] - fl1_fx * ta_z_yz_1[j] +
                                 0.5 * fl1_fx * ta_zz_y_0[j] - 0.5 * fl1_fx * ta_zz_y_1[j];

                ta_zzz_zz_0[j] = pa_z[j] * ta_zz_zz_0[j] - pc_z[j] * ta_zz_zz_1[j] + fl1_fx * ta_z_zz_0[j] - fl1_fx * ta_z_zz_1[j] +
                                 fl1_fx * ta_zz_z_0[j] - fl1_fx * ta_zz_z_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForDG(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForDG_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDG_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForDG_0_45(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_0_xxxx_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx);

            auto ta_0_xxxy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 1);

            auto ta_0_xxxz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 2);

            auto ta_0_xxyy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 3);

            auto ta_0_xxyz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 4);

            auto ta_0_xxzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 5);

            auto ta_0_xyyy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 6);

            auto ta_0_xyyz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 7);

            auto ta_0_xyzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 8);

            auto ta_0_xzzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 9);

            auto ta_0_yyyy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 10);

            auto ta_0_yyyz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 11);

            auto ta_0_yyzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 12);

            auto ta_0_yzzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 13);

            auto ta_0_zzzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 14);

            auto ta_0_xxxx_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx);

            auto ta_0_xxxy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 1);

            auto ta_0_xxxz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 2);

            auto ta_0_xxyy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 3);

            auto ta_0_xxyz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 4);

            auto ta_0_xxzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 5);

            auto ta_0_xyyy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 6);

            auto ta_0_xyyz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 7);

            auto ta_0_xyzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 8);

            auto ta_0_xzzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 9);

            auto ta_0_yyyy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 10);

            auto ta_0_yyyz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 11);

            auto ta_0_yyzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 12);

            auto ta_0_yzzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 13);

            auto ta_0_zzzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 14);

            auto ta_x_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx);

            auto ta_x_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 1);

            auto ta_x_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 2);

            auto ta_x_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 3);

            auto ta_x_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 4);

            auto ta_x_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 5);

            auto ta_x_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 6);

            auto ta_x_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 7);

            auto ta_x_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 8);

            auto ta_x_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 9);

            auto ta_y_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 10);

            auto ta_y_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 11);

            auto ta_y_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 12);

            auto ta_y_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 13);

            auto ta_y_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 14);

            auto ta_y_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 15);

            auto ta_y_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 16);

            auto ta_y_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 17);

            auto ta_y_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 18);

            auto ta_y_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 19);

            auto ta_z_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 20);

            auto ta_z_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 21);

            auto ta_z_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 22);

            auto ta_z_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 23);

            auto ta_z_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 24);

            auto ta_z_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 25);

            auto ta_z_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 26);

            auto ta_z_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 27);

            auto ta_z_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 28);

            auto ta_z_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 29);

            auto ta_x_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx);

            auto ta_x_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 1);

            auto ta_x_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 2);

            auto ta_x_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 3);

            auto ta_x_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 4);

            auto ta_x_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 5);

            auto ta_x_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 6);

            auto ta_x_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 7);

            auto ta_x_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 8);

            auto ta_x_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 9);

            auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10);

            auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11);

            auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12);

            auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13);

            auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14);

            auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15);

            auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16);

            auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17);

            auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18);

            auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19);

            auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20);

            auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21);

            auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22);

            auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23);

            auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24);

            auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25);

            auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26);

            auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27);

            auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28);

            auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29);

            // set up pointers to integrals

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

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_0_xxxx_0, ta_0_xxxx_1, ta_0_xxxy_0, ta_0_xxxy_1, \
                                         ta_0_xxxz_0, ta_0_xxxz_1, ta_0_xxyy_0, ta_0_xxyy_1, ta_0_xxyz_0, ta_0_xxyz_1, \
                                         ta_0_xxzz_0, ta_0_xxzz_1, ta_0_xyyy_0, ta_0_xyyy_1, ta_0_xyyz_0, ta_0_xyyz_1, \
                                         ta_0_xyzz_0, ta_0_xyzz_1, ta_0_xzzz_0, ta_0_xzzz_1, ta_0_yyyy_0, ta_0_yyyy_1, \
                                         ta_0_yyyz_0, ta_0_yyyz_1, ta_0_yyzz_0, ta_0_yyzz_1, ta_0_yzzz_0, ta_0_yzzz_1, \
                                         ta_0_zzzz_0, ta_0_zzzz_1, ta_x_xxx_0, ta_x_xxx_1, ta_x_xxxx_0, ta_x_xxxx_1, \
                                         ta_x_xxxy_0, ta_x_xxxy_1, ta_x_xxxz_0, ta_x_xxxz_1, ta_x_xxy_0, ta_x_xxy_1, \
                                         ta_x_xxyy_0, ta_x_xxyy_1, ta_x_xxyz_0, ta_x_xxyz_1, ta_x_xxz_0, ta_x_xxz_1, \
                                         ta_x_xxzz_0, ta_x_xxzz_1, ta_x_xyy_0, ta_x_xyy_1, ta_x_xyyy_0, ta_x_xyyy_1, \
                                         ta_x_xyyz_0, ta_x_xyyz_1, ta_x_xyz_0, ta_x_xyz_1, ta_x_xyzz_0, ta_x_xyzz_1, \
                                         ta_x_xzz_0, ta_x_xzz_1, ta_x_xzzz_0, ta_x_xzzz_1, ta_x_yyy_0, ta_x_yyy_1, \
                                         ta_x_yyyy_0, ta_x_yyyy_1, ta_x_yyyz_0, ta_x_yyyz_1, ta_x_yyz_0, ta_x_yyz_1, \
                                         ta_x_yyzz_0, ta_x_yyzz_1, ta_x_yzz_0, ta_x_yzz_1, ta_x_yzzz_0, ta_x_yzzz_1, \
                                         ta_x_zzz_0, ta_x_zzz_1, ta_x_zzzz_0, ta_x_zzzz_1, ta_xx_xxxx_0, ta_xx_xxxy_0, \
                                         ta_xx_xxxz_0, ta_xx_xxyy_0, ta_xx_xxyz_0, ta_xx_xxzz_0, ta_xx_xyyy_0, ta_xx_xyyz_0, \
                                         ta_xx_xyzz_0, ta_xx_xzzz_0, ta_xx_yyyy_0, ta_xx_yyyz_0, ta_xx_yyzz_0, ta_xx_yzzz_0, \
                                         ta_xx_zzzz_0, ta_xy_xxxx_0, ta_xy_xxxy_0, ta_xy_xxxz_0, ta_xy_xxyy_0, ta_xy_xxyz_0, \
                                         ta_xy_xxzz_0, ta_xy_xyyy_0, ta_xy_xyyz_0, ta_xy_xyzz_0, ta_xy_xzzz_0, ta_xy_yyyy_0, \
                                         ta_xy_yyyz_0, ta_xy_yyzz_0, ta_xy_yzzz_0, ta_xy_zzzz_0, ta_xz_xxxx_0, ta_xz_xxxy_0, \
                                         ta_xz_xxxz_0, ta_xz_xxyy_0, ta_xz_xxyz_0, ta_xz_xxzz_0, ta_xz_xyyy_0, ta_xz_xyyz_0, \
                                         ta_xz_xyzz_0, ta_xz_xzzz_0, ta_xz_yyyy_0, ta_xz_yyyz_0, ta_xz_yyzz_0, ta_xz_yzzz_0, \
                                         ta_xz_zzzz_0, ta_y_xxx_0, ta_y_xxx_1, ta_y_xxxx_0, ta_y_xxxx_1, ta_y_xxxy_0, \
                                         ta_y_xxxy_1, ta_y_xxxz_0, ta_y_xxxz_1, ta_y_xxy_0, ta_y_xxy_1, ta_y_xxyy_0, \
                                         ta_y_xxyy_1, ta_y_xxyz_0, ta_y_xxyz_1, ta_y_xxz_0, ta_y_xxz_1, ta_y_xxzz_0, \
                                         ta_y_xxzz_1, ta_y_xyy_0, ta_y_xyy_1, ta_y_xyyy_0, ta_y_xyyy_1, ta_y_xyyz_0, \
                                         ta_y_xyyz_1, ta_y_xyz_0, ta_y_xyz_1, ta_y_xyzz_0, ta_y_xyzz_1, ta_y_xzz_0, \
                                         ta_y_xzz_1, ta_y_xzzz_0, ta_y_xzzz_1, ta_y_yyy_0, ta_y_yyy_1, ta_y_yyyy_0, \
                                         ta_y_yyyy_1, ta_y_yyyz_0, ta_y_yyyz_1, ta_y_yyz_0, ta_y_yyz_1, ta_y_yyzz_0, \
                                         ta_y_yyzz_1, ta_y_yzz_0, ta_y_yzz_1, ta_y_yzzz_0, ta_y_yzzz_1, ta_y_zzz_0, \
                                         ta_y_zzz_1, ta_y_zzzz_0, ta_y_zzzz_1, ta_z_xxx_0, ta_z_xxx_1, ta_z_xxxx_0, \
                                         ta_z_xxxx_1, ta_z_xxxy_0, ta_z_xxxy_1, ta_z_xxxz_0, ta_z_xxxz_1, ta_z_xxy_0, \
                                         ta_z_xxy_1, ta_z_xxyy_0, ta_z_xxyy_1, ta_z_xxyz_0, ta_z_xxyz_1, ta_z_xxz_0, \
                                         ta_z_xxz_1, ta_z_xxzz_0, ta_z_xxzz_1, ta_z_xyy_0, ta_z_xyy_1, ta_z_xyyy_0, \
                                         ta_z_xyyy_1, ta_z_xyyz_0, ta_z_xyyz_1, ta_z_xyz_0, ta_z_xyz_1, ta_z_xyzz_0, \
                                         ta_z_xyzz_1, ta_z_xzz_0, ta_z_xzz_1, ta_z_xzzz_0, ta_z_xzzz_1, ta_z_yyy_0, \
                                         ta_z_yyy_1, ta_z_yyyy_0, ta_z_yyyy_1, ta_z_yyyz_0, ta_z_yyyz_1, ta_z_yyz_0, \
                                         ta_z_yyz_1, ta_z_yyzz_0, ta_z_yyzz_1, ta_z_yzz_0, ta_z_yzz_1, ta_z_yzzz_0, \
                                         ta_z_yzzz_1, ta_z_zzz_0, ta_z_zzz_1, ta_z_zzzz_0, ta_z_zzzz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xx_xxxx_0[j] = pa_x[j] * ta_x_xxxx_0[j] - pc_x[j] * ta_x_xxxx_1[j] + 0.5 * fl1_fx * ta_0_xxxx_0[j] -
                                  0.5 * fl1_fx * ta_0_xxxx_1[j] + 2.0 * fl1_fx * ta_x_xxx_0[j] - 2.0 * fl1_fx * ta_x_xxx_1[j];

                ta_xx_xxxy_0[j] = pa_x[j] * ta_x_xxxy_0[j] - pc_x[j] * ta_x_xxxy_1[j] + 0.5 * fl1_fx * ta_0_xxxy_0[j] -
                                  0.5 * fl1_fx * ta_0_xxxy_1[j] + 1.5 * fl1_fx * ta_x_xxy_0[j] - 1.5 * fl1_fx * ta_x_xxy_1[j];

                ta_xx_xxxz_0[j] = pa_x[j] * ta_x_xxxz_0[j] - pc_x[j] * ta_x_xxxz_1[j] + 0.5 * fl1_fx * ta_0_xxxz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxxz_1[j] + 1.5 * fl1_fx * ta_x_xxz_0[j] - 1.5 * fl1_fx * ta_x_xxz_1[j];

                ta_xx_xxyy_0[j] = pa_x[j] * ta_x_xxyy_0[j] - pc_x[j] * ta_x_xxyy_1[j] + 0.5 * fl1_fx * ta_0_xxyy_0[j] -
                                  0.5 * fl1_fx * ta_0_xxyy_1[j] + fl1_fx * ta_x_xyy_0[j] - fl1_fx * ta_x_xyy_1[j];

                ta_xx_xxyz_0[j] = pa_x[j] * ta_x_xxyz_0[j] - pc_x[j] * ta_x_xxyz_1[j] + 0.5 * fl1_fx * ta_0_xxyz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxyz_1[j] + fl1_fx * ta_x_xyz_0[j] - fl1_fx * ta_x_xyz_1[j];

                ta_xx_xxzz_0[j] = pa_x[j] * ta_x_xxzz_0[j] - pc_x[j] * ta_x_xxzz_1[j] + 0.5 * fl1_fx * ta_0_xxzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxzz_1[j] + fl1_fx * ta_x_xzz_0[j] - fl1_fx * ta_x_xzz_1[j];

                ta_xx_xyyy_0[j] = pa_x[j] * ta_x_xyyy_0[j] - pc_x[j] * ta_x_xyyy_1[j] + 0.5 * fl1_fx * ta_0_xyyy_0[j] -
                                  0.5 * fl1_fx * ta_0_xyyy_1[j] + 0.5 * fl1_fx * ta_x_yyy_0[j] - 0.5 * fl1_fx * ta_x_yyy_1[j];

                ta_xx_xyyz_0[j] = pa_x[j] * ta_x_xyyz_0[j] - pc_x[j] * ta_x_xyyz_1[j] + 0.5 * fl1_fx * ta_0_xyyz_0[j] -
                                  0.5 * fl1_fx * ta_0_xyyz_1[j] + 0.5 * fl1_fx * ta_x_yyz_0[j] - 0.5 * fl1_fx * ta_x_yyz_1[j];

                ta_xx_xyzz_0[j] = pa_x[j] * ta_x_xyzz_0[j] - pc_x[j] * ta_x_xyzz_1[j] + 0.5 * fl1_fx * ta_0_xyzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xyzz_1[j] + 0.5 * fl1_fx * ta_x_yzz_0[j] - 0.5 * fl1_fx * ta_x_yzz_1[j];

                ta_xx_xzzz_0[j] = pa_x[j] * ta_x_xzzz_0[j] - pc_x[j] * ta_x_xzzz_1[j] + 0.5 * fl1_fx * ta_0_xzzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xzzz_1[j] + 0.5 * fl1_fx * ta_x_zzz_0[j] - 0.5 * fl1_fx * ta_x_zzz_1[j];

                ta_xx_yyyy_0[j] = pa_x[j] * ta_x_yyyy_0[j] - pc_x[j] * ta_x_yyyy_1[j] + 0.5 * fl1_fx * ta_0_yyyy_0[j] - 0.5 * fl1_fx * ta_0_yyyy_1[j];

                ta_xx_yyyz_0[j] = pa_x[j] * ta_x_yyyz_0[j] - pc_x[j] * ta_x_yyyz_1[j] + 0.5 * fl1_fx * ta_0_yyyz_0[j] - 0.5 * fl1_fx * ta_0_yyyz_1[j];

                ta_xx_yyzz_0[j] = pa_x[j] * ta_x_yyzz_0[j] - pc_x[j] * ta_x_yyzz_1[j] + 0.5 * fl1_fx * ta_0_yyzz_0[j] - 0.5 * fl1_fx * ta_0_yyzz_1[j];

                ta_xx_yzzz_0[j] = pa_x[j] * ta_x_yzzz_0[j] - pc_x[j] * ta_x_yzzz_1[j] + 0.5 * fl1_fx * ta_0_yzzz_0[j] - 0.5 * fl1_fx * ta_0_yzzz_1[j];

                ta_xx_zzzz_0[j] = pa_x[j] * ta_x_zzzz_0[j] - pc_x[j] * ta_x_zzzz_1[j] + 0.5 * fl1_fx * ta_0_zzzz_0[j] - 0.5 * fl1_fx * ta_0_zzzz_1[j];

                ta_xy_xxxx_0[j] = pa_x[j] * ta_y_xxxx_0[j] - pc_x[j] * ta_y_xxxx_1[j] + 2.0 * fl1_fx * ta_y_xxx_0[j] - 2.0 * fl1_fx * ta_y_xxx_1[j];

                ta_xy_xxxy_0[j] = pa_x[j] * ta_y_xxxy_0[j] - pc_x[j] * ta_y_xxxy_1[j] + 1.5 * fl1_fx * ta_y_xxy_0[j] - 1.5 * fl1_fx * ta_y_xxy_1[j];

                ta_xy_xxxz_0[j] = pa_x[j] * ta_y_xxxz_0[j] - pc_x[j] * ta_y_xxxz_1[j] + 1.5 * fl1_fx * ta_y_xxz_0[j] - 1.5 * fl1_fx * ta_y_xxz_1[j];

                ta_xy_xxyy_0[j] = pa_x[j] * ta_y_xxyy_0[j] - pc_x[j] * ta_y_xxyy_1[j] + fl1_fx * ta_y_xyy_0[j] - fl1_fx * ta_y_xyy_1[j];

                ta_xy_xxyz_0[j] = pa_x[j] * ta_y_xxyz_0[j] - pc_x[j] * ta_y_xxyz_1[j] + fl1_fx * ta_y_xyz_0[j] - fl1_fx * ta_y_xyz_1[j];

                ta_xy_xxzz_0[j] = pa_x[j] * ta_y_xxzz_0[j] - pc_x[j] * ta_y_xxzz_1[j] + fl1_fx * ta_y_xzz_0[j] - fl1_fx * ta_y_xzz_1[j];

                ta_xy_xyyy_0[j] = pa_x[j] * ta_y_xyyy_0[j] - pc_x[j] * ta_y_xyyy_1[j] + 0.5 * fl1_fx * ta_y_yyy_0[j] - 0.5 * fl1_fx * ta_y_yyy_1[j];

                ta_xy_xyyz_0[j] = pa_x[j] * ta_y_xyyz_0[j] - pc_x[j] * ta_y_xyyz_1[j] + 0.5 * fl1_fx * ta_y_yyz_0[j] - 0.5 * fl1_fx * ta_y_yyz_1[j];

                ta_xy_xyzz_0[j] = pa_x[j] * ta_y_xyzz_0[j] - pc_x[j] * ta_y_xyzz_1[j] + 0.5 * fl1_fx * ta_y_yzz_0[j] - 0.5 * fl1_fx * ta_y_yzz_1[j];

                ta_xy_xzzz_0[j] = pa_x[j] * ta_y_xzzz_0[j] - pc_x[j] * ta_y_xzzz_1[j] + 0.5 * fl1_fx * ta_y_zzz_0[j] - 0.5 * fl1_fx * ta_y_zzz_1[j];

                ta_xy_yyyy_0[j] = pa_x[j] * ta_y_yyyy_0[j] - pc_x[j] * ta_y_yyyy_1[j];

                ta_xy_yyyz_0[j] = pa_x[j] * ta_y_yyyz_0[j] - pc_x[j] * ta_y_yyyz_1[j];

                ta_xy_yyzz_0[j] = pa_x[j] * ta_y_yyzz_0[j] - pc_x[j] * ta_y_yyzz_1[j];

                ta_xy_yzzz_0[j] = pa_x[j] * ta_y_yzzz_0[j] - pc_x[j] * ta_y_yzzz_1[j];

                ta_xy_zzzz_0[j] = pa_x[j] * ta_y_zzzz_0[j] - pc_x[j] * ta_y_zzzz_1[j];

                ta_xz_xxxx_0[j] = pa_x[j] * ta_z_xxxx_0[j] - pc_x[j] * ta_z_xxxx_1[j] + 2.0 * fl1_fx * ta_z_xxx_0[j] - 2.0 * fl1_fx * ta_z_xxx_1[j];

                ta_xz_xxxy_0[j] = pa_x[j] * ta_z_xxxy_0[j] - pc_x[j] * ta_z_xxxy_1[j] + 1.5 * fl1_fx * ta_z_xxy_0[j] - 1.5 * fl1_fx * ta_z_xxy_1[j];

                ta_xz_xxxz_0[j] = pa_x[j] * ta_z_xxxz_0[j] - pc_x[j] * ta_z_xxxz_1[j] + 1.5 * fl1_fx * ta_z_xxz_0[j] - 1.5 * fl1_fx * ta_z_xxz_1[j];

                ta_xz_xxyy_0[j] = pa_x[j] * ta_z_xxyy_0[j] - pc_x[j] * ta_z_xxyy_1[j] + fl1_fx * ta_z_xyy_0[j] - fl1_fx * ta_z_xyy_1[j];

                ta_xz_xxyz_0[j] = pa_x[j] * ta_z_xxyz_0[j] - pc_x[j] * ta_z_xxyz_1[j] + fl1_fx * ta_z_xyz_0[j] - fl1_fx * ta_z_xyz_1[j];

                ta_xz_xxzz_0[j] = pa_x[j] * ta_z_xxzz_0[j] - pc_x[j] * ta_z_xxzz_1[j] + fl1_fx * ta_z_xzz_0[j] - fl1_fx * ta_z_xzz_1[j];

                ta_xz_xyyy_0[j] = pa_x[j] * ta_z_xyyy_0[j] - pc_x[j] * ta_z_xyyy_1[j] + 0.5 * fl1_fx * ta_z_yyy_0[j] - 0.5 * fl1_fx * ta_z_yyy_1[j];

                ta_xz_xyyz_0[j] = pa_x[j] * ta_z_xyyz_0[j] - pc_x[j] * ta_z_xyyz_1[j] + 0.5 * fl1_fx * ta_z_yyz_0[j] - 0.5 * fl1_fx * ta_z_yyz_1[j];

                ta_xz_xyzz_0[j] = pa_x[j] * ta_z_xyzz_0[j] - pc_x[j] * ta_z_xyzz_1[j] + 0.5 * fl1_fx * ta_z_yzz_0[j] - 0.5 * fl1_fx * ta_z_yzz_1[j];

                ta_xz_xzzz_0[j] = pa_x[j] * ta_z_xzzz_0[j] - pc_x[j] * ta_z_xzzz_1[j] + 0.5 * fl1_fx * ta_z_zzz_0[j] - 0.5 * fl1_fx * ta_z_zzz_1[j];

                ta_xz_yyyy_0[j] = pa_x[j] * ta_z_yyyy_0[j] - pc_x[j] * ta_z_yyyy_1[j];

                ta_xz_yyyz_0[j] = pa_x[j] * ta_z_yyyz_0[j] - pc_x[j] * ta_z_yyyz_1[j];

                ta_xz_yyzz_0[j] = pa_x[j] * ta_z_yyzz_0[j] - pc_x[j] * ta_z_yyzz_1[j];

                ta_xz_yzzz_0[j] = pa_x[j] * ta_z_yzzz_0[j] - pc_x[j] * ta_z_yzzz_1[j];

                ta_xz_zzzz_0[j] = pa_x[j] * ta_z_zzzz_0[j] - pc_x[j] * ta_z_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForDG_45_90(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_0_xxxx_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx);

            auto ta_0_xxxy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 1);

            auto ta_0_xxxz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 2);

            auto ta_0_xxyy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 3);

            auto ta_0_xxyz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 4);

            auto ta_0_xxzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 5);

            auto ta_0_xyyy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 6);

            auto ta_0_xyyz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 7);

            auto ta_0_xyzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 8);

            auto ta_0_xzzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 9);

            auto ta_0_yyyy_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 10);

            auto ta_0_yyyz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 11);

            auto ta_0_yyzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 12);

            auto ta_0_yzzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 13);

            auto ta_0_zzzz_0 = primBuffer.data(pidx_a_0_4_m0 + 15 * idx + 14);

            auto ta_0_xxxx_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx);

            auto ta_0_xxxy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 1);

            auto ta_0_xxxz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 2);

            auto ta_0_xxyy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 3);

            auto ta_0_xxyz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 4);

            auto ta_0_xxzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 5);

            auto ta_0_xyyy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 6);

            auto ta_0_xyyz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 7);

            auto ta_0_xyzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 8);

            auto ta_0_xzzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 9);

            auto ta_0_yyyy_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 10);

            auto ta_0_yyyz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 11);

            auto ta_0_yyzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 12);

            auto ta_0_yzzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 13);

            auto ta_0_zzzz_1 = primBuffer.data(pidx_a_0_4_m1 + 15 * idx + 14);

            auto ta_y_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 10);

            auto ta_y_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 11);

            auto ta_y_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 12);

            auto ta_y_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 13);

            auto ta_y_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 14);

            auto ta_y_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 15);

            auto ta_y_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 16);

            auto ta_y_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 17);

            auto ta_y_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 18);

            auto ta_y_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 19);

            auto ta_z_xxx_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 20);

            auto ta_z_xxy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 21);

            auto ta_z_xxz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 22);

            auto ta_z_xyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 23);

            auto ta_z_xyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 24);

            auto ta_z_xzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 25);

            auto ta_z_yyy_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 26);

            auto ta_z_yyz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 27);

            auto ta_z_yzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 28);

            auto ta_z_zzz_0 = primBuffer.data(pidx_a_1_3_m0 + 30 * idx + 29);

            auto ta_y_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 10);

            auto ta_y_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 11);

            auto ta_y_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 12);

            auto ta_y_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 13);

            auto ta_y_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 14);

            auto ta_y_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 15);

            auto ta_y_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 16);

            auto ta_y_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 17);

            auto ta_y_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 18);

            auto ta_y_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 19);

            auto ta_z_xxx_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 20);

            auto ta_z_xxy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 21);

            auto ta_z_xxz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 22);

            auto ta_z_xyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 23);

            auto ta_z_xyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 24);

            auto ta_z_xzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 25);

            auto ta_z_yyy_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 26);

            auto ta_z_yyz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 27);

            auto ta_z_yzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 28);

            auto ta_z_zzz_1 = primBuffer.data(pidx_a_1_3_m1 + 30 * idx + 29);

            // set up pointers to integrals

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

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_0_xxxx_0, ta_0_xxxx_1, ta_0_xxxy_0, \
                                         ta_0_xxxy_1, ta_0_xxxz_0, ta_0_xxxz_1, ta_0_xxyy_0, ta_0_xxyy_1, ta_0_xxyz_0, \
                                         ta_0_xxyz_1, ta_0_xxzz_0, ta_0_xxzz_1, ta_0_xyyy_0, ta_0_xyyy_1, ta_0_xyyz_0, \
                                         ta_0_xyyz_1, ta_0_xyzz_0, ta_0_xyzz_1, ta_0_xzzz_0, ta_0_xzzz_1, ta_0_yyyy_0, \
                                         ta_0_yyyy_1, ta_0_yyyz_0, ta_0_yyyz_1, ta_0_yyzz_0, ta_0_yyzz_1, ta_0_yzzz_0, \
                                         ta_0_yzzz_1, ta_0_zzzz_0, ta_0_zzzz_1, ta_y_xxx_0, ta_y_xxx_1, ta_y_xxxx_0, \
                                         ta_y_xxxx_1, ta_y_xxxy_0, ta_y_xxxy_1, ta_y_xxxz_0, ta_y_xxxz_1, ta_y_xxy_0, \
                                         ta_y_xxy_1, ta_y_xxyy_0, ta_y_xxyy_1, ta_y_xxyz_0, ta_y_xxyz_1, ta_y_xxz_0, \
                                         ta_y_xxz_1, ta_y_xxzz_0, ta_y_xxzz_1, ta_y_xyy_0, ta_y_xyy_1, ta_y_xyyy_0, \
                                         ta_y_xyyy_1, ta_y_xyyz_0, ta_y_xyyz_1, ta_y_xyz_0, ta_y_xyz_1, ta_y_xyzz_0, \
                                         ta_y_xyzz_1, ta_y_xzz_0, ta_y_xzz_1, ta_y_xzzz_0, ta_y_xzzz_1, ta_y_yyy_0, \
                                         ta_y_yyy_1, ta_y_yyyy_0, ta_y_yyyy_1, ta_y_yyyz_0, ta_y_yyyz_1, ta_y_yyz_0, \
                                         ta_y_yyz_1, ta_y_yyzz_0, ta_y_yyzz_1, ta_y_yzz_0, ta_y_yzz_1, ta_y_yzzz_0, \
                                         ta_y_yzzz_1, ta_y_zzz_0, ta_y_zzz_1, ta_y_zzzz_0, ta_y_zzzz_1, ta_yy_xxxx_0, \
                                         ta_yy_xxxy_0, ta_yy_xxxz_0, ta_yy_xxyy_0, ta_yy_xxyz_0, ta_yy_xxzz_0, ta_yy_xyyy_0, \
                                         ta_yy_xyyz_0, ta_yy_xyzz_0, ta_yy_xzzz_0, ta_yy_yyyy_0, ta_yy_yyyz_0, ta_yy_yyzz_0, \
                                         ta_yy_yzzz_0, ta_yy_zzzz_0, ta_yz_xxxx_0, ta_yz_xxxy_0, ta_yz_xxxz_0, ta_yz_xxyy_0, \
                                         ta_yz_xxyz_0, ta_yz_xxzz_0, ta_yz_xyyy_0, ta_yz_xyyz_0, ta_yz_xyzz_0, ta_yz_xzzz_0, \
                                         ta_yz_yyyy_0, ta_yz_yyyz_0, ta_yz_yyzz_0, ta_yz_yzzz_0, ta_yz_zzzz_0, ta_z_xxx_0, \
                                         ta_z_xxx_1, ta_z_xxxx_0, ta_z_xxxx_1, ta_z_xxxy_0, ta_z_xxxy_1, ta_z_xxxz_0, \
                                         ta_z_xxxz_1, ta_z_xxy_0, ta_z_xxy_1, ta_z_xxyy_0, ta_z_xxyy_1, ta_z_xxyz_0, \
                                         ta_z_xxyz_1, ta_z_xxz_0, ta_z_xxz_1, ta_z_xxzz_0, ta_z_xxzz_1, ta_z_xyy_0, \
                                         ta_z_xyy_1, ta_z_xyyy_0, ta_z_xyyy_1, ta_z_xyyz_0, ta_z_xyyz_1, ta_z_xyz_0, \
                                         ta_z_xyz_1, ta_z_xyzz_0, ta_z_xyzz_1, ta_z_xzz_0, ta_z_xzz_1, ta_z_xzzz_0, \
                                         ta_z_xzzz_1, ta_z_yyy_0, ta_z_yyy_1, ta_z_yyyy_0, ta_z_yyyy_1, ta_z_yyyz_0, \
                                         ta_z_yyyz_1, ta_z_yyz_0, ta_z_yyz_1, ta_z_yyzz_0, ta_z_yyzz_1, ta_z_yzz_0, \
                                         ta_z_yzz_1, ta_z_yzzz_0, ta_z_yzzz_1, ta_z_zzz_0, ta_z_zzz_1, ta_z_zzzz_0, \
                                         ta_z_zzzz_1, ta_zz_xxxx_0, ta_zz_xxxy_0, ta_zz_xxxz_0, ta_zz_xxyy_0, ta_zz_xxyz_0, \
                                         ta_zz_xxzz_0, ta_zz_xyyy_0, ta_zz_xyyz_0, ta_zz_xyzz_0, ta_zz_xzzz_0, ta_zz_yyyy_0, \
                                         ta_zz_yyyz_0, ta_zz_yyzz_0, ta_zz_yzzz_0, ta_zz_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_yy_xxxx_0[j] = pa_y[j] * ta_y_xxxx_0[j] - pc_y[j] * ta_y_xxxx_1[j] + 0.5 * fl1_fx * ta_0_xxxx_0[j] - 0.5 * fl1_fx * ta_0_xxxx_1[j];

                ta_yy_xxxy_0[j] = pa_y[j] * ta_y_xxxy_0[j] - pc_y[j] * ta_y_xxxy_1[j] + 0.5 * fl1_fx * ta_0_xxxy_0[j] -
                                  0.5 * fl1_fx * ta_0_xxxy_1[j] + 0.5 * fl1_fx * ta_y_xxx_0[j] - 0.5 * fl1_fx * ta_y_xxx_1[j];

                ta_yy_xxxz_0[j] = pa_y[j] * ta_y_xxxz_0[j] - pc_y[j] * ta_y_xxxz_1[j] + 0.5 * fl1_fx * ta_0_xxxz_0[j] - 0.5 * fl1_fx * ta_0_xxxz_1[j];

                ta_yy_xxyy_0[j] = pa_y[j] * ta_y_xxyy_0[j] - pc_y[j] * ta_y_xxyy_1[j] + 0.5 * fl1_fx * ta_0_xxyy_0[j] -
                                  0.5 * fl1_fx * ta_0_xxyy_1[j] + fl1_fx * ta_y_xxy_0[j] - fl1_fx * ta_y_xxy_1[j];

                ta_yy_xxyz_0[j] = pa_y[j] * ta_y_xxyz_0[j] - pc_y[j] * ta_y_xxyz_1[j] + 0.5 * fl1_fx * ta_0_xxyz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxyz_1[j] + 0.5 * fl1_fx * ta_y_xxz_0[j] - 0.5 * fl1_fx * ta_y_xxz_1[j];

                ta_yy_xxzz_0[j] = pa_y[j] * ta_y_xxzz_0[j] - pc_y[j] * ta_y_xxzz_1[j] + 0.5 * fl1_fx * ta_0_xxzz_0[j] - 0.5 * fl1_fx * ta_0_xxzz_1[j];

                ta_yy_xyyy_0[j] = pa_y[j] * ta_y_xyyy_0[j] - pc_y[j] * ta_y_xyyy_1[j] + 0.5 * fl1_fx * ta_0_xyyy_0[j] -
                                  0.5 * fl1_fx * ta_0_xyyy_1[j] + 1.5 * fl1_fx * ta_y_xyy_0[j] - 1.5 * fl1_fx * ta_y_xyy_1[j];

                ta_yy_xyyz_0[j] = pa_y[j] * ta_y_xyyz_0[j] - pc_y[j] * ta_y_xyyz_1[j] + 0.5 * fl1_fx * ta_0_xyyz_0[j] -
                                  0.5 * fl1_fx * ta_0_xyyz_1[j] + fl1_fx * ta_y_xyz_0[j] - fl1_fx * ta_y_xyz_1[j];

                ta_yy_xyzz_0[j] = pa_y[j] * ta_y_xyzz_0[j] - pc_y[j] * ta_y_xyzz_1[j] + 0.5 * fl1_fx * ta_0_xyzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xyzz_1[j] + 0.5 * fl1_fx * ta_y_xzz_0[j] - 0.5 * fl1_fx * ta_y_xzz_1[j];

                ta_yy_xzzz_0[j] = pa_y[j] * ta_y_xzzz_0[j] - pc_y[j] * ta_y_xzzz_1[j] + 0.5 * fl1_fx * ta_0_xzzz_0[j] - 0.5 * fl1_fx * ta_0_xzzz_1[j];

                ta_yy_yyyy_0[j] = pa_y[j] * ta_y_yyyy_0[j] - pc_y[j] * ta_y_yyyy_1[j] + 0.5 * fl1_fx * ta_0_yyyy_0[j] -
                                  0.5 * fl1_fx * ta_0_yyyy_1[j] + 2.0 * fl1_fx * ta_y_yyy_0[j] - 2.0 * fl1_fx * ta_y_yyy_1[j];

                ta_yy_yyyz_0[j] = pa_y[j] * ta_y_yyyz_0[j] - pc_y[j] * ta_y_yyyz_1[j] + 0.5 * fl1_fx * ta_0_yyyz_0[j] -
                                  0.5 * fl1_fx * ta_0_yyyz_1[j] + 1.5 * fl1_fx * ta_y_yyz_0[j] - 1.5 * fl1_fx * ta_y_yyz_1[j];

                ta_yy_yyzz_0[j] = pa_y[j] * ta_y_yyzz_0[j] - pc_y[j] * ta_y_yyzz_1[j] + 0.5 * fl1_fx * ta_0_yyzz_0[j] -
                                  0.5 * fl1_fx * ta_0_yyzz_1[j] + fl1_fx * ta_y_yzz_0[j] - fl1_fx * ta_y_yzz_1[j];

                ta_yy_yzzz_0[j] = pa_y[j] * ta_y_yzzz_0[j] - pc_y[j] * ta_y_yzzz_1[j] + 0.5 * fl1_fx * ta_0_yzzz_0[j] -
                                  0.5 * fl1_fx * ta_0_yzzz_1[j] + 0.5 * fl1_fx * ta_y_zzz_0[j] - 0.5 * fl1_fx * ta_y_zzz_1[j];

                ta_yy_zzzz_0[j] = pa_y[j] * ta_y_zzzz_0[j] - pc_y[j] * ta_y_zzzz_1[j] + 0.5 * fl1_fx * ta_0_zzzz_0[j] - 0.5 * fl1_fx * ta_0_zzzz_1[j];

                ta_yz_xxxx_0[j] = pa_y[j] * ta_z_xxxx_0[j] - pc_y[j] * ta_z_xxxx_1[j];

                ta_yz_xxxy_0[j] = pa_y[j] * ta_z_xxxy_0[j] - pc_y[j] * ta_z_xxxy_1[j] + 0.5 * fl1_fx * ta_z_xxx_0[j] - 0.5 * fl1_fx * ta_z_xxx_1[j];

                ta_yz_xxxz_0[j] = pa_y[j] * ta_z_xxxz_0[j] - pc_y[j] * ta_z_xxxz_1[j];

                ta_yz_xxyy_0[j] = pa_y[j] * ta_z_xxyy_0[j] - pc_y[j] * ta_z_xxyy_1[j] + fl1_fx * ta_z_xxy_0[j] - fl1_fx * ta_z_xxy_1[j];

                ta_yz_xxyz_0[j] = pa_y[j] * ta_z_xxyz_0[j] - pc_y[j] * ta_z_xxyz_1[j] + 0.5 * fl1_fx * ta_z_xxz_0[j] - 0.5 * fl1_fx * ta_z_xxz_1[j];

                ta_yz_xxzz_0[j] = pa_y[j] * ta_z_xxzz_0[j] - pc_y[j] * ta_z_xxzz_1[j];

                ta_yz_xyyy_0[j] = pa_y[j] * ta_z_xyyy_0[j] - pc_y[j] * ta_z_xyyy_1[j] + 1.5 * fl1_fx * ta_z_xyy_0[j] - 1.5 * fl1_fx * ta_z_xyy_1[j];

                ta_yz_xyyz_0[j] = pa_y[j] * ta_z_xyyz_0[j] - pc_y[j] * ta_z_xyyz_1[j] + fl1_fx * ta_z_xyz_0[j] - fl1_fx * ta_z_xyz_1[j];

                ta_yz_xyzz_0[j] = pa_y[j] * ta_z_xyzz_0[j] - pc_y[j] * ta_z_xyzz_1[j] + 0.5 * fl1_fx * ta_z_xzz_0[j] - 0.5 * fl1_fx * ta_z_xzz_1[j];

                ta_yz_xzzz_0[j] = pa_y[j] * ta_z_xzzz_0[j] - pc_y[j] * ta_z_xzzz_1[j];

                ta_yz_yyyy_0[j] = pa_y[j] * ta_z_yyyy_0[j] - pc_y[j] * ta_z_yyyy_1[j] + 2.0 * fl1_fx * ta_z_yyy_0[j] - 2.0 * fl1_fx * ta_z_yyy_1[j];

                ta_yz_yyyz_0[j] = pa_y[j] * ta_z_yyyz_0[j] - pc_y[j] * ta_z_yyyz_1[j] + 1.5 * fl1_fx * ta_z_yyz_0[j] - 1.5 * fl1_fx * ta_z_yyz_1[j];

                ta_yz_yyzz_0[j] = pa_y[j] * ta_z_yyzz_0[j] - pc_y[j] * ta_z_yyzz_1[j] + fl1_fx * ta_z_yzz_0[j] - fl1_fx * ta_z_yzz_1[j];

                ta_yz_yzzz_0[j] = pa_y[j] * ta_z_yzzz_0[j] - pc_y[j] * ta_z_yzzz_1[j] + 0.5 * fl1_fx * ta_z_zzz_0[j] - 0.5 * fl1_fx * ta_z_zzz_1[j];

                ta_yz_zzzz_0[j] = pa_y[j] * ta_z_zzzz_0[j] - pc_y[j] * ta_z_zzzz_1[j];

                ta_zz_xxxx_0[j] = pa_z[j] * ta_z_xxxx_0[j] - pc_z[j] * ta_z_xxxx_1[j] + 0.5 * fl1_fx * ta_0_xxxx_0[j] - 0.5 * fl1_fx * ta_0_xxxx_1[j];

                ta_zz_xxxy_0[j] = pa_z[j] * ta_z_xxxy_0[j] - pc_z[j] * ta_z_xxxy_1[j] + 0.5 * fl1_fx * ta_0_xxxy_0[j] - 0.5 * fl1_fx * ta_0_xxxy_1[j];

                ta_zz_xxxz_0[j] = pa_z[j] * ta_z_xxxz_0[j] - pc_z[j] * ta_z_xxxz_1[j] + 0.5 * fl1_fx * ta_0_xxxz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxxz_1[j] + 0.5 * fl1_fx * ta_z_xxx_0[j] - 0.5 * fl1_fx * ta_z_xxx_1[j];

                ta_zz_xxyy_0[j] = pa_z[j] * ta_z_xxyy_0[j] - pc_z[j] * ta_z_xxyy_1[j] + 0.5 * fl1_fx * ta_0_xxyy_0[j] - 0.5 * fl1_fx * ta_0_xxyy_1[j];

                ta_zz_xxyz_0[j] = pa_z[j] * ta_z_xxyz_0[j] - pc_z[j] * ta_z_xxyz_1[j] + 0.5 * fl1_fx * ta_0_xxyz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxyz_1[j] + 0.5 * fl1_fx * ta_z_xxy_0[j] - 0.5 * fl1_fx * ta_z_xxy_1[j];

                ta_zz_xxzz_0[j] = pa_z[j] * ta_z_xxzz_0[j] - pc_z[j] * ta_z_xxzz_1[j] + 0.5 * fl1_fx * ta_0_xxzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xxzz_1[j] + fl1_fx * ta_z_xxz_0[j] - fl1_fx * ta_z_xxz_1[j];

                ta_zz_xyyy_0[j] = pa_z[j] * ta_z_xyyy_0[j] - pc_z[j] * ta_z_xyyy_1[j] + 0.5 * fl1_fx * ta_0_xyyy_0[j] - 0.5 * fl1_fx * ta_0_xyyy_1[j];

                ta_zz_xyyz_0[j] = pa_z[j] * ta_z_xyyz_0[j] - pc_z[j] * ta_z_xyyz_1[j] + 0.5 * fl1_fx * ta_0_xyyz_0[j] -
                                  0.5 * fl1_fx * ta_0_xyyz_1[j] + 0.5 * fl1_fx * ta_z_xyy_0[j] - 0.5 * fl1_fx * ta_z_xyy_1[j];

                ta_zz_xyzz_0[j] = pa_z[j] * ta_z_xyzz_0[j] - pc_z[j] * ta_z_xyzz_1[j] + 0.5 * fl1_fx * ta_0_xyzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xyzz_1[j] + fl1_fx * ta_z_xyz_0[j] - fl1_fx * ta_z_xyz_1[j];

                ta_zz_xzzz_0[j] = pa_z[j] * ta_z_xzzz_0[j] - pc_z[j] * ta_z_xzzz_1[j] + 0.5 * fl1_fx * ta_0_xzzz_0[j] -
                                  0.5 * fl1_fx * ta_0_xzzz_1[j] + 1.5 * fl1_fx * ta_z_xzz_0[j] - 1.5 * fl1_fx * ta_z_xzz_1[j];

                ta_zz_yyyy_0[j] = pa_z[j] * ta_z_yyyy_0[j] - pc_z[j] * ta_z_yyyy_1[j] + 0.5 * fl1_fx * ta_0_yyyy_0[j] - 0.5 * fl1_fx * ta_0_yyyy_1[j];

                ta_zz_yyyz_0[j] = pa_z[j] * ta_z_yyyz_0[j] - pc_z[j] * ta_z_yyyz_1[j] + 0.5 * fl1_fx * ta_0_yyyz_0[j] -
                                  0.5 * fl1_fx * ta_0_yyyz_1[j] + 0.5 * fl1_fx * ta_z_yyy_0[j] - 0.5 * fl1_fx * ta_z_yyy_1[j];

                ta_zz_yyzz_0[j] = pa_z[j] * ta_z_yyzz_0[j] - pc_z[j] * ta_z_yyzz_1[j] + 0.5 * fl1_fx * ta_0_yyzz_0[j] -
                                  0.5 * fl1_fx * ta_0_yyzz_1[j] + fl1_fx * ta_z_yyz_0[j] - fl1_fx * ta_z_yyz_1[j];

                ta_zz_yzzz_0[j] = pa_z[j] * ta_z_yzzz_0[j] - pc_z[j] * ta_z_yzzz_1[j] + 0.5 * fl1_fx * ta_0_yzzz_0[j] -
                                  0.5 * fl1_fx * ta_0_yzzz_1[j] + 1.5 * fl1_fx * ta_z_yzz_0[j] - 1.5 * fl1_fx * ta_z_yzz_1[j];

                ta_zz_zzzz_0[j] = pa_z[j] * ta_z_zzzz_0[j] - pc_z[j] * ta_z_zzzz_1[j] + 0.5 * fl1_fx * ta_0_zzzz_0[j] -
                                  0.5 * fl1_fx * ta_0_zzzz_1[j] + 2.0 * fl1_fx * ta_z_zzz_0[j] - 2.0 * fl1_fx * ta_z_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGD(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& paDistances,
                          const CMemBlock2D<double>& pcDistances,
                          const CGtoBlock&           braGtoBlock,
                          const CGtoBlock&           ketGtoBlock,
                          const int32_t              iContrGto)
{
    npotrecfunc::compNuclearPotentialForGD_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForGD_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compNuclearPotentialForGD_0_45(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
                               const CMemBlock2D<double>& pcDistances,
                               const CGtoBlock&           braGtoBlock,
                               const CGtoBlock&           ketGtoBlock,
                               const int32_t              iContrGto)
{
    // Batch of Integrals (0,45)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_xxx_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx);

            auto ta_xxx_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 1);

            auto ta_xxx_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 2);

            auto ta_xxx_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 3);

            auto ta_xxx_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 4);

            auto ta_xxx_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 5);

            auto ta_xxy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 6);

            auto ta_xxy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 7);

            auto ta_xxy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 8);

            auto ta_xxy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 9);

            auto ta_xxy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 10);

            auto ta_xxy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 11);

            auto ta_xxz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 12);

            auto ta_xxz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 13);

            auto ta_xxz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 14);

            auto ta_xxz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 15);

            auto ta_xxz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 16);

            auto ta_xxz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 17);

            auto ta_xyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 18);

            auto ta_xyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 19);

            auto ta_xyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 20);

            auto ta_xyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 21);

            auto ta_xyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 22);

            auto ta_xyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 23);

            auto ta_xyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 24);

            auto ta_xyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 25);

            auto ta_xyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 26);

            auto ta_xyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 27);

            auto ta_xyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 28);

            auto ta_xyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 29);

            auto ta_xzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 30);

            auto ta_xzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 31);

            auto ta_xzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 32);

            auto ta_xzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 33);

            auto ta_xzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 34);

            auto ta_xzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 35);

            auto ta_yyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 36);

            auto ta_yyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 37);

            auto ta_yyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 38);

            auto ta_yyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 39);

            auto ta_yyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 40);

            auto ta_yyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 41);

            auto ta_yyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 42);

            auto ta_yyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 43);

            auto ta_yyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 44);

            auto ta_xxx_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx);

            auto ta_xxx_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 1);

            auto ta_xxx_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 2);

            auto ta_xxx_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 3);

            auto ta_xxx_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 4);

            auto ta_xxx_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 5);

            auto ta_xxy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 6);

            auto ta_xxy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 7);

            auto ta_xxy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 8);

            auto ta_xxy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 9);

            auto ta_xxy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 10);

            auto ta_xxy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 11);

            auto ta_xxz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 12);

            auto ta_xxz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 13);

            auto ta_xxz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 14);

            auto ta_xxz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 15);

            auto ta_xxz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 16);

            auto ta_xxz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 17);

            auto ta_xyy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 18);

            auto ta_xyy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 19);

            auto ta_xyy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 20);

            auto ta_xyy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 21);

            auto ta_xyy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 22);

            auto ta_xyy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 23);

            auto ta_xyz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 24);

            auto ta_xyz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 25);

            auto ta_xyz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 26);

            auto ta_xyz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 27);

            auto ta_xyz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 28);

            auto ta_xyz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 29);

            auto ta_xzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 30);

            auto ta_xzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 31);

            auto ta_xzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 32);

            auto ta_xzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 33);

            auto ta_xzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 34);

            auto ta_xzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 35);

            auto ta_yyy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 36);

            auto ta_yyy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 37);

            auto ta_yyy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 38);

            auto ta_yyy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 39);

            auto ta_yyy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 40);

            auto ta_yyy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 41);

            auto ta_yyz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 42);

            auto ta_yyz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 43);

            auto ta_yyz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 44);

            auto ta_xx_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx);

            auto ta_xx_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 1);

            auto ta_xx_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 2);

            auto ta_xx_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 3);

            auto ta_xx_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 4);

            auto ta_xx_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 5);

            auto ta_xy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 6);

            auto ta_xy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 7);

            auto ta_xy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 8);

            auto ta_xy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 9);

            auto ta_xy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 10);

            auto ta_xy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 11);

            auto ta_xz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 12);

            auto ta_xz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 13);

            auto ta_xz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 14);

            auto ta_xz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 15);

            auto ta_xz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 16);

            auto ta_xz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 17);

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_zz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 30);

            auto ta_zz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 31);

            auto ta_zz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 32);

            auto ta_zz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 33);

            auto ta_zz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 34);

            auto ta_zz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 35);

            auto ta_xx_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx);

            auto ta_xx_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 1);

            auto ta_xx_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 2);

            auto ta_xx_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 3);

            auto ta_xx_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 4);

            auto ta_xx_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 5);

            auto ta_xy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 6);

            auto ta_xy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 7);

            auto ta_xy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 8);

            auto ta_xy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 9);

            auto ta_xy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 10);

            auto ta_xy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 11);

            auto ta_xz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 12);

            auto ta_xz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 13);

            auto ta_xz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 14);

            auto ta_xz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 15);

            auto ta_xz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 16);

            auto ta_xz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 17);

            auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18);

            auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19);

            auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20);

            auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21);

            auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22);

            auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23);

            auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24);

            auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25);

            auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26);

            auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27);

            auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28);

            auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29);

            auto ta_zz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 30);

            auto ta_zz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 31);

            auto ta_zz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 32);

            auto ta_zz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 33);

            auto ta_zz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 34);

            auto ta_zz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 35);

            auto ta_xxx_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx);

            auto ta_xxx_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 1);

            auto ta_xxx_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 2);

            auto ta_xxy_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 3);

            auto ta_xxy_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 4);

            auto ta_xxy_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 5);

            auto ta_xxz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 6);

            auto ta_xxz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 7);

            auto ta_xxz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 8);

            auto ta_xyy_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 9);

            auto ta_xyy_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 10);

            auto ta_xyy_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 11);

            auto ta_xyz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 12);

            auto ta_xyz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 13);

            auto ta_xyz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 14);

            auto ta_xzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 15);

            auto ta_xzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 16);

            auto ta_xzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 17);

            auto ta_yyy_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 18);

            auto ta_yyy_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 19);

            auto ta_yyy_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 20);

            auto ta_yyz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 21);

            auto ta_yyz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 22);

            auto ta_yyz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 23);

            auto ta_xxx_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx);

            auto ta_xxx_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 1);

            auto ta_xxx_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 2);

            auto ta_xxy_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 3);

            auto ta_xxy_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 4);

            auto ta_xxy_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 5);

            auto ta_xxz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 6);

            auto ta_xxz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 7);

            auto ta_xxz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 8);

            auto ta_xyy_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 9);

            auto ta_xyy_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 10);

            auto ta_xyy_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 11);

            auto ta_xyz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 12);

            auto ta_xyz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 13);

            auto ta_xyz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 14);

            auto ta_xzz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 15);

            auto ta_xzz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 16);

            auto ta_xzz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 17);

            auto ta_yyy_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 18);

            auto ta_yyy_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 19);

            auto ta_yyy_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 20);

            auto ta_yyz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 21);

            auto ta_yyz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 22);

            auto ta_yyz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 23);

            // set up pointers to integrals

            auto ta_xxxx_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx);

            auto ta_xxxx_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 1);

            auto ta_xxxx_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 2);

            auto ta_xxxx_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 3);

            auto ta_xxxx_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 4);

            auto ta_xxxx_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 5);

            auto ta_xxxy_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 6);

            auto ta_xxxy_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 7);

            auto ta_xxxy_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 8);

            auto ta_xxxy_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 9);

            auto ta_xxxy_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 10);

            auto ta_xxxy_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 11);

            auto ta_xxxz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 12);

            auto ta_xxxz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 13);

            auto ta_xxxz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 14);

            auto ta_xxxz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 15);

            auto ta_xxxz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 16);

            auto ta_xxxz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 17);

            auto ta_xxyy_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 18);

            auto ta_xxyy_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 19);

            auto ta_xxyy_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 20);

            auto ta_xxyy_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 21);

            auto ta_xxyy_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 22);

            auto ta_xxyy_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 23);

            auto ta_xxyz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 24);

            auto ta_xxyz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 25);

            auto ta_xxyz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 26);

            auto ta_xxyz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 27);

            auto ta_xxyz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 28);

            auto ta_xxyz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 29);

            auto ta_xxzz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 30);

            auto ta_xxzz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 31);

            auto ta_xxzz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 32);

            auto ta_xxzz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 33);

            auto ta_xxzz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 34);

            auto ta_xxzz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 35);

            auto ta_xyyy_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 36);

            auto ta_xyyy_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 37);

            auto ta_xyyy_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 38);

            auto ta_xyyy_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 39);

            auto ta_xyyy_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 40);

            auto ta_xyyy_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 41);

            auto ta_xyyz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 42);

            auto ta_xyyz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 43);

            auto ta_xyyz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 44);

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xx_xx_0, ta_xx_xx_1, ta_xx_xy_0, ta_xx_xy_1, ta_xx_xz_0, \
                                         ta_xx_xz_1, ta_xx_yy_0, ta_xx_yy_1, ta_xx_yz_0, ta_xx_yz_1, ta_xx_zz_0, ta_xx_zz_1, \
                                         ta_xxx_x_0, ta_xxx_x_1, ta_xxx_xx_0, ta_xxx_xx_1, ta_xxx_xy_0, ta_xxx_xy_1, \
                                         ta_xxx_xz_0, ta_xxx_xz_1, ta_xxx_y_0, ta_xxx_y_1, ta_xxx_yy_0, ta_xxx_yy_1, \
                                         ta_xxx_yz_0, ta_xxx_yz_1, ta_xxx_z_0, ta_xxx_z_1, ta_xxx_zz_0, ta_xxx_zz_1, \
                                         ta_xxxx_xx_0, ta_xxxx_xy_0, ta_xxxx_xz_0, ta_xxxx_yy_0, ta_xxxx_yz_0, ta_xxxx_zz_0, \
                                         ta_xxxy_xx_0, ta_xxxy_xy_0, ta_xxxy_xz_0, ta_xxxy_yy_0, ta_xxxy_yz_0, ta_xxxy_zz_0, \
                                         ta_xxxz_xx_0, ta_xxxz_xy_0, ta_xxxz_xz_0, ta_xxxz_yy_0, ta_xxxz_yz_0, ta_xxxz_zz_0, \
                                         ta_xxy_x_0, ta_xxy_x_1, ta_xxy_xx_0, ta_xxy_xx_1, ta_xxy_xy_0, ta_xxy_xy_1, \
                                         ta_xxy_xz_0, ta_xxy_xz_1, ta_xxy_y_0, ta_xxy_y_1, ta_xxy_yy_0, ta_xxy_yy_1, \
                                         ta_xxy_yz_0, ta_xxy_yz_1, ta_xxy_z_0, ta_xxy_z_1, ta_xxy_zz_0, ta_xxy_zz_1, \
                                         ta_xxyy_xx_0, ta_xxyy_xy_0, ta_xxyy_xz_0, ta_xxyy_yy_0, ta_xxyy_yz_0, ta_xxyy_zz_0, \
                                         ta_xxyz_xx_0, ta_xxyz_xy_0, ta_xxyz_xz_0, ta_xxyz_yy_0, ta_xxyz_yz_0, ta_xxyz_zz_0, \
                                         ta_xxz_x_0, ta_xxz_x_1, ta_xxz_xx_0, ta_xxz_xx_1, ta_xxz_xy_0, ta_xxz_xy_1, \
                                         ta_xxz_xz_0, ta_xxz_xz_1, ta_xxz_y_0, ta_xxz_y_1, ta_xxz_yy_0, ta_xxz_yy_1, \
                                         ta_xxz_yz_0, ta_xxz_yz_1, ta_xxz_z_0, ta_xxz_z_1, ta_xxz_zz_0, ta_xxz_zz_1, \
                                         ta_xxzz_xx_0, ta_xxzz_xy_0, ta_xxzz_xz_0, ta_xxzz_yy_0, ta_xxzz_yz_0, ta_xxzz_zz_0, \
                                         ta_xy_xx_0, ta_xy_xx_1, ta_xy_xy_0, ta_xy_xy_1, ta_xy_xz_0, ta_xy_xz_1, ta_xy_yy_0, \
                                         ta_xy_yy_1, ta_xy_yz_0, ta_xy_yz_1, ta_xy_zz_0, ta_xy_zz_1, ta_xyy_x_0, ta_xyy_x_1, \
                                         ta_xyy_xx_0, ta_xyy_xx_1, ta_xyy_xy_0, ta_xyy_xy_1, ta_xyy_xz_0, ta_xyy_xz_1, \
                                         ta_xyy_y_0, ta_xyy_y_1, ta_xyy_yy_0, ta_xyy_yy_1, ta_xyy_yz_0, ta_xyy_yz_1, \
                                         ta_xyy_z_0, ta_xyy_z_1, ta_xyy_zz_0, ta_xyy_zz_1, ta_xyyy_xx_0, ta_xyyy_xy_0, \
                                         ta_xyyy_xz_0, ta_xyyy_yy_0, ta_xyyy_yz_0, ta_xyyy_zz_0, ta_xyyz_xx_0, ta_xyyz_xy_0, \
                                         ta_xyyz_xz_0, ta_xyz_x_0, ta_xyz_x_1, ta_xyz_xx_0, ta_xyz_xx_1, ta_xyz_xy_0, \
                                         ta_xyz_xy_1, ta_xyz_xz_0, ta_xyz_xz_1, ta_xyz_y_0, ta_xyz_y_1, ta_xyz_yy_0, \
                                         ta_xyz_yy_1, ta_xyz_yz_0, ta_xyz_yz_1, ta_xyz_z_0, ta_xyz_z_1, ta_xyz_zz_0, \
                                         ta_xyz_zz_1, ta_xz_xx_0, ta_xz_xx_1, ta_xz_xy_0, ta_xz_xy_1, ta_xz_xz_0, ta_xz_xz_1, \
                                         ta_xz_yy_0, ta_xz_yy_1, ta_xz_yz_0, ta_xz_yz_1, ta_xz_zz_0, ta_xz_zz_1, ta_xzz_x_0, \
                                         ta_xzz_x_1, ta_xzz_xx_0, ta_xzz_xx_1, ta_xzz_xy_0, ta_xzz_xy_1, ta_xzz_xz_0, \
                                         ta_xzz_xz_1, ta_xzz_y_0, ta_xzz_y_1, ta_xzz_yy_0, ta_xzz_yy_1, ta_xzz_yz_0, \
                                         ta_xzz_yz_1, ta_xzz_z_0, ta_xzz_z_1, ta_xzz_zz_0, ta_xzz_zz_1, ta_yy_xx_0, \
                                         ta_yy_xx_1, ta_yy_xy_0, ta_yy_xy_1, ta_yy_xz_0, ta_yy_xz_1, ta_yy_yy_0, ta_yy_yy_1, \
                                         ta_yy_yz_0, ta_yy_yz_1, ta_yy_zz_0, ta_yy_zz_1, ta_yyy_x_0, ta_yyy_x_1, \
                                         ta_yyy_xx_0, ta_yyy_xx_1, ta_yyy_xy_0, ta_yyy_xy_1, ta_yyy_xz_0, ta_yyy_xz_1, \
                                         ta_yyy_y_0, ta_yyy_y_1, ta_yyy_yy_0, ta_yyy_yy_1, ta_yyy_yz_0, ta_yyy_yz_1, \
                                         ta_yyy_z_0, ta_yyy_z_1, ta_yyy_zz_0, ta_yyy_zz_1, ta_yyz_x_0, ta_yyz_x_1, \
                                         ta_yyz_xx_0, ta_yyz_xx_1, ta_yyz_xy_0, ta_yyz_xy_1, ta_yyz_xz_0, ta_yyz_xz_1, \
                                         ta_yyz_y_0, ta_yyz_y_1, ta_yyz_z_0, ta_yyz_z_1, ta_yz_xx_0, ta_yz_xx_1, ta_yz_xy_0, \
                                         ta_yz_xy_1, ta_yz_xz_0, ta_yz_xz_1, ta_yz_yy_0, ta_yz_yy_1, ta_yz_yz_0, ta_yz_yz_1, \
                                         ta_yz_zz_0, ta_yz_zz_1, ta_zz_xx_0, ta_zz_xx_1, ta_zz_xy_0, ta_zz_xy_1, ta_zz_xz_0, \
                                         ta_zz_xz_1, ta_zz_yy_0, ta_zz_yy_1, ta_zz_yz_0, ta_zz_yz_1, ta_zz_zz_0, ta_zz_zz_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxxx_xx_0[j] = pa_x[j] * ta_xxx_xx_0[j] - pc_x[j] * ta_xxx_xx_1[j] + 1.5 * fl1_fx * ta_xx_xx_0[j] - 1.5 * fl1_fx * ta_xx_xx_1[j] +
                                  fl1_fx * ta_xxx_x_0[j] - fl1_fx * ta_xxx_x_1[j];

                ta_xxxx_xy_0[j] = pa_x[j] * ta_xxx_xy_0[j] - pc_x[j] * ta_xxx_xy_1[j] + 1.5 * fl1_fx * ta_xx_xy_0[j] - 1.5 * fl1_fx * ta_xx_xy_1[j] +
                                  0.5 * fl1_fx * ta_xxx_y_0[j] - 0.5 * fl1_fx * ta_xxx_y_1[j];

                ta_xxxx_xz_0[j] = pa_x[j] * ta_xxx_xz_0[j] - pc_x[j] * ta_xxx_xz_1[j] + 1.5 * fl1_fx * ta_xx_xz_0[j] - 1.5 * fl1_fx * ta_xx_xz_1[j] +
                                  0.5 * fl1_fx * ta_xxx_z_0[j] - 0.5 * fl1_fx * ta_xxx_z_1[j];

                ta_xxxx_yy_0[j] = pa_x[j] * ta_xxx_yy_0[j] - pc_x[j] * ta_xxx_yy_1[j] + 1.5 * fl1_fx * ta_xx_yy_0[j] - 1.5 * fl1_fx * ta_xx_yy_1[j];

                ta_xxxx_yz_0[j] = pa_x[j] * ta_xxx_yz_0[j] - pc_x[j] * ta_xxx_yz_1[j] + 1.5 * fl1_fx * ta_xx_yz_0[j] - 1.5 * fl1_fx * ta_xx_yz_1[j];

                ta_xxxx_zz_0[j] = pa_x[j] * ta_xxx_zz_0[j] - pc_x[j] * ta_xxx_zz_1[j] + 1.5 * fl1_fx * ta_xx_zz_0[j] - 1.5 * fl1_fx * ta_xx_zz_1[j];

                ta_xxxy_xx_0[j] = pa_x[j] * ta_xxy_xx_0[j] - pc_x[j] * ta_xxy_xx_1[j] + fl1_fx * ta_xy_xx_0[j] - fl1_fx * ta_xy_xx_1[j] +
                                  fl1_fx * ta_xxy_x_0[j] - fl1_fx * ta_xxy_x_1[j];

                ta_xxxy_xy_0[j] = pa_x[j] * ta_xxy_xy_0[j] - pc_x[j] * ta_xxy_xy_1[j] + fl1_fx * ta_xy_xy_0[j] - fl1_fx * ta_xy_xy_1[j] +
                                  0.5 * fl1_fx * ta_xxy_y_0[j] - 0.5 * fl1_fx * ta_xxy_y_1[j];

                ta_xxxy_xz_0[j] = pa_x[j] * ta_xxy_xz_0[j] - pc_x[j] * ta_xxy_xz_1[j] + fl1_fx * ta_xy_xz_0[j] - fl1_fx * ta_xy_xz_1[j] +
                                  0.5 * fl1_fx * ta_xxy_z_0[j] - 0.5 * fl1_fx * ta_xxy_z_1[j];

                ta_xxxy_yy_0[j] = pa_x[j] * ta_xxy_yy_0[j] - pc_x[j] * ta_xxy_yy_1[j] + fl1_fx * ta_xy_yy_0[j] - fl1_fx * ta_xy_yy_1[j];

                ta_xxxy_yz_0[j] = pa_x[j] * ta_xxy_yz_0[j] - pc_x[j] * ta_xxy_yz_1[j] + fl1_fx * ta_xy_yz_0[j] - fl1_fx * ta_xy_yz_1[j];

                ta_xxxy_zz_0[j] = pa_x[j] * ta_xxy_zz_0[j] - pc_x[j] * ta_xxy_zz_1[j] + fl1_fx * ta_xy_zz_0[j] - fl1_fx * ta_xy_zz_1[j];

                ta_xxxz_xx_0[j] = pa_x[j] * ta_xxz_xx_0[j] - pc_x[j] * ta_xxz_xx_1[j] + fl1_fx * ta_xz_xx_0[j] - fl1_fx * ta_xz_xx_1[j] +
                                  fl1_fx * ta_xxz_x_0[j] - fl1_fx * ta_xxz_x_1[j];

                ta_xxxz_xy_0[j] = pa_x[j] * ta_xxz_xy_0[j] - pc_x[j] * ta_xxz_xy_1[j] + fl1_fx * ta_xz_xy_0[j] - fl1_fx * ta_xz_xy_1[j] +
                                  0.5 * fl1_fx * ta_xxz_y_0[j] - 0.5 * fl1_fx * ta_xxz_y_1[j];

                ta_xxxz_xz_0[j] = pa_x[j] * ta_xxz_xz_0[j] - pc_x[j] * ta_xxz_xz_1[j] + fl1_fx * ta_xz_xz_0[j] - fl1_fx * ta_xz_xz_1[j] +
                                  0.5 * fl1_fx * ta_xxz_z_0[j] - 0.5 * fl1_fx * ta_xxz_z_1[j];

                ta_xxxz_yy_0[j] = pa_x[j] * ta_xxz_yy_0[j] - pc_x[j] * ta_xxz_yy_1[j] + fl1_fx * ta_xz_yy_0[j] - fl1_fx * ta_xz_yy_1[j];

                ta_xxxz_yz_0[j] = pa_x[j] * ta_xxz_yz_0[j] - pc_x[j] * ta_xxz_yz_1[j] + fl1_fx * ta_xz_yz_0[j] - fl1_fx * ta_xz_yz_1[j];

                ta_xxxz_zz_0[j] = pa_x[j] * ta_xxz_zz_0[j] - pc_x[j] * ta_xxz_zz_1[j] + fl1_fx * ta_xz_zz_0[j] - fl1_fx * ta_xz_zz_1[j];

                ta_xxyy_xx_0[j] = pa_x[j] * ta_xyy_xx_0[j] - pc_x[j] * ta_xyy_xx_1[j] + 0.5 * fl1_fx * ta_yy_xx_0[j] - 0.5 * fl1_fx * ta_yy_xx_1[j] +
                                  fl1_fx * ta_xyy_x_0[j] - fl1_fx * ta_xyy_x_1[j];

                ta_xxyy_xy_0[j] = pa_x[j] * ta_xyy_xy_0[j] - pc_x[j] * ta_xyy_xy_1[j] + 0.5 * fl1_fx * ta_yy_xy_0[j] - 0.5 * fl1_fx * ta_yy_xy_1[j] +
                                  0.5 * fl1_fx * ta_xyy_y_0[j] - 0.5 * fl1_fx * ta_xyy_y_1[j];

                ta_xxyy_xz_0[j] = pa_x[j] * ta_xyy_xz_0[j] - pc_x[j] * ta_xyy_xz_1[j] + 0.5 * fl1_fx * ta_yy_xz_0[j] - 0.5 * fl1_fx * ta_yy_xz_1[j] +
                                  0.5 * fl1_fx * ta_xyy_z_0[j] - 0.5 * fl1_fx * ta_xyy_z_1[j];

                ta_xxyy_yy_0[j] = pa_x[j] * ta_xyy_yy_0[j] - pc_x[j] * ta_xyy_yy_1[j] + 0.5 * fl1_fx * ta_yy_yy_0[j] - 0.5 * fl1_fx * ta_yy_yy_1[j];

                ta_xxyy_yz_0[j] = pa_x[j] * ta_xyy_yz_0[j] - pc_x[j] * ta_xyy_yz_1[j] + 0.5 * fl1_fx * ta_yy_yz_0[j] - 0.5 * fl1_fx * ta_yy_yz_1[j];

                ta_xxyy_zz_0[j] = pa_x[j] * ta_xyy_zz_0[j] - pc_x[j] * ta_xyy_zz_1[j] + 0.5 * fl1_fx * ta_yy_zz_0[j] - 0.5 * fl1_fx * ta_yy_zz_1[j];

                ta_xxyz_xx_0[j] = pa_x[j] * ta_xyz_xx_0[j] - pc_x[j] * ta_xyz_xx_1[j] + 0.5 * fl1_fx * ta_yz_xx_0[j] - 0.5 * fl1_fx * ta_yz_xx_1[j] +
                                  fl1_fx * ta_xyz_x_0[j] - fl1_fx * ta_xyz_x_1[j];

                ta_xxyz_xy_0[j] = pa_x[j] * ta_xyz_xy_0[j] - pc_x[j] * ta_xyz_xy_1[j] + 0.5 * fl1_fx * ta_yz_xy_0[j] - 0.5 * fl1_fx * ta_yz_xy_1[j] +
                                  0.5 * fl1_fx * ta_xyz_y_0[j] - 0.5 * fl1_fx * ta_xyz_y_1[j];

                ta_xxyz_xz_0[j] = pa_x[j] * ta_xyz_xz_0[j] - pc_x[j] * ta_xyz_xz_1[j] + 0.5 * fl1_fx * ta_yz_xz_0[j] - 0.5 * fl1_fx * ta_yz_xz_1[j] +
                                  0.5 * fl1_fx * ta_xyz_z_0[j] - 0.5 * fl1_fx * ta_xyz_z_1[j];

                ta_xxyz_yy_0[j] = pa_x[j] * ta_xyz_yy_0[j] - pc_x[j] * ta_xyz_yy_1[j] + 0.5 * fl1_fx * ta_yz_yy_0[j] - 0.5 * fl1_fx * ta_yz_yy_1[j];

                ta_xxyz_yz_0[j] = pa_x[j] * ta_xyz_yz_0[j] - pc_x[j] * ta_xyz_yz_1[j] + 0.5 * fl1_fx * ta_yz_yz_0[j] - 0.5 * fl1_fx * ta_yz_yz_1[j];

                ta_xxyz_zz_0[j] = pa_x[j] * ta_xyz_zz_0[j] - pc_x[j] * ta_xyz_zz_1[j] + 0.5 * fl1_fx * ta_yz_zz_0[j] - 0.5 * fl1_fx * ta_yz_zz_1[j];

                ta_xxzz_xx_0[j] = pa_x[j] * ta_xzz_xx_0[j] - pc_x[j] * ta_xzz_xx_1[j] + 0.5 * fl1_fx * ta_zz_xx_0[j] - 0.5 * fl1_fx * ta_zz_xx_1[j] +
                                  fl1_fx * ta_xzz_x_0[j] - fl1_fx * ta_xzz_x_1[j];

                ta_xxzz_xy_0[j] = pa_x[j] * ta_xzz_xy_0[j] - pc_x[j] * ta_xzz_xy_1[j] + 0.5 * fl1_fx * ta_zz_xy_0[j] - 0.5 * fl1_fx * ta_zz_xy_1[j] +
                                  0.5 * fl1_fx * ta_xzz_y_0[j] - 0.5 * fl1_fx * ta_xzz_y_1[j];

                ta_xxzz_xz_0[j] = pa_x[j] * ta_xzz_xz_0[j] - pc_x[j] * ta_xzz_xz_1[j] + 0.5 * fl1_fx * ta_zz_xz_0[j] - 0.5 * fl1_fx * ta_zz_xz_1[j] +
                                  0.5 * fl1_fx * ta_xzz_z_0[j] - 0.5 * fl1_fx * ta_xzz_z_1[j];

                ta_xxzz_yy_0[j] = pa_x[j] * ta_xzz_yy_0[j] - pc_x[j] * ta_xzz_yy_1[j] + 0.5 * fl1_fx * ta_zz_yy_0[j] - 0.5 * fl1_fx * ta_zz_yy_1[j];

                ta_xxzz_yz_0[j] = pa_x[j] * ta_xzz_yz_0[j] - pc_x[j] * ta_xzz_yz_1[j] + 0.5 * fl1_fx * ta_zz_yz_0[j] - 0.5 * fl1_fx * ta_zz_yz_1[j];

                ta_xxzz_zz_0[j] = pa_x[j] * ta_xzz_zz_0[j] - pc_x[j] * ta_xzz_zz_1[j] + 0.5 * fl1_fx * ta_zz_zz_0[j] - 0.5 * fl1_fx * ta_zz_zz_1[j];

                ta_xyyy_xx_0[j] = pa_x[j] * ta_yyy_xx_0[j] - pc_x[j] * ta_yyy_xx_1[j] + fl1_fx * ta_yyy_x_0[j] - fl1_fx * ta_yyy_x_1[j];

                ta_xyyy_xy_0[j] = pa_x[j] * ta_yyy_xy_0[j] - pc_x[j] * ta_yyy_xy_1[j] + 0.5 * fl1_fx * ta_yyy_y_0[j] - 0.5 * fl1_fx * ta_yyy_y_1[j];

                ta_xyyy_xz_0[j] = pa_x[j] * ta_yyy_xz_0[j] - pc_x[j] * ta_yyy_xz_1[j] + 0.5 * fl1_fx * ta_yyy_z_0[j] - 0.5 * fl1_fx * ta_yyy_z_1[j];

                ta_xyyy_yy_0[j] = pa_x[j] * ta_yyy_yy_0[j] - pc_x[j] * ta_yyy_yy_1[j];

                ta_xyyy_yz_0[j] = pa_x[j] * ta_yyy_yz_0[j] - pc_x[j] * ta_yyy_yz_1[j];

                ta_xyyy_zz_0[j] = pa_x[j] * ta_yyy_zz_0[j] - pc_x[j] * ta_yyy_zz_1[j];

                ta_xyyz_xx_0[j] = pa_x[j] * ta_yyz_xx_0[j] - pc_x[j] * ta_yyz_xx_1[j] + fl1_fx * ta_yyz_x_0[j] - fl1_fx * ta_yyz_x_1[j];

                ta_xyyz_xy_0[j] = pa_x[j] * ta_yyz_xy_0[j] - pc_x[j] * ta_yyz_xy_1[j] + 0.5 * fl1_fx * ta_yyz_y_0[j] - 0.5 * fl1_fx * ta_yyz_y_1[j];

                ta_xyyz_xz_0[j] = pa_x[j] * ta_yyz_xz_0[j] - pc_x[j] * ta_yyz_xz_1[j] + 0.5 * fl1_fx * ta_yyz_z_0[j] - 0.5 * fl1_fx * ta_yyz_z_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGD_45_90(CMemBlock2D<double>&       primBuffer,
                                const CRecursionMap&       recursionMap,
                                const CMemBlock2D<double>& osFactors,
                                const CMemBlock2D<double>& paDistances,
                                const CMemBlock2D<double>& pcDistances,
                                const CGtoBlock&           braGtoBlock,
                                const CGtoBlock&           ketGtoBlock,
                                const int32_t              iContrGto)
{
    // Batch of Integrals (45,90)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto ta_yyy_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 36);

            auto ta_yyy_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 37);

            auto ta_yyy_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 38);

            auto ta_yyy_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 39);

            auto ta_yyy_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 40);

            auto ta_yyy_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 41);

            auto ta_yyz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 42);

            auto ta_yyz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 43);

            auto ta_yyz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 44);

            auto ta_yyz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 45);

            auto ta_yyz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 46);

            auto ta_yyz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 47);

            auto ta_yzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 48);

            auto ta_yzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 49);

            auto ta_yzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 50);

            auto ta_yzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 51);

            auto ta_yzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 52);

            auto ta_yzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 53);

            auto ta_zzz_xx_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 54);

            auto ta_zzz_xy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 55);

            auto ta_zzz_xz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 56);

            auto ta_zzz_yy_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 57);

            auto ta_zzz_yz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 58);

            auto ta_zzz_zz_0 = primBuffer.data(pidx_a_3_2_m0 + 60 * idx + 59);

            auto ta_yyy_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 36);

            auto ta_yyy_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 37);

            auto ta_yyy_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 38);

            auto ta_yyy_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 39);

            auto ta_yyy_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 40);

            auto ta_yyy_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 41);

            auto ta_yyz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 42);

            auto ta_yyz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 43);

            auto ta_yyz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 44);

            auto ta_yyz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 45);

            auto ta_yyz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 46);

            auto ta_yyz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 47);

            auto ta_yzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 48);

            auto ta_yzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 49);

            auto ta_yzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 50);

            auto ta_yzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 51);

            auto ta_yzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 52);

            auto ta_yzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 53);

            auto ta_zzz_xx_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 54);

            auto ta_zzz_xy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 55);

            auto ta_zzz_xz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 56);

            auto ta_zzz_yy_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 57);

            auto ta_zzz_yz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 58);

            auto ta_zzz_zz_1 = primBuffer.data(pidx_a_3_2_m1 + 60 * idx + 59);

            auto ta_yy_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 18);

            auto ta_yy_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 19);

            auto ta_yy_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 20);

            auto ta_yy_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 21);

            auto ta_yy_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 22);

            auto ta_yy_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 23);

            auto ta_yz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 24);

            auto ta_yz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 25);

            auto ta_yz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 26);

            auto ta_yz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 27);

            auto ta_yz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 28);

            auto ta_yz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 29);

            auto ta_zz_xx_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 30);

            auto ta_zz_xy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 31);

            auto ta_zz_xz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 32);

            auto ta_zz_yy_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 33);

            auto ta_zz_yz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 34);

            auto ta_zz_zz_0 = primBuffer.data(pidx_a_2_2_m0 + 36 * idx + 35);

            auto ta_yy_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 18);

            auto ta_yy_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 19);

            auto ta_yy_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 20);

            auto ta_yy_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 21);

            auto ta_yy_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 22);

            auto ta_yy_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 23);

            auto ta_yz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 24);

            auto ta_yz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 25);

            auto ta_yz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 26);

            auto ta_yz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 27);

            auto ta_yz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 28);

            auto ta_yz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 29);

            auto ta_zz_xx_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 30);

            auto ta_zz_xy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 31);

            auto ta_zz_xz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 32);

            auto ta_zz_yy_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 33);

            auto ta_zz_yz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 34);

            auto ta_zz_zz_1 = primBuffer.data(pidx_a_2_2_m1 + 36 * idx + 35);

            auto ta_yyy_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 18);

            auto ta_yyy_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 19);

            auto ta_yyy_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 20);

            auto ta_yyz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 21);

            auto ta_yyz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 22);

            auto ta_yyz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 23);

            auto ta_yzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 24);

            auto ta_yzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 25);

            auto ta_yzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 26);

            auto ta_zzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 27);

            auto ta_zzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 28);

            auto ta_zzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 29);

            auto ta_yyy_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 18);

            auto ta_yyy_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 19);

            auto ta_yyy_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 20);

            auto ta_yyz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 21);

            auto ta_yyz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 22);

            auto ta_yyz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 23);

            auto ta_yzz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 24);

            auto ta_yzz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 25);

            auto ta_yzz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 26);

            auto ta_zzz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 27);

            auto ta_zzz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 28);

            auto ta_zzz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 29);

            // set up pointers to integrals

            auto ta_xyyz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 45);

            auto ta_xyyz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 46);

            auto ta_xyyz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 47);

            auto ta_xyzz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 48);

            auto ta_xyzz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 49);

            auto ta_xyzz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 50);

            auto ta_xyzz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 51);

            auto ta_xyzz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 52);

            auto ta_xyzz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 53);

            auto ta_xzzz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 54);

            auto ta_xzzz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 55);

            auto ta_xzzz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 56);

            auto ta_xzzz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 57);

            auto ta_xzzz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 58);

            auto ta_xzzz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 59);

            auto ta_yyyy_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 60);

            auto ta_yyyy_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 61);

            auto ta_yyyy_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 62);

            auto ta_yyyy_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 63);

            auto ta_yyyy_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 64);

            auto ta_yyyy_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 65);

            auto ta_yyyz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 66);

            auto ta_yyyz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 67);

            auto ta_yyyz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 68);

            auto ta_yyyz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 69);

            auto ta_yyyz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 70);

            auto ta_yyyz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 71);

            auto ta_yyzz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 72);

            auto ta_yyzz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 73);

            auto ta_yyzz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 74);

            auto ta_yyzz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 75);

            auto ta_yyzz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 76);

            auto ta_yyzz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 77);

            auto ta_yzzz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 78);

            auto ta_yzzz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 79);

            auto ta_yzzz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 80);

            auto ta_yzzz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 81);

            auto ta_yzzz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 82);

            auto ta_yzzz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 83);

            auto ta_zzzz_xx_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 84);

            auto ta_zzzz_xy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 85);

            auto ta_zzzz_xz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 86);

            auto ta_zzzz_yy_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 87);

            auto ta_zzzz_yz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 88);

            auto ta_zzzz_zz_0 = primBuffer.data(pidx_a_4_2_m0 + 90 * idx + 89);

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xyyz_yy_0, ta_xyyz_yz_0, \
                                         ta_xyyz_zz_0, ta_xyzz_xx_0, ta_xyzz_xy_0, ta_xyzz_xz_0, ta_xyzz_yy_0, ta_xyzz_yz_0, \
                                         ta_xyzz_zz_0, ta_xzzz_xx_0, ta_xzzz_xy_0, ta_xzzz_xz_0, ta_xzzz_yy_0, ta_xzzz_yz_0, \
                                         ta_xzzz_zz_0, ta_yy_xx_0, ta_yy_xx_1, ta_yy_xy_0, ta_yy_xy_1, ta_yy_xz_0, ta_yy_xz_1, \
                                         ta_yy_yy_0, ta_yy_yy_1, ta_yy_yz_0, ta_yy_yz_1, ta_yy_zz_0, ta_yy_zz_1, ta_yyy_x_0, \
                                         ta_yyy_x_1, ta_yyy_xx_0, ta_yyy_xx_1, ta_yyy_xy_0, ta_yyy_xy_1, ta_yyy_xz_0, \
                                         ta_yyy_xz_1, ta_yyy_y_0, ta_yyy_y_1, ta_yyy_yy_0, ta_yyy_yy_1, ta_yyy_yz_0, \
                                         ta_yyy_yz_1, ta_yyy_z_0, ta_yyy_z_1, ta_yyy_zz_0, ta_yyy_zz_1, ta_yyyy_xx_0, \
                                         ta_yyyy_xy_0, ta_yyyy_xz_0, ta_yyyy_yy_0, ta_yyyy_yz_0, ta_yyyy_zz_0, ta_yyyz_xx_0, \
                                         ta_yyyz_xy_0, ta_yyyz_xz_0, ta_yyyz_yy_0, ta_yyyz_yz_0, ta_yyyz_zz_0, ta_yyz_x_0, \
                                         ta_yyz_x_1, ta_yyz_xx_0, ta_yyz_xx_1, ta_yyz_xy_0, ta_yyz_xy_1, ta_yyz_xz_0, \
                                         ta_yyz_xz_1, ta_yyz_y_0, ta_yyz_y_1, ta_yyz_yy_0, ta_yyz_yy_1, ta_yyz_yz_0, \
                                         ta_yyz_yz_1, ta_yyz_z_0, ta_yyz_z_1, ta_yyz_zz_0, ta_yyz_zz_1, ta_yyzz_xx_0, \
                                         ta_yyzz_xy_0, ta_yyzz_xz_0, ta_yyzz_yy_0, ta_yyzz_yz_0, ta_yyzz_zz_0, ta_yz_xx_0, \
                                         ta_yz_xx_1, ta_yz_xy_0, ta_yz_xy_1, ta_yz_xz_0, ta_yz_xz_1, ta_yz_yy_0, ta_yz_yy_1, \
                                         ta_yz_yz_0, ta_yz_yz_1, ta_yz_zz_0, ta_yz_zz_1, ta_yzz_x_0, ta_yzz_x_1, \
                                         ta_yzz_xx_0, ta_yzz_xx_1, ta_yzz_xy_0, ta_yzz_xy_1, ta_yzz_xz_0, ta_yzz_xz_1, \
                                         ta_yzz_y_0, ta_yzz_y_1, ta_yzz_yy_0, ta_yzz_yy_1, ta_yzz_yz_0, ta_yzz_yz_1, \
                                         ta_yzz_z_0, ta_yzz_z_1, ta_yzz_zz_0, ta_yzz_zz_1, ta_yzzz_xx_0, ta_yzzz_xy_0, \
                                         ta_yzzz_xz_0, ta_yzzz_yy_0, ta_yzzz_yz_0, ta_yzzz_zz_0, ta_zz_xx_0, ta_zz_xx_1, \
                                         ta_zz_xy_0, ta_zz_xy_1, ta_zz_xz_0, ta_zz_xz_1, ta_zz_yy_0, ta_zz_yy_1, ta_zz_yz_0, \
                                         ta_zz_yz_1, ta_zz_zz_0, ta_zz_zz_1, ta_zzz_x_0, ta_zzz_x_1, ta_zzz_xx_0, \
                                         ta_zzz_xx_1, ta_zzz_xy_0, ta_zzz_xy_1, ta_zzz_xz_0, ta_zzz_xz_1, ta_zzz_y_0, \
                                         ta_zzz_y_1, ta_zzz_yy_0, ta_zzz_yy_1, ta_zzz_yz_0, ta_zzz_yz_1, ta_zzz_z_0, \
                                         ta_zzz_z_1, ta_zzz_zz_0, ta_zzz_zz_1, ta_zzzz_xx_0, ta_zzzz_xy_0, ta_zzzz_xz_0, \
                                         ta_zzzz_yy_0, ta_zzzz_yz_0, ta_zzzz_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xyyz_yy_0[j] = pa_x[j] * ta_yyz_yy_0[j] - pc_x[j] * ta_yyz_yy_1[j];

                ta_xyyz_yz_0[j] = pa_x[j] * ta_yyz_yz_0[j] - pc_x[j] * ta_yyz_yz_1[j];

                ta_xyyz_zz_0[j] = pa_x[j] * ta_yyz_zz_0[j] - pc_x[j] * ta_yyz_zz_1[j];

                ta_xyzz_xx_0[j] = pa_x[j] * ta_yzz_xx_0[j] - pc_x[j] * ta_yzz_xx_1[j] + fl1_fx * ta_yzz_x_0[j] - fl1_fx * ta_yzz_x_1[j];

                ta_xyzz_xy_0[j] = pa_x[j] * ta_yzz_xy_0[j] - pc_x[j] * ta_yzz_xy_1[j] + 0.5 * fl1_fx * ta_yzz_y_0[j] - 0.5 * fl1_fx * ta_yzz_y_1[j];

                ta_xyzz_xz_0[j] = pa_x[j] * ta_yzz_xz_0[j] - pc_x[j] * ta_yzz_xz_1[j] + 0.5 * fl1_fx * ta_yzz_z_0[j] - 0.5 * fl1_fx * ta_yzz_z_1[j];

                ta_xyzz_yy_0[j] = pa_x[j] * ta_yzz_yy_0[j] - pc_x[j] * ta_yzz_yy_1[j];

                ta_xyzz_yz_0[j] = pa_x[j] * ta_yzz_yz_0[j] - pc_x[j] * ta_yzz_yz_1[j];

                ta_xyzz_zz_0[j] = pa_x[j] * ta_yzz_zz_0[j] - pc_x[j] * ta_yzz_zz_1[j];

                ta_xzzz_xx_0[j] = pa_x[j] * ta_zzz_xx_0[j] - pc_x[j] * ta_zzz_xx_1[j] + fl1_fx * ta_zzz_x_0[j] - fl1_fx * ta_zzz_x_1[j];

                ta_xzzz_xy_0[j] = pa_x[j] * ta_zzz_xy_0[j] - pc_x[j] * ta_zzz_xy_1[j] + 0.5 * fl1_fx * ta_zzz_y_0[j] - 0.5 * fl1_fx * ta_zzz_y_1[j];

                ta_xzzz_xz_0[j] = pa_x[j] * ta_zzz_xz_0[j] - pc_x[j] * ta_zzz_xz_1[j] + 0.5 * fl1_fx * ta_zzz_z_0[j] - 0.5 * fl1_fx * ta_zzz_z_1[j];

                ta_xzzz_yy_0[j] = pa_x[j] * ta_zzz_yy_0[j] - pc_x[j] * ta_zzz_yy_1[j];

                ta_xzzz_yz_0[j] = pa_x[j] * ta_zzz_yz_0[j] - pc_x[j] * ta_zzz_yz_1[j];

                ta_xzzz_zz_0[j] = pa_x[j] * ta_zzz_zz_0[j] - pc_x[j] * ta_zzz_zz_1[j];

                ta_yyyy_xx_0[j] = pa_y[j] * ta_yyy_xx_0[j] - pc_y[j] * ta_yyy_xx_1[j] + 1.5 * fl1_fx * ta_yy_xx_0[j] - 1.5 * fl1_fx * ta_yy_xx_1[j];

                ta_yyyy_xy_0[j] = pa_y[j] * ta_yyy_xy_0[j] - pc_y[j] * ta_yyy_xy_1[j] + 1.5 * fl1_fx * ta_yy_xy_0[j] - 1.5 * fl1_fx * ta_yy_xy_1[j] +
                                  0.5 * fl1_fx * ta_yyy_x_0[j] - 0.5 * fl1_fx * ta_yyy_x_1[j];

                ta_yyyy_xz_0[j] = pa_y[j] * ta_yyy_xz_0[j] - pc_y[j] * ta_yyy_xz_1[j] + 1.5 * fl1_fx * ta_yy_xz_0[j] - 1.5 * fl1_fx * ta_yy_xz_1[j];

                ta_yyyy_yy_0[j] = pa_y[j] * ta_yyy_yy_0[j] - pc_y[j] * ta_yyy_yy_1[j] + 1.5 * fl1_fx * ta_yy_yy_0[j] - 1.5 * fl1_fx * ta_yy_yy_1[j] +
                                  fl1_fx * ta_yyy_y_0[j] - fl1_fx * ta_yyy_y_1[j];

                ta_yyyy_yz_0[j] = pa_y[j] * ta_yyy_yz_0[j] - pc_y[j] * ta_yyy_yz_1[j] + 1.5 * fl1_fx * ta_yy_yz_0[j] - 1.5 * fl1_fx * ta_yy_yz_1[j] +
                                  0.5 * fl1_fx * ta_yyy_z_0[j] - 0.5 * fl1_fx * ta_yyy_z_1[j];

                ta_yyyy_zz_0[j] = pa_y[j] * ta_yyy_zz_0[j] - pc_y[j] * ta_yyy_zz_1[j] + 1.5 * fl1_fx * ta_yy_zz_0[j] - 1.5 * fl1_fx * ta_yy_zz_1[j];

                ta_yyyz_xx_0[j] = pa_y[j] * ta_yyz_xx_0[j] - pc_y[j] * ta_yyz_xx_1[j] + fl1_fx * ta_yz_xx_0[j] - fl1_fx * ta_yz_xx_1[j];

                ta_yyyz_xy_0[j] = pa_y[j] * ta_yyz_xy_0[j] - pc_y[j] * ta_yyz_xy_1[j] + fl1_fx * ta_yz_xy_0[j] - fl1_fx * ta_yz_xy_1[j] +
                                  0.5 * fl1_fx * ta_yyz_x_0[j] - 0.5 * fl1_fx * ta_yyz_x_1[j];

                ta_yyyz_xz_0[j] = pa_y[j] * ta_yyz_xz_0[j] - pc_y[j] * ta_yyz_xz_1[j] + fl1_fx * ta_yz_xz_0[j] - fl1_fx * ta_yz_xz_1[j];

                ta_yyyz_yy_0[j] = pa_y[j] * ta_yyz_yy_0[j] - pc_y[j] * ta_yyz_yy_1[j] + fl1_fx * ta_yz_yy_0[j] - fl1_fx * ta_yz_yy_1[j] +
                                  fl1_fx * ta_yyz_y_0[j] - fl1_fx * ta_yyz_y_1[j];

                ta_yyyz_yz_0[j] = pa_y[j] * ta_yyz_yz_0[j] - pc_y[j] * ta_yyz_yz_1[j] + fl1_fx * ta_yz_yz_0[j] - fl1_fx * ta_yz_yz_1[j] +
                                  0.5 * fl1_fx * ta_yyz_z_0[j] - 0.5 * fl1_fx * ta_yyz_z_1[j];

                ta_yyyz_zz_0[j] = pa_y[j] * ta_yyz_zz_0[j] - pc_y[j] * ta_yyz_zz_1[j] + fl1_fx * ta_yz_zz_0[j] - fl1_fx * ta_yz_zz_1[j];

                ta_yyzz_xx_0[j] = pa_y[j] * ta_yzz_xx_0[j] - pc_y[j] * ta_yzz_xx_1[j] + 0.5 * fl1_fx * ta_zz_xx_0[j] - 0.5 * fl1_fx * ta_zz_xx_1[j];

                ta_yyzz_xy_0[j] = pa_y[j] * ta_yzz_xy_0[j] - pc_y[j] * ta_yzz_xy_1[j] + 0.5 * fl1_fx * ta_zz_xy_0[j] - 0.5 * fl1_fx * ta_zz_xy_1[j] +
                                  0.5 * fl1_fx * ta_yzz_x_0[j] - 0.5 * fl1_fx * ta_yzz_x_1[j];

                ta_yyzz_xz_0[j] = pa_y[j] * ta_yzz_xz_0[j] - pc_y[j] * ta_yzz_xz_1[j] + 0.5 * fl1_fx * ta_zz_xz_0[j] - 0.5 * fl1_fx * ta_zz_xz_1[j];

                ta_yyzz_yy_0[j] = pa_y[j] * ta_yzz_yy_0[j] - pc_y[j] * ta_yzz_yy_1[j] + 0.5 * fl1_fx * ta_zz_yy_0[j] - 0.5 * fl1_fx * ta_zz_yy_1[j] +
                                  fl1_fx * ta_yzz_y_0[j] - fl1_fx * ta_yzz_y_1[j];

                ta_yyzz_yz_0[j] = pa_y[j] * ta_yzz_yz_0[j] - pc_y[j] * ta_yzz_yz_1[j] + 0.5 * fl1_fx * ta_zz_yz_0[j] - 0.5 * fl1_fx * ta_zz_yz_1[j] +
                                  0.5 * fl1_fx * ta_yzz_z_0[j] - 0.5 * fl1_fx * ta_yzz_z_1[j];

                ta_yyzz_zz_0[j] = pa_y[j] * ta_yzz_zz_0[j] - pc_y[j] * ta_yzz_zz_1[j] + 0.5 * fl1_fx * ta_zz_zz_0[j] - 0.5 * fl1_fx * ta_zz_zz_1[j];

                ta_yzzz_xx_0[j] = pa_y[j] * ta_zzz_xx_0[j] - pc_y[j] * ta_zzz_xx_1[j];

                ta_yzzz_xy_0[j] = pa_y[j] * ta_zzz_xy_0[j] - pc_y[j] * ta_zzz_xy_1[j] + 0.5 * fl1_fx * ta_zzz_x_0[j] - 0.5 * fl1_fx * ta_zzz_x_1[j];

                ta_yzzz_xz_0[j] = pa_y[j] * ta_zzz_xz_0[j] - pc_y[j] * ta_zzz_xz_1[j];

                ta_yzzz_yy_0[j] = pa_y[j] * ta_zzz_yy_0[j] - pc_y[j] * ta_zzz_yy_1[j] + fl1_fx * ta_zzz_y_0[j] - fl1_fx * ta_zzz_y_1[j];

                ta_yzzz_yz_0[j] = pa_y[j] * ta_zzz_yz_0[j] - pc_y[j] * ta_zzz_yz_1[j] + 0.5 * fl1_fx * ta_zzz_z_0[j] - 0.5 * fl1_fx * ta_zzz_z_1[j];

                ta_yzzz_zz_0[j] = pa_y[j] * ta_zzz_zz_0[j] - pc_y[j] * ta_zzz_zz_1[j];

                ta_zzzz_xx_0[j] = pa_z[j] * ta_zzz_xx_0[j] - pc_z[j] * ta_zzz_xx_1[j] + 1.5 * fl1_fx * ta_zz_xx_0[j] - 1.5 * fl1_fx * ta_zz_xx_1[j];

                ta_zzzz_xy_0[j] = pa_z[j] * ta_zzz_xy_0[j] - pc_z[j] * ta_zzz_xy_1[j] + 1.5 * fl1_fx * ta_zz_xy_0[j] - 1.5 * fl1_fx * ta_zz_xy_1[j];

                ta_zzzz_xz_0[j] = pa_z[j] * ta_zzz_xz_0[j] - pc_z[j] * ta_zzz_xz_1[j] + 1.5 * fl1_fx * ta_zz_xz_0[j] - 1.5 * fl1_fx * ta_zz_xz_1[j] +
                                  0.5 * fl1_fx * ta_zzz_x_0[j] - 0.5 * fl1_fx * ta_zzz_x_1[j];

                ta_zzzz_yy_0[j] = pa_z[j] * ta_zzz_yy_0[j] - pc_z[j] * ta_zzz_yy_1[j] + 1.5 * fl1_fx * ta_zz_yy_0[j] - 1.5 * fl1_fx * ta_zz_yy_1[j];

                ta_zzzz_yz_0[j] = pa_z[j] * ta_zzz_yz_0[j] - pc_z[j] * ta_zzz_yz_1[j] + 1.5 * fl1_fx * ta_zz_yz_0[j] - 1.5 * fl1_fx * ta_zz_yz_1[j] +
                                  0.5 * fl1_fx * ta_zzz_y_0[j] - 0.5 * fl1_fx * ta_zzz_y_1[j];

                ta_zzzz_zz_0[j] = pa_z[j] * ta_zzz_zz_0[j] - pc_z[j] * ta_zzz_zz_1[j] + 1.5 * fl1_fx * ta_zz_zz_0[j] - 1.5 * fl1_fx * ta_zz_zz_1[j] +
                                  fl1_fx * ta_zzz_z_0[j] - fl1_fx * ta_zzz_z_1[j];
            }

            idx++;
        }
    }
}

}  // namespace npotrecfunc
