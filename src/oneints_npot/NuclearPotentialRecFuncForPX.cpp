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

#include "NuclearPotentialRecFuncForPX.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

void
compNuclearPotentialForPP(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_1_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_1_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_0_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_0_x_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx);

            auto ta_0_y_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 1);

            auto ta_0_z_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 2);

            auto ta_0_x_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx);

            auto ta_0_y_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 1);

            auto ta_0_z_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 2);

            auto ta_0_0_0 = primBuffer.data(pidx_a_0_0_m0 + idx);

            auto ta_0_0_1 = primBuffer.data(pidx_a_0_0_m1 + idx);

            // set up pointers to integrals

            auto ta_x_x_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx);

            auto ta_x_y_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 1);

            auto ta_x_z_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 2);

            auto ta_y_x_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 3);

            auto ta_y_y_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 4);

            auto ta_y_z_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 5);

            auto ta_z_x_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 6);

            auto ta_z_y_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 7);

            auto ta_z_z_0 = primBuffer.data(pidx_a_1_1_m0 + 9 * idx + 8);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_0_0, ta_0_0_1, ta_0_x_0, ta_0_x_1, \
                                         ta_0_y_0, ta_0_y_1, ta_0_z_0, ta_0_z_1, ta_x_x_0, ta_x_y_0, ta_x_z_0, ta_y_x_0, \
                                         ta_y_y_0, ta_y_z_0, ta_z_x_0, ta_z_y_0, ta_z_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_x_x_0[j] = pa_x[j] * ta_0_x_0[j] - pc_x[j] * ta_0_x_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];

                ta_x_y_0[j] = pa_x[j] * ta_0_y_0[j] - pc_x[j] * ta_0_y_1[j];

                ta_x_z_0[j] = pa_x[j] * ta_0_z_0[j] - pc_x[j] * ta_0_z_1[j];

                ta_y_x_0[j] = pa_y[j] * ta_0_x_0[j] - pc_y[j] * ta_0_x_1[j];

                ta_y_y_0[j] = pa_y[j] * ta_0_y_0[j] - pc_y[j] * ta_0_y_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];

                ta_y_z_0[j] = pa_y[j] * ta_0_z_0[j] - pc_y[j] * ta_0_z_1[j];

                ta_z_x_0[j] = pa_z[j] * ta_0_x_0[j] - pc_z[j] * ta_0_x_1[j];

                ta_z_y_0[j] = pa_z[j] * ta_0_y_0[j] - pc_z[j] * ta_0_y_1[j];

                ta_z_z_0[j] = pa_z[j] * ta_0_z_0[j] - pc_z[j] * ta_0_z_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForPD(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_1_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_1_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_0_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_0_x_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx);

            auto ta_0_y_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 1);

            auto ta_0_z_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 2);

            auto ta_0_x_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx);

            auto ta_0_y_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 1);

            auto ta_0_z_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 2);

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_x_0, ta_0_x_1, ta_0_xx_0, ta_0_xx_1, \
                                         ta_0_xy_0, ta_0_xy_1, ta_0_xz_0, ta_0_xz_1, ta_0_y_0, ta_0_y_1, ta_0_yy_0, \
                                         ta_0_yy_1, ta_0_yz_0, ta_0_yz_1, ta_0_z_0, ta_0_z_1, ta_0_zz_0, ta_0_zz_1, \
                                         ta_x_xx_0, ta_x_xy_0, ta_x_xz_0, ta_x_yy_0, ta_x_yz_0, ta_x_zz_0, ta_y_xx_0, \
                                         ta_y_xy_0, ta_y_xz_0, ta_y_yy_0, ta_y_yz_0, ta_y_zz_0, ta_z_xx_0, ta_z_xy_0, \
                                         ta_z_xz_0, ta_z_yy_0, ta_z_yz_0, ta_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_x_xx_0[j] = pa_x[j] * ta_0_xx_0[j] - pc_x[j] * ta_0_xx_1[j] + fl1_fx * ta_0_x_0[j] - fl1_fx * ta_0_x_1[j];

                ta_x_xy_0[j] = pa_x[j] * ta_0_xy_0[j] - pc_x[j] * ta_0_xy_1[j] + 0.5 * fl1_fx * ta_0_y_0[j] - 0.5 * fl1_fx * ta_0_y_1[j];

                ta_x_xz_0[j] = pa_x[j] * ta_0_xz_0[j] - pc_x[j] * ta_0_xz_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j];

                ta_x_yy_0[j] = pa_x[j] * ta_0_yy_0[j] - pc_x[j] * ta_0_yy_1[j];

                ta_x_yz_0[j] = pa_x[j] * ta_0_yz_0[j] - pc_x[j] * ta_0_yz_1[j];

                ta_x_zz_0[j] = pa_x[j] * ta_0_zz_0[j] - pc_x[j] * ta_0_zz_1[j];

                ta_y_xx_0[j] = pa_y[j] * ta_0_xx_0[j] - pc_y[j] * ta_0_xx_1[j];

                ta_y_xy_0[j] = pa_y[j] * ta_0_xy_0[j] - pc_y[j] * ta_0_xy_1[j] + 0.5 * fl1_fx * ta_0_x_0[j] - 0.5 * fl1_fx * ta_0_x_1[j];

                ta_y_xz_0[j] = pa_y[j] * ta_0_xz_0[j] - pc_y[j] * ta_0_xz_1[j];

                ta_y_yy_0[j] = pa_y[j] * ta_0_yy_0[j] - pc_y[j] * ta_0_yy_1[j] + fl1_fx * ta_0_y_0[j] - fl1_fx * ta_0_y_1[j];

                ta_y_yz_0[j] = pa_y[j] * ta_0_yz_0[j] - pc_y[j] * ta_0_yz_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j];

                ta_y_zz_0[j] = pa_y[j] * ta_0_zz_0[j] - pc_y[j] * ta_0_zz_1[j];

                ta_z_xx_0[j] = pa_z[j] * ta_0_xx_0[j] - pc_z[j] * ta_0_xx_1[j];

                ta_z_xy_0[j] = pa_z[j] * ta_0_xy_0[j] - pc_z[j] * ta_0_xy_1[j];

                ta_z_xz_0[j] = pa_z[j] * ta_0_xz_0[j] - pc_z[j] * ta_0_xz_1[j] + 0.5 * fl1_fx * ta_0_x_0[j] - 0.5 * fl1_fx * ta_0_x_1[j];

                ta_z_yy_0[j] = pa_z[j] * ta_0_yy_0[j] - pc_z[j] * ta_0_yy_1[j];

                ta_z_yz_0[j] = pa_z[j] * ta_0_yz_0[j] - pc_z[j] * ta_0_yz_1[j] + 0.5 * fl1_fx * ta_0_y_0[j] - 0.5 * fl1_fx * ta_0_y_1[j];

                ta_z_zz_0[j] = pa_z[j] * ta_0_zz_0[j] - pc_z[j] * ta_0_zz_1[j] + fl1_fx * ta_0_z_0[j] - fl1_fx * ta_0_z_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForDP(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_0_x_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx);

            auto ta_0_y_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 1);

            auto ta_0_z_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 2);

            auto ta_0_x_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx);

            auto ta_0_y_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 1);

            auto ta_0_z_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 2);

            auto ta_x_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx);

            auto ta_y_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 1);

            auto ta_z_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 2);

            auto ta_x_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx);

            auto ta_y_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 1);

            auto ta_z_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 2);

            // set up pointers to integrals

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

            auto ta_zz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 15);

            auto ta_zz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 16);

            auto ta_zz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 17);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_x_0, ta_0_x_1, ta_0_y_0, ta_0_y_1, \
                                         ta_0_z_0, ta_0_z_1, ta_x_0_0, ta_x_0_1, ta_x_x_0, ta_x_x_1, ta_x_y_0, ta_x_y_1, \
                                         ta_x_z_0, ta_x_z_1, ta_xx_x_0, ta_xx_y_0, ta_xx_z_0, ta_xy_x_0, ta_xy_y_0, \
                                         ta_xy_z_0, ta_xz_x_0, ta_xz_y_0, ta_xz_z_0, ta_y_0_0, ta_y_0_1, ta_y_x_0, ta_y_x_1, \
                                         ta_y_y_0, ta_y_y_1, ta_y_z_0, ta_y_z_1, ta_yy_x_0, ta_yy_y_0, ta_yy_z_0, ta_yz_x_0, \
                                         ta_yz_y_0, ta_yz_z_0, ta_z_0_0, ta_z_0_1, ta_z_x_0, ta_z_x_1, ta_z_y_0, ta_z_y_1, \
                                         ta_z_z_0, ta_z_z_1, ta_zz_x_0, ta_zz_y_0, ta_zz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xx_x_0[j] = pa_x[j] * ta_x_x_0[j] - pc_x[j] * ta_x_x_1[j] + 0.5 * fl1_fx * ta_0_x_0[j] - 0.5 * fl1_fx * ta_0_x_1[j] +
                               0.5 * fl1_fx * ta_x_0_0[j] - 0.5 * fl1_fx * ta_x_0_1[j];

                ta_xx_y_0[j] = pa_x[j] * ta_x_y_0[j] - pc_x[j] * ta_x_y_1[j] + 0.5 * fl1_fx * ta_0_y_0[j] - 0.5 * fl1_fx * ta_0_y_1[j];

                ta_xx_z_0[j] = pa_x[j] * ta_x_z_0[j] - pc_x[j] * ta_x_z_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j];

                ta_xy_x_0[j] = pa_x[j] * ta_y_x_0[j] - pc_x[j] * ta_y_x_1[j] + 0.5 * fl1_fx * ta_y_0_0[j] - 0.5 * fl1_fx * ta_y_0_1[j];

                ta_xy_y_0[j] = pa_x[j] * ta_y_y_0[j] - pc_x[j] * ta_y_y_1[j];

                ta_xy_z_0[j] = pa_x[j] * ta_y_z_0[j] - pc_x[j] * ta_y_z_1[j];

                ta_xz_x_0[j] = pa_x[j] * ta_z_x_0[j] - pc_x[j] * ta_z_x_1[j] + 0.5 * fl1_fx * ta_z_0_0[j] - 0.5 * fl1_fx * ta_z_0_1[j];

                ta_xz_y_0[j] = pa_x[j] * ta_z_y_0[j] - pc_x[j] * ta_z_y_1[j];

                ta_xz_z_0[j] = pa_x[j] * ta_z_z_0[j] - pc_x[j] * ta_z_z_1[j];

                ta_yy_x_0[j] = pa_y[j] * ta_y_x_0[j] - pc_y[j] * ta_y_x_1[j] + 0.5 * fl1_fx * ta_0_x_0[j] - 0.5 * fl1_fx * ta_0_x_1[j];

                ta_yy_y_0[j] = pa_y[j] * ta_y_y_0[j] - pc_y[j] * ta_y_y_1[j] + 0.5 * fl1_fx * ta_0_y_0[j] - 0.5 * fl1_fx * ta_0_y_1[j] +
                               0.5 * fl1_fx * ta_y_0_0[j] - 0.5 * fl1_fx * ta_y_0_1[j];

                ta_yy_z_0[j] = pa_y[j] * ta_y_z_0[j] - pc_y[j] * ta_y_z_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j];

                ta_yz_x_0[j] = pa_y[j] * ta_z_x_0[j] - pc_y[j] * ta_z_x_1[j];

                ta_yz_y_0[j] = pa_y[j] * ta_z_y_0[j] - pc_y[j] * ta_z_y_1[j] + 0.5 * fl1_fx * ta_z_0_0[j] - 0.5 * fl1_fx * ta_z_0_1[j];

                ta_yz_z_0[j] = pa_y[j] * ta_z_z_0[j] - pc_y[j] * ta_z_z_1[j];

                ta_zz_x_0[j] = pa_z[j] * ta_z_x_0[j] - pc_z[j] * ta_z_x_1[j] + 0.5 * fl1_fx * ta_0_x_0[j] - 0.5 * fl1_fx * ta_0_x_1[j];

                ta_zz_y_0[j] = pa_z[j] * ta_z_y_0[j] - pc_z[j] * ta_z_y_1[j] + 0.5 * fl1_fx * ta_0_y_0[j] - 0.5 * fl1_fx * ta_0_y_1[j];

                ta_zz_z_0[j] = pa_z[j] * ta_z_z_0[j] - pc_z[j] * ta_z_z_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j] +
                               0.5 * fl1_fx * ta_z_0_0[j] - 0.5 * fl1_fx * ta_z_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForPF(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_1_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_1_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_0_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_xx_0, ta_0_xx_1, ta_0_xxx_0, \
                                         ta_0_xxx_1, ta_0_xxy_0, ta_0_xxy_1, ta_0_xxz_0, ta_0_xxz_1, ta_0_xy_0, ta_0_xy_1, \
                                         ta_0_xyy_0, ta_0_xyy_1, ta_0_xyz_0, ta_0_xyz_1, ta_0_xz_0, ta_0_xz_1, ta_0_xzz_0, \
                                         ta_0_xzz_1, ta_0_yy_0, ta_0_yy_1, ta_0_yyy_0, ta_0_yyy_1, ta_0_yyz_0, ta_0_yyz_1, \
                                         ta_0_yz_0, ta_0_yz_1, ta_0_yzz_0, ta_0_yzz_1, ta_0_zz_0, ta_0_zz_1, ta_0_zzz_0, \
                                         ta_0_zzz_1, ta_x_xxx_0, ta_x_xxy_0, ta_x_xxz_0, ta_x_xyy_0, ta_x_xyz_0, ta_x_xzz_0, \
                                         ta_x_yyy_0, ta_x_yyz_0, ta_x_yzz_0, ta_x_zzz_0, ta_y_xxx_0, ta_y_xxy_0, ta_y_xxz_0, \
                                         ta_y_xyy_0, ta_y_xyz_0, ta_y_xzz_0, ta_y_yyy_0, ta_y_yyz_0, ta_y_yzz_0, ta_y_zzz_0, \
                                         ta_z_xxx_0, ta_z_xxy_0, ta_z_xxz_0, ta_z_xyy_0, ta_z_xyz_0, ta_z_xzz_0, ta_z_yyy_0, \
                                         ta_z_yyz_0, ta_z_yzz_0, ta_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_x_xxx_0[j] = pa_x[j] * ta_0_xxx_0[j] - pc_x[j] * ta_0_xxx_1[j] + 1.5 * fl1_fx * ta_0_xx_0[j] - 1.5 * fl1_fx * ta_0_xx_1[j];

                ta_x_xxy_0[j] = pa_x[j] * ta_0_xxy_0[j] - pc_x[j] * ta_0_xxy_1[j] + fl1_fx * ta_0_xy_0[j] - fl1_fx * ta_0_xy_1[j];

                ta_x_xxz_0[j] = pa_x[j] * ta_0_xxz_0[j] - pc_x[j] * ta_0_xxz_1[j] + fl1_fx * ta_0_xz_0[j] - fl1_fx * ta_0_xz_1[j];

                ta_x_xyy_0[j] = pa_x[j] * ta_0_xyy_0[j] - pc_x[j] * ta_0_xyy_1[j] + 0.5 * fl1_fx * ta_0_yy_0[j] - 0.5 * fl1_fx * ta_0_yy_1[j];

                ta_x_xyz_0[j] = pa_x[j] * ta_0_xyz_0[j] - pc_x[j] * ta_0_xyz_1[j] + 0.5 * fl1_fx * ta_0_yz_0[j] - 0.5 * fl1_fx * ta_0_yz_1[j];

                ta_x_xzz_0[j] = pa_x[j] * ta_0_xzz_0[j] - pc_x[j] * ta_0_xzz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j];

                ta_x_yyy_0[j] = pa_x[j] * ta_0_yyy_0[j] - pc_x[j] * ta_0_yyy_1[j];

                ta_x_yyz_0[j] = pa_x[j] * ta_0_yyz_0[j] - pc_x[j] * ta_0_yyz_1[j];

                ta_x_yzz_0[j] = pa_x[j] * ta_0_yzz_0[j] - pc_x[j] * ta_0_yzz_1[j];

                ta_x_zzz_0[j] = pa_x[j] * ta_0_zzz_0[j] - pc_x[j] * ta_0_zzz_1[j];

                ta_y_xxx_0[j] = pa_y[j] * ta_0_xxx_0[j] - pc_y[j] * ta_0_xxx_1[j];

                ta_y_xxy_0[j] = pa_y[j] * ta_0_xxy_0[j] - pc_y[j] * ta_0_xxy_1[j] + 0.5 * fl1_fx * ta_0_xx_0[j] - 0.5 * fl1_fx * ta_0_xx_1[j];

                ta_y_xxz_0[j] = pa_y[j] * ta_0_xxz_0[j] - pc_y[j] * ta_0_xxz_1[j];

                ta_y_xyy_0[j] = pa_y[j] * ta_0_xyy_0[j] - pc_y[j] * ta_0_xyy_1[j] + fl1_fx * ta_0_xy_0[j] - fl1_fx * ta_0_xy_1[j];

                ta_y_xyz_0[j] = pa_y[j] * ta_0_xyz_0[j] - pc_y[j] * ta_0_xyz_1[j] + 0.5 * fl1_fx * ta_0_xz_0[j] - 0.5 * fl1_fx * ta_0_xz_1[j];

                ta_y_xzz_0[j] = pa_y[j] * ta_0_xzz_0[j] - pc_y[j] * ta_0_xzz_1[j];

                ta_y_yyy_0[j] = pa_y[j] * ta_0_yyy_0[j] - pc_y[j] * ta_0_yyy_1[j] + 1.5 * fl1_fx * ta_0_yy_0[j] - 1.5 * fl1_fx * ta_0_yy_1[j];

                ta_y_yyz_0[j] = pa_y[j] * ta_0_yyz_0[j] - pc_y[j] * ta_0_yyz_1[j] + fl1_fx * ta_0_yz_0[j] - fl1_fx * ta_0_yz_1[j];

                ta_y_yzz_0[j] = pa_y[j] * ta_0_yzz_0[j] - pc_y[j] * ta_0_yzz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j];

                ta_y_zzz_0[j] = pa_y[j] * ta_0_zzz_0[j] - pc_y[j] * ta_0_zzz_1[j];

                ta_z_xxx_0[j] = pa_z[j] * ta_0_xxx_0[j] - pc_z[j] * ta_0_xxx_1[j];

                ta_z_xxy_0[j] = pa_z[j] * ta_0_xxy_0[j] - pc_z[j] * ta_0_xxy_1[j];

                ta_z_xxz_0[j] = pa_z[j] * ta_0_xxz_0[j] - pc_z[j] * ta_0_xxz_1[j] + 0.5 * fl1_fx * ta_0_xx_0[j] - 0.5 * fl1_fx * ta_0_xx_1[j];

                ta_z_xyy_0[j] = pa_z[j] * ta_0_xyy_0[j] - pc_z[j] * ta_0_xyy_1[j];

                ta_z_xyz_0[j] = pa_z[j] * ta_0_xyz_0[j] - pc_z[j] * ta_0_xyz_1[j] + 0.5 * fl1_fx * ta_0_xy_0[j] - 0.5 * fl1_fx * ta_0_xy_1[j];

                ta_z_xzz_0[j] = pa_z[j] * ta_0_xzz_0[j] - pc_z[j] * ta_0_xzz_1[j] + fl1_fx * ta_0_xz_0[j] - fl1_fx * ta_0_xz_1[j];

                ta_z_yyy_0[j] = pa_z[j] * ta_0_yyy_0[j] - pc_z[j] * ta_0_yyy_1[j];

                ta_z_yyz_0[j] = pa_z[j] * ta_0_yyz_0[j] - pc_z[j] * ta_0_yyz_1[j] + 0.5 * fl1_fx * ta_0_yy_0[j] - 0.5 * fl1_fx * ta_0_yy_1[j];

                ta_z_yzz_0[j] = pa_z[j] * ta_0_yzz_0[j] - pc_z[j] * ta_0_yzz_1[j] + fl1_fx * ta_0_yz_0[j] - fl1_fx * ta_0_yz_1[j];

                ta_z_zzz_0[j] = pa_z[j] * ta_0_zzz_0[j] - pc_z[j] * ta_0_zzz_1[j] + 1.5 * fl1_fx * ta_0_zz_0[j] - 1.5 * fl1_fx * ta_0_zz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFP(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_zz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 15);

            auto ta_zz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 16);

            auto ta_zz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 17);

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

            auto ta_zz_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 15);

            auto ta_zz_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 16);

            auto ta_zz_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 17);

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

            auto ta_xx_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx);

            auto ta_xy_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 1);

            auto ta_xz_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 2);

            auto ta_yy_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 3);

            auto ta_yz_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 4);

            auto ta_zz_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 5);

            auto ta_xx_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx);

            auto ta_xy_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 1);

            auto ta_xz_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 2);

            auto ta_yy_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 3);

            auto ta_yz_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 4);

            auto ta_zz_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 5);

            // set up pointers to integrals

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

            auto ta_yzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 24);

            auto ta_yzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 25);

            auto ta_yzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 26);

            auto ta_zzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 27);

            auto ta_zzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 28);

            auto ta_zzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 29);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_x_x_0, ta_x_x_1, ta_x_y_0, ta_x_y_1, \
                                         ta_x_z_0, ta_x_z_1, ta_xx_0_0, ta_xx_0_1, ta_xx_x_0, ta_xx_x_1, ta_xx_y_0, \
                                         ta_xx_y_1, ta_xx_z_0, ta_xx_z_1, ta_xxx_x_0, ta_xxx_y_0, ta_xxx_z_0, ta_xxy_x_0, \
                                         ta_xxy_y_0, ta_xxy_z_0, ta_xxz_x_0, ta_xxz_y_0, ta_xxz_z_0, ta_xy_0_0, ta_xy_0_1, \
                                         ta_xy_x_0, ta_xy_x_1, ta_xy_y_0, ta_xy_y_1, ta_xy_z_0, ta_xy_z_1, ta_xyy_x_0, \
                                         ta_xyy_y_0, ta_xyy_z_0, ta_xyz_x_0, ta_xyz_y_0, ta_xyz_z_0, ta_xz_0_0, ta_xz_0_1, \
                                         ta_xz_x_0, ta_xz_x_1, ta_xz_y_0, ta_xz_y_1, ta_xz_z_0, ta_xz_z_1, ta_xzz_x_0, \
                                         ta_xzz_y_0, ta_xzz_z_0, ta_y_x_0, ta_y_x_1, ta_y_y_0, ta_y_y_1, ta_y_z_0, ta_y_z_1, \
                                         ta_yy_0_0, ta_yy_0_1, ta_yy_x_0, ta_yy_x_1, ta_yy_y_0, ta_yy_y_1, ta_yy_z_0, \
                                         ta_yy_z_1, ta_yyy_x_0, ta_yyy_y_0, ta_yyy_z_0, ta_yyz_x_0, ta_yyz_y_0, ta_yyz_z_0, \
                                         ta_yz_0_0, ta_yz_0_1, ta_yz_x_0, ta_yz_x_1, ta_yz_y_0, ta_yz_y_1, ta_yz_z_0, \
                                         ta_yz_z_1, ta_yzz_x_0, ta_yzz_y_0, ta_yzz_z_0, ta_z_x_0, ta_z_x_1, ta_z_y_0, \
                                         ta_z_y_1, ta_z_z_0, ta_z_z_1, ta_zz_0_0, ta_zz_0_1, ta_zz_x_0, ta_zz_x_1, \
                                         ta_zz_y_0, ta_zz_y_1, ta_zz_z_0, ta_zz_z_1, ta_zzz_x_0, ta_zzz_y_0, ta_zzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxx_x_0[j] = pa_x[j] * ta_xx_x_0[j] - pc_x[j] * ta_xx_x_1[j] + fl1_fx * ta_x_x_0[j] - fl1_fx * ta_x_x_1[j] +
                                0.5 * fl1_fx * ta_xx_0_0[j] - 0.5 * fl1_fx * ta_xx_0_1[j];

                ta_xxx_y_0[j] = pa_x[j] * ta_xx_y_0[j] - pc_x[j] * ta_xx_y_1[j] + fl1_fx * ta_x_y_0[j] - fl1_fx * ta_x_y_1[j];

                ta_xxx_z_0[j] = pa_x[j] * ta_xx_z_0[j] - pc_x[j] * ta_xx_z_1[j] + fl1_fx * ta_x_z_0[j] - fl1_fx * ta_x_z_1[j];

                ta_xxy_x_0[j] = pa_x[j] * ta_xy_x_0[j] - pc_x[j] * ta_xy_x_1[j] + 0.5 * fl1_fx * ta_y_x_0[j] - 0.5 * fl1_fx * ta_y_x_1[j] +
                                0.5 * fl1_fx * ta_xy_0_0[j] - 0.5 * fl1_fx * ta_xy_0_1[j];

                ta_xxy_y_0[j] = pa_x[j] * ta_xy_y_0[j] - pc_x[j] * ta_xy_y_1[j] + 0.5 * fl1_fx * ta_y_y_0[j] - 0.5 * fl1_fx * ta_y_y_1[j];

                ta_xxy_z_0[j] = pa_x[j] * ta_xy_z_0[j] - pc_x[j] * ta_xy_z_1[j] + 0.5 * fl1_fx * ta_y_z_0[j] - 0.5 * fl1_fx * ta_y_z_1[j];

                ta_xxz_x_0[j] = pa_x[j] * ta_xz_x_0[j] - pc_x[j] * ta_xz_x_1[j] + 0.5 * fl1_fx * ta_z_x_0[j] - 0.5 * fl1_fx * ta_z_x_1[j] +
                                0.5 * fl1_fx * ta_xz_0_0[j] - 0.5 * fl1_fx * ta_xz_0_1[j];

                ta_xxz_y_0[j] = pa_x[j] * ta_xz_y_0[j] - pc_x[j] * ta_xz_y_1[j] + 0.5 * fl1_fx * ta_z_y_0[j] - 0.5 * fl1_fx * ta_z_y_1[j];

                ta_xxz_z_0[j] = pa_x[j] * ta_xz_z_0[j] - pc_x[j] * ta_xz_z_1[j] + 0.5 * fl1_fx * ta_z_z_0[j] - 0.5 * fl1_fx * ta_z_z_1[j];

                ta_xyy_x_0[j] = pa_x[j] * ta_yy_x_0[j] - pc_x[j] * ta_yy_x_1[j] + 0.5 * fl1_fx * ta_yy_0_0[j] - 0.5 * fl1_fx * ta_yy_0_1[j];

                ta_xyy_y_0[j] = pa_x[j] * ta_yy_y_0[j] - pc_x[j] * ta_yy_y_1[j];

                ta_xyy_z_0[j] = pa_x[j] * ta_yy_z_0[j] - pc_x[j] * ta_yy_z_1[j];

                ta_xyz_x_0[j] = pa_x[j] * ta_yz_x_0[j] - pc_x[j] * ta_yz_x_1[j] + 0.5 * fl1_fx * ta_yz_0_0[j] - 0.5 * fl1_fx * ta_yz_0_1[j];

                ta_xyz_y_0[j] = pa_x[j] * ta_yz_y_0[j] - pc_x[j] * ta_yz_y_1[j];

                ta_xyz_z_0[j] = pa_x[j] * ta_yz_z_0[j] - pc_x[j] * ta_yz_z_1[j];

                ta_xzz_x_0[j] = pa_x[j] * ta_zz_x_0[j] - pc_x[j] * ta_zz_x_1[j] + 0.5 * fl1_fx * ta_zz_0_0[j] - 0.5 * fl1_fx * ta_zz_0_1[j];

                ta_xzz_y_0[j] = pa_x[j] * ta_zz_y_0[j] - pc_x[j] * ta_zz_y_1[j];

                ta_xzz_z_0[j] = pa_x[j] * ta_zz_z_0[j] - pc_x[j] * ta_zz_z_1[j];

                ta_yyy_x_0[j] = pa_y[j] * ta_yy_x_0[j] - pc_y[j] * ta_yy_x_1[j] + fl1_fx * ta_y_x_0[j] - fl1_fx * ta_y_x_1[j];

                ta_yyy_y_0[j] = pa_y[j] * ta_yy_y_0[j] - pc_y[j] * ta_yy_y_1[j] + fl1_fx * ta_y_y_0[j] - fl1_fx * ta_y_y_1[j] +
                                0.5 * fl1_fx * ta_yy_0_0[j] - 0.5 * fl1_fx * ta_yy_0_1[j];

                ta_yyy_z_0[j] = pa_y[j] * ta_yy_z_0[j] - pc_y[j] * ta_yy_z_1[j] + fl1_fx * ta_y_z_0[j] - fl1_fx * ta_y_z_1[j];

                ta_yyz_x_0[j] = pa_y[j] * ta_yz_x_0[j] - pc_y[j] * ta_yz_x_1[j] + 0.5 * fl1_fx * ta_z_x_0[j] - 0.5 * fl1_fx * ta_z_x_1[j];

                ta_yyz_y_0[j] = pa_y[j] * ta_yz_y_0[j] - pc_y[j] * ta_yz_y_1[j] + 0.5 * fl1_fx * ta_z_y_0[j] - 0.5 * fl1_fx * ta_z_y_1[j] +
                                0.5 * fl1_fx * ta_yz_0_0[j] - 0.5 * fl1_fx * ta_yz_0_1[j];

                ta_yyz_z_0[j] = pa_y[j] * ta_yz_z_0[j] - pc_y[j] * ta_yz_z_1[j] + 0.5 * fl1_fx * ta_z_z_0[j] - 0.5 * fl1_fx * ta_z_z_1[j];

                ta_yzz_x_0[j] = pa_y[j] * ta_zz_x_0[j] - pc_y[j] * ta_zz_x_1[j];

                ta_yzz_y_0[j] = pa_y[j] * ta_zz_y_0[j] - pc_y[j] * ta_zz_y_1[j] + 0.5 * fl1_fx * ta_zz_0_0[j] - 0.5 * fl1_fx * ta_zz_0_1[j];

                ta_yzz_z_0[j] = pa_y[j] * ta_zz_z_0[j] - pc_y[j] * ta_zz_z_1[j];

                ta_zzz_x_0[j] = pa_z[j] * ta_zz_x_0[j] - pc_z[j] * ta_zz_x_1[j] + fl1_fx * ta_z_x_0[j] - fl1_fx * ta_z_x_1[j];

                ta_zzz_y_0[j] = pa_z[j] * ta_zz_y_0[j] - pc_z[j] * ta_zz_y_1[j] + fl1_fx * ta_z_y_0[j] - fl1_fx * ta_z_y_1[j];

                ta_zzz_z_0[j] = pa_z[j] * ta_zz_z_0[j] - pc_z[j] * ta_zz_z_1[j] + fl1_fx * ta_z_z_0[j] - fl1_fx * ta_z_z_1[j] +
                                0.5 * fl1_fx * ta_zz_0_0[j] - 0.5 * fl1_fx * ta_zz_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForPG(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_1_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_1_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_0_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

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

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_xxx_0, ta_0_xxx_1, ta_0_xxxx_0, \
                                         ta_0_xxxx_1, ta_0_xxxy_0, ta_0_xxxy_1, ta_0_xxxz_0, ta_0_xxxz_1, ta_0_xxy_0, \
                                         ta_0_xxy_1, ta_0_xxyy_0, ta_0_xxyy_1, ta_0_xxyz_0, ta_0_xxyz_1, ta_0_xxz_0, \
                                         ta_0_xxz_1, ta_0_xxzz_0, ta_0_xxzz_1, ta_0_xyy_0, ta_0_xyy_1, ta_0_xyyy_0, \
                                         ta_0_xyyy_1, ta_0_xyyz_0, ta_0_xyyz_1, ta_0_xyz_0, ta_0_xyz_1, ta_0_xyzz_0, \
                                         ta_0_xyzz_1, ta_0_xzz_0, ta_0_xzz_1, ta_0_xzzz_0, ta_0_xzzz_1, ta_0_yyy_0, \
                                         ta_0_yyy_1, ta_0_yyyy_0, ta_0_yyyy_1, ta_0_yyyz_0, ta_0_yyyz_1, ta_0_yyz_0, \
                                         ta_0_yyz_1, ta_0_yyzz_0, ta_0_yyzz_1, ta_0_yzz_0, ta_0_yzz_1, ta_0_yzzz_0, \
                                         ta_0_yzzz_1, ta_0_zzz_0, ta_0_zzz_1, ta_0_zzzz_0, ta_0_zzzz_1, ta_x_xxxx_0, \
                                         ta_x_xxxy_0, ta_x_xxxz_0, ta_x_xxyy_0, ta_x_xxyz_0, ta_x_xxzz_0, ta_x_xyyy_0, \
                                         ta_x_xyyz_0, ta_x_xyzz_0, ta_x_xzzz_0, ta_x_yyyy_0, ta_x_yyyz_0, ta_x_yyzz_0, \
                                         ta_x_yzzz_0, ta_x_zzzz_0, ta_y_xxxx_0, ta_y_xxxy_0, ta_y_xxxz_0, ta_y_xxyy_0, \
                                         ta_y_xxyz_0, ta_y_xxzz_0, ta_y_xyyy_0, ta_y_xyyz_0, ta_y_xyzz_0, ta_y_xzzz_0, \
                                         ta_y_yyyy_0, ta_y_yyyz_0, ta_y_yyzz_0, ta_y_yzzz_0, ta_y_zzzz_0, ta_z_xxxx_0, \
                                         ta_z_xxxy_0, ta_z_xxxz_0, ta_z_xxyy_0, ta_z_xxyz_0, ta_z_xxzz_0, ta_z_xyyy_0, \
                                         ta_z_xyyz_0, ta_z_xyzz_0, ta_z_xzzz_0, ta_z_yyyy_0, ta_z_yyyz_0, ta_z_yyzz_0, \
                                         ta_z_yzzz_0, ta_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_x_xxxx_0[j] = pa_x[j] * ta_0_xxxx_0[j] - pc_x[j] * ta_0_xxxx_1[j] + 2.0 * fl1_fx * ta_0_xxx_0[j] - 2.0 * fl1_fx * ta_0_xxx_1[j];

                ta_x_xxxy_0[j] = pa_x[j] * ta_0_xxxy_0[j] - pc_x[j] * ta_0_xxxy_1[j] + 1.5 * fl1_fx * ta_0_xxy_0[j] - 1.5 * fl1_fx * ta_0_xxy_1[j];

                ta_x_xxxz_0[j] = pa_x[j] * ta_0_xxxz_0[j] - pc_x[j] * ta_0_xxxz_1[j] + 1.5 * fl1_fx * ta_0_xxz_0[j] - 1.5 * fl1_fx * ta_0_xxz_1[j];

                ta_x_xxyy_0[j] = pa_x[j] * ta_0_xxyy_0[j] - pc_x[j] * ta_0_xxyy_1[j] + fl1_fx * ta_0_xyy_0[j] - fl1_fx * ta_0_xyy_1[j];

                ta_x_xxyz_0[j] = pa_x[j] * ta_0_xxyz_0[j] - pc_x[j] * ta_0_xxyz_1[j] + fl1_fx * ta_0_xyz_0[j] - fl1_fx * ta_0_xyz_1[j];

                ta_x_xxzz_0[j] = pa_x[j] * ta_0_xxzz_0[j] - pc_x[j] * ta_0_xxzz_1[j] + fl1_fx * ta_0_xzz_0[j] - fl1_fx * ta_0_xzz_1[j];

                ta_x_xyyy_0[j] = pa_x[j] * ta_0_xyyy_0[j] - pc_x[j] * ta_0_xyyy_1[j] + 0.5 * fl1_fx * ta_0_yyy_0[j] - 0.5 * fl1_fx * ta_0_yyy_1[j];

                ta_x_xyyz_0[j] = pa_x[j] * ta_0_xyyz_0[j] - pc_x[j] * ta_0_xyyz_1[j] + 0.5 * fl1_fx * ta_0_yyz_0[j] - 0.5 * fl1_fx * ta_0_yyz_1[j];

                ta_x_xyzz_0[j] = pa_x[j] * ta_0_xyzz_0[j] - pc_x[j] * ta_0_xyzz_1[j] + 0.5 * fl1_fx * ta_0_yzz_0[j] - 0.5 * fl1_fx * ta_0_yzz_1[j];

                ta_x_xzzz_0[j] = pa_x[j] * ta_0_xzzz_0[j] - pc_x[j] * ta_0_xzzz_1[j] + 0.5 * fl1_fx * ta_0_zzz_0[j] - 0.5 * fl1_fx * ta_0_zzz_1[j];

                ta_x_yyyy_0[j] = pa_x[j] * ta_0_yyyy_0[j] - pc_x[j] * ta_0_yyyy_1[j];

                ta_x_yyyz_0[j] = pa_x[j] * ta_0_yyyz_0[j] - pc_x[j] * ta_0_yyyz_1[j];

                ta_x_yyzz_0[j] = pa_x[j] * ta_0_yyzz_0[j] - pc_x[j] * ta_0_yyzz_1[j];

                ta_x_yzzz_0[j] = pa_x[j] * ta_0_yzzz_0[j] - pc_x[j] * ta_0_yzzz_1[j];

                ta_x_zzzz_0[j] = pa_x[j] * ta_0_zzzz_0[j] - pc_x[j] * ta_0_zzzz_1[j];

                ta_y_xxxx_0[j] = pa_y[j] * ta_0_xxxx_0[j] - pc_y[j] * ta_0_xxxx_1[j];

                ta_y_xxxy_0[j] = pa_y[j] * ta_0_xxxy_0[j] - pc_y[j] * ta_0_xxxy_1[j] + 0.5 * fl1_fx * ta_0_xxx_0[j] - 0.5 * fl1_fx * ta_0_xxx_1[j];

                ta_y_xxxz_0[j] = pa_y[j] * ta_0_xxxz_0[j] - pc_y[j] * ta_0_xxxz_1[j];

                ta_y_xxyy_0[j] = pa_y[j] * ta_0_xxyy_0[j] - pc_y[j] * ta_0_xxyy_1[j] + fl1_fx * ta_0_xxy_0[j] - fl1_fx * ta_0_xxy_1[j];

                ta_y_xxyz_0[j] = pa_y[j] * ta_0_xxyz_0[j] - pc_y[j] * ta_0_xxyz_1[j] + 0.5 * fl1_fx * ta_0_xxz_0[j] - 0.5 * fl1_fx * ta_0_xxz_1[j];

                ta_y_xxzz_0[j] = pa_y[j] * ta_0_xxzz_0[j] - pc_y[j] * ta_0_xxzz_1[j];

                ta_y_xyyy_0[j] = pa_y[j] * ta_0_xyyy_0[j] - pc_y[j] * ta_0_xyyy_1[j] + 1.5 * fl1_fx * ta_0_xyy_0[j] - 1.5 * fl1_fx * ta_0_xyy_1[j];

                ta_y_xyyz_0[j] = pa_y[j] * ta_0_xyyz_0[j] - pc_y[j] * ta_0_xyyz_1[j] + fl1_fx * ta_0_xyz_0[j] - fl1_fx * ta_0_xyz_1[j];

                ta_y_xyzz_0[j] = pa_y[j] * ta_0_xyzz_0[j] - pc_y[j] * ta_0_xyzz_1[j] + 0.5 * fl1_fx * ta_0_xzz_0[j] - 0.5 * fl1_fx * ta_0_xzz_1[j];

                ta_y_xzzz_0[j] = pa_y[j] * ta_0_xzzz_0[j] - pc_y[j] * ta_0_xzzz_1[j];

                ta_y_yyyy_0[j] = pa_y[j] * ta_0_yyyy_0[j] - pc_y[j] * ta_0_yyyy_1[j] + 2.0 * fl1_fx * ta_0_yyy_0[j] - 2.0 * fl1_fx * ta_0_yyy_1[j];

                ta_y_yyyz_0[j] = pa_y[j] * ta_0_yyyz_0[j] - pc_y[j] * ta_0_yyyz_1[j] + 1.5 * fl1_fx * ta_0_yyz_0[j] - 1.5 * fl1_fx * ta_0_yyz_1[j];

                ta_y_yyzz_0[j] = pa_y[j] * ta_0_yyzz_0[j] - pc_y[j] * ta_0_yyzz_1[j] + fl1_fx * ta_0_yzz_0[j] - fl1_fx * ta_0_yzz_1[j];

                ta_y_yzzz_0[j] = pa_y[j] * ta_0_yzzz_0[j] - pc_y[j] * ta_0_yzzz_1[j] + 0.5 * fl1_fx * ta_0_zzz_0[j] - 0.5 * fl1_fx * ta_0_zzz_1[j];

                ta_y_zzzz_0[j] = pa_y[j] * ta_0_zzzz_0[j] - pc_y[j] * ta_0_zzzz_1[j];

                ta_z_xxxx_0[j] = pa_z[j] * ta_0_xxxx_0[j] - pc_z[j] * ta_0_xxxx_1[j];

                ta_z_xxxy_0[j] = pa_z[j] * ta_0_xxxy_0[j] - pc_z[j] * ta_0_xxxy_1[j];

                ta_z_xxxz_0[j] = pa_z[j] * ta_0_xxxz_0[j] - pc_z[j] * ta_0_xxxz_1[j] + 0.5 * fl1_fx * ta_0_xxx_0[j] - 0.5 * fl1_fx * ta_0_xxx_1[j];

                ta_z_xxyy_0[j] = pa_z[j] * ta_0_xxyy_0[j] - pc_z[j] * ta_0_xxyy_1[j];

                ta_z_xxyz_0[j] = pa_z[j] * ta_0_xxyz_0[j] - pc_z[j] * ta_0_xxyz_1[j] + 0.5 * fl1_fx * ta_0_xxy_0[j] - 0.5 * fl1_fx * ta_0_xxy_1[j];

                ta_z_xxzz_0[j] = pa_z[j] * ta_0_xxzz_0[j] - pc_z[j] * ta_0_xxzz_1[j] + fl1_fx * ta_0_xxz_0[j] - fl1_fx * ta_0_xxz_1[j];

                ta_z_xyyy_0[j] = pa_z[j] * ta_0_xyyy_0[j] - pc_z[j] * ta_0_xyyy_1[j];

                ta_z_xyyz_0[j] = pa_z[j] * ta_0_xyyz_0[j] - pc_z[j] * ta_0_xyyz_1[j] + 0.5 * fl1_fx * ta_0_xyy_0[j] - 0.5 * fl1_fx * ta_0_xyy_1[j];

                ta_z_xyzz_0[j] = pa_z[j] * ta_0_xyzz_0[j] - pc_z[j] * ta_0_xyzz_1[j] + fl1_fx * ta_0_xyz_0[j] - fl1_fx * ta_0_xyz_1[j];

                ta_z_xzzz_0[j] = pa_z[j] * ta_0_xzzz_0[j] - pc_z[j] * ta_0_xzzz_1[j] + 1.5 * fl1_fx * ta_0_xzz_0[j] - 1.5 * fl1_fx * ta_0_xzz_1[j];

                ta_z_yyyy_0[j] = pa_z[j] * ta_0_yyyy_0[j] - pc_z[j] * ta_0_yyyy_1[j];

                ta_z_yyyz_0[j] = pa_z[j] * ta_0_yyyz_0[j] - pc_z[j] * ta_0_yyyz_1[j] + 0.5 * fl1_fx * ta_0_yyy_0[j] - 0.5 * fl1_fx * ta_0_yyy_1[j];

                ta_z_yyzz_0[j] = pa_z[j] * ta_0_yyzz_0[j] - pc_z[j] * ta_0_yyzz_1[j] + fl1_fx * ta_0_yyz_0[j] - fl1_fx * ta_0_yyz_1[j];

                ta_z_yzzz_0[j] = pa_z[j] * ta_0_yzzz_0[j] - pc_z[j] * ta_0_yzzz_1[j] + 1.5 * fl1_fx * ta_0_yzz_0[j] - 1.5 * fl1_fx * ta_0_yzz_1[j];

                ta_z_zzzz_0[j] = pa_z[j] * ta_0_zzzz_0[j] - pc_z[j] * ta_0_zzzz_1[j] + 2.0 * fl1_fx * ta_0_zzz_0[j] - 2.0 * fl1_fx * ta_0_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGP(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_yzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 24);

            auto ta_yzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 25);

            auto ta_yzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 26);

            auto ta_zzz_x_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 27);

            auto ta_zzz_y_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 28);

            auto ta_zzz_z_0 = primBuffer.data(pidx_a_3_1_m0 + 30 * idx + 29);

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

            auto ta_yzz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 24);

            auto ta_yzz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 25);

            auto ta_yzz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 26);

            auto ta_zzz_x_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 27);

            auto ta_zzz_y_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 28);

            auto ta_zzz_z_1 = primBuffer.data(pidx_a_3_1_m1 + 30 * idx + 29);

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

            auto ta_zz_x_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 15);

            auto ta_zz_y_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 16);

            auto ta_zz_z_0 = primBuffer.data(pidx_a_2_1_m0 + 18 * idx + 17);

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

            auto ta_zz_x_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 15);

            auto ta_zz_y_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 16);

            auto ta_zz_z_1 = primBuffer.data(pidx_a_2_1_m1 + 18 * idx + 17);

            auto ta_xxx_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx);

            auto ta_xxy_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 1);

            auto ta_xxz_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 2);

            auto ta_xyy_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 3);

            auto ta_xyz_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 4);

            auto ta_xzz_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 5);

            auto ta_yyy_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 6);

            auto ta_yyz_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 7);

            auto ta_yzz_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 8);

            auto ta_zzz_0_0 = primBuffer.data(pidx_a_3_0_m0 + 10 * idx + 9);

            auto ta_xxx_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx);

            auto ta_xxy_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 1);

            auto ta_xxz_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 2);

            auto ta_xyy_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 3);

            auto ta_xyz_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 4);

            auto ta_xzz_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 5);

            auto ta_yyy_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 6);

            auto ta_yyz_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 7);

            auto ta_yzz_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 8);

            auto ta_zzz_0_1 = primBuffer.data(pidx_a_3_0_m1 + 10 * idx + 9);

            // set up pointers to integrals

            auto ta_xxxx_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx);

            auto ta_xxxx_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 1);

            auto ta_xxxx_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 2);

            auto ta_xxxy_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 3);

            auto ta_xxxy_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 4);

            auto ta_xxxy_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 5);

            auto ta_xxxz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 6);

            auto ta_xxxz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 7);

            auto ta_xxxz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 8);

            auto ta_xxyy_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 9);

            auto ta_xxyy_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 10);

            auto ta_xxyy_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 11);

            auto ta_xxyz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 12);

            auto ta_xxyz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 13);

            auto ta_xxyz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 14);

            auto ta_xxzz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 15);

            auto ta_xxzz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 16);

            auto ta_xxzz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 17);

            auto ta_xyyy_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 18);

            auto ta_xyyy_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 19);

            auto ta_xyyy_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 20);

            auto ta_xyyz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 21);

            auto ta_xyyz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 22);

            auto ta_xyyz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 23);

            auto ta_xyzz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 24);

            auto ta_xyzz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 25);

            auto ta_xyzz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 26);

            auto ta_xzzz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 27);

            auto ta_xzzz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 28);

            auto ta_xzzz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 29);

            auto ta_yyyy_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 30);

            auto ta_yyyy_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 31);

            auto ta_yyyy_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 32);

            auto ta_yyyz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 33);

            auto ta_yyyz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 34);

            auto ta_yyyz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 35);

            auto ta_yyzz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 36);

            auto ta_yyzz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 37);

            auto ta_yyzz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 38);

            auto ta_yzzz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 39);

            auto ta_yzzz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 40);

            auto ta_yzzz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 41);

            auto ta_zzzz_x_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 42);

            auto ta_zzzz_y_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 43);

            auto ta_zzzz_z_0 = primBuffer.data(pidx_a_4_1_m0 + 45 * idx + 44);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xx_x_0, ta_xx_x_1, ta_xx_y_0, \
                                         ta_xx_y_1, ta_xx_z_0, ta_xx_z_1, ta_xxx_0_0, ta_xxx_0_1, ta_xxx_x_0, ta_xxx_x_1, \
                                         ta_xxx_y_0, ta_xxx_y_1, ta_xxx_z_0, ta_xxx_z_1, ta_xxxx_x_0, ta_xxxx_y_0, \
                                         ta_xxxx_z_0, ta_xxxy_x_0, ta_xxxy_y_0, ta_xxxy_z_0, ta_xxxz_x_0, ta_xxxz_y_0, \
                                         ta_xxxz_z_0, ta_xxy_0_0, ta_xxy_0_1, ta_xxy_x_0, ta_xxy_x_1, ta_xxy_y_0, ta_xxy_y_1, \
                                         ta_xxy_z_0, ta_xxy_z_1, ta_xxyy_x_0, ta_xxyy_y_0, ta_xxyy_z_0, ta_xxyz_x_0, \
                                         ta_xxyz_y_0, ta_xxyz_z_0, ta_xxz_0_0, ta_xxz_0_1, ta_xxz_x_0, ta_xxz_x_1, \
                                         ta_xxz_y_0, ta_xxz_y_1, ta_xxz_z_0, ta_xxz_z_1, ta_xxzz_x_0, ta_xxzz_y_0, \
                                         ta_xxzz_z_0, ta_xy_x_0, ta_xy_x_1, ta_xy_y_0, ta_xy_y_1, ta_xy_z_0, ta_xy_z_1, \
                                         ta_xyy_0_0, ta_xyy_0_1, ta_xyy_x_0, ta_xyy_x_1, ta_xyy_y_0, ta_xyy_y_1, ta_xyy_z_0, \
                                         ta_xyy_z_1, ta_xyyy_x_0, ta_xyyy_y_0, ta_xyyy_z_0, ta_xyyz_x_0, ta_xyyz_y_0, \
                                         ta_xyyz_z_0, ta_xyz_0_0, ta_xyz_0_1, ta_xyz_x_0, ta_xyz_x_1, ta_xyz_y_0, ta_xyz_y_1, \
                                         ta_xyz_z_0, ta_xyz_z_1, ta_xyzz_x_0, ta_xyzz_y_0, ta_xyzz_z_0, ta_xz_x_0, \
                                         ta_xz_x_1, ta_xz_y_0, ta_xz_y_1, ta_xz_z_0, ta_xz_z_1, ta_xzz_0_0, ta_xzz_0_1, \
                                         ta_xzz_x_0, ta_xzz_x_1, ta_xzz_y_0, ta_xzz_y_1, ta_xzz_z_0, ta_xzz_z_1, \
                                         ta_xzzz_x_0, ta_xzzz_y_0, ta_xzzz_z_0, ta_yy_x_0, ta_yy_x_1, ta_yy_y_0, ta_yy_y_1, \
                                         ta_yy_z_0, ta_yy_z_1, ta_yyy_0_0, ta_yyy_0_1, ta_yyy_x_0, ta_yyy_x_1, ta_yyy_y_0, \
                                         ta_yyy_y_1, ta_yyy_z_0, ta_yyy_z_1, ta_yyyy_x_0, ta_yyyy_y_0, ta_yyyy_z_0, \
                                         ta_yyyz_x_0, ta_yyyz_y_0, ta_yyyz_z_0, ta_yyz_0_0, ta_yyz_0_1, ta_yyz_x_0, \
                                         ta_yyz_x_1, ta_yyz_y_0, ta_yyz_y_1, ta_yyz_z_0, ta_yyz_z_1, ta_yyzz_x_0, \
                                         ta_yyzz_y_0, ta_yyzz_z_0, ta_yz_x_0, ta_yz_x_1, ta_yz_y_0, ta_yz_y_1, ta_yz_z_0, \
                                         ta_yz_z_1, ta_yzz_0_0, ta_yzz_0_1, ta_yzz_x_0, ta_yzz_x_1, ta_yzz_y_0, ta_yzz_y_1, \
                                         ta_yzz_z_0, ta_yzz_z_1, ta_yzzz_x_0, ta_yzzz_y_0, ta_yzzz_z_0, ta_zz_x_0, \
                                         ta_zz_x_1, ta_zz_y_0, ta_zz_y_1, ta_zz_z_0, ta_zz_z_1, ta_zzz_0_0, ta_zzz_0_1, \
                                         ta_zzz_x_0, ta_zzz_x_1, ta_zzz_y_0, ta_zzz_y_1, ta_zzz_z_0, ta_zzz_z_1, \
                                         ta_zzzz_x_0, ta_zzzz_y_0, ta_zzzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxxx_x_0[j] = pa_x[j] * ta_xxx_x_0[j] - pc_x[j] * ta_xxx_x_1[j] + 1.5 * fl1_fx * ta_xx_x_0[j] - 1.5 * fl1_fx * ta_xx_x_1[j] +
                                 0.5 * fl1_fx * ta_xxx_0_0[j] - 0.5 * fl1_fx * ta_xxx_0_1[j];

                ta_xxxx_y_0[j] = pa_x[j] * ta_xxx_y_0[j] - pc_x[j] * ta_xxx_y_1[j] + 1.5 * fl1_fx * ta_xx_y_0[j] - 1.5 * fl1_fx * ta_xx_y_1[j];

                ta_xxxx_z_0[j] = pa_x[j] * ta_xxx_z_0[j] - pc_x[j] * ta_xxx_z_1[j] + 1.5 * fl1_fx * ta_xx_z_0[j] - 1.5 * fl1_fx * ta_xx_z_1[j];

                ta_xxxy_x_0[j] = pa_x[j] * ta_xxy_x_0[j] - pc_x[j] * ta_xxy_x_1[j] + fl1_fx * ta_xy_x_0[j] - fl1_fx * ta_xy_x_1[j] +
                                 0.5 * fl1_fx * ta_xxy_0_0[j] - 0.5 * fl1_fx * ta_xxy_0_1[j];

                ta_xxxy_y_0[j] = pa_x[j] * ta_xxy_y_0[j] - pc_x[j] * ta_xxy_y_1[j] + fl1_fx * ta_xy_y_0[j] - fl1_fx * ta_xy_y_1[j];

                ta_xxxy_z_0[j] = pa_x[j] * ta_xxy_z_0[j] - pc_x[j] * ta_xxy_z_1[j] + fl1_fx * ta_xy_z_0[j] - fl1_fx * ta_xy_z_1[j];

                ta_xxxz_x_0[j] = pa_x[j] * ta_xxz_x_0[j] - pc_x[j] * ta_xxz_x_1[j] + fl1_fx * ta_xz_x_0[j] - fl1_fx * ta_xz_x_1[j] +
                                 0.5 * fl1_fx * ta_xxz_0_0[j] - 0.5 * fl1_fx * ta_xxz_0_1[j];

                ta_xxxz_y_0[j] = pa_x[j] * ta_xxz_y_0[j] - pc_x[j] * ta_xxz_y_1[j] + fl1_fx * ta_xz_y_0[j] - fl1_fx * ta_xz_y_1[j];

                ta_xxxz_z_0[j] = pa_x[j] * ta_xxz_z_0[j] - pc_x[j] * ta_xxz_z_1[j] + fl1_fx * ta_xz_z_0[j] - fl1_fx * ta_xz_z_1[j];

                ta_xxyy_x_0[j] = pa_x[j] * ta_xyy_x_0[j] - pc_x[j] * ta_xyy_x_1[j] + 0.5 * fl1_fx * ta_yy_x_0[j] - 0.5 * fl1_fx * ta_yy_x_1[j] +
                                 0.5 * fl1_fx * ta_xyy_0_0[j] - 0.5 * fl1_fx * ta_xyy_0_1[j];

                ta_xxyy_y_0[j] = pa_x[j] * ta_xyy_y_0[j] - pc_x[j] * ta_xyy_y_1[j] + 0.5 * fl1_fx * ta_yy_y_0[j] - 0.5 * fl1_fx * ta_yy_y_1[j];

                ta_xxyy_z_0[j] = pa_x[j] * ta_xyy_z_0[j] - pc_x[j] * ta_xyy_z_1[j] + 0.5 * fl1_fx * ta_yy_z_0[j] - 0.5 * fl1_fx * ta_yy_z_1[j];

                ta_xxyz_x_0[j] = pa_x[j] * ta_xyz_x_0[j] - pc_x[j] * ta_xyz_x_1[j] + 0.5 * fl1_fx * ta_yz_x_0[j] - 0.5 * fl1_fx * ta_yz_x_1[j] +
                                 0.5 * fl1_fx * ta_xyz_0_0[j] - 0.5 * fl1_fx * ta_xyz_0_1[j];

                ta_xxyz_y_0[j] = pa_x[j] * ta_xyz_y_0[j] - pc_x[j] * ta_xyz_y_1[j] + 0.5 * fl1_fx * ta_yz_y_0[j] - 0.5 * fl1_fx * ta_yz_y_1[j];

                ta_xxyz_z_0[j] = pa_x[j] * ta_xyz_z_0[j] - pc_x[j] * ta_xyz_z_1[j] + 0.5 * fl1_fx * ta_yz_z_0[j] - 0.5 * fl1_fx * ta_yz_z_1[j];

                ta_xxzz_x_0[j] = pa_x[j] * ta_xzz_x_0[j] - pc_x[j] * ta_xzz_x_1[j] + 0.5 * fl1_fx * ta_zz_x_0[j] - 0.5 * fl1_fx * ta_zz_x_1[j] +
                                 0.5 * fl1_fx * ta_xzz_0_0[j] - 0.5 * fl1_fx * ta_xzz_0_1[j];

                ta_xxzz_y_0[j] = pa_x[j] * ta_xzz_y_0[j] - pc_x[j] * ta_xzz_y_1[j] + 0.5 * fl1_fx * ta_zz_y_0[j] - 0.5 * fl1_fx * ta_zz_y_1[j];

                ta_xxzz_z_0[j] = pa_x[j] * ta_xzz_z_0[j] - pc_x[j] * ta_xzz_z_1[j] + 0.5 * fl1_fx * ta_zz_z_0[j] - 0.5 * fl1_fx * ta_zz_z_1[j];

                ta_xyyy_x_0[j] = pa_x[j] * ta_yyy_x_0[j] - pc_x[j] * ta_yyy_x_1[j] + 0.5 * fl1_fx * ta_yyy_0_0[j] - 0.5 * fl1_fx * ta_yyy_0_1[j];

                ta_xyyy_y_0[j] = pa_x[j] * ta_yyy_y_0[j] - pc_x[j] * ta_yyy_y_1[j];

                ta_xyyy_z_0[j] = pa_x[j] * ta_yyy_z_0[j] - pc_x[j] * ta_yyy_z_1[j];

                ta_xyyz_x_0[j] = pa_x[j] * ta_yyz_x_0[j] - pc_x[j] * ta_yyz_x_1[j] + 0.5 * fl1_fx * ta_yyz_0_0[j] - 0.5 * fl1_fx * ta_yyz_0_1[j];

                ta_xyyz_y_0[j] = pa_x[j] * ta_yyz_y_0[j] - pc_x[j] * ta_yyz_y_1[j];

                ta_xyyz_z_0[j] = pa_x[j] * ta_yyz_z_0[j] - pc_x[j] * ta_yyz_z_1[j];

                ta_xyzz_x_0[j] = pa_x[j] * ta_yzz_x_0[j] - pc_x[j] * ta_yzz_x_1[j] + 0.5 * fl1_fx * ta_yzz_0_0[j] - 0.5 * fl1_fx * ta_yzz_0_1[j];

                ta_xyzz_y_0[j] = pa_x[j] * ta_yzz_y_0[j] - pc_x[j] * ta_yzz_y_1[j];

                ta_xyzz_z_0[j] = pa_x[j] * ta_yzz_z_0[j] - pc_x[j] * ta_yzz_z_1[j];

                ta_xzzz_x_0[j] = pa_x[j] * ta_zzz_x_0[j] - pc_x[j] * ta_zzz_x_1[j] + 0.5 * fl1_fx * ta_zzz_0_0[j] - 0.5 * fl1_fx * ta_zzz_0_1[j];

                ta_xzzz_y_0[j] = pa_x[j] * ta_zzz_y_0[j] - pc_x[j] * ta_zzz_y_1[j];

                ta_xzzz_z_0[j] = pa_x[j] * ta_zzz_z_0[j] - pc_x[j] * ta_zzz_z_1[j];

                ta_yyyy_x_0[j] = pa_y[j] * ta_yyy_x_0[j] - pc_y[j] * ta_yyy_x_1[j] + 1.5 * fl1_fx * ta_yy_x_0[j] - 1.5 * fl1_fx * ta_yy_x_1[j];

                ta_yyyy_y_0[j] = pa_y[j] * ta_yyy_y_0[j] - pc_y[j] * ta_yyy_y_1[j] + 1.5 * fl1_fx * ta_yy_y_0[j] - 1.5 * fl1_fx * ta_yy_y_1[j] +
                                 0.5 * fl1_fx * ta_yyy_0_0[j] - 0.5 * fl1_fx * ta_yyy_0_1[j];

                ta_yyyy_z_0[j] = pa_y[j] * ta_yyy_z_0[j] - pc_y[j] * ta_yyy_z_1[j] + 1.5 * fl1_fx * ta_yy_z_0[j] - 1.5 * fl1_fx * ta_yy_z_1[j];

                ta_yyyz_x_0[j] = pa_y[j] * ta_yyz_x_0[j] - pc_y[j] * ta_yyz_x_1[j] + fl1_fx * ta_yz_x_0[j] - fl1_fx * ta_yz_x_1[j];

                ta_yyyz_y_0[j] = pa_y[j] * ta_yyz_y_0[j] - pc_y[j] * ta_yyz_y_1[j] + fl1_fx * ta_yz_y_0[j] - fl1_fx * ta_yz_y_1[j] +
                                 0.5 * fl1_fx * ta_yyz_0_0[j] - 0.5 * fl1_fx * ta_yyz_0_1[j];

                ta_yyyz_z_0[j] = pa_y[j] * ta_yyz_z_0[j] - pc_y[j] * ta_yyz_z_1[j] + fl1_fx * ta_yz_z_0[j] - fl1_fx * ta_yz_z_1[j];

                ta_yyzz_x_0[j] = pa_y[j] * ta_yzz_x_0[j] - pc_y[j] * ta_yzz_x_1[j] + 0.5 * fl1_fx * ta_zz_x_0[j] - 0.5 * fl1_fx * ta_zz_x_1[j];

                ta_yyzz_y_0[j] = pa_y[j] * ta_yzz_y_0[j] - pc_y[j] * ta_yzz_y_1[j] + 0.5 * fl1_fx * ta_zz_y_0[j] - 0.5 * fl1_fx * ta_zz_y_1[j] +
                                 0.5 * fl1_fx * ta_yzz_0_0[j] - 0.5 * fl1_fx * ta_yzz_0_1[j];

                ta_yyzz_z_0[j] = pa_y[j] * ta_yzz_z_0[j] - pc_y[j] * ta_yzz_z_1[j] + 0.5 * fl1_fx * ta_zz_z_0[j] - 0.5 * fl1_fx * ta_zz_z_1[j];

                ta_yzzz_x_0[j] = pa_y[j] * ta_zzz_x_0[j] - pc_y[j] * ta_zzz_x_1[j];

                ta_yzzz_y_0[j] = pa_y[j] * ta_zzz_y_0[j] - pc_y[j] * ta_zzz_y_1[j] + 0.5 * fl1_fx * ta_zzz_0_0[j] - 0.5 * fl1_fx * ta_zzz_0_1[j];

                ta_yzzz_z_0[j] = pa_y[j] * ta_zzz_z_0[j] - pc_y[j] * ta_zzz_z_1[j];

                ta_zzzz_x_0[j] = pa_z[j] * ta_zzz_x_0[j] - pc_z[j] * ta_zzz_x_1[j] + 1.5 * fl1_fx * ta_zz_x_0[j] - 1.5 * fl1_fx * ta_zz_x_1[j];

                ta_zzzz_y_0[j] = pa_z[j] * ta_zzz_y_0[j] - pc_z[j] * ta_zzz_y_1[j] + 1.5 * fl1_fx * ta_zz_y_0[j] - 1.5 * fl1_fx * ta_zz_y_1[j];

                ta_zzzz_z_0[j] = pa_z[j] * ta_zzz_z_0[j] - pc_z[j] * ta_zzz_z_1[j] + 1.5 * fl1_fx * ta_zz_z_0[j] - 1.5 * fl1_fx * ta_zz_z_1[j] +
                                 0.5 * fl1_fx * ta_zzz_0_0[j] - 0.5 * fl1_fx * ta_zzz_0_1[j];
            }

            idx++;
        }
    }
}

}  // namespace npotrecfunc
