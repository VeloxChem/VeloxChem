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

#include "ElectricFieldRecFuncForPX.hpp"

namespace efieldrecfunc {  // efieldrecfunc namespace

void
compElectricFieldForPP(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx);

            auto tey_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx);

            auto tez_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx);

            auto tex_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 1);

            auto tey_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 2);

            auto tey_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx);

            auto tey_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx);

            auto tez_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx);

            auto tex_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 1);

            auto tey_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 2);

            auto tey_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 2);

            auto tex_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + idx);

            auto tey_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + bdim + idx);

            auto tez_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + 2 * bdim + idx);

            auto tex_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + idx);

            auto tey_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + bdim + idx);

            auto tez_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + 2 * bdim + idx);

            auto ta_0_x_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx);

            auto ta_0_y_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 1);

            auto ta_0_z_1 = primBuffer.data(pidx_a_0_1_m1 + 3 * idx + 2);

            // set up pointers to integrals

            auto tex_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx);

            auto tey_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx);

            auto tez_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx);

            auto tex_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 1);

            auto tey_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 1);

            auto tez_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 1);

            auto tex_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 2);

            auto tey_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 2);

            auto tez_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 2);

            auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3);

            auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4);

            auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5);

            auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6);

            auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7);

            auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8);

            auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_x_1, ta_0_y_1, ta_0_z_1, tex_0_0_0, \
                                         tex_0_0_1, tex_0_x_0, tex_0_x_1, tex_0_y_0, tex_0_y_1, tex_0_z_0, tex_0_z_1, \
                                         tex_x_x_0, tex_x_y_0, tex_x_z_0, tex_y_x_0, tex_y_y_0, tex_y_z_0, tex_z_x_0, \
                                         tex_z_y_0, tex_z_z_0, tey_0_0_0, tey_0_0_1, tey_0_x_0, tey_0_x_1, tey_0_y_0, \
                                         tey_0_y_1, tey_0_z_0, tey_0_z_1, tey_x_x_0, tey_x_y_0, tey_x_z_0, tey_y_x_0, \
                                         tey_y_y_0, tey_y_z_0, tey_z_x_0, tey_z_y_0, tey_z_z_0, tez_0_0_0, tez_0_0_1, \
                                         tez_0_x_0, tez_0_x_1, tez_0_y_0, tez_0_y_1, tez_0_z_0, tez_0_z_1, tez_x_x_0, \
                                         tez_x_y_0, tez_x_z_0, tez_y_x_0, tez_y_y_0, tez_y_z_0, tez_z_x_0, tez_z_y_0, \
                                         tez_z_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_x_x_0[j] =
                    pa_x[j] * tex_0_x_0[j] - pc_x[j] * tex_0_x_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j] + ta_0_x_1[j];

                tey_x_x_0[j] = pa_x[j] * tey_0_x_0[j] - pc_x[j] * tey_0_x_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j];

                tez_x_x_0[j] = pa_x[j] * tez_0_x_0[j] - pc_x[j] * tez_0_x_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j];

                tex_x_y_0[j] = pa_x[j] * tex_0_y_0[j] - pc_x[j] * tex_0_y_1[j] + ta_0_y_1[j];

                tey_x_y_0[j] = pa_x[j] * tey_0_y_0[j] - pc_x[j] * tey_0_y_1[j];

                tez_x_y_0[j] = pa_x[j] * tez_0_y_0[j] - pc_x[j] * tez_0_y_1[j];

                tex_x_z_0[j] = pa_x[j] * tex_0_z_0[j] - pc_x[j] * tex_0_z_1[j] + ta_0_z_1[j];

                tey_x_z_0[j] = pa_x[j] * tey_0_z_0[j] - pc_x[j] * tey_0_z_1[j];

                tez_x_z_0[j] = pa_x[j] * tez_0_z_0[j] - pc_x[j] * tez_0_z_1[j];

                tex_y_x_0[j] = pa_y[j] * tex_0_x_0[j] - pc_y[j] * tex_0_x_1[j];

                tey_y_x_0[j] = pa_y[j] * tey_0_x_0[j] - pc_y[j] * tey_0_x_1[j] + ta_0_x_1[j];

                tez_y_x_0[j] = pa_y[j] * tez_0_x_0[j] - pc_y[j] * tez_0_x_1[j];

                tex_y_y_0[j] = pa_y[j] * tex_0_y_0[j] - pc_y[j] * tex_0_y_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j];

                tey_y_y_0[j] =
                    pa_y[j] * tey_0_y_0[j] - pc_y[j] * tey_0_y_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j] + ta_0_y_1[j];

                tez_y_y_0[j] = pa_y[j] * tez_0_y_0[j] - pc_y[j] * tez_0_y_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j];

                tex_y_z_0[j] = pa_y[j] * tex_0_z_0[j] - pc_y[j] * tex_0_z_1[j];

                tey_y_z_0[j] = pa_y[j] * tey_0_z_0[j] - pc_y[j] * tey_0_z_1[j] + ta_0_z_1[j];

                tez_y_z_0[j] = pa_y[j] * tez_0_z_0[j] - pc_y[j] * tez_0_z_1[j];

                tex_z_x_0[j] = pa_z[j] * tex_0_x_0[j] - pc_z[j] * tex_0_x_1[j];

                tey_z_x_0[j] = pa_z[j] * tey_0_x_0[j] - pc_z[j] * tey_0_x_1[j];

                tez_z_x_0[j] = pa_z[j] * tez_0_x_0[j] - pc_z[j] * tez_0_x_1[j] + ta_0_x_1[j];

                tex_z_y_0[j] = pa_z[j] * tex_0_y_0[j] - pc_z[j] * tex_0_y_1[j];

                tey_z_y_0[j] = pa_z[j] * tey_0_y_0[j] - pc_z[j] * tey_0_y_1[j];

                tez_z_y_0[j] = pa_z[j] * tez_0_y_0[j] - pc_z[j] * tez_0_y_1[j] + ta_0_y_1[j];

                tex_z_z_0[j] = pa_z[j] * tex_0_z_0[j] - pc_z[j] * tex_0_z_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j];

                tey_z_z_0[j] = pa_z[j] * tey_0_z_0[j] - pc_z[j] * tey_0_z_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j];

                tez_z_z_0[j] =
                    pa_z[j] * tez_0_z_0[j] - pc_z[j] * tez_0_z_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j] + ta_0_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPD(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForPD_0_27(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPD_27_54(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForPD_0_27(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CMemBlock2D<double>& pcDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
{
    // Batch of Integrals (0,27)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

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

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx);

            auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx);

            auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx);

            auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1);

            auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2);

            auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3);

            auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4);

            auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5);

            auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5);

            auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx);

            auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx);

            auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx);

            auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1);

            auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2);

            auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3);

            auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4);

            auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5);

            auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5);

            auto tex_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx);

            auto tey_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx);

            auto tez_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx);

            auto tex_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 1);

            auto tey_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 2);

            auto tey_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx);

            auto tey_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx);

            auto tez_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx);

            auto tex_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 1);

            auto tey_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 2);

            auto tey_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 2);

            auto ta_0_xx_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx);

            auto ta_0_xy_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 1);

            auto ta_0_xz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 2);

            auto ta_0_yy_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 3);

            auto ta_0_yz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 4);

            auto ta_0_zz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 5);

            // set up pointers to integrals

            auto tex_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx);

            auto tey_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx);

            auto tez_x_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx);

            auto tex_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 1);

            auto tey_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 1);

            auto tez_x_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 1);

            auto tex_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 2);

            auto tey_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 2);

            auto tez_x_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 2);

            auto tex_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 3);

            auto tey_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 3);

            auto tez_x_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 3);

            auto tex_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 4);

            auto tey_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 4);

            auto tez_x_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 4);

            auto tex_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 5);

            auto tey_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 5);

            auto tez_x_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 5);

            auto tex_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 6);

            auto tey_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 6);

            auto tez_y_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 6);

            auto tex_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 7);

            auto tey_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 7);

            auto tez_y_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 7);

            auto tex_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 8);

            auto tey_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 8);

            auto tez_y_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 8);

            // Batch of Integrals (0,27)

            #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_0_xx_1, ta_0_xy_1, ta_0_xz_1, ta_0_yy_1, \
                                         ta_0_yz_1, ta_0_zz_1, tex_0_x_0, tex_0_x_1, tex_0_xx_0, tex_0_xx_1, tex_0_xy_0, \
                                         tex_0_xy_1, tex_0_xz_0, tex_0_xz_1, tex_0_y_0, tex_0_y_1, tex_0_yy_0, tex_0_yy_1, \
                                         tex_0_yz_0, tex_0_yz_1, tex_0_z_0, tex_0_z_1, tex_0_zz_0, tex_0_zz_1, tex_x_xx_0, \
                                         tex_x_xy_0, tex_x_xz_0, tex_x_yy_0, tex_x_yz_0, tex_x_zz_0, tex_y_xx_0, tex_y_xy_0, \
                                         tex_y_xz_0, tey_0_x_0, tey_0_x_1, tey_0_xx_0, tey_0_xx_1, tey_0_xy_0, tey_0_xy_1, \
                                         tey_0_xz_0, tey_0_xz_1, tey_0_y_0, tey_0_y_1, tey_0_yy_0, tey_0_yy_1, tey_0_yz_0, \
                                         tey_0_yz_1, tey_0_z_0, tey_0_z_1, tey_0_zz_0, tey_0_zz_1, tey_x_xx_0, tey_x_xy_0, \
                                         tey_x_xz_0, tey_x_yy_0, tey_x_yz_0, tey_x_zz_0, tey_y_xx_0, tey_y_xy_0, tey_y_xz_0, \
                                         tez_0_x_0, tez_0_x_1, tez_0_xx_0, tez_0_xx_1, tez_0_xy_0, tez_0_xy_1, tez_0_xz_0, \
                                         tez_0_xz_1, tez_0_y_0, tez_0_y_1, tez_0_yy_0, tez_0_yy_1, tez_0_yz_0, tez_0_yz_1, \
                                         tez_0_z_0, tez_0_z_1, tez_0_zz_0, tez_0_zz_1, tez_x_xx_0, tez_x_xy_0, tez_x_xz_0, \
                                         tez_x_yy_0, tez_x_yz_0, tez_x_zz_0, tez_y_xx_0, tez_y_xy_0, tez_y_xz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_x_xx_0[j] = pa_x[j] * tex_0_xx_0[j] - pc_x[j] * tex_0_xx_1[j] + fl1_fx * tex_0_x_0[j] - fl1_fx * tex_0_x_1[j] + ta_0_xx_1[j];

                tey_x_xx_0[j] = pa_x[j] * tey_0_xx_0[j] - pc_x[j] * tey_0_xx_1[j] + fl1_fx * tey_0_x_0[j] - fl1_fx * tey_0_x_1[j];

                tez_x_xx_0[j] = pa_x[j] * tez_0_xx_0[j] - pc_x[j] * tez_0_xx_1[j] + fl1_fx * tez_0_x_0[j] - fl1_fx * tez_0_x_1[j];

                tex_x_xy_0[j] =
                    pa_x[j] * tex_0_xy_0[j] - pc_x[j] * tex_0_xy_1[j] + 0.5 * fl1_fx * tex_0_y_0[j] - 0.5 * fl1_fx * tex_0_y_1[j] + ta_0_xy_1[j];

                tey_x_xy_0[j] = pa_x[j] * tey_0_xy_0[j] - pc_x[j] * tey_0_xy_1[j] + 0.5 * fl1_fx * tey_0_y_0[j] - 0.5 * fl1_fx * tey_0_y_1[j];

                tez_x_xy_0[j] = pa_x[j] * tez_0_xy_0[j] - pc_x[j] * tez_0_xy_1[j] + 0.5 * fl1_fx * tez_0_y_0[j] - 0.5 * fl1_fx * tez_0_y_1[j];

                tex_x_xz_0[j] =
                    pa_x[j] * tex_0_xz_0[j] - pc_x[j] * tex_0_xz_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j] + ta_0_xz_1[j];

                tey_x_xz_0[j] = pa_x[j] * tey_0_xz_0[j] - pc_x[j] * tey_0_xz_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j];

                tez_x_xz_0[j] = pa_x[j] * tez_0_xz_0[j] - pc_x[j] * tez_0_xz_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j];

                tex_x_yy_0[j] = pa_x[j] * tex_0_yy_0[j] - pc_x[j] * tex_0_yy_1[j] + ta_0_yy_1[j];

                tey_x_yy_0[j] = pa_x[j] * tey_0_yy_0[j] - pc_x[j] * tey_0_yy_1[j];

                tez_x_yy_0[j] = pa_x[j] * tez_0_yy_0[j] - pc_x[j] * tez_0_yy_1[j];

                tex_x_yz_0[j] = pa_x[j] * tex_0_yz_0[j] - pc_x[j] * tex_0_yz_1[j] + ta_0_yz_1[j];

                tey_x_yz_0[j] = pa_x[j] * tey_0_yz_0[j] - pc_x[j] * tey_0_yz_1[j];

                tez_x_yz_0[j] = pa_x[j] * tez_0_yz_0[j] - pc_x[j] * tez_0_yz_1[j];

                tex_x_zz_0[j] = pa_x[j] * tex_0_zz_0[j] - pc_x[j] * tex_0_zz_1[j] + ta_0_zz_1[j];

                tey_x_zz_0[j] = pa_x[j] * tey_0_zz_0[j] - pc_x[j] * tey_0_zz_1[j];

                tez_x_zz_0[j] = pa_x[j] * tez_0_zz_0[j] - pc_x[j] * tez_0_zz_1[j];

                tex_y_xx_0[j] = pa_y[j] * tex_0_xx_0[j] - pc_y[j] * tex_0_xx_1[j];

                tey_y_xx_0[j] = pa_y[j] * tey_0_xx_0[j] - pc_y[j] * tey_0_xx_1[j] + ta_0_xx_1[j];

                tez_y_xx_0[j] = pa_y[j] * tez_0_xx_0[j] - pc_y[j] * tez_0_xx_1[j];

                tex_y_xy_0[j] = pa_y[j] * tex_0_xy_0[j] - pc_y[j] * tex_0_xy_1[j] + 0.5 * fl1_fx * tex_0_x_0[j] - 0.5 * fl1_fx * tex_0_x_1[j];

                tey_y_xy_0[j] =
                    pa_y[j] * tey_0_xy_0[j] - pc_y[j] * tey_0_xy_1[j] + 0.5 * fl1_fx * tey_0_x_0[j] - 0.5 * fl1_fx * tey_0_x_1[j] + ta_0_xy_1[j];

                tez_y_xy_0[j] = pa_y[j] * tez_0_xy_0[j] - pc_y[j] * tez_0_xy_1[j] + 0.5 * fl1_fx * tez_0_x_0[j] - 0.5 * fl1_fx * tez_0_x_1[j];

                tex_y_xz_0[j] = pa_y[j] * tex_0_xz_0[j] - pc_y[j] * tex_0_xz_1[j];

                tey_y_xz_0[j] = pa_y[j] * tey_0_xz_0[j] - pc_y[j] * tey_0_xz_1[j] + ta_0_xz_1[j];

                tez_y_xz_0[j] = pa_y[j] * tez_0_xz_0[j] - pc_y[j] * tez_0_xz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPD_27_54(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pcDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
{
    // Batch of Integrals (27,54)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_2_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx);

            auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx);

            auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx);

            auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1);

            auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2);

            auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3);

            auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4);

            auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5);

            auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5);

            auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx);

            auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx);

            auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx);

            auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1);

            auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2);

            auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3);

            auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4);

            auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5);

            auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5);

            auto tex_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx);

            auto tey_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx);

            auto tez_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx);

            auto tex_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 1);

            auto tey_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 2);

            auto tey_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx);

            auto tey_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx);

            auto tez_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx);

            auto tex_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 1);

            auto tey_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 2);

            auto tey_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 2);

            auto ta_0_xx_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx);

            auto ta_0_xy_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 1);

            auto ta_0_xz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 2);

            auto ta_0_yy_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 3);

            auto ta_0_yz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 4);

            auto ta_0_zz_1 = primBuffer.data(pidx_a_0_2_m1 + 6 * idx + 5);

            // set up pointers to integrals

            auto tex_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 9);

            auto tey_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 9);

            auto tez_y_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 9);

            auto tex_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 10);

            auto tey_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 10);

            auto tez_y_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 10);

            auto tex_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 11);

            auto tey_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 11);

            auto tez_y_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 11);

            auto tex_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 12);

            auto tey_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 12);

            auto tez_z_xx_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 12);

            auto tex_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 13);

            auto tey_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 13);

            auto tez_z_xy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 13);

            auto tex_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 14);

            auto tey_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 14);

            auto tez_z_xz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 14);

            auto tex_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 15);

            auto tey_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 15);

            auto tez_z_yy_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 15);

            auto tex_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 16);

            auto tey_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 16);

            auto tez_z_yz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 16);

            auto tex_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * idx + 17);

            auto tey_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 18 * bdim + 18 * idx + 17);

            auto tez_z_zz_0 = primBuffer.data(pidx_e_1_2_m0 + 36 * bdim + 18 * idx + 17);

            // Batch of Integrals (27,54)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_0_xx_1, ta_0_xy_1, ta_0_xz_1, ta_0_yy_1, \
                                         ta_0_yz_1, ta_0_zz_1, tex_0_x_0, tex_0_x_1, tex_0_xx_0, tex_0_xx_1, tex_0_xy_0, \
                                         tex_0_xy_1, tex_0_xz_0, tex_0_xz_1, tex_0_y_0, tex_0_y_1, tex_0_yy_0, tex_0_yy_1, \
                                         tex_0_yz_0, tex_0_yz_1, tex_0_z_0, tex_0_z_1, tex_0_zz_0, tex_0_zz_1, tex_y_yy_0, \
                                         tex_y_yz_0, tex_y_zz_0, tex_z_xx_0, tex_z_xy_0, tex_z_xz_0, tex_z_yy_0, tex_z_yz_0, \
                                         tex_z_zz_0, tey_0_x_0, tey_0_x_1, tey_0_xx_0, tey_0_xx_1, tey_0_xy_0, tey_0_xy_1, \
                                         tey_0_xz_0, tey_0_xz_1, tey_0_y_0, tey_0_y_1, tey_0_yy_0, tey_0_yy_1, tey_0_yz_0, \
                                         tey_0_yz_1, tey_0_z_0, tey_0_z_1, tey_0_zz_0, tey_0_zz_1, tey_y_yy_0, tey_y_yz_0, \
                                         tey_y_zz_0, tey_z_xx_0, tey_z_xy_0, tey_z_xz_0, tey_z_yy_0, tey_z_yz_0, tey_z_zz_0, \
                                         tez_0_x_0, tez_0_x_1, tez_0_xx_0, tez_0_xx_1, tez_0_xy_0, tez_0_xy_1, tez_0_xz_0, \
                                         tez_0_xz_1, tez_0_y_0, tez_0_y_1, tez_0_yy_0, tez_0_yy_1, tez_0_yz_0, tez_0_yz_1, \
                                         tez_0_z_0, tez_0_z_1, tez_0_zz_0, tez_0_zz_1, tez_y_yy_0, tez_y_yz_0, tez_y_zz_0, \
                                         tez_z_xx_0, tez_z_xy_0, tez_z_xz_0, tez_z_yy_0, tez_z_yz_0, tez_z_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_y_yy_0[j] = pa_y[j] * tex_0_yy_0[j] - pc_y[j] * tex_0_yy_1[j] + fl1_fx * tex_0_y_0[j] - fl1_fx * tex_0_y_1[j];

                tey_y_yy_0[j] = pa_y[j] * tey_0_yy_0[j] - pc_y[j] * tey_0_yy_1[j] + fl1_fx * tey_0_y_0[j] - fl1_fx * tey_0_y_1[j] + ta_0_yy_1[j];

                tez_y_yy_0[j] = pa_y[j] * tez_0_yy_0[j] - pc_y[j] * tez_0_yy_1[j] + fl1_fx * tez_0_y_0[j] - fl1_fx * tez_0_y_1[j];

                tex_y_yz_0[j] = pa_y[j] * tex_0_yz_0[j] - pc_y[j] * tex_0_yz_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j];

                tey_y_yz_0[j] =
                    pa_y[j] * tey_0_yz_0[j] - pc_y[j] * tey_0_yz_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j] + ta_0_yz_1[j];

                tez_y_yz_0[j] = pa_y[j] * tez_0_yz_0[j] - pc_y[j] * tez_0_yz_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j];

                tex_y_zz_0[j] = pa_y[j] * tex_0_zz_0[j] - pc_y[j] * tex_0_zz_1[j];

                tey_y_zz_0[j] = pa_y[j] * tey_0_zz_0[j] - pc_y[j] * tey_0_zz_1[j] + ta_0_zz_1[j];

                tez_y_zz_0[j] = pa_y[j] * tez_0_zz_0[j] - pc_y[j] * tez_0_zz_1[j];

                tex_z_xx_0[j] = pa_z[j] * tex_0_xx_0[j] - pc_z[j] * tex_0_xx_1[j];

                tey_z_xx_0[j] = pa_z[j] * tey_0_xx_0[j] - pc_z[j] * tey_0_xx_1[j];

                tez_z_xx_0[j] = pa_z[j] * tez_0_xx_0[j] - pc_z[j] * tez_0_xx_1[j] + ta_0_xx_1[j];

                tex_z_xy_0[j] = pa_z[j] * tex_0_xy_0[j] - pc_z[j] * tex_0_xy_1[j];

                tey_z_xy_0[j] = pa_z[j] * tey_0_xy_0[j] - pc_z[j] * tey_0_xy_1[j];

                tez_z_xy_0[j] = pa_z[j] * tez_0_xy_0[j] - pc_z[j] * tez_0_xy_1[j] + ta_0_xy_1[j];

                tex_z_xz_0[j] = pa_z[j] * tex_0_xz_0[j] - pc_z[j] * tex_0_xz_1[j] + 0.5 * fl1_fx * tex_0_x_0[j] - 0.5 * fl1_fx * tex_0_x_1[j];

                tey_z_xz_0[j] = pa_z[j] * tey_0_xz_0[j] - pc_z[j] * tey_0_xz_1[j] + 0.5 * fl1_fx * tey_0_x_0[j] - 0.5 * fl1_fx * tey_0_x_1[j];

                tez_z_xz_0[j] =
                    pa_z[j] * tez_0_xz_0[j] - pc_z[j] * tez_0_xz_1[j] + 0.5 * fl1_fx * tez_0_x_0[j] - 0.5 * fl1_fx * tez_0_x_1[j] + ta_0_xz_1[j];

                tex_z_yy_0[j] = pa_z[j] * tex_0_yy_0[j] - pc_z[j] * tex_0_yy_1[j];

                tey_z_yy_0[j] = pa_z[j] * tey_0_yy_0[j] - pc_z[j] * tey_0_yy_1[j];

                tez_z_yy_0[j] = pa_z[j] * tez_0_yy_0[j] - pc_z[j] * tez_0_yy_1[j] + ta_0_yy_1[j];

                tex_z_yz_0[j] = pa_z[j] * tex_0_yz_0[j] - pc_z[j] * tex_0_yz_1[j] + 0.5 * fl1_fx * tex_0_y_0[j] - 0.5 * fl1_fx * tex_0_y_1[j];

                tey_z_yz_0[j] = pa_z[j] * tey_0_yz_0[j] - pc_z[j] * tey_0_yz_1[j] + 0.5 * fl1_fx * tey_0_y_0[j] - 0.5 * fl1_fx * tey_0_y_1[j];

                tez_z_yz_0[j] =
                    pa_z[j] * tez_0_yz_0[j] - pc_z[j] * tez_0_yz_1[j] + 0.5 * fl1_fx * tez_0_y_0[j] - 0.5 * fl1_fx * tez_0_y_1[j] + ta_0_yz_1[j];

                tex_z_zz_0[j] = pa_z[j] * tex_0_zz_0[j] - pc_z[j] * tex_0_zz_1[j] + fl1_fx * tex_0_z_0[j] - fl1_fx * tex_0_z_1[j];

                tey_z_zz_0[j] = pa_z[j] * tey_0_zz_0[j] - pc_z[j] * tey_0_zz_1[j] + fl1_fx * tey_0_z_0[j] - fl1_fx * tey_0_z_1[j];

                tez_z_zz_0[j] = pa_z[j] * tez_0_zz_0[j] - pc_z[j] * tez_0_zz_1[j] + fl1_fx * tez_0_z_0[j] - fl1_fx * tez_0_z_1[j] + ta_0_zz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForDP(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForDP_0_27(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForDP_27_54(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForDP_0_27(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
                            const CMemBlock2D<double>& pcDistances,
                            const CGtoBlock&           braGtoBlock,
                            const CGtoBlock&           ketGtoBlock,
                            const int32_t              iContrGto)
{
    // Batch of Integrals (0,27)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_2_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            // set up pointers to auxilary integrals

            auto tex_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx);

            auto tey_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx);

            auto tez_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx);

            auto tex_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 1);

            auto tey_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 1);

            auto tez_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 1);

            auto tex_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 2);

            auto tey_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 2);

            auto tez_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 2);

            auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3);

            auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4);

            auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5);

            auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6);

            auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7);

            auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8);

            auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8);

            auto tex_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx);

            auto tey_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx);

            auto tez_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx);

            auto tex_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 1);

            auto tey_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 1);

            auto tez_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 1);

            auto tex_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 2);

            auto tey_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 2);

            auto tez_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 2);

            auto tex_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 3);

            auto tey_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 4);

            auto tey_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 5);

            auto tey_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 6);

            auto tey_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 7);

            auto tey_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 8);

            auto tey_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 8);

            auto tex_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx);

            auto tey_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx);

            auto tez_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx);

            auto tex_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 1);

            auto tey_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 2);

            auto tey_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx);

            auto tey_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx);

            auto tez_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx);

            auto tex_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 1);

            auto tey_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 2);

            auto tey_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 2);

            auto tex_x_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx);

            auto tey_x_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx);

            auto tez_x_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx);

            auto tex_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx + 1);

            auto tey_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx + 2);

            auto tey_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_x_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * idx);

            auto tey_x_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * bdim + 3 * idx);

            auto tez_x_0_1 = primBuffer.data(pidx_e_1_0_m1 + 6 * bdim + 3 * idx);

            auto tex_y_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * idx + 1);

            auto tey_y_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_y_0_1 = primBuffer.data(pidx_e_1_0_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_z_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * idx + 2);

            auto tey_z_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_z_0_1 = primBuffer.data(pidx_e_1_0_m1 + 6 * bdim + 3 * idx + 2);

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

            auto tex_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx);

            auto tey_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx);

            auto tez_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx);

            auto tex_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 1);

            auto tey_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 1);

            auto tez_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 1);

            auto tex_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 2);

            auto tey_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 2);

            auto tez_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 2);

            auto tex_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 3);

            auto tey_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 3);

            auto tez_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 3);

            auto tex_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 4);

            auto tey_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 4);

            auto tez_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 4);

            auto tex_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 5);

            auto tey_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 5);

            auto tez_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 5);

            auto tex_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 6);

            auto tey_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 6);

            auto tez_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 6);

            auto tex_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 7);

            auto tey_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 7);

            auto tez_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 7);

            auto tex_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 8);

            auto tey_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 8);

            auto tez_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 8);

            // Batch of Integrals (0,27)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_x_x_1, ta_x_y_1, ta_x_z_1, ta_y_x_1, ta_y_y_1, ta_y_z_1, \
                                         ta_z_x_1, ta_z_y_1, ta_z_z_1, tex_0_x_0, tex_0_x_1, tex_0_y_0, tex_0_y_1, \
                                         tex_0_z_0, tex_0_z_1, tex_x_0_0, tex_x_0_1, tex_x_x_0, tex_x_x_1, tex_x_y_0, \
                                         tex_x_y_1, tex_x_z_0, tex_x_z_1, tex_xx_x_0, tex_xx_y_0, tex_xx_z_0, tex_xy_x_0, \
                                         tex_xy_y_0, tex_xy_z_0, tex_xz_x_0, tex_xz_y_0, tex_xz_z_0, tex_y_0_0, tex_y_0_1, \
                                         tex_y_x_0, tex_y_x_1, tex_y_y_0, tex_y_y_1, tex_y_z_0, tex_y_z_1, tex_z_0_0, \
                                         tex_z_0_1, tex_z_x_0, tex_z_x_1, tex_z_y_0, tex_z_y_1, tex_z_z_0, tex_z_z_1, \
                                         tey_0_x_0, tey_0_x_1, tey_0_y_0, tey_0_y_1, tey_0_z_0, tey_0_z_1, tey_x_0_0, \
                                         tey_x_0_1, tey_x_x_0, tey_x_x_1, tey_x_y_0, tey_x_y_1, tey_x_z_0, tey_x_z_1, \
                                         tey_xx_x_0, tey_xx_y_0, tey_xx_z_0, tey_xy_x_0, tey_xy_y_0, tey_xy_z_0, tey_xz_x_0, \
                                         tey_xz_y_0, tey_xz_z_0, tey_y_0_0, tey_y_0_1, tey_y_x_0, tey_y_x_1, tey_y_y_0, \
                                         tey_y_y_1, tey_y_z_0, tey_y_z_1, tey_z_0_0, tey_z_0_1, tey_z_x_0, tey_z_x_1, \
                                         tey_z_y_0, tey_z_y_1, tey_z_z_0, tey_z_z_1, tez_0_x_0, tez_0_x_1, tez_0_y_0, \
                                         tez_0_y_1, tez_0_z_0, tez_0_z_1, tez_x_0_0, tez_x_0_1, tez_x_x_0, tez_x_x_1, \
                                         tez_x_y_0, tez_x_y_1, tez_x_z_0, tez_x_z_1, tez_xx_x_0, tez_xx_y_0, tez_xx_z_0, \
                                         tez_xy_x_0, tez_xy_y_0, tez_xy_z_0, tez_xz_x_0, tez_xz_y_0, tez_xz_z_0, tez_y_0_0, \
                                         tez_y_0_1, tez_y_x_0, tez_y_x_1, tez_y_y_0, tez_y_y_1, tez_y_z_0, tez_y_z_1, \
                                         tez_z_0_0, tez_z_0_1, tez_z_x_0, tez_z_x_1, tez_z_y_0, tez_z_y_1, tez_z_z_0, \
                                         tez_z_z_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xx_x_0[j] = pa_x[j] * tex_x_x_0[j] - pc_x[j] * tex_x_x_1[j] + 0.5 * fl1_fx * tex_0_x_0[j] - 0.5 * fl1_fx * tex_0_x_1[j] +
                                0.5 * fl1_fx * tex_x_0_0[j] - 0.5 * fl1_fx * tex_x_0_1[j] + ta_x_x_1[j];

                tey_xx_x_0[j] = pa_x[j] * tey_x_x_0[j] - pc_x[j] * tey_x_x_1[j] + 0.5 * fl1_fx * tey_0_x_0[j] - 0.5 * fl1_fx * tey_0_x_1[j] +
                                0.5 * fl1_fx * tey_x_0_0[j] - 0.5 * fl1_fx * tey_x_0_1[j];

                tez_xx_x_0[j] = pa_x[j] * tez_x_x_0[j] - pc_x[j] * tez_x_x_1[j] + 0.5 * fl1_fx * tez_0_x_0[j] - 0.5 * fl1_fx * tez_0_x_1[j] +
                                0.5 * fl1_fx * tez_x_0_0[j] - 0.5 * fl1_fx * tez_x_0_1[j];

                tex_xx_y_0[j] =
                    pa_x[j] * tex_x_y_0[j] - pc_x[j] * tex_x_y_1[j] + 0.5 * fl1_fx * tex_0_y_0[j] - 0.5 * fl1_fx * tex_0_y_1[j] + ta_x_y_1[j];

                tey_xx_y_0[j] = pa_x[j] * tey_x_y_0[j] - pc_x[j] * tey_x_y_1[j] + 0.5 * fl1_fx * tey_0_y_0[j] - 0.5 * fl1_fx * tey_0_y_1[j];

                tez_xx_y_0[j] = pa_x[j] * tez_x_y_0[j] - pc_x[j] * tez_x_y_1[j] + 0.5 * fl1_fx * tez_0_y_0[j] - 0.5 * fl1_fx * tez_0_y_1[j];

                tex_xx_z_0[j] =
                    pa_x[j] * tex_x_z_0[j] - pc_x[j] * tex_x_z_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j] + ta_x_z_1[j];

                tey_xx_z_0[j] = pa_x[j] * tey_x_z_0[j] - pc_x[j] * tey_x_z_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j];

                tez_xx_z_0[j] = pa_x[j] * tez_x_z_0[j] - pc_x[j] * tez_x_z_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j];

                tex_xy_x_0[j] =
                    pa_x[j] * tex_y_x_0[j] - pc_x[j] * tex_y_x_1[j] + 0.5 * fl1_fx * tex_y_0_0[j] - 0.5 * fl1_fx * tex_y_0_1[j] + ta_y_x_1[j];

                tey_xy_x_0[j] = pa_x[j] * tey_y_x_0[j] - pc_x[j] * tey_y_x_1[j] + 0.5 * fl1_fx * tey_y_0_0[j] - 0.5 * fl1_fx * tey_y_0_1[j];

                tez_xy_x_0[j] = pa_x[j] * tez_y_x_0[j] - pc_x[j] * tez_y_x_1[j] + 0.5 * fl1_fx * tez_y_0_0[j] - 0.5 * fl1_fx * tez_y_0_1[j];

                tex_xy_y_0[j] = pa_x[j] * tex_y_y_0[j] - pc_x[j] * tex_y_y_1[j] + ta_y_y_1[j];

                tey_xy_y_0[j] = pa_x[j] * tey_y_y_0[j] - pc_x[j] * tey_y_y_1[j];

                tez_xy_y_0[j] = pa_x[j] * tez_y_y_0[j] - pc_x[j] * tez_y_y_1[j];

                tex_xy_z_0[j] = pa_x[j] * tex_y_z_0[j] - pc_x[j] * tex_y_z_1[j] + ta_y_z_1[j];

                tey_xy_z_0[j] = pa_x[j] * tey_y_z_0[j] - pc_x[j] * tey_y_z_1[j];

                tez_xy_z_0[j] = pa_x[j] * tez_y_z_0[j] - pc_x[j] * tez_y_z_1[j];

                tex_xz_x_0[j] =
                    pa_x[j] * tex_z_x_0[j] - pc_x[j] * tex_z_x_1[j] + 0.5 * fl1_fx * tex_z_0_0[j] - 0.5 * fl1_fx * tex_z_0_1[j] + ta_z_x_1[j];

                tey_xz_x_0[j] = pa_x[j] * tey_z_x_0[j] - pc_x[j] * tey_z_x_1[j] + 0.5 * fl1_fx * tey_z_0_0[j] - 0.5 * fl1_fx * tey_z_0_1[j];

                tez_xz_x_0[j] = pa_x[j] * tez_z_x_0[j] - pc_x[j] * tez_z_x_1[j] + 0.5 * fl1_fx * tez_z_0_0[j] - 0.5 * fl1_fx * tez_z_0_1[j];

                tex_xz_y_0[j] = pa_x[j] * tex_z_y_0[j] - pc_x[j] * tex_z_y_1[j] + ta_z_y_1[j];

                tey_xz_y_0[j] = pa_x[j] * tey_z_y_0[j] - pc_x[j] * tey_z_y_1[j];

                tez_xz_y_0[j] = pa_x[j] * tez_z_y_0[j] - pc_x[j] * tez_z_y_1[j];

                tex_xz_z_0[j] = pa_x[j] * tex_z_z_0[j] - pc_x[j] * tex_z_z_1[j] + ta_z_z_1[j];

                tey_xz_z_0[j] = pa_x[j] * tey_z_z_0[j] - pc_x[j] * tey_z_z_1[j];

                tez_xz_z_0[j] = pa_x[j] * tez_z_z_0[j] - pc_x[j] * tez_z_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForDP_27_54(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
                             const CMemBlock2D<double>& pcDistances,
                             const CGtoBlock&           braGtoBlock,
                             const CGtoBlock&           ketGtoBlock,
                             const int32_t              iContrGto)
{
    // Batch of Integrals (27,54)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_2_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3);

            auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4);

            auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5);

            auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6);

            auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7);

            auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8);

            auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8);

            auto tex_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 3);

            auto tey_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 4);

            auto tey_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 5);

            auto tey_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 6);

            auto tey_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 7);

            auto tey_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 8);

            auto tey_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 8);

            auto tex_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx);

            auto tey_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx);

            auto tez_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx);

            auto tex_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 1);

            auto tey_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 2);

            auto tey_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx);

            auto tey_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx);

            auto tez_0_x_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx);

            auto tex_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 1);

            auto tey_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * idx + 2);

            auto tey_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_1 = primBuffer.data(pidx_e_0_1_m1 + 6 * bdim + 3 * idx + 2);

            auto tex_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx + 1);

            auto tey_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx + 2);

            auto tey_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx + 2);

            auto tex_y_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * idx + 1);

            auto tey_y_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * bdim + 3 * idx + 1);

            auto tez_y_0_1 = primBuffer.data(pidx_e_1_0_m1 + 6 * bdim + 3 * idx + 1);

            auto tex_z_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * idx + 2);

            auto tey_z_0_1 = primBuffer.data(pidx_e_1_0_m1 + 3 * bdim + 3 * idx + 2);

            auto tez_z_0_1 = primBuffer.data(pidx_e_1_0_m1 + 6 * bdim + 3 * idx + 2);

            auto ta_y_x_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 3);

            auto ta_y_y_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 4);

            auto ta_y_z_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 5);

            auto ta_z_x_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 6);

            auto ta_z_y_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 7);

            auto ta_z_z_1 = primBuffer.data(pidx_a_1_1_m1 + 9 * idx + 8);

            // set up pointers to integrals

            auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9);

            auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10);

            auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11);

            auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12);

            auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13);

            auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14);

            auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14);

            auto tex_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 15);

            auto tey_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 16);

            auto tey_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 17);

            auto tey_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 17);

            // Batch of Integrals (27,54)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_y_x_1, ta_y_y_1, ta_y_z_1, ta_z_x_1, ta_z_y_1, \
                                         ta_z_z_1, tex_0_x_0, tex_0_x_1, tex_0_y_0, tex_0_y_1, tex_0_z_0, tex_0_z_1, \
                                         tex_y_0_0, tex_y_0_1, tex_y_x_0, tex_y_x_1, tex_y_y_0, tex_y_y_1, tex_y_z_0, \
                                         tex_y_z_1, tex_yy_x_0, tex_yy_y_0, tex_yy_z_0, tex_yz_x_0, tex_yz_y_0, tex_yz_z_0, \
                                         tex_z_0_0, tex_z_0_1, tex_z_x_0, tex_z_x_1, tex_z_y_0, tex_z_y_1, tex_z_z_0, \
                                         tex_z_z_1, tex_zz_x_0, tex_zz_y_0, tex_zz_z_0, tey_0_x_0, tey_0_x_1, tey_0_y_0, \
                                         tey_0_y_1, tey_0_z_0, tey_0_z_1, tey_y_0_0, tey_y_0_1, tey_y_x_0, tey_y_x_1, \
                                         tey_y_y_0, tey_y_y_1, tey_y_z_0, tey_y_z_1, tey_yy_x_0, tey_yy_y_0, tey_yy_z_0, \
                                         tey_yz_x_0, tey_yz_y_0, tey_yz_z_0, tey_z_0_0, tey_z_0_1, tey_z_x_0, tey_z_x_1, \
                                         tey_z_y_0, tey_z_y_1, tey_z_z_0, tey_z_z_1, tey_zz_x_0, tey_zz_y_0, tey_zz_z_0, \
                                         tez_0_x_0, tez_0_x_1, tez_0_y_0, tez_0_y_1, tez_0_z_0, tez_0_z_1, tez_y_0_0, \
                                         tez_y_0_1, tez_y_x_0, tez_y_x_1, tez_y_y_0, tez_y_y_1, tez_y_z_0, tez_y_z_1, \
                                         tez_yy_x_0, tez_yy_y_0, tez_yy_z_0, tez_yz_x_0, tez_yz_y_0, tez_yz_z_0, tez_z_0_0, \
                                         tez_z_0_1, tez_z_x_0, tez_z_x_1, tez_z_y_0, tez_z_y_1, tez_z_z_0, tez_z_z_1, \
                                         tez_zz_x_0, tez_zz_y_0, tez_zz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yy_x_0[j] = pa_y[j] * tex_y_x_0[j] - pc_y[j] * tex_y_x_1[j] + 0.5 * fl1_fx * tex_0_x_0[j] - 0.5 * fl1_fx * tex_0_x_1[j];

                tey_yy_x_0[j] =
                    pa_y[j] * tey_y_x_0[j] - pc_y[j] * tey_y_x_1[j] + 0.5 * fl1_fx * tey_0_x_0[j] - 0.5 * fl1_fx * tey_0_x_1[j] + ta_y_x_1[j];

                tez_yy_x_0[j] = pa_y[j] * tez_y_x_0[j] - pc_y[j] * tez_y_x_1[j] + 0.5 * fl1_fx * tez_0_x_0[j] - 0.5 * fl1_fx * tez_0_x_1[j];

                tex_yy_y_0[j] = pa_y[j] * tex_y_y_0[j] - pc_y[j] * tex_y_y_1[j] + 0.5 * fl1_fx * tex_0_y_0[j] - 0.5 * fl1_fx * tex_0_y_1[j] +
                                0.5 * fl1_fx * tex_y_0_0[j] - 0.5 * fl1_fx * tex_y_0_1[j];

                tey_yy_y_0[j] = pa_y[j] * tey_y_y_0[j] - pc_y[j] * tey_y_y_1[j] + 0.5 * fl1_fx * tey_0_y_0[j] - 0.5 * fl1_fx * tey_0_y_1[j] +
                                0.5 * fl1_fx * tey_y_0_0[j] - 0.5 * fl1_fx * tey_y_0_1[j] + ta_y_y_1[j];

                tez_yy_y_0[j] = pa_y[j] * tez_y_y_0[j] - pc_y[j] * tez_y_y_1[j] + 0.5 * fl1_fx * tez_0_y_0[j] - 0.5 * fl1_fx * tez_0_y_1[j] +
                                0.5 * fl1_fx * tez_y_0_0[j] - 0.5 * fl1_fx * tez_y_0_1[j];

                tex_yy_z_0[j] = pa_y[j] * tex_y_z_0[j] - pc_y[j] * tex_y_z_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j];

                tey_yy_z_0[j] =
                    pa_y[j] * tey_y_z_0[j] - pc_y[j] * tey_y_z_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j] + ta_y_z_1[j];

                tez_yy_z_0[j] = pa_y[j] * tez_y_z_0[j] - pc_y[j] * tez_y_z_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j];

                tex_yz_x_0[j] = pa_y[j] * tex_z_x_0[j] - pc_y[j] * tex_z_x_1[j];

                tey_yz_x_0[j] = pa_y[j] * tey_z_x_0[j] - pc_y[j] * tey_z_x_1[j] + ta_z_x_1[j];

                tez_yz_x_0[j] = pa_y[j] * tez_z_x_0[j] - pc_y[j] * tez_z_x_1[j];

                tex_yz_y_0[j] = pa_y[j] * tex_z_y_0[j] - pc_y[j] * tex_z_y_1[j] + 0.5 * fl1_fx * tex_z_0_0[j] - 0.5 * fl1_fx * tex_z_0_1[j];

                tey_yz_y_0[j] =
                    pa_y[j] * tey_z_y_0[j] - pc_y[j] * tey_z_y_1[j] + 0.5 * fl1_fx * tey_z_0_0[j] - 0.5 * fl1_fx * tey_z_0_1[j] + ta_z_y_1[j];

                tez_yz_y_0[j] = pa_y[j] * tez_z_y_0[j] - pc_y[j] * tez_z_y_1[j] + 0.5 * fl1_fx * tez_z_0_0[j] - 0.5 * fl1_fx * tez_z_0_1[j];

                tex_yz_z_0[j] = pa_y[j] * tex_z_z_0[j] - pc_y[j] * tex_z_z_1[j];

                tey_yz_z_0[j] = pa_y[j] * tey_z_z_0[j] - pc_y[j] * tey_z_z_1[j] + ta_z_z_1[j];

                tez_yz_z_0[j] = pa_y[j] * tez_z_z_0[j] - pc_y[j] * tez_z_z_1[j];

                tex_zz_x_0[j] = pa_z[j] * tex_z_x_0[j] - pc_z[j] * tex_z_x_1[j] + 0.5 * fl1_fx * tex_0_x_0[j] - 0.5 * fl1_fx * tex_0_x_1[j];

                tey_zz_x_0[j] = pa_z[j] * tey_z_x_0[j] - pc_z[j] * tey_z_x_1[j] + 0.5 * fl1_fx * tey_0_x_0[j] - 0.5 * fl1_fx * tey_0_x_1[j];

                tez_zz_x_0[j] =
                    pa_z[j] * tez_z_x_0[j] - pc_z[j] * tez_z_x_1[j] + 0.5 * fl1_fx * tez_0_x_0[j] - 0.5 * fl1_fx * tez_0_x_1[j] + ta_z_x_1[j];

                tex_zz_y_0[j] = pa_z[j] * tex_z_y_0[j] - pc_z[j] * tex_z_y_1[j] + 0.5 * fl1_fx * tex_0_y_0[j] - 0.5 * fl1_fx * tex_0_y_1[j];

                tey_zz_y_0[j] = pa_z[j] * tey_z_y_0[j] - pc_z[j] * tey_z_y_1[j] + 0.5 * fl1_fx * tey_0_y_0[j] - 0.5 * fl1_fx * tey_0_y_1[j];

                tez_zz_y_0[j] =
                    pa_z[j] * tez_z_y_0[j] - pc_z[j] * tez_z_y_1[j] + 0.5 * fl1_fx * tez_0_y_0[j] - 0.5 * fl1_fx * tez_0_y_1[j] + ta_z_y_1[j];

                tex_zz_z_0[j] = pa_z[j] * tex_z_z_0[j] - pc_z[j] * tex_z_z_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j] +
                                0.5 * fl1_fx * tex_z_0_0[j] - 0.5 * fl1_fx * tex_z_0_1[j];

                tey_zz_z_0[j] = pa_z[j] * tey_z_z_0[j] - pc_z[j] * tey_z_z_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j] +
                                0.5 * fl1_fx * tey_z_0_0[j] - 0.5 * fl1_fx * tey_z_0_1[j];

                tez_zz_z_0[j] = pa_z[j] * tez_z_z_0[j] - pc_z[j] * tez_z_z_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j] +
                                0.5 * fl1_fx * tez_z_0_0[j] - 0.5 * fl1_fx * tez_z_0_1[j] + ta_z_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPF(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForPF_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPF_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForPF_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx);

            auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1);

            auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2);

            auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3);

            auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4);

            auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5);

            auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6);

            auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7);

            auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8);

            auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9);

            auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx);

            auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1);

            auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2);

            auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3);

            auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4);

            auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5);

            auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6);

            auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7);

            auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8);

            auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9);

            auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx);

            auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx);

            auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx);

            auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1);

            auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2);

            auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3);

            auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4);

            auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5);

            auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5);

            auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx);

            auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx);

            auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx);

            auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1);

            auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2);

            auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3);

            auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4);

            auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5);

            auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5);

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

            auto tex_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx);

            auto tey_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx);

            auto tez_x_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx);

            auto tex_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 1);

            auto tey_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 1);

            auto tez_x_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 1);

            auto tex_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 2);

            auto tey_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 2);

            auto tez_x_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 2);

            auto tex_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 3);

            auto tey_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 3);

            auto tez_x_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 3);

            auto tex_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 4);

            auto tey_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 4);

            auto tez_x_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 4);

            auto tex_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 5);

            auto tey_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 5);

            auto tez_x_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 5);

            auto tex_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 6);

            auto tey_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 6);

            auto tez_x_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 6);

            auto tex_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 7);

            auto tey_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 7);

            auto tez_x_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 7);

            auto tex_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 8);

            auto tey_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 8);

            auto tez_x_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 8);

            auto tex_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 9);

            auto tey_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 9);

            auto tez_x_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 9);

            auto tex_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 10);

            auto tey_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 10);

            auto tez_y_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 10);

            auto tex_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 11);

            auto tey_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 11);

            auto tez_y_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 11);

            auto tex_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 12);

            auto tey_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 12);

            auto tez_y_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 12);

            auto tex_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 13);

            auto tey_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 13);

            auto tez_y_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 13);

            auto tex_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 14);

            auto tey_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 14);

            auto tez_y_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 14);

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pa_y, pc_x, pc_y, ta_0_xxx_1, ta_0_xxy_1, ta_0_xxz_1, ta_0_xyy_1, \
                                         ta_0_xyz_1, ta_0_xzz_1, ta_0_yyy_1, ta_0_yyz_1, ta_0_yzz_1, ta_0_zzz_1, tex_0_xx_0, \
                                         tex_0_xx_1, tex_0_xxx_0, tex_0_xxx_1, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxz_0, \
                                         tex_0_xxz_1, tex_0_xy_0, tex_0_xy_1, tex_0_xyy_0, tex_0_xyy_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xz_0, tex_0_xz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_yy_0, \
                                         tex_0_yy_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyz_0, tex_0_yyz_1, tex_0_yz_0, \
                                         tex_0_yz_1, tex_0_yzz_0, tex_0_yzz_1, tex_0_zz_0, tex_0_zz_1, tex_0_zzz_0, \
                                         tex_0_zzz_1, tex_x_xxx_0, tex_x_xxy_0, tex_x_xxz_0, tex_x_xyy_0, tex_x_xyz_0, \
                                         tex_x_xzz_0, tex_x_yyy_0, tex_x_yyz_0, tex_x_yzz_0, tex_x_zzz_0, tex_y_xxx_0, \
                                         tex_y_xxy_0, tex_y_xxz_0, tex_y_xyy_0, tex_y_xyz_0, tey_0_xx_0, tey_0_xx_1, \
                                         tey_0_xxx_0, tey_0_xxx_1, tey_0_xxy_0, tey_0_xxy_1, tey_0_xxz_0, tey_0_xxz_1, \
                                         tey_0_xy_0, tey_0_xy_1, tey_0_xyy_0, tey_0_xyy_1, tey_0_xyz_0, tey_0_xyz_1, \
                                         tey_0_xz_0, tey_0_xz_1, tey_0_xzz_0, tey_0_xzz_1, tey_0_yy_0, tey_0_yy_1, \
                                         tey_0_yyy_0, tey_0_yyy_1, tey_0_yyz_0, tey_0_yyz_1, tey_0_yz_0, tey_0_yz_1, \
                                         tey_0_yzz_0, tey_0_yzz_1, tey_0_zz_0, tey_0_zz_1, tey_0_zzz_0, tey_0_zzz_1, \
                                         tey_x_xxx_0, tey_x_xxy_0, tey_x_xxz_0, tey_x_xyy_0, tey_x_xyz_0, tey_x_xzz_0, \
                                         tey_x_yyy_0, tey_x_yyz_0, tey_x_yzz_0, tey_x_zzz_0, tey_y_xxx_0, tey_y_xxy_0, \
                                         tey_y_xxz_0, tey_y_xyy_0, tey_y_xyz_0, tez_0_xx_0, tez_0_xx_1, tez_0_xxx_0, \
                                         tez_0_xxx_1, tez_0_xxy_0, tez_0_xxy_1, tez_0_xxz_0, tez_0_xxz_1, tez_0_xy_0, \
                                         tez_0_xy_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyz_0, tez_0_xyz_1, tez_0_xz_0, \
                                         tez_0_xz_1, tez_0_xzz_0, tez_0_xzz_1, tez_0_yy_0, tez_0_yy_1, tez_0_yyy_0, \
                                         tez_0_yyy_1, tez_0_yyz_0, tez_0_yyz_1, tez_0_yz_0, tez_0_yz_1, tez_0_yzz_0, \
                                         tez_0_yzz_1, tez_0_zz_0, tez_0_zz_1, tez_0_zzz_0, tez_0_zzz_1, tez_x_xxx_0, \
                                         tez_x_xxy_0, tez_x_xxz_0, tez_x_xyy_0, tez_x_xyz_0, tez_x_xzz_0, tez_x_yyy_0, \
                                         tez_x_yyz_0, tez_x_yzz_0, tez_x_zzz_0, tez_y_xxx_0, tez_y_xxy_0, tez_y_xxz_0, \
                                         tez_y_xyy_0, tez_y_xyz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_x_xxx_0[j] =
                    pa_x[j] * tex_0_xxx_0[j] - pc_x[j] * tex_0_xxx_1[j] + 1.5 * fl1_fx * tex_0_xx_0[j] - 1.5 * fl1_fx * tex_0_xx_1[j] + ta_0_xxx_1[j];

                tey_x_xxx_0[j] = pa_x[j] * tey_0_xxx_0[j] - pc_x[j] * tey_0_xxx_1[j] + 1.5 * fl1_fx * tey_0_xx_0[j] - 1.5 * fl1_fx * tey_0_xx_1[j];

                tez_x_xxx_0[j] = pa_x[j] * tez_0_xxx_0[j] - pc_x[j] * tez_0_xxx_1[j] + 1.5 * fl1_fx * tez_0_xx_0[j] - 1.5 * fl1_fx * tez_0_xx_1[j];

                tex_x_xxy_0[j] =
                    pa_x[j] * tex_0_xxy_0[j] - pc_x[j] * tex_0_xxy_1[j] + fl1_fx * tex_0_xy_0[j] - fl1_fx * tex_0_xy_1[j] + ta_0_xxy_1[j];

                tey_x_xxy_0[j] = pa_x[j] * tey_0_xxy_0[j] - pc_x[j] * tey_0_xxy_1[j] + fl1_fx * tey_0_xy_0[j] - fl1_fx * tey_0_xy_1[j];

                tez_x_xxy_0[j] = pa_x[j] * tez_0_xxy_0[j] - pc_x[j] * tez_0_xxy_1[j] + fl1_fx * tez_0_xy_0[j] - fl1_fx * tez_0_xy_1[j];

                tex_x_xxz_0[j] =
                    pa_x[j] * tex_0_xxz_0[j] - pc_x[j] * tex_0_xxz_1[j] + fl1_fx * tex_0_xz_0[j] - fl1_fx * tex_0_xz_1[j] + ta_0_xxz_1[j];

                tey_x_xxz_0[j] = pa_x[j] * tey_0_xxz_0[j] - pc_x[j] * tey_0_xxz_1[j] + fl1_fx * tey_0_xz_0[j] - fl1_fx * tey_0_xz_1[j];

                tez_x_xxz_0[j] = pa_x[j] * tez_0_xxz_0[j] - pc_x[j] * tez_0_xxz_1[j] + fl1_fx * tez_0_xz_0[j] - fl1_fx * tez_0_xz_1[j];

                tex_x_xyy_0[j] =
                    pa_x[j] * tex_0_xyy_0[j] - pc_x[j] * tex_0_xyy_1[j] + 0.5 * fl1_fx * tex_0_yy_0[j] - 0.5 * fl1_fx * tex_0_yy_1[j] + ta_0_xyy_1[j];

                tey_x_xyy_0[j] = pa_x[j] * tey_0_xyy_0[j] - pc_x[j] * tey_0_xyy_1[j] + 0.5 * fl1_fx * tey_0_yy_0[j] - 0.5 * fl1_fx * tey_0_yy_1[j];

                tez_x_xyy_0[j] = pa_x[j] * tez_0_xyy_0[j] - pc_x[j] * tez_0_xyy_1[j] + 0.5 * fl1_fx * tez_0_yy_0[j] - 0.5 * fl1_fx * tez_0_yy_1[j];

                tex_x_xyz_0[j] =
                    pa_x[j] * tex_0_xyz_0[j] - pc_x[j] * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_0_yz_0[j] - 0.5 * fl1_fx * tex_0_yz_1[j] + ta_0_xyz_1[j];

                tey_x_xyz_0[j] = pa_x[j] * tey_0_xyz_0[j] - pc_x[j] * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_0_yz_0[j] - 0.5 * fl1_fx * tey_0_yz_1[j];

                tez_x_xyz_0[j] = pa_x[j] * tez_0_xyz_0[j] - pc_x[j] * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_0_yz_0[j] - 0.5 * fl1_fx * tez_0_yz_1[j];

                tex_x_xzz_0[j] =
                    pa_x[j] * tex_0_xzz_0[j] - pc_x[j] * tex_0_xzz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j] + ta_0_xzz_1[j];

                tey_x_xzz_0[j] = pa_x[j] * tey_0_xzz_0[j] - pc_x[j] * tey_0_xzz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j];

                tez_x_xzz_0[j] = pa_x[j] * tez_0_xzz_0[j] - pc_x[j] * tez_0_xzz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j];

                tex_x_yyy_0[j] = pa_x[j] * tex_0_yyy_0[j] - pc_x[j] * tex_0_yyy_1[j] + ta_0_yyy_1[j];

                tey_x_yyy_0[j] = pa_x[j] * tey_0_yyy_0[j] - pc_x[j] * tey_0_yyy_1[j];

                tez_x_yyy_0[j] = pa_x[j] * tez_0_yyy_0[j] - pc_x[j] * tez_0_yyy_1[j];

                tex_x_yyz_0[j] = pa_x[j] * tex_0_yyz_0[j] - pc_x[j] * tex_0_yyz_1[j] + ta_0_yyz_1[j];

                tey_x_yyz_0[j] = pa_x[j] * tey_0_yyz_0[j] - pc_x[j] * tey_0_yyz_1[j];

                tez_x_yyz_0[j] = pa_x[j] * tez_0_yyz_0[j] - pc_x[j] * tez_0_yyz_1[j];

                tex_x_yzz_0[j] = pa_x[j] * tex_0_yzz_0[j] - pc_x[j] * tex_0_yzz_1[j] + ta_0_yzz_1[j];

                tey_x_yzz_0[j] = pa_x[j] * tey_0_yzz_0[j] - pc_x[j] * tey_0_yzz_1[j];

                tez_x_yzz_0[j] = pa_x[j] * tez_0_yzz_0[j] - pc_x[j] * tez_0_yzz_1[j];

                tex_x_zzz_0[j] = pa_x[j] * tex_0_zzz_0[j] - pc_x[j] * tex_0_zzz_1[j] + ta_0_zzz_1[j];

                tey_x_zzz_0[j] = pa_x[j] * tey_0_zzz_0[j] - pc_x[j] * tey_0_zzz_1[j];

                tez_x_zzz_0[j] = pa_x[j] * tez_0_zzz_0[j] - pc_x[j] * tez_0_zzz_1[j];

                tex_y_xxx_0[j] = pa_y[j] * tex_0_xxx_0[j] - pc_y[j] * tex_0_xxx_1[j];

                tey_y_xxx_0[j] = pa_y[j] * tey_0_xxx_0[j] - pc_y[j] * tey_0_xxx_1[j] + ta_0_xxx_1[j];

                tez_y_xxx_0[j] = pa_y[j] * tez_0_xxx_0[j] - pc_y[j] * tez_0_xxx_1[j];

                tex_y_xxy_0[j] = pa_y[j] * tex_0_xxy_0[j] - pc_y[j] * tex_0_xxy_1[j] + 0.5 * fl1_fx * tex_0_xx_0[j] - 0.5 * fl1_fx * tex_0_xx_1[j];

                tey_y_xxy_0[j] =
                    pa_y[j] * tey_0_xxy_0[j] - pc_y[j] * tey_0_xxy_1[j] + 0.5 * fl1_fx * tey_0_xx_0[j] - 0.5 * fl1_fx * tey_0_xx_1[j] + ta_0_xxy_1[j];

                tez_y_xxy_0[j] = pa_y[j] * tez_0_xxy_0[j] - pc_y[j] * tez_0_xxy_1[j] + 0.5 * fl1_fx * tez_0_xx_0[j] - 0.5 * fl1_fx * tez_0_xx_1[j];

                tex_y_xxz_0[j] = pa_y[j] * tex_0_xxz_0[j] - pc_y[j] * tex_0_xxz_1[j];

                tey_y_xxz_0[j] = pa_y[j] * tey_0_xxz_0[j] - pc_y[j] * tey_0_xxz_1[j] + ta_0_xxz_1[j];

                tez_y_xxz_0[j] = pa_y[j] * tez_0_xxz_0[j] - pc_y[j] * tez_0_xxz_1[j];

                tex_y_xyy_0[j] = pa_y[j] * tex_0_xyy_0[j] - pc_y[j] * tex_0_xyy_1[j] + fl1_fx * tex_0_xy_0[j] - fl1_fx * tex_0_xy_1[j];

                tey_y_xyy_0[j] =
                    pa_y[j] * tey_0_xyy_0[j] - pc_y[j] * tey_0_xyy_1[j] + fl1_fx * tey_0_xy_0[j] - fl1_fx * tey_0_xy_1[j] + ta_0_xyy_1[j];

                tez_y_xyy_0[j] = pa_y[j] * tez_0_xyy_0[j] - pc_y[j] * tez_0_xyy_1[j] + fl1_fx * tez_0_xy_0[j] - fl1_fx * tez_0_xy_1[j];

                tex_y_xyz_0[j] = pa_y[j] * tex_0_xyz_0[j] - pc_y[j] * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_0_xz_0[j] - 0.5 * fl1_fx * tex_0_xz_1[j];

                tey_y_xyz_0[j] =
                    pa_y[j] * tey_0_xyz_0[j] - pc_y[j] * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_0_xz_0[j] - 0.5 * fl1_fx * tey_0_xz_1[j] + ta_0_xyz_1[j];

                tez_y_xyz_0[j] = pa_y[j] * tez_0_xyz_0[j] - pc_y[j] * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_0_xz_0[j] - 0.5 * fl1_fx * tez_0_xz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPF_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_3_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_2_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx);

            auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1);

            auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2);

            auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3);

            auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4);

            auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5);

            auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6);

            auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7);

            auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8);

            auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9);

            auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx);

            auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1);

            auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2);

            auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3);

            auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4);

            auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5);

            auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6);

            auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7);

            auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8);

            auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9);

            auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx);

            auto tey_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx);

            auto tez_0_xx_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx);

            auto tex_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 1);

            auto tey_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 2);

            auto tey_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 3);

            auto tey_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 4);

            auto tey_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * idx + 5);

            auto tey_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_0 = primBuffer.data(pidx_e_0_2_m0 + 12 * bdim + 6 * idx + 5);

            auto tex_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx);

            auto tey_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx);

            auto tez_0_xx_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx);

            auto tex_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 1);

            auto tey_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 1);

            auto tez_0_xy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 1);

            auto tex_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 2);

            auto tey_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 2);

            auto tez_0_xz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 2);

            auto tex_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 3);

            auto tey_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 3);

            auto tez_0_yy_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 3);

            auto tex_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 4);

            auto tey_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 4);

            auto tez_0_yz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 4);

            auto tex_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * idx + 5);

            auto tey_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_0_zz_1 = primBuffer.data(pidx_e_0_2_m1 + 12 * bdim + 6 * idx + 5);

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

            auto tex_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 15);

            auto tey_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 15);

            auto tez_y_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 15);

            auto tex_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 16);

            auto tey_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 16);

            auto tez_y_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 16);

            auto tex_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 17);

            auto tey_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 17);

            auto tez_y_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 17);

            auto tex_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 18);

            auto tey_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 18);

            auto tez_y_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 18);

            auto tex_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 19);

            auto tey_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 19);

            auto tez_y_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 19);

            auto tex_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 20);

            auto tey_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 20);

            auto tez_z_xxx_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 20);

            auto tex_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 21);

            auto tey_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 21);

            auto tez_z_xxy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 21);

            auto tex_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 22);

            auto tey_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 22);

            auto tez_z_xxz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 22);

            auto tex_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 23);

            auto tey_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 23);

            auto tez_z_xyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 23);

            auto tex_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 24);

            auto tey_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 24);

            auto tez_z_xyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 24);

            auto tex_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 25);

            auto tey_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 25);

            auto tez_z_xzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 25);

            auto tex_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 26);

            auto tey_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 26);

            auto tez_z_yyy_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 26);

            auto tex_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 27);

            auto tey_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 27);

            auto tez_z_yyz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 27);

            auto tex_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 28);

            auto tey_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 28);

            auto tez_z_yzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 28);

            auto tex_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * idx + 29);

            auto tey_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 30 * bdim + 30 * idx + 29);

            auto tez_z_zzz_0 = primBuffer.data(pidx_e_1_3_m0 + 60 * bdim + 30 * idx + 29);

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_0_xxx_1, ta_0_xxy_1, ta_0_xxz_1, ta_0_xyy_1, \
                                         ta_0_xyz_1, ta_0_xzz_1, ta_0_yyy_1, ta_0_yyz_1, ta_0_yzz_1, ta_0_zzz_1, tex_0_xx_0, \
                                         tex_0_xx_1, tex_0_xxx_0, tex_0_xxx_1, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxz_0, \
                                         tex_0_xxz_1, tex_0_xy_0, tex_0_xy_1, tex_0_xyy_0, tex_0_xyy_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xz_0, tex_0_xz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_yy_0, \
                                         tex_0_yy_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyz_0, tex_0_yyz_1, tex_0_yz_0, \
                                         tex_0_yz_1, tex_0_yzz_0, tex_0_yzz_1, tex_0_zz_0, tex_0_zz_1, tex_0_zzz_0, \
                                         tex_0_zzz_1, tex_y_xzz_0, tex_y_yyy_0, tex_y_yyz_0, tex_y_yzz_0, tex_y_zzz_0, \
                                         tex_z_xxx_0, tex_z_xxy_0, tex_z_xxz_0, tex_z_xyy_0, tex_z_xyz_0, tex_z_xzz_0, \
                                         tex_z_yyy_0, tex_z_yyz_0, tex_z_yzz_0, tex_z_zzz_0, tey_0_xx_0, tey_0_xx_1, \
                                         tey_0_xxx_0, tey_0_xxx_1, tey_0_xxy_0, tey_0_xxy_1, tey_0_xxz_0, tey_0_xxz_1, \
                                         tey_0_xy_0, tey_0_xy_1, tey_0_xyy_0, tey_0_xyy_1, tey_0_xyz_0, tey_0_xyz_1, \
                                         tey_0_xz_0, tey_0_xz_1, tey_0_xzz_0, tey_0_xzz_1, tey_0_yy_0, tey_0_yy_1, \
                                         tey_0_yyy_0, tey_0_yyy_1, tey_0_yyz_0, tey_0_yyz_1, tey_0_yz_0, tey_0_yz_1, \
                                         tey_0_yzz_0, tey_0_yzz_1, tey_0_zz_0, tey_0_zz_1, tey_0_zzz_0, tey_0_zzz_1, \
                                         tey_y_xzz_0, tey_y_yyy_0, tey_y_yyz_0, tey_y_yzz_0, tey_y_zzz_0, tey_z_xxx_0, \
                                         tey_z_xxy_0, tey_z_xxz_0, tey_z_xyy_0, tey_z_xyz_0, tey_z_xzz_0, tey_z_yyy_0, \
                                         tey_z_yyz_0, tey_z_yzz_0, tey_z_zzz_0, tez_0_xx_0, tez_0_xx_1, tez_0_xxx_0, \
                                         tez_0_xxx_1, tez_0_xxy_0, tez_0_xxy_1, tez_0_xxz_0, tez_0_xxz_1, tez_0_xy_0, \
                                         tez_0_xy_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyz_0, tez_0_xyz_1, tez_0_xz_0, \
                                         tez_0_xz_1, tez_0_xzz_0, tez_0_xzz_1, tez_0_yy_0, tez_0_yy_1, tez_0_yyy_0, \
                                         tez_0_yyy_1, tez_0_yyz_0, tez_0_yyz_1, tez_0_yz_0, tez_0_yz_1, tez_0_yzz_0, \
                                         tez_0_yzz_1, tez_0_zz_0, tez_0_zz_1, tez_0_zzz_0, tez_0_zzz_1, tez_y_xzz_0, \
                                         tez_y_yyy_0, tez_y_yyz_0, tez_y_yzz_0, tez_y_zzz_0, tez_z_xxx_0, tez_z_xxy_0, \
                                         tez_z_xxz_0, tez_z_xyy_0, tez_z_xyz_0, tez_z_xzz_0, tez_z_yyy_0, tez_z_yyz_0, \
                                         tez_z_yzz_0, tez_z_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_y_xzz_0[j] = pa_y[j] * tex_0_xzz_0[j] - pc_y[j] * tex_0_xzz_1[j];

                tey_y_xzz_0[j] = pa_y[j] * tey_0_xzz_0[j] - pc_y[j] * tey_0_xzz_1[j] + ta_0_xzz_1[j];

                tez_y_xzz_0[j] = pa_y[j] * tez_0_xzz_0[j] - pc_y[j] * tez_0_xzz_1[j];

                tex_y_yyy_0[j] = pa_y[j] * tex_0_yyy_0[j] - pc_y[j] * tex_0_yyy_1[j] + 1.5 * fl1_fx * tex_0_yy_0[j] - 1.5 * fl1_fx * tex_0_yy_1[j];

                tey_y_yyy_0[j] =
                    pa_y[j] * tey_0_yyy_0[j] - pc_y[j] * tey_0_yyy_1[j] + 1.5 * fl1_fx * tey_0_yy_0[j] - 1.5 * fl1_fx * tey_0_yy_1[j] + ta_0_yyy_1[j];

                tez_y_yyy_0[j] = pa_y[j] * tez_0_yyy_0[j] - pc_y[j] * tez_0_yyy_1[j] + 1.5 * fl1_fx * tez_0_yy_0[j] - 1.5 * fl1_fx * tez_0_yy_1[j];

                tex_y_yyz_0[j] = pa_y[j] * tex_0_yyz_0[j] - pc_y[j] * tex_0_yyz_1[j] + fl1_fx * tex_0_yz_0[j] - fl1_fx * tex_0_yz_1[j];

                tey_y_yyz_0[j] =
                    pa_y[j] * tey_0_yyz_0[j] - pc_y[j] * tey_0_yyz_1[j] + fl1_fx * tey_0_yz_0[j] - fl1_fx * tey_0_yz_1[j] + ta_0_yyz_1[j];

                tez_y_yyz_0[j] = pa_y[j] * tez_0_yyz_0[j] - pc_y[j] * tez_0_yyz_1[j] + fl1_fx * tez_0_yz_0[j] - fl1_fx * tez_0_yz_1[j];

                tex_y_yzz_0[j] = pa_y[j] * tex_0_yzz_0[j] - pc_y[j] * tex_0_yzz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j];

                tey_y_yzz_0[j] =
                    pa_y[j] * tey_0_yzz_0[j] - pc_y[j] * tey_0_yzz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j] + ta_0_yzz_1[j];

                tez_y_yzz_0[j] = pa_y[j] * tez_0_yzz_0[j] - pc_y[j] * tez_0_yzz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j];

                tex_y_zzz_0[j] = pa_y[j] * tex_0_zzz_0[j] - pc_y[j] * tex_0_zzz_1[j];

                tey_y_zzz_0[j] = pa_y[j] * tey_0_zzz_0[j] - pc_y[j] * tey_0_zzz_1[j] + ta_0_zzz_1[j];

                tez_y_zzz_0[j] = pa_y[j] * tez_0_zzz_0[j] - pc_y[j] * tez_0_zzz_1[j];

                tex_z_xxx_0[j] = pa_z[j] * tex_0_xxx_0[j] - pc_z[j] * tex_0_xxx_1[j];

                tey_z_xxx_0[j] = pa_z[j] * tey_0_xxx_0[j] - pc_z[j] * tey_0_xxx_1[j];

                tez_z_xxx_0[j] = pa_z[j] * tez_0_xxx_0[j] - pc_z[j] * tez_0_xxx_1[j] + ta_0_xxx_1[j];

                tex_z_xxy_0[j] = pa_z[j] * tex_0_xxy_0[j] - pc_z[j] * tex_0_xxy_1[j];

                tey_z_xxy_0[j] = pa_z[j] * tey_0_xxy_0[j] - pc_z[j] * tey_0_xxy_1[j];

                tez_z_xxy_0[j] = pa_z[j] * tez_0_xxy_0[j] - pc_z[j] * tez_0_xxy_1[j] + ta_0_xxy_1[j];

                tex_z_xxz_0[j] = pa_z[j] * tex_0_xxz_0[j] - pc_z[j] * tex_0_xxz_1[j] + 0.5 * fl1_fx * tex_0_xx_0[j] - 0.5 * fl1_fx * tex_0_xx_1[j];

                tey_z_xxz_0[j] = pa_z[j] * tey_0_xxz_0[j] - pc_z[j] * tey_0_xxz_1[j] + 0.5 * fl1_fx * tey_0_xx_0[j] - 0.5 * fl1_fx * tey_0_xx_1[j];

                tez_z_xxz_0[j] =
                    pa_z[j] * tez_0_xxz_0[j] - pc_z[j] * tez_0_xxz_1[j] + 0.5 * fl1_fx * tez_0_xx_0[j] - 0.5 * fl1_fx * tez_0_xx_1[j] + ta_0_xxz_1[j];

                tex_z_xyy_0[j] = pa_z[j] * tex_0_xyy_0[j] - pc_z[j] * tex_0_xyy_1[j];

                tey_z_xyy_0[j] = pa_z[j] * tey_0_xyy_0[j] - pc_z[j] * tey_0_xyy_1[j];

                tez_z_xyy_0[j] = pa_z[j] * tez_0_xyy_0[j] - pc_z[j] * tez_0_xyy_1[j] + ta_0_xyy_1[j];

                tex_z_xyz_0[j] = pa_z[j] * tex_0_xyz_0[j] - pc_z[j] * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_0_xy_0[j] - 0.5 * fl1_fx * tex_0_xy_1[j];

                tey_z_xyz_0[j] = pa_z[j] * tey_0_xyz_0[j] - pc_z[j] * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_0_xy_0[j] - 0.5 * fl1_fx * tey_0_xy_1[j];

                tez_z_xyz_0[j] =
                    pa_z[j] * tez_0_xyz_0[j] - pc_z[j] * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_0_xy_0[j] - 0.5 * fl1_fx * tez_0_xy_1[j] + ta_0_xyz_1[j];

                tex_z_xzz_0[j] = pa_z[j] * tex_0_xzz_0[j] - pc_z[j] * tex_0_xzz_1[j] + fl1_fx * tex_0_xz_0[j] - fl1_fx * tex_0_xz_1[j];

                tey_z_xzz_0[j] = pa_z[j] * tey_0_xzz_0[j] - pc_z[j] * tey_0_xzz_1[j] + fl1_fx * tey_0_xz_0[j] - fl1_fx * tey_0_xz_1[j];

                tez_z_xzz_0[j] =
                    pa_z[j] * tez_0_xzz_0[j] - pc_z[j] * tez_0_xzz_1[j] + fl1_fx * tez_0_xz_0[j] - fl1_fx * tez_0_xz_1[j] + ta_0_xzz_1[j];

                tex_z_yyy_0[j] = pa_z[j] * tex_0_yyy_0[j] - pc_z[j] * tex_0_yyy_1[j];

                tey_z_yyy_0[j] = pa_z[j] * tey_0_yyy_0[j] - pc_z[j] * tey_0_yyy_1[j];

                tez_z_yyy_0[j] = pa_z[j] * tez_0_yyy_0[j] - pc_z[j] * tez_0_yyy_1[j] + ta_0_yyy_1[j];

                tex_z_yyz_0[j] = pa_z[j] * tex_0_yyz_0[j] - pc_z[j] * tex_0_yyz_1[j] + 0.5 * fl1_fx * tex_0_yy_0[j] - 0.5 * fl1_fx * tex_0_yy_1[j];

                tey_z_yyz_0[j] = pa_z[j] * tey_0_yyz_0[j] - pc_z[j] * tey_0_yyz_1[j] + 0.5 * fl1_fx * tey_0_yy_0[j] - 0.5 * fl1_fx * tey_0_yy_1[j];

                tez_z_yyz_0[j] =
                    pa_z[j] * tez_0_yyz_0[j] - pc_z[j] * tez_0_yyz_1[j] + 0.5 * fl1_fx * tez_0_yy_0[j] - 0.5 * fl1_fx * tez_0_yy_1[j] + ta_0_yyz_1[j];

                tex_z_yzz_0[j] = pa_z[j] * tex_0_yzz_0[j] - pc_z[j] * tex_0_yzz_1[j] + fl1_fx * tex_0_yz_0[j] - fl1_fx * tex_0_yz_1[j];

                tey_z_yzz_0[j] = pa_z[j] * tey_0_yzz_0[j] - pc_z[j] * tey_0_yzz_1[j] + fl1_fx * tey_0_yz_0[j] - fl1_fx * tey_0_yz_1[j];

                tez_z_yzz_0[j] =
                    pa_z[j] * tez_0_yzz_0[j] - pc_z[j] * tez_0_yzz_1[j] + fl1_fx * tez_0_yz_0[j] - fl1_fx * tez_0_yz_1[j] + ta_0_yzz_1[j];

                tex_z_zzz_0[j] = pa_z[j] * tex_0_zzz_0[j] - pc_z[j] * tex_0_zzz_1[j] + 1.5 * fl1_fx * tex_0_zz_0[j] - 1.5 * fl1_fx * tex_0_zz_1[j];

                tey_z_zzz_0[j] = pa_z[j] * tey_0_zzz_0[j] - pc_z[j] * tey_0_zzz_1[j] + 1.5 * fl1_fx * tey_0_zz_0[j] - 1.5 * fl1_fx * tey_0_zz_1[j];

                tez_z_zzz_0[j] =
                    pa_z[j] * tez_0_zzz_0[j] - pc_z[j] * tez_0_zzz_1[j] + 1.5 * fl1_fx * tez_0_zz_0[j] - 1.5 * fl1_fx * tez_0_zz_1[j] + ta_0_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForFP(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForFP_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForFP_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForFP_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_3_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx);

            auto tey_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx);

            auto tez_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx);

            auto tex_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 1);

            auto tey_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 1);

            auto tez_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 1);

            auto tex_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 2);

            auto tey_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 2);

            auto tez_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 2);

            auto tex_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 3);

            auto tey_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 3);

            auto tez_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 3);

            auto tex_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 4);

            auto tey_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 4);

            auto tez_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 4);

            auto tex_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 5);

            auto tey_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 5);

            auto tez_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 5);

            auto tex_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 6);

            auto tey_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 6);

            auto tez_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 6);

            auto tex_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 7);

            auto tey_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 7);

            auto tez_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 7);

            auto tex_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 8);

            auto tey_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 8);

            auto tez_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 8);

            auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9);

            auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10);

            auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11);

            auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12);

            auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13);

            auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14);

            auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14);

            auto tex_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx);

            auto tey_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx);

            auto tez_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx);

            auto tex_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 1);

            auto tey_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 1);

            auto tez_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 1);

            auto tex_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 2);

            auto tey_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 2);

            auto tez_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 2);

            auto tex_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 3);

            auto tey_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 3);

            auto tez_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 3);

            auto tex_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 4);

            auto tey_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 4);

            auto tez_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 4);

            auto tex_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 5);

            auto tey_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 5);

            auto tez_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 5);

            auto tex_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 6);

            auto tey_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 6);

            auto tez_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 6);

            auto tex_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 7);

            auto tey_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 7);

            auto tez_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 7);

            auto tex_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 8);

            auto tey_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 8);

            auto tez_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 8);

            auto tex_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 9);

            auto tey_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 10);

            auto tey_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 11);

            auto tey_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 12);

            auto tey_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 13);

            auto tey_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 14);

            auto tey_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 14);

            auto tex_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx);

            auto tey_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx);

            auto tez_x_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx);

            auto tex_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 1);

            auto tey_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 1);

            auto tez_x_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 1);

            auto tex_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 2);

            auto tey_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 2);

            auto tez_x_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 2);

            auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3);

            auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4);

            auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5);

            auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6);

            auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7);

            auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8);

            auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8);

            auto tex_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx);

            auto tey_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx);

            auto tez_x_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx);

            auto tex_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 1);

            auto tey_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 1);

            auto tez_x_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 1);

            auto tex_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 2);

            auto tey_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 2);

            auto tez_x_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 2);

            auto tex_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 3);

            auto tey_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 4);

            auto tey_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 5);

            auto tey_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 6);

            auto tey_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 7);

            auto tey_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 8);

            auto tey_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 8);

            auto tex_xx_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx);

            auto tey_xx_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx);

            auto tez_xx_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx);

            auto tex_xy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 1);

            auto tey_xy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 1);

            auto tez_xy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 1);

            auto tex_xz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 2);

            auto tey_xz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 2);

            auto tez_xz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 2);

            auto tex_yy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 3);

            auto tey_yy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 3);

            auto tez_yy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 3);

            auto tex_yz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 4);

            auto tey_yz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 4);

            auto tez_yz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 4);

            auto tex_xx_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx);

            auto tey_xx_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx);

            auto tez_xx_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx);

            auto tex_xy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 1);

            auto tey_xy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 1);

            auto tez_xy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 1);

            auto tex_xz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 2);

            auto tey_xz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 2);

            auto tez_xz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 2);

            auto tex_yy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 3);

            auto tey_yy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 3);

            auto tez_yy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 3);

            auto tex_yz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 4);

            auto tey_yz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 4);

            auto tez_yz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 4);

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

            auto tex_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx);

            auto tey_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx);

            auto tez_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx);

            auto tex_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 1);

            auto tey_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 1);

            auto tez_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 1);

            auto tex_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 2);

            auto tey_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 2);

            auto tez_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 2);

            auto tex_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 3);

            auto tey_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 3);

            auto tez_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 3);

            auto tex_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 4);

            auto tey_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 4);

            auto tez_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 4);

            auto tex_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 5);

            auto tey_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 5);

            auto tez_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 5);

            auto tex_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 6);

            auto tey_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 6);

            auto tez_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 6);

            auto tex_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 7);

            auto tey_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 7);

            auto tez_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 7);

            auto tex_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 8);

            auto tey_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 8);

            auto tez_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 8);

            auto tex_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 9);

            auto tey_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 9);

            auto tez_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 9);

            auto tex_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 10);

            auto tey_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 10);

            auto tez_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 10);

            auto tex_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 11);

            auto tey_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 11);

            auto tez_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 11);

            auto tex_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 12);

            auto tey_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 12);

            auto tez_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 12);

            auto tex_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 13);

            auto tey_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 13);

            auto tez_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 13);

            auto tex_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 14);

            auto tey_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 14);

            auto tez_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 14);

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xx_x_1, ta_xx_y_1, ta_xx_z_1, ta_xy_x_1, ta_xy_y_1, \
                                         ta_xy_z_1, ta_xz_x_1, ta_xz_y_1, ta_xz_z_1, ta_yy_x_1, ta_yy_y_1, ta_yy_z_1, \
                                         ta_yz_x_1, ta_yz_y_1, ta_yz_z_1, tex_x_x_0, tex_x_x_1, tex_x_y_0, tex_x_y_1, \
                                         tex_x_z_0, tex_x_z_1, tex_xx_0_0, tex_xx_0_1, tex_xx_x_0, tex_xx_x_1, tex_xx_y_0, \
                                         tex_xx_y_1, tex_xx_z_0, tex_xx_z_1, tex_xxx_x_0, tex_xxx_y_0, tex_xxx_z_0, \
                                         tex_xxy_x_0, tex_xxy_y_0, tex_xxy_z_0, tex_xxz_x_0, tex_xxz_y_0, tex_xxz_z_0, \
                                         tex_xy_0_0, tex_xy_0_1, tex_xy_x_0, tex_xy_x_1, tex_xy_y_0, tex_xy_y_1, tex_xy_z_0, \
                                         tex_xy_z_1, tex_xyy_x_0, tex_xyy_y_0, tex_xyy_z_0, tex_xyz_x_0, tex_xyz_y_0, \
                                         tex_xyz_z_0, tex_xz_0_0, tex_xz_0_1, tex_xz_x_0, tex_xz_x_1, tex_xz_y_0, tex_xz_y_1, \
                                         tex_xz_z_0, tex_xz_z_1, tex_y_x_0, tex_y_x_1, tex_y_y_0, tex_y_y_1, tex_y_z_0, \
                                         tex_y_z_1, tex_yy_0_0, tex_yy_0_1, tex_yy_x_0, tex_yy_x_1, tex_yy_y_0, tex_yy_y_1, \
                                         tex_yy_z_0, tex_yy_z_1, tex_yz_0_0, tex_yz_0_1, tex_yz_x_0, tex_yz_x_1, tex_yz_y_0, \
                                         tex_yz_y_1, tex_yz_z_0, tex_yz_z_1, tex_z_x_0, tex_z_x_1, tex_z_y_0, tex_z_y_1, \
                                         tex_z_z_0, tex_z_z_1, tey_x_x_0, tey_x_x_1, tey_x_y_0, tey_x_y_1, tey_x_z_0, \
                                         tey_x_z_1, tey_xx_0_0, tey_xx_0_1, tey_xx_x_0, tey_xx_x_1, tey_xx_y_0, tey_xx_y_1, \
                                         tey_xx_z_0, tey_xx_z_1, tey_xxx_x_0, tey_xxx_y_0, tey_xxx_z_0, tey_xxy_x_0, \
                                         tey_xxy_y_0, tey_xxy_z_0, tey_xxz_x_0, tey_xxz_y_0, tey_xxz_z_0, tey_xy_0_0, \
                                         tey_xy_0_1, tey_xy_x_0, tey_xy_x_1, tey_xy_y_0, tey_xy_y_1, tey_xy_z_0, tey_xy_z_1, \
                                         tey_xyy_x_0, tey_xyy_y_0, tey_xyy_z_0, tey_xyz_x_0, tey_xyz_y_0, tey_xyz_z_0, \
                                         tey_xz_0_0, tey_xz_0_1, tey_xz_x_0, tey_xz_x_1, tey_xz_y_0, tey_xz_y_1, tey_xz_z_0, \
                                         tey_xz_z_1, tey_y_x_0, tey_y_x_1, tey_y_y_0, tey_y_y_1, tey_y_z_0, tey_y_z_1, \
                                         tey_yy_0_0, tey_yy_0_1, tey_yy_x_0, tey_yy_x_1, tey_yy_y_0, tey_yy_y_1, tey_yy_z_0, \
                                         tey_yy_z_1, tey_yz_0_0, tey_yz_0_1, tey_yz_x_0, tey_yz_x_1, tey_yz_y_0, tey_yz_y_1, \
                                         tey_yz_z_0, tey_yz_z_1, tey_z_x_0, tey_z_x_1, tey_z_y_0, tey_z_y_1, tey_z_z_0, \
                                         tey_z_z_1, tez_x_x_0, tez_x_x_1, tez_x_y_0, tez_x_y_1, tez_x_z_0, tez_x_z_1, \
                                         tez_xx_0_0, tez_xx_0_1, tez_xx_x_0, tez_xx_x_1, tez_xx_y_0, tez_xx_y_1, tez_xx_z_0, \
                                         tez_xx_z_1, tez_xxx_x_0, tez_xxx_y_0, tez_xxx_z_0, tez_xxy_x_0, tez_xxy_y_0, \
                                         tez_xxy_z_0, tez_xxz_x_0, tez_xxz_y_0, tez_xxz_z_0, tez_xy_0_0, tez_xy_0_1, \
                                         tez_xy_x_0, tez_xy_x_1, tez_xy_y_0, tez_xy_y_1, tez_xy_z_0, tez_xy_z_1, \
                                         tez_xyy_x_0, tez_xyy_y_0, tez_xyy_z_0, tez_xyz_x_0, tez_xyz_y_0, tez_xyz_z_0, \
                                         tez_xz_0_0, tez_xz_0_1, tez_xz_x_0, tez_xz_x_1, tez_xz_y_0, tez_xz_y_1, tez_xz_z_0, \
                                         tez_xz_z_1, tez_y_x_0, tez_y_x_1, tez_y_y_0, tez_y_y_1, tez_y_z_0, tez_y_z_1, \
                                         tez_yy_0_0, tez_yy_0_1, tez_yy_x_0, tez_yy_x_1, tez_yy_y_0, tez_yy_y_1, tez_yy_z_0, \
                                         tez_yy_z_1, tez_yz_0_0, tez_yz_0_1, tez_yz_x_0, tez_yz_x_1, tez_yz_y_0, tez_yz_y_1, \
                                         tez_yz_z_0, tez_yz_z_1, tez_z_x_0, tez_z_x_1, tez_z_y_0, tez_z_y_1, tez_z_z_0, \
                                         tez_z_z_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxx_x_0[j] = pa_x[j] * tex_xx_x_0[j] - pc_x[j] * tex_xx_x_1[j] + fl1_fx * tex_x_x_0[j] - fl1_fx * tex_x_x_1[j] +
                                 0.5 * fl1_fx * tex_xx_0_0[j] - 0.5 * fl1_fx * tex_xx_0_1[j] + ta_xx_x_1[j];

                tey_xxx_x_0[j] = pa_x[j] * tey_xx_x_0[j] - pc_x[j] * tey_xx_x_1[j] + fl1_fx * tey_x_x_0[j] - fl1_fx * tey_x_x_1[j] +
                                 0.5 * fl1_fx * tey_xx_0_0[j] - 0.5 * fl1_fx * tey_xx_0_1[j];

                tez_xxx_x_0[j] = pa_x[j] * tez_xx_x_0[j] - pc_x[j] * tez_xx_x_1[j] + fl1_fx * tez_x_x_0[j] - fl1_fx * tez_x_x_1[j] +
                                 0.5 * fl1_fx * tez_xx_0_0[j] - 0.5 * fl1_fx * tez_xx_0_1[j];

                tex_xxx_y_0[j] = pa_x[j] * tex_xx_y_0[j] - pc_x[j] * tex_xx_y_1[j] + fl1_fx * tex_x_y_0[j] - fl1_fx * tex_x_y_1[j] + ta_xx_y_1[j];

                tey_xxx_y_0[j] = pa_x[j] * tey_xx_y_0[j] - pc_x[j] * tey_xx_y_1[j] + fl1_fx * tey_x_y_0[j] - fl1_fx * tey_x_y_1[j];

                tez_xxx_y_0[j] = pa_x[j] * tez_xx_y_0[j] - pc_x[j] * tez_xx_y_1[j] + fl1_fx * tez_x_y_0[j] - fl1_fx * tez_x_y_1[j];

                tex_xxx_z_0[j] = pa_x[j] * tex_xx_z_0[j] - pc_x[j] * tex_xx_z_1[j] + fl1_fx * tex_x_z_0[j] - fl1_fx * tex_x_z_1[j] + ta_xx_z_1[j];

                tey_xxx_z_0[j] = pa_x[j] * tey_xx_z_0[j] - pc_x[j] * tey_xx_z_1[j] + fl1_fx * tey_x_z_0[j] - fl1_fx * tey_x_z_1[j];

                tez_xxx_z_0[j] = pa_x[j] * tez_xx_z_0[j] - pc_x[j] * tez_xx_z_1[j] + fl1_fx * tez_x_z_0[j] - fl1_fx * tez_x_z_1[j];

                tex_xxy_x_0[j] = pa_x[j] * tex_xy_x_0[j] - pc_x[j] * tex_xy_x_1[j] + 0.5 * fl1_fx * tex_y_x_0[j] - 0.5 * fl1_fx * tex_y_x_1[j] +
                                 0.5 * fl1_fx * tex_xy_0_0[j] - 0.5 * fl1_fx * tex_xy_0_1[j] + ta_xy_x_1[j];

                tey_xxy_x_0[j] = pa_x[j] * tey_xy_x_0[j] - pc_x[j] * tey_xy_x_1[j] + 0.5 * fl1_fx * tey_y_x_0[j] - 0.5 * fl1_fx * tey_y_x_1[j] +
                                 0.5 * fl1_fx * tey_xy_0_0[j] - 0.5 * fl1_fx * tey_xy_0_1[j];

                tez_xxy_x_0[j] = pa_x[j] * tez_xy_x_0[j] - pc_x[j] * tez_xy_x_1[j] + 0.5 * fl1_fx * tez_y_x_0[j] - 0.5 * fl1_fx * tez_y_x_1[j] +
                                 0.5 * fl1_fx * tez_xy_0_0[j] - 0.5 * fl1_fx * tez_xy_0_1[j];

                tex_xxy_y_0[j] =
                    pa_x[j] * tex_xy_y_0[j] - pc_x[j] * tex_xy_y_1[j] + 0.5 * fl1_fx * tex_y_y_0[j] - 0.5 * fl1_fx * tex_y_y_1[j] + ta_xy_y_1[j];

                tey_xxy_y_0[j] = pa_x[j] * tey_xy_y_0[j] - pc_x[j] * tey_xy_y_1[j] + 0.5 * fl1_fx * tey_y_y_0[j] - 0.5 * fl1_fx * tey_y_y_1[j];

                tez_xxy_y_0[j] = pa_x[j] * tez_xy_y_0[j] - pc_x[j] * tez_xy_y_1[j] + 0.5 * fl1_fx * tez_y_y_0[j] - 0.5 * fl1_fx * tez_y_y_1[j];

                tex_xxy_z_0[j] =
                    pa_x[j] * tex_xy_z_0[j] - pc_x[j] * tex_xy_z_1[j] + 0.5 * fl1_fx * tex_y_z_0[j] - 0.5 * fl1_fx * tex_y_z_1[j] + ta_xy_z_1[j];

                tey_xxy_z_0[j] = pa_x[j] * tey_xy_z_0[j] - pc_x[j] * tey_xy_z_1[j] + 0.5 * fl1_fx * tey_y_z_0[j] - 0.5 * fl1_fx * tey_y_z_1[j];

                tez_xxy_z_0[j] = pa_x[j] * tez_xy_z_0[j] - pc_x[j] * tez_xy_z_1[j] + 0.5 * fl1_fx * tez_y_z_0[j] - 0.5 * fl1_fx * tez_y_z_1[j];

                tex_xxz_x_0[j] = pa_x[j] * tex_xz_x_0[j] - pc_x[j] * tex_xz_x_1[j] + 0.5 * fl1_fx * tex_z_x_0[j] - 0.5 * fl1_fx * tex_z_x_1[j] +
                                 0.5 * fl1_fx * tex_xz_0_0[j] - 0.5 * fl1_fx * tex_xz_0_1[j] + ta_xz_x_1[j];

                tey_xxz_x_0[j] = pa_x[j] * tey_xz_x_0[j] - pc_x[j] * tey_xz_x_1[j] + 0.5 * fl1_fx * tey_z_x_0[j] - 0.5 * fl1_fx * tey_z_x_1[j] +
                                 0.5 * fl1_fx * tey_xz_0_0[j] - 0.5 * fl1_fx * tey_xz_0_1[j];

                tez_xxz_x_0[j] = pa_x[j] * tez_xz_x_0[j] - pc_x[j] * tez_xz_x_1[j] + 0.5 * fl1_fx * tez_z_x_0[j] - 0.5 * fl1_fx * tez_z_x_1[j] +
                                 0.5 * fl1_fx * tez_xz_0_0[j] - 0.5 * fl1_fx * tez_xz_0_1[j];

                tex_xxz_y_0[j] =
                    pa_x[j] * tex_xz_y_0[j] - pc_x[j] * tex_xz_y_1[j] + 0.5 * fl1_fx * tex_z_y_0[j] - 0.5 * fl1_fx * tex_z_y_1[j] + ta_xz_y_1[j];

                tey_xxz_y_0[j] = pa_x[j] * tey_xz_y_0[j] - pc_x[j] * tey_xz_y_1[j] + 0.5 * fl1_fx * tey_z_y_0[j] - 0.5 * fl1_fx * tey_z_y_1[j];

                tez_xxz_y_0[j] = pa_x[j] * tez_xz_y_0[j] - pc_x[j] * tez_xz_y_1[j] + 0.5 * fl1_fx * tez_z_y_0[j] - 0.5 * fl1_fx * tez_z_y_1[j];

                tex_xxz_z_0[j] =
                    pa_x[j] * tex_xz_z_0[j] - pc_x[j] * tex_xz_z_1[j] + 0.5 * fl1_fx * tex_z_z_0[j] - 0.5 * fl1_fx * tex_z_z_1[j] + ta_xz_z_1[j];

                tey_xxz_z_0[j] = pa_x[j] * tey_xz_z_0[j] - pc_x[j] * tey_xz_z_1[j] + 0.5 * fl1_fx * tey_z_z_0[j] - 0.5 * fl1_fx * tey_z_z_1[j];

                tez_xxz_z_0[j] = pa_x[j] * tez_xz_z_0[j] - pc_x[j] * tez_xz_z_1[j] + 0.5 * fl1_fx * tez_z_z_0[j] - 0.5 * fl1_fx * tez_z_z_1[j];

                tex_xyy_x_0[j] =
                    pa_x[j] * tex_yy_x_0[j] - pc_x[j] * tex_yy_x_1[j] + 0.5 * fl1_fx * tex_yy_0_0[j] - 0.5 * fl1_fx * tex_yy_0_1[j] + ta_yy_x_1[j];

                tey_xyy_x_0[j] = pa_x[j] * tey_yy_x_0[j] - pc_x[j] * tey_yy_x_1[j] + 0.5 * fl1_fx * tey_yy_0_0[j] - 0.5 * fl1_fx * tey_yy_0_1[j];

                tez_xyy_x_0[j] = pa_x[j] * tez_yy_x_0[j] - pc_x[j] * tez_yy_x_1[j] + 0.5 * fl1_fx * tez_yy_0_0[j] - 0.5 * fl1_fx * tez_yy_0_1[j];

                tex_xyy_y_0[j] = pa_x[j] * tex_yy_y_0[j] - pc_x[j] * tex_yy_y_1[j] + ta_yy_y_1[j];

                tey_xyy_y_0[j] = pa_x[j] * tey_yy_y_0[j] - pc_x[j] * tey_yy_y_1[j];

                tez_xyy_y_0[j] = pa_x[j] * tez_yy_y_0[j] - pc_x[j] * tez_yy_y_1[j];

                tex_xyy_z_0[j] = pa_x[j] * tex_yy_z_0[j] - pc_x[j] * tex_yy_z_1[j] + ta_yy_z_1[j];

                tey_xyy_z_0[j] = pa_x[j] * tey_yy_z_0[j] - pc_x[j] * tey_yy_z_1[j];

                tez_xyy_z_0[j] = pa_x[j] * tez_yy_z_0[j] - pc_x[j] * tez_yy_z_1[j];

                tex_xyz_x_0[j] =
                    pa_x[j] * tex_yz_x_0[j] - pc_x[j] * tex_yz_x_1[j] + 0.5 * fl1_fx * tex_yz_0_0[j] - 0.5 * fl1_fx * tex_yz_0_1[j] + ta_yz_x_1[j];

                tey_xyz_x_0[j] = pa_x[j] * tey_yz_x_0[j] - pc_x[j] * tey_yz_x_1[j] + 0.5 * fl1_fx * tey_yz_0_0[j] - 0.5 * fl1_fx * tey_yz_0_1[j];

                tez_xyz_x_0[j] = pa_x[j] * tez_yz_x_0[j] - pc_x[j] * tez_yz_x_1[j] + 0.5 * fl1_fx * tez_yz_0_0[j] - 0.5 * fl1_fx * tez_yz_0_1[j];

                tex_xyz_y_0[j] = pa_x[j] * tex_yz_y_0[j] - pc_x[j] * tex_yz_y_1[j] + ta_yz_y_1[j];

                tey_xyz_y_0[j] = pa_x[j] * tey_yz_y_0[j] - pc_x[j] * tey_yz_y_1[j];

                tez_xyz_y_0[j] = pa_x[j] * tez_yz_y_0[j] - pc_x[j] * tez_yz_y_1[j];

                tex_xyz_z_0[j] = pa_x[j] * tex_yz_z_0[j] - pc_x[j] * tex_yz_z_1[j] + ta_yz_z_1[j];

                tey_xyz_z_0[j] = pa_x[j] * tey_yz_z_0[j] - pc_x[j] * tey_yz_z_1[j];

                tez_xyz_z_0[j] = pa_x[j] * tez_yz_z_0[j] - pc_x[j] * tez_yz_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForFP_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_3_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9);

            auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10);

            auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11);

            auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12);

            auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13);

            auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14);

            auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14);

            auto tex_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 15);

            auto tey_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 16);

            auto tey_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 17);

            auto tey_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 17);

            auto tex_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 9);

            auto tey_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 10);

            auto tey_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 11);

            auto tey_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 12);

            auto tey_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 13);

            auto tey_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 14);

            auto tey_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 14);

            auto tex_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 15);

            auto tey_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 16);

            auto tey_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 17);

            auto tey_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 17);

            auto tex_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 3);

            auto tey_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 4);

            auto tey_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 5);

            auto tey_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 6);

            auto tey_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 7);

            auto tey_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * idx + 8);

            auto tey_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_0 = primBuffer.data(pidx_e_1_1_m0 + 18 * bdim + 9 * idx + 8);

            auto tex_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 3);

            auto tey_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 3);

            auto tez_y_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 3);

            auto tex_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 4);

            auto tey_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 4);

            auto tez_y_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 4);

            auto tex_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 5);

            auto tey_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 5);

            auto tez_y_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 5);

            auto tex_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 6);

            auto tey_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 6);

            auto tez_z_x_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 6);

            auto tex_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 7);

            auto tey_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 7);

            auto tez_z_y_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 7);

            auto tex_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * idx + 8);

            auto tey_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 9 * bdim + 9 * idx + 8);

            auto tez_z_z_1 = primBuffer.data(pidx_e_1_1_m1 + 18 * bdim + 9 * idx + 8);

            auto tex_yy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 3);

            auto tey_yy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 3);

            auto tez_yy_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 3);

            auto tex_yz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 4);

            auto tey_yz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 4);

            auto tez_yz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 4);

            auto tex_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 5);

            auto tey_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 5);

            auto tex_yy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 3);

            auto tey_yy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 3);

            auto tez_yy_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 3);

            auto tex_yz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 4);

            auto tey_yz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 4);

            auto tez_yz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 4);

            auto tex_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 5);

            auto tey_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 5);

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

            auto tex_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 15);

            auto tey_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 15);

            auto tez_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 15);

            auto tex_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 16);

            auto tey_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 16);

            auto tez_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 16);

            auto tex_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 17);

            auto tey_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 17);

            auto tez_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 17);

            auto tex_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 18);

            auto tey_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 18);

            auto tez_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 18);

            auto tex_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 19);

            auto tey_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 19);

            auto tez_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 19);

            auto tex_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 20);

            auto tey_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 20);

            auto tez_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 20);

            auto tex_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 21);

            auto tey_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 21);

            auto tez_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 21);

            auto tex_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 22);

            auto tey_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 22);

            auto tez_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 22);

            auto tex_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 23);

            auto tey_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 23);

            auto tez_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 23);

            auto tex_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 24);

            auto tey_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 24);

            auto tez_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 24);

            auto tex_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 25);

            auto tey_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 25);

            auto tez_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 25);

            auto tex_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 26);

            auto tey_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 26);

            auto tez_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 26);

            auto tex_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 27);

            auto tey_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 27);

            auto tez_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 27);

            auto tex_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 28);

            auto tey_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 28);

            auto tez_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 28);

            auto tex_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 29);

            auto tey_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 29);

            auto tez_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 29);

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_yy_x_1, ta_yy_y_1, ta_yy_z_1, \
                                         ta_yz_x_1, ta_yz_y_1, ta_yz_z_1, ta_zz_x_1, ta_zz_y_1, ta_zz_z_1, tex_xzz_x_0, \
                                         tex_xzz_y_0, tex_xzz_z_0, tex_y_x_0, tex_y_x_1, tex_y_y_0, tex_y_y_1, tex_y_z_0, \
                                         tex_y_z_1, tex_yy_0_0, tex_yy_0_1, tex_yy_x_0, tex_yy_x_1, tex_yy_y_0, tex_yy_y_1, \
                                         tex_yy_z_0, tex_yy_z_1, tex_yyy_x_0, tex_yyy_y_0, tex_yyy_z_0, tex_yyz_x_0, \
                                         tex_yyz_y_0, tex_yyz_z_0, tex_yz_0_0, tex_yz_0_1, tex_yz_x_0, tex_yz_x_1, \
                                         tex_yz_y_0, tex_yz_y_1, tex_yz_z_0, tex_yz_z_1, tex_yzz_x_0, tex_yzz_y_0, \
                                         tex_yzz_z_0, tex_z_x_0, tex_z_x_1, tex_z_y_0, tex_z_y_1, tex_z_z_0, tex_z_z_1, \
                                         tex_zz_0_0, tex_zz_0_1, tex_zz_x_0, tex_zz_x_1, tex_zz_y_0, tex_zz_y_1, tex_zz_z_0, \
                                         tex_zz_z_1, tex_zzz_x_0, tex_zzz_y_0, tex_zzz_z_0, tey_xzz_x_0, tey_xzz_y_0, \
                                         tey_xzz_z_0, tey_y_x_0, tey_y_x_1, tey_y_y_0, tey_y_y_1, tey_y_z_0, tey_y_z_1, \
                                         tey_yy_0_0, tey_yy_0_1, tey_yy_x_0, tey_yy_x_1, tey_yy_y_0, tey_yy_y_1, tey_yy_z_0, \
                                         tey_yy_z_1, tey_yyy_x_0, tey_yyy_y_0, tey_yyy_z_0, tey_yyz_x_0, tey_yyz_y_0, \
                                         tey_yyz_z_0, tey_yz_0_0, tey_yz_0_1, tey_yz_x_0, tey_yz_x_1, tey_yz_y_0, tey_yz_y_1, \
                                         tey_yz_z_0, tey_yz_z_1, tey_yzz_x_0, tey_yzz_y_0, tey_yzz_z_0, tey_z_x_0, \
                                         tey_z_x_1, tey_z_y_0, tey_z_y_1, tey_z_z_0, tey_z_z_1, tey_zz_0_0, tey_zz_0_1, \
                                         tey_zz_x_0, tey_zz_x_1, tey_zz_y_0, tey_zz_y_1, tey_zz_z_0, tey_zz_z_1, \
                                         tey_zzz_x_0, tey_zzz_y_0, tey_zzz_z_0, tez_xzz_x_0, tez_xzz_y_0, tez_xzz_z_0, \
                                         tez_y_x_0, tez_y_x_1, tez_y_y_0, tez_y_y_1, tez_y_z_0, tez_y_z_1, tez_yy_0_0, \
                                         tez_yy_0_1, tez_yy_x_0, tez_yy_x_1, tez_yy_y_0, tez_yy_y_1, tez_yy_z_0, tez_yy_z_1, \
                                         tez_yyy_x_0, tez_yyy_y_0, tez_yyy_z_0, tez_yyz_x_0, tez_yyz_y_0, tez_yyz_z_0, \
                                         tez_yz_0_0, tez_yz_0_1, tez_yz_x_0, tez_yz_x_1, tez_yz_y_0, tez_yz_y_1, tez_yz_z_0, \
                                         tez_yz_z_1, tez_yzz_x_0, tez_yzz_y_0, tez_yzz_z_0, tez_z_x_0, tez_z_x_1, tez_z_y_0, \
                                         tez_z_y_1, tez_z_z_0, tez_z_z_1, tez_zz_0_0, tez_zz_0_1, tez_zz_x_0, tez_zz_x_1, \
                                         tez_zz_y_0, tez_zz_y_1, tez_zz_z_0, tez_zz_z_1, tez_zzz_x_0, tez_zzz_y_0, \
                                         tez_zzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xzz_x_0[j] =
                    pa_x[j] * tex_zz_x_0[j] - pc_x[j] * tex_zz_x_1[j] + 0.5 * fl1_fx * tex_zz_0_0[j] - 0.5 * fl1_fx * tex_zz_0_1[j] + ta_zz_x_1[j];

                tey_xzz_x_0[j] = pa_x[j] * tey_zz_x_0[j] - pc_x[j] * tey_zz_x_1[j] + 0.5 * fl1_fx * tey_zz_0_0[j] - 0.5 * fl1_fx * tey_zz_0_1[j];

                tez_xzz_x_0[j] = pa_x[j] * tez_zz_x_0[j] - pc_x[j] * tez_zz_x_1[j] + 0.5 * fl1_fx * tez_zz_0_0[j] - 0.5 * fl1_fx * tez_zz_0_1[j];

                tex_xzz_y_0[j] = pa_x[j] * tex_zz_y_0[j] - pc_x[j] * tex_zz_y_1[j] + ta_zz_y_1[j];

                tey_xzz_y_0[j] = pa_x[j] * tey_zz_y_0[j] - pc_x[j] * tey_zz_y_1[j];

                tez_xzz_y_0[j] = pa_x[j] * tez_zz_y_0[j] - pc_x[j] * tez_zz_y_1[j];

                tex_xzz_z_0[j] = pa_x[j] * tex_zz_z_0[j] - pc_x[j] * tex_zz_z_1[j] + ta_zz_z_1[j];

                tey_xzz_z_0[j] = pa_x[j] * tey_zz_z_0[j] - pc_x[j] * tey_zz_z_1[j];

                tez_xzz_z_0[j] = pa_x[j] * tez_zz_z_0[j] - pc_x[j] * tez_zz_z_1[j];

                tex_yyy_x_0[j] = pa_y[j] * tex_yy_x_0[j] - pc_y[j] * tex_yy_x_1[j] + fl1_fx * tex_y_x_0[j] - fl1_fx * tex_y_x_1[j];

                tey_yyy_x_0[j] = pa_y[j] * tey_yy_x_0[j] - pc_y[j] * tey_yy_x_1[j] + fl1_fx * tey_y_x_0[j] - fl1_fx * tey_y_x_1[j] + ta_yy_x_1[j];

                tez_yyy_x_0[j] = pa_y[j] * tez_yy_x_0[j] - pc_y[j] * tez_yy_x_1[j] + fl1_fx * tez_y_x_0[j] - fl1_fx * tez_y_x_1[j];

                tex_yyy_y_0[j] = pa_y[j] * tex_yy_y_0[j] - pc_y[j] * tex_yy_y_1[j] + fl1_fx * tex_y_y_0[j] - fl1_fx * tex_y_y_1[j] +
                                 0.5 * fl1_fx * tex_yy_0_0[j] - 0.5 * fl1_fx * tex_yy_0_1[j];

                tey_yyy_y_0[j] = pa_y[j] * tey_yy_y_0[j] - pc_y[j] * tey_yy_y_1[j] + fl1_fx * tey_y_y_0[j] - fl1_fx * tey_y_y_1[j] +
                                 0.5 * fl1_fx * tey_yy_0_0[j] - 0.5 * fl1_fx * tey_yy_0_1[j] + ta_yy_y_1[j];

                tez_yyy_y_0[j] = pa_y[j] * tez_yy_y_0[j] - pc_y[j] * tez_yy_y_1[j] + fl1_fx * tez_y_y_0[j] - fl1_fx * tez_y_y_1[j] +
                                 0.5 * fl1_fx * tez_yy_0_0[j] - 0.5 * fl1_fx * tez_yy_0_1[j];

                tex_yyy_z_0[j] = pa_y[j] * tex_yy_z_0[j] - pc_y[j] * tex_yy_z_1[j] + fl1_fx * tex_y_z_0[j] - fl1_fx * tex_y_z_1[j];

                tey_yyy_z_0[j] = pa_y[j] * tey_yy_z_0[j] - pc_y[j] * tey_yy_z_1[j] + fl1_fx * tey_y_z_0[j] - fl1_fx * tey_y_z_1[j] + ta_yy_z_1[j];

                tez_yyy_z_0[j] = pa_y[j] * tez_yy_z_0[j] - pc_y[j] * tez_yy_z_1[j] + fl1_fx * tez_y_z_0[j] - fl1_fx * tez_y_z_1[j];

                tex_yyz_x_0[j] = pa_y[j] * tex_yz_x_0[j] - pc_y[j] * tex_yz_x_1[j] + 0.5 * fl1_fx * tex_z_x_0[j] - 0.5 * fl1_fx * tex_z_x_1[j];

                tey_yyz_x_0[j] =
                    pa_y[j] * tey_yz_x_0[j] - pc_y[j] * tey_yz_x_1[j] + 0.5 * fl1_fx * tey_z_x_0[j] - 0.5 * fl1_fx * tey_z_x_1[j] + ta_yz_x_1[j];

                tez_yyz_x_0[j] = pa_y[j] * tez_yz_x_0[j] - pc_y[j] * tez_yz_x_1[j] + 0.5 * fl1_fx * tez_z_x_0[j] - 0.5 * fl1_fx * tez_z_x_1[j];

                tex_yyz_y_0[j] = pa_y[j] * tex_yz_y_0[j] - pc_y[j] * tex_yz_y_1[j] + 0.5 * fl1_fx * tex_z_y_0[j] - 0.5 * fl1_fx * tex_z_y_1[j] +
                                 0.5 * fl1_fx * tex_yz_0_0[j] - 0.5 * fl1_fx * tex_yz_0_1[j];

                tey_yyz_y_0[j] = pa_y[j] * tey_yz_y_0[j] - pc_y[j] * tey_yz_y_1[j] + 0.5 * fl1_fx * tey_z_y_0[j] - 0.5 * fl1_fx * tey_z_y_1[j] +
                                 0.5 * fl1_fx * tey_yz_0_0[j] - 0.5 * fl1_fx * tey_yz_0_1[j] + ta_yz_y_1[j];

                tez_yyz_y_0[j] = pa_y[j] * tez_yz_y_0[j] - pc_y[j] * tez_yz_y_1[j] + 0.5 * fl1_fx * tez_z_y_0[j] - 0.5 * fl1_fx * tez_z_y_1[j] +
                                 0.5 * fl1_fx * tez_yz_0_0[j] - 0.5 * fl1_fx * tez_yz_0_1[j];

                tex_yyz_z_0[j] = pa_y[j] * tex_yz_z_0[j] - pc_y[j] * tex_yz_z_1[j] + 0.5 * fl1_fx * tex_z_z_0[j] - 0.5 * fl1_fx * tex_z_z_1[j];

                tey_yyz_z_0[j] =
                    pa_y[j] * tey_yz_z_0[j] - pc_y[j] * tey_yz_z_1[j] + 0.5 * fl1_fx * tey_z_z_0[j] - 0.5 * fl1_fx * tey_z_z_1[j] + ta_yz_z_1[j];

                tez_yyz_z_0[j] = pa_y[j] * tez_yz_z_0[j] - pc_y[j] * tez_yz_z_1[j] + 0.5 * fl1_fx * tez_z_z_0[j] - 0.5 * fl1_fx * tez_z_z_1[j];

                tex_yzz_x_0[j] = pa_y[j] * tex_zz_x_0[j] - pc_y[j] * tex_zz_x_1[j];

                tey_yzz_x_0[j] = pa_y[j] * tey_zz_x_0[j] - pc_y[j] * tey_zz_x_1[j] + ta_zz_x_1[j];

                tez_yzz_x_0[j] = pa_y[j] * tez_zz_x_0[j] - pc_y[j] * tez_zz_x_1[j];

                tex_yzz_y_0[j] = pa_y[j] * tex_zz_y_0[j] - pc_y[j] * tex_zz_y_1[j] + 0.5 * fl1_fx * tex_zz_0_0[j] - 0.5 * fl1_fx * tex_zz_0_1[j];

                tey_yzz_y_0[j] =
                    pa_y[j] * tey_zz_y_0[j] - pc_y[j] * tey_zz_y_1[j] + 0.5 * fl1_fx * tey_zz_0_0[j] - 0.5 * fl1_fx * tey_zz_0_1[j] + ta_zz_y_1[j];

                tez_yzz_y_0[j] = pa_y[j] * tez_zz_y_0[j] - pc_y[j] * tez_zz_y_1[j] + 0.5 * fl1_fx * tez_zz_0_0[j] - 0.5 * fl1_fx * tez_zz_0_1[j];

                tex_yzz_z_0[j] = pa_y[j] * tex_zz_z_0[j] - pc_y[j] * tex_zz_z_1[j];

                tey_yzz_z_0[j] = pa_y[j] * tey_zz_z_0[j] - pc_y[j] * tey_zz_z_1[j] + ta_zz_z_1[j];

                tez_yzz_z_0[j] = pa_y[j] * tez_zz_z_0[j] - pc_y[j] * tez_zz_z_1[j];

                tex_zzz_x_0[j] = pa_z[j] * tex_zz_x_0[j] - pc_z[j] * tex_zz_x_1[j] + fl1_fx * tex_z_x_0[j] - fl1_fx * tex_z_x_1[j];

                tey_zzz_x_0[j] = pa_z[j] * tey_zz_x_0[j] - pc_z[j] * tey_zz_x_1[j] + fl1_fx * tey_z_x_0[j] - fl1_fx * tey_z_x_1[j];

                tez_zzz_x_0[j] = pa_z[j] * tez_zz_x_0[j] - pc_z[j] * tez_zz_x_1[j] + fl1_fx * tez_z_x_0[j] - fl1_fx * tez_z_x_1[j] + ta_zz_x_1[j];

                tex_zzz_y_0[j] = pa_z[j] * tex_zz_y_0[j] - pc_z[j] * tex_zz_y_1[j] + fl1_fx * tex_z_y_0[j] - fl1_fx * tex_z_y_1[j];

                tey_zzz_y_0[j] = pa_z[j] * tey_zz_y_0[j] - pc_z[j] * tey_zz_y_1[j] + fl1_fx * tey_z_y_0[j] - fl1_fx * tey_z_y_1[j];

                tez_zzz_y_0[j] = pa_z[j] * tez_zz_y_0[j] - pc_z[j] * tez_zz_y_1[j] + fl1_fx * tez_z_y_0[j] - fl1_fx * tez_z_y_1[j] + ta_zz_y_1[j];

                tex_zzz_z_0[j] = pa_z[j] * tex_zz_z_0[j] - pc_z[j] * tex_zz_z_1[j] + fl1_fx * tex_z_z_0[j] - fl1_fx * tex_z_z_1[j] +
                                 0.5 * fl1_fx * tex_zz_0_0[j] - 0.5 * fl1_fx * tex_zz_0_1[j];

                tey_zzz_z_0[j] = pa_z[j] * tey_zz_z_0[j] - pc_z[j] * tey_zz_z_1[j] + fl1_fx * tey_z_z_0[j] - fl1_fx * tey_z_z_1[j] +
                                 0.5 * fl1_fx * tey_zz_0_0[j] - 0.5 * fl1_fx * tey_zz_0_1[j];

                tez_zzz_z_0[j] = pa_z[j] * tez_zz_z_0[j] - pc_z[j] * tez_zz_z_1[j] + fl1_fx * tez_z_z_0[j] - fl1_fx * tez_z_z_1[j] +
                                 0.5 * fl1_fx * tez_zz_0_0[j] - 0.5 * fl1_fx * tez_zz_0_1[j] + ta_zz_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPG(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForPG_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPG_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPG_90_135(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForPG_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx);

            auto tey_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx);

            auto tez_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx);

            auto tex_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 1);

            auto tey_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 1);

            auto tez_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 1);

            auto tex_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 2);

            auto tey_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 2);

            auto tez_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 2);

            auto tex_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 3);

            auto tey_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 3);

            auto tez_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 3);

            auto tex_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 4);

            auto tey_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 4);

            auto tez_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 4);

            auto tex_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 5);

            auto tey_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 5);

            auto tez_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 5);

            auto tex_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 6);

            auto tey_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 6);

            auto tez_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 6);

            auto tex_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 7);

            auto tey_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 7);

            auto tez_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 7);

            auto tex_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 8);

            auto tey_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 8);

            auto tez_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 8);

            auto tex_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 9);

            auto tey_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 9);

            auto tez_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 9);

            auto tex_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 10);

            auto tey_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 10);

            auto tez_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 10);

            auto tex_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 11);

            auto tey_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 11);

            auto tez_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 11);

            auto tex_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 12);

            auto tey_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 12);

            auto tez_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 12);

            auto tex_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 13);

            auto tey_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 13);

            auto tez_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 13);

            auto tex_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 14);

            auto tey_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 14);

            auto tez_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 14);

            auto tex_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx);

            auto tey_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx);

            auto tez_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx);

            auto tex_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 1);

            auto tey_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 1);

            auto tez_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 1);

            auto tex_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 2);

            auto tey_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 2);

            auto tez_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 2);

            auto tex_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 3);

            auto tey_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 3);

            auto tez_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 3);

            auto tex_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 4);

            auto tey_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 4);

            auto tez_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 4);

            auto tex_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 5);

            auto tey_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 5);

            auto tez_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 5);

            auto tex_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 6);

            auto tey_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 6);

            auto tez_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 6);

            auto tex_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 7);

            auto tey_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 7);

            auto tez_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 7);

            auto tex_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 8);

            auto tey_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 8);

            auto tez_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 8);

            auto tex_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 9);

            auto tey_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 9);

            auto tez_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 9);

            auto tex_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 10);

            auto tey_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 10);

            auto tez_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 10);

            auto tex_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 11);

            auto tey_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 11);

            auto tez_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 11);

            auto tex_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 12);

            auto tey_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 12);

            auto tez_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 12);

            auto tex_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 13);

            auto tey_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 13);

            auto tez_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 13);

            auto tex_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 14);

            auto tey_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 14);

            auto tez_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 14);

            auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx);

            auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1);

            auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2);

            auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3);

            auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4);

            auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5);

            auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6);

            auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7);

            auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8);

            auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9);

            auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx);

            auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1);

            auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2);

            auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3);

            auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4);

            auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5);

            auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6);

            auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7);

            auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8);

            auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9);

            auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9);

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

            // set up pointers to integrals

            auto tex_x_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx);

            auto tey_x_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx);

            auto tez_x_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx);

            auto tex_x_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 1);

            auto tey_x_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 1);

            auto tez_x_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 1);

            auto tex_x_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 2);

            auto tey_x_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 2);

            auto tez_x_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 2);

            auto tex_x_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 3);

            auto tey_x_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 3);

            auto tez_x_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 3);

            auto tex_x_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 4);

            auto tey_x_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 4);

            auto tez_x_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 4);

            auto tex_x_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 5);

            auto tey_x_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 5);

            auto tez_x_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 5);

            auto tex_x_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 6);

            auto tey_x_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 6);

            auto tez_x_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 6);

            auto tex_x_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 7);

            auto tey_x_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 7);

            auto tez_x_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 7);

            auto tex_x_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 8);

            auto tey_x_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 8);

            auto tez_x_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 8);

            auto tex_x_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 9);

            auto tey_x_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 9);

            auto tez_x_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 9);

            auto tex_x_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 10);

            auto tey_x_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 10);

            auto tez_x_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 10);

            auto tex_x_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 11);

            auto tey_x_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 11);

            auto tez_x_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 11);

            auto tex_x_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 12);

            auto tey_x_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 12);

            auto tez_x_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 12);

            auto tex_x_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 13);

            auto tey_x_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 13);

            auto tez_x_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 13);

            auto tex_x_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 14);

            auto tey_x_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 14);

            auto tez_x_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 14);

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_0_xxxx_1, ta_0_xxxy_1, ta_0_xxxz_1, ta_0_xxyy_1, \
                                         ta_0_xxyz_1, ta_0_xxzz_1, ta_0_xyyy_1, ta_0_xyyz_1, ta_0_xyzz_1, ta_0_xzzz_1, \
                                         ta_0_yyyy_1, ta_0_yyyz_1, ta_0_yyzz_1, ta_0_yzzz_1, ta_0_zzzz_1, tex_0_xxx_0, \
                                         tex_0_xxx_1, tex_0_xxxx_0, tex_0_xxxx_1, tex_0_xxxy_0, tex_0_xxxy_1, tex_0_xxxz_0, \
                                         tex_0_xxxz_1, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxyy_0, tex_0_xxyy_1, tex_0_xxyz_0, \
                                         tex_0_xxyz_1, tex_0_xxz_0, tex_0_xxz_1, tex_0_xxzz_0, tex_0_xxzz_1, tex_0_xyy_0, \
                                         tex_0_xyy_1, tex_0_xyyy_0, tex_0_xyyy_1, tex_0_xyyz_0, tex_0_xyyz_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xyzz_0, tex_0_xyzz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_xzzz_0, \
                                         tex_0_xzzz_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyyy_0, tex_0_yyyy_1, tex_0_yyyz_0, \
                                         tex_0_yyyz_1, tex_0_yyz_0, tex_0_yyz_1, tex_0_yyzz_0, tex_0_yyzz_1, tex_0_yzz_0, \
                                         tex_0_yzz_1, tex_0_yzzz_0, tex_0_yzzz_1, tex_0_zzz_0, tex_0_zzz_1, tex_0_zzzz_0, \
                                         tex_0_zzzz_1, tex_x_xxxx_0, tex_x_xxxy_0, tex_x_xxxz_0, tex_x_xxyy_0, tex_x_xxyz_0, \
                                         tex_x_xxzz_0, tex_x_xyyy_0, tex_x_xyyz_0, tex_x_xyzz_0, tex_x_xzzz_0, tex_x_yyyy_0, \
                                         tex_x_yyyz_0, tex_x_yyzz_0, tex_x_yzzz_0, tex_x_zzzz_0, tey_0_xxx_0, tey_0_xxx_1, \
                                         tey_0_xxxx_0, tey_0_xxxx_1, tey_0_xxxy_0, tey_0_xxxy_1, tey_0_xxxz_0, tey_0_xxxz_1, \
                                         tey_0_xxy_0, tey_0_xxy_1, tey_0_xxyy_0, tey_0_xxyy_1, tey_0_xxyz_0, tey_0_xxyz_1, \
                                         tey_0_xxz_0, tey_0_xxz_1, tey_0_xxzz_0, tey_0_xxzz_1, tey_0_xyy_0, tey_0_xyy_1, \
                                         tey_0_xyyy_0, tey_0_xyyy_1, tey_0_xyyz_0, tey_0_xyyz_1, tey_0_xyz_0, tey_0_xyz_1, \
                                         tey_0_xyzz_0, tey_0_xyzz_1, tey_0_xzz_0, tey_0_xzz_1, tey_0_xzzz_0, tey_0_xzzz_1, \
                                         tey_0_yyy_0, tey_0_yyy_1, tey_0_yyyy_0, tey_0_yyyy_1, tey_0_yyyz_0, tey_0_yyyz_1, \
                                         tey_0_yyz_0, tey_0_yyz_1, tey_0_yyzz_0, tey_0_yyzz_1, tey_0_yzz_0, tey_0_yzz_1, \
                                         tey_0_yzzz_0, tey_0_yzzz_1, tey_0_zzz_0, tey_0_zzz_1, tey_0_zzzz_0, tey_0_zzzz_1, \
                                         tey_x_xxxx_0, tey_x_xxxy_0, tey_x_xxxz_0, tey_x_xxyy_0, tey_x_xxyz_0, tey_x_xxzz_0, \
                                         tey_x_xyyy_0, tey_x_xyyz_0, tey_x_xyzz_0, tey_x_xzzz_0, tey_x_yyyy_0, tey_x_yyyz_0, \
                                         tey_x_yyzz_0, tey_x_yzzz_0, tey_x_zzzz_0, tez_0_xxx_0, tez_0_xxx_1, tez_0_xxxx_0, \
                                         tez_0_xxxx_1, tez_0_xxxy_0, tez_0_xxxy_1, tez_0_xxxz_0, tez_0_xxxz_1, tez_0_xxy_0, \
                                         tez_0_xxy_1, tez_0_xxyy_0, tez_0_xxyy_1, tez_0_xxyz_0, tez_0_xxyz_1, tez_0_xxz_0, \
                                         tez_0_xxz_1, tez_0_xxzz_0, tez_0_xxzz_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyyy_0, \
                                         tez_0_xyyy_1, tez_0_xyyz_0, tez_0_xyyz_1, tez_0_xyz_0, tez_0_xyz_1, tez_0_xyzz_0, \
                                         tez_0_xyzz_1, tez_0_xzz_0, tez_0_xzz_1, tez_0_xzzz_0, tez_0_xzzz_1, tez_0_yyy_0, \
                                         tez_0_yyy_1, tez_0_yyyy_0, tez_0_yyyy_1, tez_0_yyyz_0, tez_0_yyyz_1, tez_0_yyz_0, \
                                         tez_0_yyz_1, tez_0_yyzz_0, tez_0_yyzz_1, tez_0_yzz_0, tez_0_yzz_1, tez_0_yzzz_0, \
                                         tez_0_yzzz_1, tez_0_zzz_0, tez_0_zzz_1, tez_0_zzzz_0, tez_0_zzzz_1, tez_x_xxxx_0, \
                                         tez_x_xxxy_0, tez_x_xxxz_0, tez_x_xxyy_0, tez_x_xxyz_0, tez_x_xxzz_0, tez_x_xyyy_0, \
                                         tez_x_xyyz_0, tez_x_xyzz_0, tez_x_xzzz_0, tez_x_yyyy_0, tez_x_yyyz_0, tez_x_yyzz_0, \
                                         tez_x_yzzz_0, tez_x_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_x_xxxx_0[j] = pa_x[j] * tex_0_xxxx_0[j] - pc_x[j] * tex_0_xxxx_1[j] + 2.0 * fl1_fx * tex_0_xxx_0[j] -
                                  2.0 * fl1_fx * tex_0_xxx_1[j] + ta_0_xxxx_1[j];

                tey_x_xxxx_0[j] =
                    pa_x[j] * tey_0_xxxx_0[j] - pc_x[j] * tey_0_xxxx_1[j] + 2.0 * fl1_fx * tey_0_xxx_0[j] - 2.0 * fl1_fx * tey_0_xxx_1[j];

                tez_x_xxxx_0[j] =
                    pa_x[j] * tez_0_xxxx_0[j] - pc_x[j] * tez_0_xxxx_1[j] + 2.0 * fl1_fx * tez_0_xxx_0[j] - 2.0 * fl1_fx * tez_0_xxx_1[j];

                tex_x_xxxy_0[j] = pa_x[j] * tex_0_xxxy_0[j] - pc_x[j] * tex_0_xxxy_1[j] + 1.5 * fl1_fx * tex_0_xxy_0[j] -
                                  1.5 * fl1_fx * tex_0_xxy_1[j] + ta_0_xxxy_1[j];

                tey_x_xxxy_0[j] =
                    pa_x[j] * tey_0_xxxy_0[j] - pc_x[j] * tey_0_xxxy_1[j] + 1.5 * fl1_fx * tey_0_xxy_0[j] - 1.5 * fl1_fx * tey_0_xxy_1[j];

                tez_x_xxxy_0[j] =
                    pa_x[j] * tez_0_xxxy_0[j] - pc_x[j] * tez_0_xxxy_1[j] + 1.5 * fl1_fx * tez_0_xxy_0[j] - 1.5 * fl1_fx * tez_0_xxy_1[j];

                tex_x_xxxz_0[j] = pa_x[j] * tex_0_xxxz_0[j] - pc_x[j] * tex_0_xxxz_1[j] + 1.5 * fl1_fx * tex_0_xxz_0[j] -
                                  1.5 * fl1_fx * tex_0_xxz_1[j] + ta_0_xxxz_1[j];

                tey_x_xxxz_0[j] =
                    pa_x[j] * tey_0_xxxz_0[j] - pc_x[j] * tey_0_xxxz_1[j] + 1.5 * fl1_fx * tey_0_xxz_0[j] - 1.5 * fl1_fx * tey_0_xxz_1[j];

                tez_x_xxxz_0[j] =
                    pa_x[j] * tez_0_xxxz_0[j] - pc_x[j] * tez_0_xxxz_1[j] + 1.5 * fl1_fx * tez_0_xxz_0[j] - 1.5 * fl1_fx * tez_0_xxz_1[j];

                tex_x_xxyy_0[j] =
                    pa_x[j] * tex_0_xxyy_0[j] - pc_x[j] * tex_0_xxyy_1[j] + fl1_fx * tex_0_xyy_0[j] - fl1_fx * tex_0_xyy_1[j] + ta_0_xxyy_1[j];

                tey_x_xxyy_0[j] = pa_x[j] * tey_0_xxyy_0[j] - pc_x[j] * tey_0_xxyy_1[j] + fl1_fx * tey_0_xyy_0[j] - fl1_fx * tey_0_xyy_1[j];

                tez_x_xxyy_0[j] = pa_x[j] * tez_0_xxyy_0[j] - pc_x[j] * tez_0_xxyy_1[j] + fl1_fx * tez_0_xyy_0[j] - fl1_fx * tez_0_xyy_1[j];

                tex_x_xxyz_0[j] =
                    pa_x[j] * tex_0_xxyz_0[j] - pc_x[j] * tex_0_xxyz_1[j] + fl1_fx * tex_0_xyz_0[j] - fl1_fx * tex_0_xyz_1[j] + ta_0_xxyz_1[j];

                tey_x_xxyz_0[j] = pa_x[j] * tey_0_xxyz_0[j] - pc_x[j] * tey_0_xxyz_1[j] + fl1_fx * tey_0_xyz_0[j] - fl1_fx * tey_0_xyz_1[j];

                tez_x_xxyz_0[j] = pa_x[j] * tez_0_xxyz_0[j] - pc_x[j] * tez_0_xxyz_1[j] + fl1_fx * tez_0_xyz_0[j] - fl1_fx * tez_0_xyz_1[j];

                tex_x_xxzz_0[j] =
                    pa_x[j] * tex_0_xxzz_0[j] - pc_x[j] * tex_0_xxzz_1[j] + fl1_fx * tex_0_xzz_0[j] - fl1_fx * tex_0_xzz_1[j] + ta_0_xxzz_1[j];

                tey_x_xxzz_0[j] = pa_x[j] * tey_0_xxzz_0[j] - pc_x[j] * tey_0_xxzz_1[j] + fl1_fx * tey_0_xzz_0[j] - fl1_fx * tey_0_xzz_1[j];

                tez_x_xxzz_0[j] = pa_x[j] * tez_0_xxzz_0[j] - pc_x[j] * tez_0_xxzz_1[j] + fl1_fx * tez_0_xzz_0[j] - fl1_fx * tez_0_xzz_1[j];

                tex_x_xyyy_0[j] = pa_x[j] * tex_0_xyyy_0[j] - pc_x[j] * tex_0_xyyy_1[j] + 0.5 * fl1_fx * tex_0_yyy_0[j] -
                                  0.5 * fl1_fx * tex_0_yyy_1[j] + ta_0_xyyy_1[j];

                tey_x_xyyy_0[j] =
                    pa_x[j] * tey_0_xyyy_0[j] - pc_x[j] * tey_0_xyyy_1[j] + 0.5 * fl1_fx * tey_0_yyy_0[j] - 0.5 * fl1_fx * tey_0_yyy_1[j];

                tez_x_xyyy_0[j] =
                    pa_x[j] * tez_0_xyyy_0[j] - pc_x[j] * tez_0_xyyy_1[j] + 0.5 * fl1_fx * tez_0_yyy_0[j] - 0.5 * fl1_fx * tez_0_yyy_1[j];

                tex_x_xyyz_0[j] = pa_x[j] * tex_0_xyyz_0[j] - pc_x[j] * tex_0_xyyz_1[j] + 0.5 * fl1_fx * tex_0_yyz_0[j] -
                                  0.5 * fl1_fx * tex_0_yyz_1[j] + ta_0_xyyz_1[j];

                tey_x_xyyz_0[j] =
                    pa_x[j] * tey_0_xyyz_0[j] - pc_x[j] * tey_0_xyyz_1[j] + 0.5 * fl1_fx * tey_0_yyz_0[j] - 0.5 * fl1_fx * tey_0_yyz_1[j];

                tez_x_xyyz_0[j] =
                    pa_x[j] * tez_0_xyyz_0[j] - pc_x[j] * tez_0_xyyz_1[j] + 0.5 * fl1_fx * tez_0_yyz_0[j] - 0.5 * fl1_fx * tez_0_yyz_1[j];

                tex_x_xyzz_0[j] = pa_x[j] * tex_0_xyzz_0[j] - pc_x[j] * tex_0_xyzz_1[j] + 0.5 * fl1_fx * tex_0_yzz_0[j] -
                                  0.5 * fl1_fx * tex_0_yzz_1[j] + ta_0_xyzz_1[j];

                tey_x_xyzz_0[j] =
                    pa_x[j] * tey_0_xyzz_0[j] - pc_x[j] * tey_0_xyzz_1[j] + 0.5 * fl1_fx * tey_0_yzz_0[j] - 0.5 * fl1_fx * tey_0_yzz_1[j];

                tez_x_xyzz_0[j] =
                    pa_x[j] * tez_0_xyzz_0[j] - pc_x[j] * tez_0_xyzz_1[j] + 0.5 * fl1_fx * tez_0_yzz_0[j] - 0.5 * fl1_fx * tez_0_yzz_1[j];

                tex_x_xzzz_0[j] = pa_x[j] * tex_0_xzzz_0[j] - pc_x[j] * tex_0_xzzz_1[j] + 0.5 * fl1_fx * tex_0_zzz_0[j] -
                                  0.5 * fl1_fx * tex_0_zzz_1[j] + ta_0_xzzz_1[j];

                tey_x_xzzz_0[j] =
                    pa_x[j] * tey_0_xzzz_0[j] - pc_x[j] * tey_0_xzzz_1[j] + 0.5 * fl1_fx * tey_0_zzz_0[j] - 0.5 * fl1_fx * tey_0_zzz_1[j];

                tez_x_xzzz_0[j] =
                    pa_x[j] * tez_0_xzzz_0[j] - pc_x[j] * tez_0_xzzz_1[j] + 0.5 * fl1_fx * tez_0_zzz_0[j] - 0.5 * fl1_fx * tez_0_zzz_1[j];

                tex_x_yyyy_0[j] = pa_x[j] * tex_0_yyyy_0[j] - pc_x[j] * tex_0_yyyy_1[j] + ta_0_yyyy_1[j];

                tey_x_yyyy_0[j] = pa_x[j] * tey_0_yyyy_0[j] - pc_x[j] * tey_0_yyyy_1[j];

                tez_x_yyyy_0[j] = pa_x[j] * tez_0_yyyy_0[j] - pc_x[j] * tez_0_yyyy_1[j];

                tex_x_yyyz_0[j] = pa_x[j] * tex_0_yyyz_0[j] - pc_x[j] * tex_0_yyyz_1[j] + ta_0_yyyz_1[j];

                tey_x_yyyz_0[j] = pa_x[j] * tey_0_yyyz_0[j] - pc_x[j] * tey_0_yyyz_1[j];

                tez_x_yyyz_0[j] = pa_x[j] * tez_0_yyyz_0[j] - pc_x[j] * tez_0_yyyz_1[j];

                tex_x_yyzz_0[j] = pa_x[j] * tex_0_yyzz_0[j] - pc_x[j] * tex_0_yyzz_1[j] + ta_0_yyzz_1[j];

                tey_x_yyzz_0[j] = pa_x[j] * tey_0_yyzz_0[j] - pc_x[j] * tey_0_yyzz_1[j];

                tez_x_yyzz_0[j] = pa_x[j] * tez_0_yyzz_0[j] - pc_x[j] * tez_0_yyzz_1[j];

                tex_x_yzzz_0[j] = pa_x[j] * tex_0_yzzz_0[j] - pc_x[j] * tex_0_yzzz_1[j] + ta_0_yzzz_1[j];

                tey_x_yzzz_0[j] = pa_x[j] * tey_0_yzzz_0[j] - pc_x[j] * tey_0_yzzz_1[j];

                tez_x_yzzz_0[j] = pa_x[j] * tez_0_yzzz_0[j] - pc_x[j] * tez_0_yzzz_1[j];

                tex_x_zzzz_0[j] = pa_x[j] * tex_0_zzzz_0[j] - pc_x[j] * tex_0_zzzz_1[j] + ta_0_zzzz_1[j];

                tey_x_zzzz_0[j] = pa_x[j] * tey_0_zzzz_0[j] - pc_x[j] * tey_0_zzzz_1[j];

                tez_x_zzzz_0[j] = pa_x[j] * tez_0_zzzz_0[j] - pc_x[j] * tez_0_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPG_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_y = paDistances.data(3 * idx + 1);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_y = pcDistances.data(3 * idx + 1);

            // set up pointers to auxilary integrals

            auto tex_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx);

            auto tey_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx);

            auto tez_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx);

            auto tex_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 1);

            auto tey_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 1);

            auto tez_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 1);

            auto tex_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 2);

            auto tey_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 2);

            auto tez_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 2);

            auto tex_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 3);

            auto tey_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 3);

            auto tez_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 3);

            auto tex_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 4);

            auto tey_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 4);

            auto tez_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 4);

            auto tex_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 5);

            auto tey_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 5);

            auto tez_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 5);

            auto tex_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 6);

            auto tey_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 6);

            auto tez_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 6);

            auto tex_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 7);

            auto tey_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 7);

            auto tez_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 7);

            auto tex_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 8);

            auto tey_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 8);

            auto tez_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 8);

            auto tex_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 9);

            auto tey_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 9);

            auto tez_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 9);

            auto tex_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 10);

            auto tey_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 10);

            auto tez_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 10);

            auto tex_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 11);

            auto tey_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 11);

            auto tez_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 11);

            auto tex_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 12);

            auto tey_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 12);

            auto tez_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 12);

            auto tex_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 13);

            auto tey_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 13);

            auto tez_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 13);

            auto tex_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 14);

            auto tey_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 14);

            auto tez_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 14);

            auto tex_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx);

            auto tey_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx);

            auto tez_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx);

            auto tex_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 1);

            auto tey_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 1);

            auto tez_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 1);

            auto tex_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 2);

            auto tey_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 2);

            auto tez_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 2);

            auto tex_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 3);

            auto tey_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 3);

            auto tez_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 3);

            auto tex_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 4);

            auto tey_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 4);

            auto tez_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 4);

            auto tex_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 5);

            auto tey_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 5);

            auto tez_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 5);

            auto tex_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 6);

            auto tey_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 6);

            auto tez_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 6);

            auto tex_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 7);

            auto tey_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 7);

            auto tez_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 7);

            auto tex_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 8);

            auto tey_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 8);

            auto tez_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 8);

            auto tex_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 9);

            auto tey_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 9);

            auto tez_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 9);

            auto tex_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 10);

            auto tey_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 10);

            auto tez_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 10);

            auto tex_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 11);

            auto tey_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 11);

            auto tez_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 11);

            auto tex_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 12);

            auto tey_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 12);

            auto tez_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 12);

            auto tex_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 13);

            auto tey_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 13);

            auto tez_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 13);

            auto tex_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 14);

            auto tey_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 14);

            auto tez_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 14);

            auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx);

            auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1);

            auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2);

            auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3);

            auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4);

            auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5);

            auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6);

            auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7);

            auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8);

            auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9);

            auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx);

            auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1);

            auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2);

            auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3);

            auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4);

            auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5);

            auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6);

            auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7);

            auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8);

            auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9);

            auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9);

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

            // set up pointers to integrals

            auto tex_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 15);

            auto tey_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 15);

            auto tez_y_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 15);

            auto tex_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 16);

            auto tey_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 16);

            auto tez_y_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 16);

            auto tex_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 17);

            auto tey_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 17);

            auto tez_y_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 17);

            auto tex_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 18);

            auto tey_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 18);

            auto tez_y_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 18);

            auto tex_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 19);

            auto tey_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 19);

            auto tez_y_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 19);

            auto tex_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 20);

            auto tey_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 20);

            auto tez_y_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 20);

            auto tex_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 21);

            auto tey_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 21);

            auto tez_y_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 21);

            auto tex_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 22);

            auto tey_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 22);

            auto tez_y_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 22);

            auto tex_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 23);

            auto tey_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 23);

            auto tez_y_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 23);

            auto tex_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 24);

            auto tey_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 24);

            auto tez_y_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 24);

            auto tex_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 25);

            auto tey_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 25);

            auto tez_y_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 25);

            auto tex_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 26);

            auto tey_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 26);

            auto tez_y_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 26);

            auto tex_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 27);

            auto tey_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 27);

            auto tez_y_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 27);

            auto tex_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 28);

            auto tey_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 28);

            auto tez_y_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 28);

            auto tex_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 29);

            auto tey_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 29);

            auto tez_y_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 29);

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_y, pc_y, ta_0_xxxx_1, ta_0_xxxy_1, ta_0_xxxz_1, ta_0_xxyy_1, \
                                         ta_0_xxyz_1, ta_0_xxzz_1, ta_0_xyyy_1, ta_0_xyyz_1, ta_0_xyzz_1, ta_0_xzzz_1, \
                                         ta_0_yyyy_1, ta_0_yyyz_1, ta_0_yyzz_1, ta_0_yzzz_1, ta_0_zzzz_1, tex_0_xxx_0, \
                                         tex_0_xxx_1, tex_0_xxxx_0, tex_0_xxxx_1, tex_0_xxxy_0, tex_0_xxxy_1, tex_0_xxxz_0, \
                                         tex_0_xxxz_1, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxyy_0, tex_0_xxyy_1, tex_0_xxyz_0, \
                                         tex_0_xxyz_1, tex_0_xxz_0, tex_0_xxz_1, tex_0_xxzz_0, tex_0_xxzz_1, tex_0_xyy_0, \
                                         tex_0_xyy_1, tex_0_xyyy_0, tex_0_xyyy_1, tex_0_xyyz_0, tex_0_xyyz_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xyzz_0, tex_0_xyzz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_xzzz_0, \
                                         tex_0_xzzz_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyyy_0, tex_0_yyyy_1, tex_0_yyyz_0, \
                                         tex_0_yyyz_1, tex_0_yyz_0, tex_0_yyz_1, tex_0_yyzz_0, tex_0_yyzz_1, tex_0_yzz_0, \
                                         tex_0_yzz_1, tex_0_yzzz_0, tex_0_yzzz_1, tex_0_zzz_0, tex_0_zzz_1, tex_0_zzzz_0, \
                                         tex_0_zzzz_1, tex_y_xxxx_0, tex_y_xxxy_0, tex_y_xxxz_0, tex_y_xxyy_0, tex_y_xxyz_0, \
                                         tex_y_xxzz_0, tex_y_xyyy_0, tex_y_xyyz_0, tex_y_xyzz_0, tex_y_xzzz_0, tex_y_yyyy_0, \
                                         tex_y_yyyz_0, tex_y_yyzz_0, tex_y_yzzz_0, tex_y_zzzz_0, tey_0_xxx_0, tey_0_xxx_1, \
                                         tey_0_xxxx_0, tey_0_xxxx_1, tey_0_xxxy_0, tey_0_xxxy_1, tey_0_xxxz_0, tey_0_xxxz_1, \
                                         tey_0_xxy_0, tey_0_xxy_1, tey_0_xxyy_0, tey_0_xxyy_1, tey_0_xxyz_0, tey_0_xxyz_1, \
                                         tey_0_xxz_0, tey_0_xxz_1, tey_0_xxzz_0, tey_0_xxzz_1, tey_0_xyy_0, tey_0_xyy_1, \
                                         tey_0_xyyy_0, tey_0_xyyy_1, tey_0_xyyz_0, tey_0_xyyz_1, tey_0_xyz_0, tey_0_xyz_1, \
                                         tey_0_xyzz_0, tey_0_xyzz_1, tey_0_xzz_0, tey_0_xzz_1, tey_0_xzzz_0, tey_0_xzzz_1, \
                                         tey_0_yyy_0, tey_0_yyy_1, tey_0_yyyy_0, tey_0_yyyy_1, tey_0_yyyz_0, tey_0_yyyz_1, \
                                         tey_0_yyz_0, tey_0_yyz_1, tey_0_yyzz_0, tey_0_yyzz_1, tey_0_yzz_0, tey_0_yzz_1, \
                                         tey_0_yzzz_0, tey_0_yzzz_1, tey_0_zzz_0, tey_0_zzz_1, tey_0_zzzz_0, tey_0_zzzz_1, \
                                         tey_y_xxxx_0, tey_y_xxxy_0, tey_y_xxxz_0, tey_y_xxyy_0, tey_y_xxyz_0, tey_y_xxzz_0, \
                                         tey_y_xyyy_0, tey_y_xyyz_0, tey_y_xyzz_0, tey_y_xzzz_0, tey_y_yyyy_0, tey_y_yyyz_0, \
                                         tey_y_yyzz_0, tey_y_yzzz_0, tey_y_zzzz_0, tez_0_xxx_0, tez_0_xxx_1, tez_0_xxxx_0, \
                                         tez_0_xxxx_1, tez_0_xxxy_0, tez_0_xxxy_1, tez_0_xxxz_0, tez_0_xxxz_1, tez_0_xxy_0, \
                                         tez_0_xxy_1, tez_0_xxyy_0, tez_0_xxyy_1, tez_0_xxyz_0, tez_0_xxyz_1, tez_0_xxz_0, \
                                         tez_0_xxz_1, tez_0_xxzz_0, tez_0_xxzz_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyyy_0, \
                                         tez_0_xyyy_1, tez_0_xyyz_0, tez_0_xyyz_1, tez_0_xyz_0, tez_0_xyz_1, tez_0_xyzz_0, \
                                         tez_0_xyzz_1, tez_0_xzz_0, tez_0_xzz_1, tez_0_xzzz_0, tez_0_xzzz_1, tez_0_yyy_0, \
                                         tez_0_yyy_1, tez_0_yyyy_0, tez_0_yyyy_1, tez_0_yyyz_0, tez_0_yyyz_1, tez_0_yyz_0, \
                                         tez_0_yyz_1, tez_0_yyzz_0, tez_0_yyzz_1, tez_0_yzz_0, tez_0_yzz_1, tez_0_yzzz_0, \
                                         tez_0_yzzz_1, tez_0_zzz_0, tez_0_zzz_1, tez_0_zzzz_0, tez_0_zzzz_1, tez_y_xxxx_0, \
                                         tez_y_xxxy_0, tez_y_xxxz_0, tez_y_xxyy_0, tez_y_xxyz_0, tez_y_xxzz_0, tez_y_xyyy_0, \
                                         tez_y_xyyz_0, tez_y_xyzz_0, tez_y_xzzz_0, tez_y_yyyy_0, tez_y_yyyz_0, tez_y_yyzz_0, \
                                         tez_y_yzzz_0, tez_y_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_y_xxxx_0[j] = pa_y[j] * tex_0_xxxx_0[j] - pc_y[j] * tex_0_xxxx_1[j];

                tey_y_xxxx_0[j] = pa_y[j] * tey_0_xxxx_0[j] - pc_y[j] * tey_0_xxxx_1[j] + ta_0_xxxx_1[j];

                tez_y_xxxx_0[j] = pa_y[j] * tez_0_xxxx_0[j] - pc_y[j] * tez_0_xxxx_1[j];

                tex_y_xxxy_0[j] =
                    pa_y[j] * tex_0_xxxy_0[j] - pc_y[j] * tex_0_xxxy_1[j] + 0.5 * fl1_fx * tex_0_xxx_0[j] - 0.5 * fl1_fx * tex_0_xxx_1[j];

                tey_y_xxxy_0[j] = pa_y[j] * tey_0_xxxy_0[j] - pc_y[j] * tey_0_xxxy_1[j] + 0.5 * fl1_fx * tey_0_xxx_0[j] -
                                  0.5 * fl1_fx * tey_0_xxx_1[j] + ta_0_xxxy_1[j];

                tez_y_xxxy_0[j] =
                    pa_y[j] * tez_0_xxxy_0[j] - pc_y[j] * tez_0_xxxy_1[j] + 0.5 * fl1_fx * tez_0_xxx_0[j] - 0.5 * fl1_fx * tez_0_xxx_1[j];

                tex_y_xxxz_0[j] = pa_y[j] * tex_0_xxxz_0[j] - pc_y[j] * tex_0_xxxz_1[j];

                tey_y_xxxz_0[j] = pa_y[j] * tey_0_xxxz_0[j] - pc_y[j] * tey_0_xxxz_1[j] + ta_0_xxxz_1[j];

                tez_y_xxxz_0[j] = pa_y[j] * tez_0_xxxz_0[j] - pc_y[j] * tez_0_xxxz_1[j];

                tex_y_xxyy_0[j] = pa_y[j] * tex_0_xxyy_0[j] - pc_y[j] * tex_0_xxyy_1[j] + fl1_fx * tex_0_xxy_0[j] - fl1_fx * tex_0_xxy_1[j];

                tey_y_xxyy_0[j] =
                    pa_y[j] * tey_0_xxyy_0[j] - pc_y[j] * tey_0_xxyy_1[j] + fl1_fx * tey_0_xxy_0[j] - fl1_fx * tey_0_xxy_1[j] + ta_0_xxyy_1[j];

                tez_y_xxyy_0[j] = pa_y[j] * tez_0_xxyy_0[j] - pc_y[j] * tez_0_xxyy_1[j] + fl1_fx * tez_0_xxy_0[j] - fl1_fx * tez_0_xxy_1[j];

                tex_y_xxyz_0[j] =
                    pa_y[j] * tex_0_xxyz_0[j] - pc_y[j] * tex_0_xxyz_1[j] + 0.5 * fl1_fx * tex_0_xxz_0[j] - 0.5 * fl1_fx * tex_0_xxz_1[j];

                tey_y_xxyz_0[j] = pa_y[j] * tey_0_xxyz_0[j] - pc_y[j] * tey_0_xxyz_1[j] + 0.5 * fl1_fx * tey_0_xxz_0[j] -
                                  0.5 * fl1_fx * tey_0_xxz_1[j] + ta_0_xxyz_1[j];

                tez_y_xxyz_0[j] =
                    pa_y[j] * tez_0_xxyz_0[j] - pc_y[j] * tez_0_xxyz_1[j] + 0.5 * fl1_fx * tez_0_xxz_0[j] - 0.5 * fl1_fx * tez_0_xxz_1[j];

                tex_y_xxzz_0[j] = pa_y[j] * tex_0_xxzz_0[j] - pc_y[j] * tex_0_xxzz_1[j];

                tey_y_xxzz_0[j] = pa_y[j] * tey_0_xxzz_0[j] - pc_y[j] * tey_0_xxzz_1[j] + ta_0_xxzz_1[j];

                tez_y_xxzz_0[j] = pa_y[j] * tez_0_xxzz_0[j] - pc_y[j] * tez_0_xxzz_1[j];

                tex_y_xyyy_0[j] =
                    pa_y[j] * tex_0_xyyy_0[j] - pc_y[j] * tex_0_xyyy_1[j] + 1.5 * fl1_fx * tex_0_xyy_0[j] - 1.5 * fl1_fx * tex_0_xyy_1[j];

                tey_y_xyyy_0[j] = pa_y[j] * tey_0_xyyy_0[j] - pc_y[j] * tey_0_xyyy_1[j] + 1.5 * fl1_fx * tey_0_xyy_0[j] -
                                  1.5 * fl1_fx * tey_0_xyy_1[j] + ta_0_xyyy_1[j];

                tez_y_xyyy_0[j] =
                    pa_y[j] * tez_0_xyyy_0[j] - pc_y[j] * tez_0_xyyy_1[j] + 1.5 * fl1_fx * tez_0_xyy_0[j] - 1.5 * fl1_fx * tez_0_xyy_1[j];

                tex_y_xyyz_0[j] = pa_y[j] * tex_0_xyyz_0[j] - pc_y[j] * tex_0_xyyz_1[j] + fl1_fx * tex_0_xyz_0[j] - fl1_fx * tex_0_xyz_1[j];

                tey_y_xyyz_0[j] =
                    pa_y[j] * tey_0_xyyz_0[j] - pc_y[j] * tey_0_xyyz_1[j] + fl1_fx * tey_0_xyz_0[j] - fl1_fx * tey_0_xyz_1[j] + ta_0_xyyz_1[j];

                tez_y_xyyz_0[j] = pa_y[j] * tez_0_xyyz_0[j] - pc_y[j] * tez_0_xyyz_1[j] + fl1_fx * tez_0_xyz_0[j] - fl1_fx * tez_0_xyz_1[j];

                tex_y_xyzz_0[j] =
                    pa_y[j] * tex_0_xyzz_0[j] - pc_y[j] * tex_0_xyzz_1[j] + 0.5 * fl1_fx * tex_0_xzz_0[j] - 0.5 * fl1_fx * tex_0_xzz_1[j];

                tey_y_xyzz_0[j] = pa_y[j] * tey_0_xyzz_0[j] - pc_y[j] * tey_0_xyzz_1[j] + 0.5 * fl1_fx * tey_0_xzz_0[j] -
                                  0.5 * fl1_fx * tey_0_xzz_1[j] + ta_0_xyzz_1[j];

                tez_y_xyzz_0[j] =
                    pa_y[j] * tez_0_xyzz_0[j] - pc_y[j] * tez_0_xyzz_1[j] + 0.5 * fl1_fx * tez_0_xzz_0[j] - 0.5 * fl1_fx * tez_0_xzz_1[j];

                tex_y_xzzz_0[j] = pa_y[j] * tex_0_xzzz_0[j] - pc_y[j] * tex_0_xzzz_1[j];

                tey_y_xzzz_0[j] = pa_y[j] * tey_0_xzzz_0[j] - pc_y[j] * tey_0_xzzz_1[j] + ta_0_xzzz_1[j];

                tez_y_xzzz_0[j] = pa_y[j] * tez_0_xzzz_0[j] - pc_y[j] * tez_0_xzzz_1[j];

                tex_y_yyyy_0[j] =
                    pa_y[j] * tex_0_yyyy_0[j] - pc_y[j] * tex_0_yyyy_1[j] + 2.0 * fl1_fx * tex_0_yyy_0[j] - 2.0 * fl1_fx * tex_0_yyy_1[j];

                tey_y_yyyy_0[j] = pa_y[j] * tey_0_yyyy_0[j] - pc_y[j] * tey_0_yyyy_1[j] + 2.0 * fl1_fx * tey_0_yyy_0[j] -
                                  2.0 * fl1_fx * tey_0_yyy_1[j] + ta_0_yyyy_1[j];

                tez_y_yyyy_0[j] =
                    pa_y[j] * tez_0_yyyy_0[j] - pc_y[j] * tez_0_yyyy_1[j] + 2.0 * fl1_fx * tez_0_yyy_0[j] - 2.0 * fl1_fx * tez_0_yyy_1[j];

                tex_y_yyyz_0[j] =
                    pa_y[j] * tex_0_yyyz_0[j] - pc_y[j] * tex_0_yyyz_1[j] + 1.5 * fl1_fx * tex_0_yyz_0[j] - 1.5 * fl1_fx * tex_0_yyz_1[j];

                tey_y_yyyz_0[j] = pa_y[j] * tey_0_yyyz_0[j] - pc_y[j] * tey_0_yyyz_1[j] + 1.5 * fl1_fx * tey_0_yyz_0[j] -
                                  1.5 * fl1_fx * tey_0_yyz_1[j] + ta_0_yyyz_1[j];

                tez_y_yyyz_0[j] =
                    pa_y[j] * tez_0_yyyz_0[j] - pc_y[j] * tez_0_yyyz_1[j] + 1.5 * fl1_fx * tez_0_yyz_0[j] - 1.5 * fl1_fx * tez_0_yyz_1[j];

                tex_y_yyzz_0[j] = pa_y[j] * tex_0_yyzz_0[j] - pc_y[j] * tex_0_yyzz_1[j] + fl1_fx * tex_0_yzz_0[j] - fl1_fx * tex_0_yzz_1[j];

                tey_y_yyzz_0[j] =
                    pa_y[j] * tey_0_yyzz_0[j] - pc_y[j] * tey_0_yyzz_1[j] + fl1_fx * tey_0_yzz_0[j] - fl1_fx * tey_0_yzz_1[j] + ta_0_yyzz_1[j];

                tez_y_yyzz_0[j] = pa_y[j] * tez_0_yyzz_0[j] - pc_y[j] * tez_0_yyzz_1[j] + fl1_fx * tez_0_yzz_0[j] - fl1_fx * tez_0_yzz_1[j];

                tex_y_yzzz_0[j] =
                    pa_y[j] * tex_0_yzzz_0[j] - pc_y[j] * tex_0_yzzz_1[j] + 0.5 * fl1_fx * tex_0_zzz_0[j] - 0.5 * fl1_fx * tex_0_zzz_1[j];

                tey_y_yzzz_0[j] = pa_y[j] * tey_0_yzzz_0[j] - pc_y[j] * tey_0_yzzz_1[j] + 0.5 * fl1_fx * tey_0_zzz_0[j] -
                                  0.5 * fl1_fx * tey_0_zzz_1[j] + ta_0_yzzz_1[j];

                tez_y_yzzz_0[j] =
                    pa_y[j] * tez_0_yzzz_0[j] - pc_y[j] * tez_0_yzzz_1[j] + 0.5 * fl1_fx * tez_0_zzz_0[j] - 0.5 * fl1_fx * tez_0_zzz_1[j];

                tex_y_zzzz_0[j] = pa_y[j] * tex_0_zzzz_0[j] - pc_y[j] * tex_0_zzzz_1[j];

                tey_y_zzzz_0[j] = pa_y[j] * tey_0_zzzz_0[j] - pc_y[j] * tey_0_zzzz_1[j] + ta_0_zzzz_1[j];

                tez_y_zzzz_0[j] = pa_y[j] * tez_0_zzzz_0[j] - pc_y[j] * tez_0_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPG_90_135(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (90,135)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_4_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_4_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_3_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to Obara-Saika factors

            auto fx = osFactors.data(3 * idx);

            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tex_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx);

            auto tey_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx);

            auto tez_0_xxxx_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx);

            auto tex_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 1);

            auto tey_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 1);

            auto tez_0_xxxy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 1);

            auto tex_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 2);

            auto tey_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 2);

            auto tez_0_xxxz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 2);

            auto tex_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 3);

            auto tey_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 3);

            auto tez_0_xxyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 3);

            auto tex_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 4);

            auto tey_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 4);

            auto tez_0_xxyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 4);

            auto tex_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 5);

            auto tey_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 5);

            auto tez_0_xxzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 5);

            auto tex_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 6);

            auto tey_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 6);

            auto tez_0_xyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 6);

            auto tex_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 7);

            auto tey_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 7);

            auto tez_0_xyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 7);

            auto tex_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 8);

            auto tey_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 8);

            auto tez_0_xyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 8);

            auto tex_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 9);

            auto tey_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 9);

            auto tez_0_xzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 9);

            auto tex_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 10);

            auto tey_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 10);

            auto tez_0_yyyy_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 10);

            auto tex_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 11);

            auto tey_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 11);

            auto tez_0_yyyz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 11);

            auto tex_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 12);

            auto tey_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 12);

            auto tez_0_yyzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 12);

            auto tex_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 13);

            auto tey_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 13);

            auto tez_0_yzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 13);

            auto tex_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * idx + 14);

            auto tey_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 15 * bdim + 15 * idx + 14);

            auto tez_0_zzzz_0 = primBuffer.data(pidx_e_0_4_m0 + 30 * bdim + 15 * idx + 14);

            auto tex_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx);

            auto tey_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx);

            auto tez_0_xxxx_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx);

            auto tex_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 1);

            auto tey_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 1);

            auto tez_0_xxxy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 1);

            auto tex_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 2);

            auto tey_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 2);

            auto tez_0_xxxz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 2);

            auto tex_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 3);

            auto tey_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 3);

            auto tez_0_xxyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 3);

            auto tex_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 4);

            auto tey_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 4);

            auto tez_0_xxyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 4);

            auto tex_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 5);

            auto tey_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 5);

            auto tez_0_xxzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 5);

            auto tex_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 6);

            auto tey_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 6);

            auto tez_0_xyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 6);

            auto tex_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 7);

            auto tey_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 7);

            auto tez_0_xyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 7);

            auto tex_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 8);

            auto tey_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 8);

            auto tez_0_xyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 8);

            auto tex_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 9);

            auto tey_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 9);

            auto tez_0_xzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 9);

            auto tex_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 10);

            auto tey_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 10);

            auto tez_0_yyyy_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 10);

            auto tex_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 11);

            auto tey_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 11);

            auto tez_0_yyyz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 11);

            auto tex_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 12);

            auto tey_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 12);

            auto tez_0_yyzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 12);

            auto tex_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 13);

            auto tey_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 13);

            auto tez_0_yzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 13);

            auto tex_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * idx + 14);

            auto tey_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 15 * bdim + 15 * idx + 14);

            auto tez_0_zzzz_1 = primBuffer.data(pidx_e_0_4_m1 + 30 * bdim + 15 * idx + 14);

            auto tex_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx);

            auto tey_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 1);

            auto tey_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 2);

            auto tey_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 3);

            auto tey_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 4);

            auto tey_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 5);

            auto tey_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 6);

            auto tey_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 7);

            auto tey_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 8);

            auto tey_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * idx + 9);

            auto tey_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_0 = primBuffer.data(pidx_e_0_3_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx);

            auto tey_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx);

            auto tez_0_xxx_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx);

            auto tex_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 1);

            auto tey_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 1);

            auto tez_0_xxy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 1);

            auto tex_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 2);

            auto tey_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 2);

            auto tez_0_xxz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 2);

            auto tex_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 3);

            auto tey_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 3);

            auto tez_0_xyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 3);

            auto tex_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 4);

            auto tey_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 4);

            auto tez_0_xyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 4);

            auto tex_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 5);

            auto tey_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 5);

            auto tez_0_xzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 5);

            auto tex_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 6);

            auto tey_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_0_yyy_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 7);

            auto tey_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_0_yyz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 8);

            auto tey_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_0_yzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * idx + 9);

            auto tey_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_0_zzz_1 = primBuffer.data(pidx_e_0_3_m1 + 20 * bdim + 10 * idx + 9);

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

            // set up pointers to integrals

            auto tex_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 30);

            auto tey_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 30);

            auto tez_z_xxxx_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 30);

            auto tex_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 31);

            auto tey_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 31);

            auto tez_z_xxxy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 31);

            auto tex_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 32);

            auto tey_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 32);

            auto tez_z_xxxz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 32);

            auto tex_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 33);

            auto tey_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 33);

            auto tez_z_xxyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 33);

            auto tex_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 34);

            auto tey_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 34);

            auto tez_z_xxyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 34);

            auto tex_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 35);

            auto tey_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 35);

            auto tez_z_xxzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 35);

            auto tex_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 36);

            auto tey_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 36);

            auto tez_z_xyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 36);

            auto tex_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 37);

            auto tey_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 37);

            auto tez_z_xyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 37);

            auto tex_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 38);

            auto tey_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 38);

            auto tez_z_xyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 38);

            auto tex_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 39);

            auto tey_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 39);

            auto tez_z_xzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 39);

            auto tex_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 40);

            auto tey_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 40);

            auto tez_z_yyyy_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 40);

            auto tex_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 41);

            auto tey_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 41);

            auto tez_z_yyyz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 41);

            auto tex_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 42);

            auto tey_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 42);

            auto tez_z_yyzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 42);

            auto tex_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 43);

            auto tey_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 43);

            auto tez_z_yzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 43);

            auto tex_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * idx + 44);

            auto tey_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 45 * bdim + 45 * idx + 44);

            auto tez_z_zzzz_0 = primBuffer.data(pidx_e_1_4_m0 + 90 * bdim + 45 * idx + 44);

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_z, pc_z, ta_0_xxxx_1, ta_0_xxxy_1, ta_0_xxxz_1, ta_0_xxyy_1, \
                                         ta_0_xxyz_1, ta_0_xxzz_1, ta_0_xyyy_1, ta_0_xyyz_1, ta_0_xyzz_1, ta_0_xzzz_1, \
                                         ta_0_yyyy_1, ta_0_yyyz_1, ta_0_yyzz_1, ta_0_yzzz_1, ta_0_zzzz_1, tex_0_xxx_0, \
                                         tex_0_xxx_1, tex_0_xxxx_0, tex_0_xxxx_1, tex_0_xxxy_0, tex_0_xxxy_1, tex_0_xxxz_0, \
                                         tex_0_xxxz_1, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxyy_0, tex_0_xxyy_1, tex_0_xxyz_0, \
                                         tex_0_xxyz_1, tex_0_xxz_0, tex_0_xxz_1, tex_0_xxzz_0, tex_0_xxzz_1, tex_0_xyy_0, \
                                         tex_0_xyy_1, tex_0_xyyy_0, tex_0_xyyy_1, tex_0_xyyz_0, tex_0_xyyz_1, tex_0_xyz_0, \
                                         tex_0_xyz_1, tex_0_xyzz_0, tex_0_xyzz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_xzzz_0, \
                                         tex_0_xzzz_1, tex_0_yyy_0, tex_0_yyy_1, tex_0_yyyy_0, tex_0_yyyy_1, tex_0_yyyz_0, \
                                         tex_0_yyyz_1, tex_0_yyz_0, tex_0_yyz_1, tex_0_yyzz_0, tex_0_yyzz_1, tex_0_yzz_0, \
                                         tex_0_yzz_1, tex_0_yzzz_0, tex_0_yzzz_1, tex_0_zzz_0, tex_0_zzz_1, tex_0_zzzz_0, \
                                         tex_0_zzzz_1, tex_z_xxxx_0, tex_z_xxxy_0, tex_z_xxxz_0, tex_z_xxyy_0, tex_z_xxyz_0, \
                                         tex_z_xxzz_0, tex_z_xyyy_0, tex_z_xyyz_0, tex_z_xyzz_0, tex_z_xzzz_0, tex_z_yyyy_0, \
                                         tex_z_yyyz_0, tex_z_yyzz_0, tex_z_yzzz_0, tex_z_zzzz_0, tey_0_xxx_0, tey_0_xxx_1, \
                                         tey_0_xxxx_0, tey_0_xxxx_1, tey_0_xxxy_0, tey_0_xxxy_1, tey_0_xxxz_0, tey_0_xxxz_1, \
                                         tey_0_xxy_0, tey_0_xxy_1, tey_0_xxyy_0, tey_0_xxyy_1, tey_0_xxyz_0, tey_0_xxyz_1, \
                                         tey_0_xxz_0, tey_0_xxz_1, tey_0_xxzz_0, tey_0_xxzz_1, tey_0_xyy_0, tey_0_xyy_1, \
                                         tey_0_xyyy_0, tey_0_xyyy_1, tey_0_xyyz_0, tey_0_xyyz_1, tey_0_xyz_0, tey_0_xyz_1, \
                                         tey_0_xyzz_0, tey_0_xyzz_1, tey_0_xzz_0, tey_0_xzz_1, tey_0_xzzz_0, tey_0_xzzz_1, \
                                         tey_0_yyy_0, tey_0_yyy_1, tey_0_yyyy_0, tey_0_yyyy_1, tey_0_yyyz_0, tey_0_yyyz_1, \
                                         tey_0_yyz_0, tey_0_yyz_1, tey_0_yyzz_0, tey_0_yyzz_1, tey_0_yzz_0, tey_0_yzz_1, \
                                         tey_0_yzzz_0, tey_0_yzzz_1, tey_0_zzz_0, tey_0_zzz_1, tey_0_zzzz_0, tey_0_zzzz_1, \
                                         tey_z_xxxx_0, tey_z_xxxy_0, tey_z_xxxz_0, tey_z_xxyy_0, tey_z_xxyz_0, tey_z_xxzz_0, \
                                         tey_z_xyyy_0, tey_z_xyyz_0, tey_z_xyzz_0, tey_z_xzzz_0, tey_z_yyyy_0, tey_z_yyyz_0, \
                                         tey_z_yyzz_0, tey_z_yzzz_0, tey_z_zzzz_0, tez_0_xxx_0, tez_0_xxx_1, tez_0_xxxx_0, \
                                         tez_0_xxxx_1, tez_0_xxxy_0, tez_0_xxxy_1, tez_0_xxxz_0, tez_0_xxxz_1, tez_0_xxy_0, \
                                         tez_0_xxy_1, tez_0_xxyy_0, tez_0_xxyy_1, tez_0_xxyz_0, tez_0_xxyz_1, tez_0_xxz_0, \
                                         tez_0_xxz_1, tez_0_xxzz_0, tez_0_xxzz_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyyy_0, \
                                         tez_0_xyyy_1, tez_0_xyyz_0, tez_0_xyyz_1, tez_0_xyz_0, tez_0_xyz_1, tez_0_xyzz_0, \
                                         tez_0_xyzz_1, tez_0_xzz_0, tez_0_xzz_1, tez_0_xzzz_0, tez_0_xzzz_1, tez_0_yyy_0, \
                                         tez_0_yyy_1, tez_0_yyyy_0, tez_0_yyyy_1, tez_0_yyyz_0, tez_0_yyyz_1, tez_0_yyz_0, \
                                         tez_0_yyz_1, tez_0_yyzz_0, tez_0_yyzz_1, tez_0_yzz_0, tez_0_yzz_1, tez_0_yzzz_0, \
                                         tez_0_yzzz_1, tez_0_zzz_0, tez_0_zzz_1, tez_0_zzzz_0, tez_0_zzzz_1, tez_z_xxxx_0, \
                                         tez_z_xxxy_0, tez_z_xxxz_0, tez_z_xxyy_0, tez_z_xxyz_0, tez_z_xxzz_0, tez_z_xyyy_0, \
                                         tez_z_xyyz_0, tez_z_xyzz_0, tez_z_xzzz_0, tez_z_yyyy_0, tez_z_yyyz_0, tez_z_yyzz_0, \
                                         tez_z_yzzz_0, tez_z_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_z_xxxx_0[j] = pa_z[j] * tex_0_xxxx_0[j] - pc_z[j] * tex_0_xxxx_1[j];

                tey_z_xxxx_0[j] = pa_z[j] * tey_0_xxxx_0[j] - pc_z[j] * tey_0_xxxx_1[j];

                tez_z_xxxx_0[j] = pa_z[j] * tez_0_xxxx_0[j] - pc_z[j] * tez_0_xxxx_1[j] + ta_0_xxxx_1[j];

                tex_z_xxxy_0[j] = pa_z[j] * tex_0_xxxy_0[j] - pc_z[j] * tex_0_xxxy_1[j];

                tey_z_xxxy_0[j] = pa_z[j] * tey_0_xxxy_0[j] - pc_z[j] * tey_0_xxxy_1[j];

                tez_z_xxxy_0[j] = pa_z[j] * tez_0_xxxy_0[j] - pc_z[j] * tez_0_xxxy_1[j] + ta_0_xxxy_1[j];

                tex_z_xxxz_0[j] =
                    pa_z[j] * tex_0_xxxz_0[j] - pc_z[j] * tex_0_xxxz_1[j] + 0.5 * fl1_fx * tex_0_xxx_0[j] - 0.5 * fl1_fx * tex_0_xxx_1[j];

                tey_z_xxxz_0[j] =
                    pa_z[j] * tey_0_xxxz_0[j] - pc_z[j] * tey_0_xxxz_1[j] + 0.5 * fl1_fx * tey_0_xxx_0[j] - 0.5 * fl1_fx * tey_0_xxx_1[j];

                tez_z_xxxz_0[j] = pa_z[j] * tez_0_xxxz_0[j] - pc_z[j] * tez_0_xxxz_1[j] + 0.5 * fl1_fx * tez_0_xxx_0[j] -
                                  0.5 * fl1_fx * tez_0_xxx_1[j] + ta_0_xxxz_1[j];

                tex_z_xxyy_0[j] = pa_z[j] * tex_0_xxyy_0[j] - pc_z[j] * tex_0_xxyy_1[j];

                tey_z_xxyy_0[j] = pa_z[j] * tey_0_xxyy_0[j] - pc_z[j] * tey_0_xxyy_1[j];

                tez_z_xxyy_0[j] = pa_z[j] * tez_0_xxyy_0[j] - pc_z[j] * tez_0_xxyy_1[j] + ta_0_xxyy_1[j];

                tex_z_xxyz_0[j] =
                    pa_z[j] * tex_0_xxyz_0[j] - pc_z[j] * tex_0_xxyz_1[j] + 0.5 * fl1_fx * tex_0_xxy_0[j] - 0.5 * fl1_fx * tex_0_xxy_1[j];

                tey_z_xxyz_0[j] =
                    pa_z[j] * tey_0_xxyz_0[j] - pc_z[j] * tey_0_xxyz_1[j] + 0.5 * fl1_fx * tey_0_xxy_0[j] - 0.5 * fl1_fx * tey_0_xxy_1[j];

                tez_z_xxyz_0[j] = pa_z[j] * tez_0_xxyz_0[j] - pc_z[j] * tez_0_xxyz_1[j] + 0.5 * fl1_fx * tez_0_xxy_0[j] -
                                  0.5 * fl1_fx * tez_0_xxy_1[j] + ta_0_xxyz_1[j];

                tex_z_xxzz_0[j] = pa_z[j] * tex_0_xxzz_0[j] - pc_z[j] * tex_0_xxzz_1[j] + fl1_fx * tex_0_xxz_0[j] - fl1_fx * tex_0_xxz_1[j];

                tey_z_xxzz_0[j] = pa_z[j] * tey_0_xxzz_0[j] - pc_z[j] * tey_0_xxzz_1[j] + fl1_fx * tey_0_xxz_0[j] - fl1_fx * tey_0_xxz_1[j];

                tez_z_xxzz_0[j] =
                    pa_z[j] * tez_0_xxzz_0[j] - pc_z[j] * tez_0_xxzz_1[j] + fl1_fx * tez_0_xxz_0[j] - fl1_fx * tez_0_xxz_1[j] + ta_0_xxzz_1[j];

                tex_z_xyyy_0[j] = pa_z[j] * tex_0_xyyy_0[j] - pc_z[j] * tex_0_xyyy_1[j];

                tey_z_xyyy_0[j] = pa_z[j] * tey_0_xyyy_0[j] - pc_z[j] * tey_0_xyyy_1[j];

                tez_z_xyyy_0[j] = pa_z[j] * tez_0_xyyy_0[j] - pc_z[j] * tez_0_xyyy_1[j] + ta_0_xyyy_1[j];

                tex_z_xyyz_0[j] =
                    pa_z[j] * tex_0_xyyz_0[j] - pc_z[j] * tex_0_xyyz_1[j] + 0.5 * fl1_fx * tex_0_xyy_0[j] - 0.5 * fl1_fx * tex_0_xyy_1[j];

                tey_z_xyyz_0[j] =
                    pa_z[j] * tey_0_xyyz_0[j] - pc_z[j] * tey_0_xyyz_1[j] + 0.5 * fl1_fx * tey_0_xyy_0[j] - 0.5 * fl1_fx * tey_0_xyy_1[j];

                tez_z_xyyz_0[j] = pa_z[j] * tez_0_xyyz_0[j] - pc_z[j] * tez_0_xyyz_1[j] + 0.5 * fl1_fx * tez_0_xyy_0[j] -
                                  0.5 * fl1_fx * tez_0_xyy_1[j] + ta_0_xyyz_1[j];

                tex_z_xyzz_0[j] = pa_z[j] * tex_0_xyzz_0[j] - pc_z[j] * tex_0_xyzz_1[j] + fl1_fx * tex_0_xyz_0[j] - fl1_fx * tex_0_xyz_1[j];

                tey_z_xyzz_0[j] = pa_z[j] * tey_0_xyzz_0[j] - pc_z[j] * tey_0_xyzz_1[j] + fl1_fx * tey_0_xyz_0[j] - fl1_fx * tey_0_xyz_1[j];

                tez_z_xyzz_0[j] =
                    pa_z[j] * tez_0_xyzz_0[j] - pc_z[j] * tez_0_xyzz_1[j] + fl1_fx * tez_0_xyz_0[j] - fl1_fx * tez_0_xyz_1[j] + ta_0_xyzz_1[j];

                tex_z_xzzz_0[j] =
                    pa_z[j] * tex_0_xzzz_0[j] - pc_z[j] * tex_0_xzzz_1[j] + 1.5 * fl1_fx * tex_0_xzz_0[j] - 1.5 * fl1_fx * tex_0_xzz_1[j];

                tey_z_xzzz_0[j] =
                    pa_z[j] * tey_0_xzzz_0[j] - pc_z[j] * tey_0_xzzz_1[j] + 1.5 * fl1_fx * tey_0_xzz_0[j] - 1.5 * fl1_fx * tey_0_xzz_1[j];

                tez_z_xzzz_0[j] = pa_z[j] * tez_0_xzzz_0[j] - pc_z[j] * tez_0_xzzz_1[j] + 1.5 * fl1_fx * tez_0_xzz_0[j] -
                                  1.5 * fl1_fx * tez_0_xzz_1[j] + ta_0_xzzz_1[j];

                tex_z_yyyy_0[j] = pa_z[j] * tex_0_yyyy_0[j] - pc_z[j] * tex_0_yyyy_1[j];

                tey_z_yyyy_0[j] = pa_z[j] * tey_0_yyyy_0[j] - pc_z[j] * tey_0_yyyy_1[j];

                tez_z_yyyy_0[j] = pa_z[j] * tez_0_yyyy_0[j] - pc_z[j] * tez_0_yyyy_1[j] + ta_0_yyyy_1[j];

                tex_z_yyyz_0[j] =
                    pa_z[j] * tex_0_yyyz_0[j] - pc_z[j] * tex_0_yyyz_1[j] + 0.5 * fl1_fx * tex_0_yyy_0[j] - 0.5 * fl1_fx * tex_0_yyy_1[j];

                tey_z_yyyz_0[j] =
                    pa_z[j] * tey_0_yyyz_0[j] - pc_z[j] * tey_0_yyyz_1[j] + 0.5 * fl1_fx * tey_0_yyy_0[j] - 0.5 * fl1_fx * tey_0_yyy_1[j];

                tez_z_yyyz_0[j] = pa_z[j] * tez_0_yyyz_0[j] - pc_z[j] * tez_0_yyyz_1[j] + 0.5 * fl1_fx * tez_0_yyy_0[j] -
                                  0.5 * fl1_fx * tez_0_yyy_1[j] + ta_0_yyyz_1[j];

                tex_z_yyzz_0[j] = pa_z[j] * tex_0_yyzz_0[j] - pc_z[j] * tex_0_yyzz_1[j] + fl1_fx * tex_0_yyz_0[j] - fl1_fx * tex_0_yyz_1[j];

                tey_z_yyzz_0[j] = pa_z[j] * tey_0_yyzz_0[j] - pc_z[j] * tey_0_yyzz_1[j] + fl1_fx * tey_0_yyz_0[j] - fl1_fx * tey_0_yyz_1[j];

                tez_z_yyzz_0[j] =
                    pa_z[j] * tez_0_yyzz_0[j] - pc_z[j] * tez_0_yyzz_1[j] + fl1_fx * tez_0_yyz_0[j] - fl1_fx * tez_0_yyz_1[j] + ta_0_yyzz_1[j];

                tex_z_yzzz_0[j] =
                    pa_z[j] * tex_0_yzzz_0[j] - pc_z[j] * tex_0_yzzz_1[j] + 1.5 * fl1_fx * tex_0_yzz_0[j] - 1.5 * fl1_fx * tex_0_yzz_1[j];

                tey_z_yzzz_0[j] =
                    pa_z[j] * tey_0_yzzz_0[j] - pc_z[j] * tey_0_yzzz_1[j] + 1.5 * fl1_fx * tey_0_yzz_0[j] - 1.5 * fl1_fx * tey_0_yzz_1[j];

                tez_z_yzzz_0[j] = pa_z[j] * tez_0_yzzz_0[j] - pc_z[j] * tez_0_yzzz_1[j] + 1.5 * fl1_fx * tez_0_yzz_0[j] -
                                  1.5 * fl1_fx * tez_0_yzz_1[j] + ta_0_yzzz_1[j];

                tex_z_zzzz_0[j] =
                    pa_z[j] * tex_0_zzzz_0[j] - pc_z[j] * tex_0_zzzz_1[j] + 2.0 * fl1_fx * tex_0_zzz_0[j] - 2.0 * fl1_fx * tex_0_zzz_1[j];

                tey_z_zzzz_0[j] =
                    pa_z[j] * tey_0_zzzz_0[j] - pc_z[j] * tey_0_zzzz_1[j] + 2.0 * fl1_fx * tey_0_zzz_0[j] - 2.0 * fl1_fx * tey_0_zzz_1[j];

                tez_z_zzzz_0[j] = pa_z[j] * tez_0_zzzz_0[j] - pc_z[j] * tez_0_zzzz_1[j] + 2.0 * fl1_fx * tez_0_zzz_0[j] -
                                  2.0 * fl1_fx * tez_0_zzz_1[j] + ta_0_zzzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGP(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CMemBlock2D<double>& pcDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    efieldrecfunc::compElectricFieldForGP_0_45(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGP_45_90(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGP_90_135(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compElectricFieldForGP_0_45(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx);

            auto tey_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx);

            auto tez_xxx_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx);

            auto tex_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 1);

            auto tey_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 1);

            auto tez_xxx_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 1);

            auto tex_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 2);

            auto tey_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 2);

            auto tez_xxx_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 2);

            auto tex_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 3);

            auto tey_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 3);

            auto tez_xxy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 3);

            auto tex_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 4);

            auto tey_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 4);

            auto tez_xxy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 4);

            auto tex_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 5);

            auto tey_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 5);

            auto tez_xxy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 5);

            auto tex_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 6);

            auto tey_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 6);

            auto tez_xxz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 6);

            auto tex_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 7);

            auto tey_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 7);

            auto tez_xxz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 7);

            auto tex_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 8);

            auto tey_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 8);

            auto tez_xxz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 8);

            auto tex_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 9);

            auto tey_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 9);

            auto tez_xyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 9);

            auto tex_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 10);

            auto tey_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 10);

            auto tez_xyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 10);

            auto tex_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 11);

            auto tey_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 11);

            auto tez_xyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 11);

            auto tex_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 12);

            auto tey_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 12);

            auto tez_xyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 12);

            auto tex_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 13);

            auto tey_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 13);

            auto tez_xyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 13);

            auto tex_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 14);

            auto tey_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 14);

            auto tez_xyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 14);

            auto tex_xxx_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx);

            auto tey_xxx_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx);

            auto tez_xxx_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx);

            auto tex_xxx_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 1);

            auto tey_xxx_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 1);

            auto tez_xxx_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 1);

            auto tex_xxx_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 2);

            auto tey_xxx_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 2);

            auto tez_xxx_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 2);

            auto tex_xxy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 3);

            auto tey_xxy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 3);

            auto tez_xxy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 3);

            auto tex_xxy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 4);

            auto tey_xxy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 4);

            auto tez_xxy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 4);

            auto tex_xxy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 5);

            auto tey_xxy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 5);

            auto tez_xxy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 5);

            auto tex_xxz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 6);

            auto tey_xxz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 6);

            auto tez_xxz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 6);

            auto tex_xxz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 7);

            auto tey_xxz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 7);

            auto tez_xxz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 7);

            auto tex_xxz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 8);

            auto tey_xxz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 8);

            auto tez_xxz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 8);

            auto tex_xyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 9);

            auto tey_xyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 9);

            auto tez_xyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 9);

            auto tex_xyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 10);

            auto tey_xyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 10);

            auto tez_xyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 10);

            auto tex_xyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 11);

            auto tey_xyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 11);

            auto tez_xyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 11);

            auto tex_xyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 12);

            auto tey_xyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 12);

            auto tez_xyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 12);

            auto tex_xyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 13);

            auto tey_xyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 13);

            auto tez_xyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 13);

            auto tex_xyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 14);

            auto tey_xyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 14);

            auto tez_xyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 14);

            auto tex_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx);

            auto tey_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx);

            auto tez_xx_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx);

            auto tex_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 1);

            auto tey_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 1);

            auto tez_xx_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 1);

            auto tex_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 2);

            auto tey_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 2);

            auto tez_xx_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 2);

            auto tex_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 3);

            auto tey_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 3);

            auto tez_xy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 3);

            auto tex_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 4);

            auto tey_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 4);

            auto tez_xy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 4);

            auto tex_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 5);

            auto tey_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 5);

            auto tez_xy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 5);

            auto tex_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 6);

            auto tey_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 6);

            auto tez_xz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 6);

            auto tex_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 7);

            auto tey_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 7);

            auto tez_xz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 7);

            auto tex_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 8);

            auto tey_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 8);

            auto tez_xz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 8);

            auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9);

            auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10);

            auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11);

            auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12);

            auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13);

            auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14);

            auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14);

            auto tex_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx);

            auto tey_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx);

            auto tez_xx_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx);

            auto tex_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 1);

            auto tey_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 1);

            auto tez_xx_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 1);

            auto tex_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 2);

            auto tey_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 2);

            auto tez_xx_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 2);

            auto tex_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 3);

            auto tey_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 3);

            auto tez_xy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 3);

            auto tex_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 4);

            auto tey_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 4);

            auto tez_xy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 4);

            auto tex_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 5);

            auto tey_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 5);

            auto tez_xy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 5);

            auto tex_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 6);

            auto tey_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 6);

            auto tez_xz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 6);

            auto tex_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 7);

            auto tey_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 7);

            auto tez_xz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 7);

            auto tex_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 8);

            auto tey_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 8);

            auto tez_xz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 8);

            auto tex_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 9);

            auto tey_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 10);

            auto tey_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 11);

            auto tey_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 12);

            auto tey_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 13);

            auto tey_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 14);

            auto tey_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 14);

            auto tex_xxx_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx);

            auto tey_xxx_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx);

            auto tez_xxx_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx);

            auto tex_xxy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 1);

            auto tey_xxy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 1);

            auto tez_xxy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 1);

            auto tex_xxz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 2);

            auto tey_xxz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 2);

            auto tez_xxz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 2);

            auto tex_xyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 3);

            auto tey_xyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 3);

            auto tez_xyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 3);

            auto tex_xyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 4);

            auto tey_xyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 4);

            auto tez_xyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 4);

            auto tex_xxx_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx);

            auto tey_xxx_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx);

            auto tez_xxx_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx);

            auto tex_xxy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 1);

            auto tey_xxy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 1);

            auto tez_xxy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 1);

            auto tex_xxz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 2);

            auto tey_xxz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 2);

            auto tez_xxz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 2);

            auto tex_xyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 3);

            auto tey_xyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 3);

            auto tez_xyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 3);

            auto tex_xyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 4);

            auto tey_xyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 4);

            auto tez_xyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 4);

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

            // set up pointers to integrals

            auto tex_xxxx_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx);

            auto tey_xxxx_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx);

            auto tez_xxxx_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx);

            auto tex_xxxx_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 1);

            auto tey_xxxx_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 1);

            auto tez_xxxx_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 1);

            auto tex_xxxx_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 2);

            auto tey_xxxx_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 2);

            auto tez_xxxx_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 2);

            auto tex_xxxy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 3);

            auto tey_xxxy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 3);

            auto tez_xxxy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 3);

            auto tex_xxxy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 4);

            auto tey_xxxy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 4);

            auto tez_xxxy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 4);

            auto tex_xxxy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 5);

            auto tey_xxxy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 5);

            auto tez_xxxy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 5);

            auto tex_xxxz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 6);

            auto tey_xxxz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 6);

            auto tez_xxxz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 6);

            auto tex_xxxz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 7);

            auto tey_xxxz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 7);

            auto tez_xxxz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 7);

            auto tex_xxxz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 8);

            auto tey_xxxz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 8);

            auto tez_xxxz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 8);

            auto tex_xxyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 9);

            auto tey_xxyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 9);

            auto tez_xxyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 9);

            auto tex_xxyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 10);

            auto tey_xxyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 10);

            auto tez_xxyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 10);

            auto tex_xxyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 11);

            auto tey_xxyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 11);

            auto tez_xxyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 11);

            auto tex_xxyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 12);

            auto tey_xxyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 12);

            auto tez_xxyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 12);

            auto tex_xxyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 13);

            auto tey_xxyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 13);

            auto tez_xxyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 13);

            auto tex_xxyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 14);

            auto tey_xxyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 14);

            auto tez_xxyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 14);

            // Batch of Integrals (0,45)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xxx_x_1, ta_xxx_y_1, ta_xxx_z_1, ta_xxy_x_1, ta_xxy_y_1, \
                                         ta_xxy_z_1, ta_xxz_x_1, ta_xxz_y_1, ta_xxz_z_1, ta_xyy_x_1, ta_xyy_y_1, ta_xyy_z_1, \
                                         ta_xyz_x_1, ta_xyz_y_1, ta_xyz_z_1, tex_xx_x_0, tex_xx_x_1, tex_xx_y_0, tex_xx_y_1, \
                                         tex_xx_z_0, tex_xx_z_1, tex_xxx_0_0, tex_xxx_0_1, tex_xxx_x_0, tex_xxx_x_1, \
                                         tex_xxx_y_0, tex_xxx_y_1, tex_xxx_z_0, tex_xxx_z_1, tex_xxxx_x_0, tex_xxxx_y_0, \
                                         tex_xxxx_z_0, tex_xxxy_x_0, tex_xxxy_y_0, tex_xxxy_z_0, tex_xxxz_x_0, tex_xxxz_y_0, \
                                         tex_xxxz_z_0, tex_xxy_0_0, tex_xxy_0_1, tex_xxy_x_0, tex_xxy_x_1, tex_xxy_y_0, \
                                         tex_xxy_y_1, tex_xxy_z_0, tex_xxy_z_1, tex_xxyy_x_0, tex_xxyy_y_0, tex_xxyy_z_0, \
                                         tex_xxyz_x_0, tex_xxyz_y_0, tex_xxyz_z_0, tex_xxz_0_0, tex_xxz_0_1, tex_xxz_x_0, \
                                         tex_xxz_x_1, tex_xxz_y_0, tex_xxz_y_1, tex_xxz_z_0, tex_xxz_z_1, tex_xy_x_0, \
                                         tex_xy_x_1, tex_xy_y_0, tex_xy_y_1, tex_xy_z_0, tex_xy_z_1, tex_xyy_0_0, \
                                         tex_xyy_0_1, tex_xyy_x_0, tex_xyy_x_1, tex_xyy_y_0, tex_xyy_y_1, tex_xyy_z_0, \
                                         tex_xyy_z_1, tex_xyz_0_0, tex_xyz_0_1, tex_xyz_x_0, tex_xyz_x_1, tex_xyz_y_0, \
                                         tex_xyz_y_1, tex_xyz_z_0, tex_xyz_z_1, tex_xz_x_0, tex_xz_x_1, tex_xz_y_0, \
                                         tex_xz_y_1, tex_xz_z_0, tex_xz_z_1, tex_yy_x_0, tex_yy_x_1, tex_yy_y_0, tex_yy_y_1, \
                                         tex_yy_z_0, tex_yy_z_1, tex_yz_x_0, tex_yz_x_1, tex_yz_y_0, tex_yz_y_1, tex_yz_z_0, \
                                         tex_yz_z_1, tey_xx_x_0, tey_xx_x_1, tey_xx_y_0, tey_xx_y_1, tey_xx_z_0, tey_xx_z_1, \
                                         tey_xxx_0_0, tey_xxx_0_1, tey_xxx_x_0, tey_xxx_x_1, tey_xxx_y_0, tey_xxx_y_1, \
                                         tey_xxx_z_0, tey_xxx_z_1, tey_xxxx_x_0, tey_xxxx_y_0, tey_xxxx_z_0, tey_xxxy_x_0, \
                                         tey_xxxy_y_0, tey_xxxy_z_0, tey_xxxz_x_0, tey_xxxz_y_0, tey_xxxz_z_0, tey_xxy_0_0, \
                                         tey_xxy_0_1, tey_xxy_x_0, tey_xxy_x_1, tey_xxy_y_0, tey_xxy_y_1, tey_xxy_z_0, \
                                         tey_xxy_z_1, tey_xxyy_x_0, tey_xxyy_y_0, tey_xxyy_z_0, tey_xxyz_x_0, tey_xxyz_y_0, \
                                         tey_xxyz_z_0, tey_xxz_0_0, tey_xxz_0_1, tey_xxz_x_0, tey_xxz_x_1, tey_xxz_y_0, \
                                         tey_xxz_y_1, tey_xxz_z_0, tey_xxz_z_1, tey_xy_x_0, tey_xy_x_1, tey_xy_y_0, \
                                         tey_xy_y_1, tey_xy_z_0, tey_xy_z_1, tey_xyy_0_0, tey_xyy_0_1, tey_xyy_x_0, \
                                         tey_xyy_x_1, tey_xyy_y_0, tey_xyy_y_1, tey_xyy_z_0, tey_xyy_z_1, tey_xyz_0_0, \
                                         tey_xyz_0_1, tey_xyz_x_0, tey_xyz_x_1, tey_xyz_y_0, tey_xyz_y_1, tey_xyz_z_0, \
                                         tey_xyz_z_1, tey_xz_x_0, tey_xz_x_1, tey_xz_y_0, tey_xz_y_1, tey_xz_z_0, tey_xz_z_1, \
                                         tey_yy_x_0, tey_yy_x_1, tey_yy_y_0, tey_yy_y_1, tey_yy_z_0, tey_yy_z_1, tey_yz_x_0, \
                                         tey_yz_x_1, tey_yz_y_0, tey_yz_y_1, tey_yz_z_0, tey_yz_z_1, tez_xx_x_0, tez_xx_x_1, \
                                         tez_xx_y_0, tez_xx_y_1, tez_xx_z_0, tez_xx_z_1, tez_xxx_0_0, tez_xxx_0_1, \
                                         tez_xxx_x_0, tez_xxx_x_1, tez_xxx_y_0, tez_xxx_y_1, tez_xxx_z_0, tez_xxx_z_1, \
                                         tez_xxxx_x_0, tez_xxxx_y_0, tez_xxxx_z_0, tez_xxxy_x_0, tez_xxxy_y_0, tez_xxxy_z_0, \
                                         tez_xxxz_x_0, tez_xxxz_y_0, tez_xxxz_z_0, tez_xxy_0_0, tez_xxy_0_1, tez_xxy_x_0, \
                                         tez_xxy_x_1, tez_xxy_y_0, tez_xxy_y_1, tez_xxy_z_0, tez_xxy_z_1, tez_xxyy_x_0, \
                                         tez_xxyy_y_0, tez_xxyy_z_0, tez_xxyz_x_0, tez_xxyz_y_0, tez_xxyz_z_0, tez_xxz_0_0, \
                                         tez_xxz_0_1, tez_xxz_x_0, tez_xxz_x_1, tez_xxz_y_0, tez_xxz_y_1, tez_xxz_z_0, \
                                         tez_xxz_z_1, tez_xy_x_0, tez_xy_x_1, tez_xy_y_0, tez_xy_y_1, tez_xy_z_0, tez_xy_z_1, \
                                         tez_xyy_0_0, tez_xyy_0_1, tez_xyy_x_0, tez_xyy_x_1, tez_xyy_y_0, tez_xyy_y_1, \
                                         tez_xyy_z_0, tez_xyy_z_1, tez_xyz_0_0, tez_xyz_0_1, tez_xyz_x_0, tez_xyz_x_1, \
                                         tez_xyz_y_0, tez_xyz_y_1, tez_xyz_z_0, tez_xyz_z_1, tez_xz_x_0, tez_xz_x_1, \
                                         tez_xz_y_0, tez_xz_y_1, tez_xz_z_0, tez_xz_z_1, tez_yy_x_0, tez_yy_x_1, tez_yy_y_0, \
                                         tez_yy_y_1, tez_yy_z_0, tez_yy_z_1, tez_yz_x_0, tez_yz_x_1, tez_yz_y_0, tez_yz_y_1, \
                                         tez_yz_z_0, tez_yz_z_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxxx_x_0[j] = pa_x[j] * tex_xxx_x_0[j] - pc_x[j] * tex_xxx_x_1[j] + 1.5 * fl1_fx * tex_xx_x_0[j] - 1.5 * fl1_fx * tex_xx_x_1[j] +
                                  0.5 * fl1_fx * tex_xxx_0_0[j] - 0.5 * fl1_fx * tex_xxx_0_1[j] + ta_xxx_x_1[j];

                tey_xxxx_x_0[j] = pa_x[j] * tey_xxx_x_0[j] - pc_x[j] * tey_xxx_x_1[j] + 1.5 * fl1_fx * tey_xx_x_0[j] - 1.5 * fl1_fx * tey_xx_x_1[j] +
                                  0.5 * fl1_fx * tey_xxx_0_0[j] - 0.5 * fl1_fx * tey_xxx_0_1[j];

                tez_xxxx_x_0[j] = pa_x[j] * tez_xxx_x_0[j] - pc_x[j] * tez_xxx_x_1[j] + 1.5 * fl1_fx * tez_xx_x_0[j] - 1.5 * fl1_fx * tez_xx_x_1[j] +
                                  0.5 * fl1_fx * tez_xxx_0_0[j] - 0.5 * fl1_fx * tez_xxx_0_1[j];

                tex_xxxx_y_0[j] =
                    pa_x[j] * tex_xxx_y_0[j] - pc_x[j] * tex_xxx_y_1[j] + 1.5 * fl1_fx * tex_xx_y_0[j] - 1.5 * fl1_fx * tex_xx_y_1[j] + ta_xxx_y_1[j];

                tey_xxxx_y_0[j] = pa_x[j] * tey_xxx_y_0[j] - pc_x[j] * tey_xxx_y_1[j] + 1.5 * fl1_fx * tey_xx_y_0[j] - 1.5 * fl1_fx * tey_xx_y_1[j];

                tez_xxxx_y_0[j] = pa_x[j] * tez_xxx_y_0[j] - pc_x[j] * tez_xxx_y_1[j] + 1.5 * fl1_fx * tez_xx_y_0[j] - 1.5 * fl1_fx * tez_xx_y_1[j];

                tex_xxxx_z_0[j] =
                    pa_x[j] * tex_xxx_z_0[j] - pc_x[j] * tex_xxx_z_1[j] + 1.5 * fl1_fx * tex_xx_z_0[j] - 1.5 * fl1_fx * tex_xx_z_1[j] + ta_xxx_z_1[j];

                tey_xxxx_z_0[j] = pa_x[j] * tey_xxx_z_0[j] - pc_x[j] * tey_xxx_z_1[j] + 1.5 * fl1_fx * tey_xx_z_0[j] - 1.5 * fl1_fx * tey_xx_z_1[j];

                tez_xxxx_z_0[j] = pa_x[j] * tez_xxx_z_0[j] - pc_x[j] * tez_xxx_z_1[j] + 1.5 * fl1_fx * tez_xx_z_0[j] - 1.5 * fl1_fx * tez_xx_z_1[j];

                tex_xxxy_x_0[j] = pa_x[j] * tex_xxy_x_0[j] - pc_x[j] * tex_xxy_x_1[j] + fl1_fx * tex_xy_x_0[j] - fl1_fx * tex_xy_x_1[j] +
                                  0.5 * fl1_fx * tex_xxy_0_0[j] - 0.5 * fl1_fx * tex_xxy_0_1[j] + ta_xxy_x_1[j];

                tey_xxxy_x_0[j] = pa_x[j] * tey_xxy_x_0[j] - pc_x[j] * tey_xxy_x_1[j] + fl1_fx * tey_xy_x_0[j] - fl1_fx * tey_xy_x_1[j] +
                                  0.5 * fl1_fx * tey_xxy_0_0[j] - 0.5 * fl1_fx * tey_xxy_0_1[j];

                tez_xxxy_x_0[j] = pa_x[j] * tez_xxy_x_0[j] - pc_x[j] * tez_xxy_x_1[j] + fl1_fx * tez_xy_x_0[j] - fl1_fx * tez_xy_x_1[j] +
                                  0.5 * fl1_fx * tez_xxy_0_0[j] - 0.5 * fl1_fx * tez_xxy_0_1[j];

                tex_xxxy_y_0[j] =
                    pa_x[j] * tex_xxy_y_0[j] - pc_x[j] * tex_xxy_y_1[j] + fl1_fx * tex_xy_y_0[j] - fl1_fx * tex_xy_y_1[j] + ta_xxy_y_1[j];

                tey_xxxy_y_0[j] = pa_x[j] * tey_xxy_y_0[j] - pc_x[j] * tey_xxy_y_1[j] + fl1_fx * tey_xy_y_0[j] - fl1_fx * tey_xy_y_1[j];

                tez_xxxy_y_0[j] = pa_x[j] * tez_xxy_y_0[j] - pc_x[j] * tez_xxy_y_1[j] + fl1_fx * tez_xy_y_0[j] - fl1_fx * tez_xy_y_1[j];

                tex_xxxy_z_0[j] =
                    pa_x[j] * tex_xxy_z_0[j] - pc_x[j] * tex_xxy_z_1[j] + fl1_fx * tex_xy_z_0[j] - fl1_fx * tex_xy_z_1[j] + ta_xxy_z_1[j];

                tey_xxxy_z_0[j] = pa_x[j] * tey_xxy_z_0[j] - pc_x[j] * tey_xxy_z_1[j] + fl1_fx * tey_xy_z_0[j] - fl1_fx * tey_xy_z_1[j];

                tez_xxxy_z_0[j] = pa_x[j] * tez_xxy_z_0[j] - pc_x[j] * tez_xxy_z_1[j] + fl1_fx * tez_xy_z_0[j] - fl1_fx * tez_xy_z_1[j];

                tex_xxxz_x_0[j] = pa_x[j] * tex_xxz_x_0[j] - pc_x[j] * tex_xxz_x_1[j] + fl1_fx * tex_xz_x_0[j] - fl1_fx * tex_xz_x_1[j] +
                                  0.5 * fl1_fx * tex_xxz_0_0[j] - 0.5 * fl1_fx * tex_xxz_0_1[j] + ta_xxz_x_1[j];

                tey_xxxz_x_0[j] = pa_x[j] * tey_xxz_x_0[j] - pc_x[j] * tey_xxz_x_1[j] + fl1_fx * tey_xz_x_0[j] - fl1_fx * tey_xz_x_1[j] +
                                  0.5 * fl1_fx * tey_xxz_0_0[j] - 0.5 * fl1_fx * tey_xxz_0_1[j];

                tez_xxxz_x_0[j] = pa_x[j] * tez_xxz_x_0[j] - pc_x[j] * tez_xxz_x_1[j] + fl1_fx * tez_xz_x_0[j] - fl1_fx * tez_xz_x_1[j] +
                                  0.5 * fl1_fx * tez_xxz_0_0[j] - 0.5 * fl1_fx * tez_xxz_0_1[j];

                tex_xxxz_y_0[j] =
                    pa_x[j] * tex_xxz_y_0[j] - pc_x[j] * tex_xxz_y_1[j] + fl1_fx * tex_xz_y_0[j] - fl1_fx * tex_xz_y_1[j] + ta_xxz_y_1[j];

                tey_xxxz_y_0[j] = pa_x[j] * tey_xxz_y_0[j] - pc_x[j] * tey_xxz_y_1[j] + fl1_fx * tey_xz_y_0[j] - fl1_fx * tey_xz_y_1[j];

                tez_xxxz_y_0[j] = pa_x[j] * tez_xxz_y_0[j] - pc_x[j] * tez_xxz_y_1[j] + fl1_fx * tez_xz_y_0[j] - fl1_fx * tez_xz_y_1[j];

                tex_xxxz_z_0[j] =
                    pa_x[j] * tex_xxz_z_0[j] - pc_x[j] * tex_xxz_z_1[j] + fl1_fx * tex_xz_z_0[j] - fl1_fx * tex_xz_z_1[j] + ta_xxz_z_1[j];

                tey_xxxz_z_0[j] = pa_x[j] * tey_xxz_z_0[j] - pc_x[j] * tey_xxz_z_1[j] + fl1_fx * tey_xz_z_0[j] - fl1_fx * tey_xz_z_1[j];

                tez_xxxz_z_0[j] = pa_x[j] * tez_xxz_z_0[j] - pc_x[j] * tez_xxz_z_1[j] + fl1_fx * tez_xz_z_0[j] - fl1_fx * tez_xz_z_1[j];

                tex_xxyy_x_0[j] = pa_x[j] * tex_xyy_x_0[j] - pc_x[j] * tex_xyy_x_1[j] + 0.5 * fl1_fx * tex_yy_x_0[j] - 0.5 * fl1_fx * tex_yy_x_1[j] +
                                  0.5 * fl1_fx * tex_xyy_0_0[j] - 0.5 * fl1_fx * tex_xyy_0_1[j] + ta_xyy_x_1[j];

                tey_xxyy_x_0[j] = pa_x[j] * tey_xyy_x_0[j] - pc_x[j] * tey_xyy_x_1[j] + 0.5 * fl1_fx * tey_yy_x_0[j] - 0.5 * fl1_fx * tey_yy_x_1[j] +
                                  0.5 * fl1_fx * tey_xyy_0_0[j] - 0.5 * fl1_fx * tey_xyy_0_1[j];

                tez_xxyy_x_0[j] = pa_x[j] * tez_xyy_x_0[j] - pc_x[j] * tez_xyy_x_1[j] + 0.5 * fl1_fx * tez_yy_x_0[j] - 0.5 * fl1_fx * tez_yy_x_1[j] +
                                  0.5 * fl1_fx * tez_xyy_0_0[j] - 0.5 * fl1_fx * tez_xyy_0_1[j];

                tex_xxyy_y_0[j] =
                    pa_x[j] * tex_xyy_y_0[j] - pc_x[j] * tex_xyy_y_1[j] + 0.5 * fl1_fx * tex_yy_y_0[j] - 0.5 * fl1_fx * tex_yy_y_1[j] + ta_xyy_y_1[j];

                tey_xxyy_y_0[j] = pa_x[j] * tey_xyy_y_0[j] - pc_x[j] * tey_xyy_y_1[j] + 0.5 * fl1_fx * tey_yy_y_0[j] - 0.5 * fl1_fx * tey_yy_y_1[j];

                tez_xxyy_y_0[j] = pa_x[j] * tez_xyy_y_0[j] - pc_x[j] * tez_xyy_y_1[j] + 0.5 * fl1_fx * tez_yy_y_0[j] - 0.5 * fl1_fx * tez_yy_y_1[j];

                tex_xxyy_z_0[j] =
                    pa_x[j] * tex_xyy_z_0[j] - pc_x[j] * tex_xyy_z_1[j] + 0.5 * fl1_fx * tex_yy_z_0[j] - 0.5 * fl1_fx * tex_yy_z_1[j] + ta_xyy_z_1[j];

                tey_xxyy_z_0[j] = pa_x[j] * tey_xyy_z_0[j] - pc_x[j] * tey_xyy_z_1[j] + 0.5 * fl1_fx * tey_yy_z_0[j] - 0.5 * fl1_fx * tey_yy_z_1[j];

                tez_xxyy_z_0[j] = pa_x[j] * tez_xyy_z_0[j] - pc_x[j] * tez_xyy_z_1[j] + 0.5 * fl1_fx * tez_yy_z_0[j] - 0.5 * fl1_fx * tez_yy_z_1[j];

                tex_xxyz_x_0[j] = pa_x[j] * tex_xyz_x_0[j] - pc_x[j] * tex_xyz_x_1[j] + 0.5 * fl1_fx * tex_yz_x_0[j] - 0.5 * fl1_fx * tex_yz_x_1[j] +
                                  0.5 * fl1_fx * tex_xyz_0_0[j] - 0.5 * fl1_fx * tex_xyz_0_1[j] + ta_xyz_x_1[j];

                tey_xxyz_x_0[j] = pa_x[j] * tey_xyz_x_0[j] - pc_x[j] * tey_xyz_x_1[j] + 0.5 * fl1_fx * tey_yz_x_0[j] - 0.5 * fl1_fx * tey_yz_x_1[j] +
                                  0.5 * fl1_fx * tey_xyz_0_0[j] - 0.5 * fl1_fx * tey_xyz_0_1[j];

                tez_xxyz_x_0[j] = pa_x[j] * tez_xyz_x_0[j] - pc_x[j] * tez_xyz_x_1[j] + 0.5 * fl1_fx * tez_yz_x_0[j] - 0.5 * fl1_fx * tez_yz_x_1[j] +
                                  0.5 * fl1_fx * tez_xyz_0_0[j] - 0.5 * fl1_fx * tez_xyz_0_1[j];

                tex_xxyz_y_0[j] =
                    pa_x[j] * tex_xyz_y_0[j] - pc_x[j] * tex_xyz_y_1[j] + 0.5 * fl1_fx * tex_yz_y_0[j] - 0.5 * fl1_fx * tex_yz_y_1[j] + ta_xyz_y_1[j];

                tey_xxyz_y_0[j] = pa_x[j] * tey_xyz_y_0[j] - pc_x[j] * tey_xyz_y_1[j] + 0.5 * fl1_fx * tey_yz_y_0[j] - 0.5 * fl1_fx * tey_yz_y_1[j];

                tez_xxyz_y_0[j] = pa_x[j] * tez_xyz_y_0[j] - pc_x[j] * tez_xyz_y_1[j] + 0.5 * fl1_fx * tez_yz_y_0[j] - 0.5 * fl1_fx * tez_yz_y_1[j];

                tex_xxyz_z_0[j] =
                    pa_x[j] * tex_xyz_z_0[j] - pc_x[j] * tex_xyz_z_1[j] + 0.5 * fl1_fx * tex_yz_z_0[j] - 0.5 * fl1_fx * tex_yz_z_1[j] + ta_xyz_z_1[j];

                tey_xxyz_z_0[j] = pa_x[j] * tey_xyz_z_0[j] - pc_x[j] * tey_xyz_z_1[j] + 0.5 * fl1_fx * tey_yz_z_0[j] - 0.5 * fl1_fx * tey_yz_z_1[j];

                tez_xxyz_z_0[j] = pa_x[j] * tez_xyz_z_0[j] - pc_x[j] * tez_xyz_z_1[j] + 0.5 * fl1_fx * tez_yz_z_0[j] - 0.5 * fl1_fx * tez_yz_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGP_45_90(CMemBlock2D<double>&       primBuffer,
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

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 15);

            auto tey_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 15);

            auto tez_xzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 15);

            auto tex_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 16);

            auto tey_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 16);

            auto tez_xzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 16);

            auto tex_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 17);

            auto tey_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 17);

            auto tez_xzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 17);

            auto tex_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 18);

            auto tey_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 18);

            auto tez_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 18);

            auto tex_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 19);

            auto tey_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 19);

            auto tez_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 19);

            auto tex_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 20);

            auto tey_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 20);

            auto tez_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 20);

            auto tex_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 21);

            auto tey_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 21);

            auto tez_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 21);

            auto tex_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 22);

            auto tey_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 22);

            auto tez_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 22);

            auto tex_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 23);

            auto tey_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 23);

            auto tez_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 23);

            auto tex_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 24);

            auto tey_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 24);

            auto tez_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 24);

            auto tex_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 25);

            auto tey_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 25);

            auto tez_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 25);

            auto tex_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 26);

            auto tey_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 26);

            auto tez_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 26);

            auto tex_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 27);

            auto tey_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 27);

            auto tez_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 27);

            auto tex_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 28);

            auto tey_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 28);

            auto tez_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 28);

            auto tex_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 29);

            auto tey_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 29);

            auto tez_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 29);

            auto tex_xzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 15);

            auto tey_xzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 15);

            auto tez_xzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 15);

            auto tex_xzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 16);

            auto tey_xzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 16);

            auto tez_xzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 16);

            auto tex_xzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 17);

            auto tey_xzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 17);

            auto tez_xzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 17);

            auto tex_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 18);

            auto tey_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 18);

            auto tez_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 18);

            auto tex_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 19);

            auto tey_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 19);

            auto tez_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 19);

            auto tex_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 20);

            auto tey_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 20);

            auto tez_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 20);

            auto tex_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 21);

            auto tey_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 21);

            auto tez_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 21);

            auto tex_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 22);

            auto tey_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 22);

            auto tez_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 22);

            auto tex_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 23);

            auto tey_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 23);

            auto tez_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 23);

            auto tex_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 24);

            auto tey_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 24);

            auto tez_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 24);

            auto tex_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 25);

            auto tey_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 25);

            auto tez_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 25);

            auto tex_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 26);

            auto tey_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 26);

            auto tez_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 26);

            auto tex_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 27);

            auto tey_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 27);

            auto tez_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 27);

            auto tex_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 28);

            auto tey_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 28);

            auto tez_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 28);

            auto tex_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 29);

            auto tey_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 29);

            auto tez_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 29);

            auto tex_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 15);

            auto tey_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 16);

            auto tey_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 17);

            auto tey_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 17);

            auto tex_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 15);

            auto tey_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 16);

            auto tey_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 17);

            auto tey_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 17);

            auto tex_xzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 5);

            auto tey_xzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 5);

            auto tez_xzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 5);

            auto tex_yyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 6);

            auto tey_yyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_yyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_yyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 7);

            auto tey_yyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_yyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_yzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 8);

            auto tey_yzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_yzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_zzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 9);

            auto tey_zzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_zzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_xzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 5);

            auto tey_xzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 5);

            auto tez_xzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 5);

            auto tex_yyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 6);

            auto tey_yyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_yyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_yyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 7);

            auto tey_yyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_yyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_yzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 8);

            auto tey_yzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_yzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_zzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 9);

            auto tey_zzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_zzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 9);

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

            // set up pointers to integrals

            auto tex_xxzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 15);

            auto tey_xxzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 15);

            auto tez_xxzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 15);

            auto tex_xxzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 16);

            auto tey_xxzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 16);

            auto tez_xxzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 16);

            auto tex_xxzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 17);

            auto tey_xxzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 17);

            auto tez_xxzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 17);

            auto tex_xyyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 18);

            auto tey_xyyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 18);

            auto tez_xyyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 18);

            auto tex_xyyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 19);

            auto tey_xyyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 19);

            auto tez_xyyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 19);

            auto tex_xyyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 20);

            auto tey_xyyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 20);

            auto tez_xyyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 20);

            auto tex_xyyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 21);

            auto tey_xyyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 21);

            auto tez_xyyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 21);

            auto tex_xyyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 22);

            auto tey_xyyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 22);

            auto tez_xyyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 22);

            auto tex_xyyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 23);

            auto tey_xyyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 23);

            auto tez_xyyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 23);

            auto tex_xyzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 24);

            auto tey_xyzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 24);

            auto tez_xyzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 24);

            auto tex_xyzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 25);

            auto tey_xyzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 25);

            auto tez_xyzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 25);

            auto tex_xyzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 26);

            auto tey_xyzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 26);

            auto tez_xyzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 26);

            auto tex_xzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 27);

            auto tey_xzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 27);

            auto tez_xzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 27);

            auto tex_xzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 28);

            auto tey_xzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 28);

            auto tez_xzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 28);

            auto tex_xzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 29);

            auto tey_xzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 29);

            auto tez_xzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 29);

            // Batch of Integrals (45,90)

            #pragma omp simd aligned(fx, pa_x, pc_x, ta_xzz_x_1, ta_xzz_y_1, ta_xzz_z_1, ta_yyy_x_1, ta_yyy_y_1, \
                                         ta_yyy_z_1, ta_yyz_x_1, ta_yyz_y_1, ta_yyz_z_1, ta_yzz_x_1, ta_yzz_y_1, ta_yzz_z_1, \
                                         ta_zzz_x_1, ta_zzz_y_1, ta_zzz_z_1, tex_xxzz_x_0, tex_xxzz_y_0, tex_xxzz_z_0, \
                                         tex_xyyy_x_0, tex_xyyy_y_0, tex_xyyy_z_0, tex_xyyz_x_0, tex_xyyz_y_0, tex_xyyz_z_0, \
                                         tex_xyzz_x_0, tex_xyzz_y_0, tex_xyzz_z_0, tex_xzz_0_0, tex_xzz_0_1, tex_xzz_x_0, \
                                         tex_xzz_x_1, tex_xzz_y_0, tex_xzz_y_1, tex_xzz_z_0, tex_xzz_z_1, tex_xzzz_x_0, \
                                         tex_xzzz_y_0, tex_xzzz_z_0, tex_yyy_0_0, tex_yyy_0_1, tex_yyy_x_0, tex_yyy_x_1, \
                                         tex_yyy_y_0, tex_yyy_y_1, tex_yyy_z_0, tex_yyy_z_1, tex_yyz_0_0, tex_yyz_0_1, \
                                         tex_yyz_x_0, tex_yyz_x_1, tex_yyz_y_0, tex_yyz_y_1, tex_yyz_z_0, tex_yyz_z_1, \
                                         tex_yzz_0_0, tex_yzz_0_1, tex_yzz_x_0, tex_yzz_x_1, tex_yzz_y_0, tex_yzz_y_1, \
                                         tex_yzz_z_0, tex_yzz_z_1, tex_zz_x_0, tex_zz_x_1, tex_zz_y_0, tex_zz_y_1, \
                                         tex_zz_z_0, tex_zz_z_1, tex_zzz_0_0, tex_zzz_0_1, tex_zzz_x_0, tex_zzz_x_1, \
                                         tex_zzz_y_0, tex_zzz_y_1, tex_zzz_z_0, tex_zzz_z_1, tey_xxzz_x_0, tey_xxzz_y_0, \
                                         tey_xxzz_z_0, tey_xyyy_x_0, tey_xyyy_y_0, tey_xyyy_z_0, tey_xyyz_x_0, tey_xyyz_y_0, \
                                         tey_xyyz_z_0, tey_xyzz_x_0, tey_xyzz_y_0, tey_xyzz_z_0, tey_xzz_0_0, tey_xzz_0_1, \
                                         tey_xzz_x_0, tey_xzz_x_1, tey_xzz_y_0, tey_xzz_y_1, tey_xzz_z_0, tey_xzz_z_1, \
                                         tey_xzzz_x_0, tey_xzzz_y_0, tey_xzzz_z_0, tey_yyy_0_0, tey_yyy_0_1, tey_yyy_x_0, \
                                         tey_yyy_x_1, tey_yyy_y_0, tey_yyy_y_1, tey_yyy_z_0, tey_yyy_z_1, tey_yyz_0_0, \
                                         tey_yyz_0_1, tey_yyz_x_0, tey_yyz_x_1, tey_yyz_y_0, tey_yyz_y_1, tey_yyz_z_0, \
                                         tey_yyz_z_1, tey_yzz_0_0, tey_yzz_0_1, tey_yzz_x_0, tey_yzz_x_1, tey_yzz_y_0, \
                                         tey_yzz_y_1, tey_yzz_z_0, tey_yzz_z_1, tey_zz_x_0, tey_zz_x_1, tey_zz_y_0, \
                                         tey_zz_y_1, tey_zz_z_0, tey_zz_z_1, tey_zzz_0_0, tey_zzz_0_1, tey_zzz_x_0, \
                                         tey_zzz_x_1, tey_zzz_y_0, tey_zzz_y_1, tey_zzz_z_0, tey_zzz_z_1, tez_xxzz_x_0, \
                                         tez_xxzz_y_0, tez_xxzz_z_0, tez_xyyy_x_0, tez_xyyy_y_0, tez_xyyy_z_0, tez_xyyz_x_0, \
                                         tez_xyyz_y_0, tez_xyyz_z_0, tez_xyzz_x_0, tez_xyzz_y_0, tez_xyzz_z_0, tez_xzz_0_0, \
                                         tez_xzz_0_1, tez_xzz_x_0, tez_xzz_x_1, tez_xzz_y_0, tez_xzz_y_1, tez_xzz_z_0, \
                                         tez_xzz_z_1, tez_xzzz_x_0, tez_xzzz_y_0, tez_xzzz_z_0, tez_yyy_0_0, tez_yyy_0_1, \
                                         tez_yyy_x_0, tez_yyy_x_1, tez_yyy_y_0, tez_yyy_y_1, tez_yyy_z_0, tez_yyy_z_1, \
                                         tez_yyz_0_0, tez_yyz_0_1, tez_yyz_x_0, tez_yyz_x_1, tez_yyz_y_0, tez_yyz_y_1, \
                                         tez_yyz_z_0, tez_yyz_z_1, tez_yzz_0_0, tez_yzz_0_1, tez_yzz_x_0, tez_yzz_x_1, \
                                         tez_yzz_y_0, tez_yzz_y_1, tez_yzz_z_0, tez_yzz_z_1, tez_zz_x_0, tez_zz_x_1, \
                                         tez_zz_y_0, tez_zz_y_1, tez_zz_z_0, tez_zz_z_1, tez_zzz_0_0, tez_zzz_0_1, \
                                         tez_zzz_x_0, tez_zzz_x_1, tez_zzz_y_0, tez_zzz_y_1, tez_zzz_z_0, tez_zzz_z_1: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxzz_x_0[j] = pa_x[j] * tex_xzz_x_0[j] - pc_x[j] * tex_xzz_x_1[j] + 0.5 * fl1_fx * tex_zz_x_0[j] - 0.5 * fl1_fx * tex_zz_x_1[j] +
                                  0.5 * fl1_fx * tex_xzz_0_0[j] - 0.5 * fl1_fx * tex_xzz_0_1[j] + ta_xzz_x_1[j];

                tey_xxzz_x_0[j] = pa_x[j] * tey_xzz_x_0[j] - pc_x[j] * tey_xzz_x_1[j] + 0.5 * fl1_fx * tey_zz_x_0[j] - 0.5 * fl1_fx * tey_zz_x_1[j] +
                                  0.5 * fl1_fx * tey_xzz_0_0[j] - 0.5 * fl1_fx * tey_xzz_0_1[j];

                tez_xxzz_x_0[j] = pa_x[j] * tez_xzz_x_0[j] - pc_x[j] * tez_xzz_x_1[j] + 0.5 * fl1_fx * tez_zz_x_0[j] - 0.5 * fl1_fx * tez_zz_x_1[j] +
                                  0.5 * fl1_fx * tez_xzz_0_0[j] - 0.5 * fl1_fx * tez_xzz_0_1[j];

                tex_xxzz_y_0[j] =
                    pa_x[j] * tex_xzz_y_0[j] - pc_x[j] * tex_xzz_y_1[j] + 0.5 * fl1_fx * tex_zz_y_0[j] - 0.5 * fl1_fx * tex_zz_y_1[j] + ta_xzz_y_1[j];

                tey_xxzz_y_0[j] = pa_x[j] * tey_xzz_y_0[j] - pc_x[j] * tey_xzz_y_1[j] + 0.5 * fl1_fx * tey_zz_y_0[j] - 0.5 * fl1_fx * tey_zz_y_1[j];

                tez_xxzz_y_0[j] = pa_x[j] * tez_xzz_y_0[j] - pc_x[j] * tez_xzz_y_1[j] + 0.5 * fl1_fx * tez_zz_y_0[j] - 0.5 * fl1_fx * tez_zz_y_1[j];

                tex_xxzz_z_0[j] =
                    pa_x[j] * tex_xzz_z_0[j] - pc_x[j] * tex_xzz_z_1[j] + 0.5 * fl1_fx * tex_zz_z_0[j] - 0.5 * fl1_fx * tex_zz_z_1[j] + ta_xzz_z_1[j];

                tey_xxzz_z_0[j] = pa_x[j] * tey_xzz_z_0[j] - pc_x[j] * tey_xzz_z_1[j] + 0.5 * fl1_fx * tey_zz_z_0[j] - 0.5 * fl1_fx * tey_zz_z_1[j];

                tez_xxzz_z_0[j] = pa_x[j] * tez_xzz_z_0[j] - pc_x[j] * tez_xzz_z_1[j] + 0.5 * fl1_fx * tez_zz_z_0[j] - 0.5 * fl1_fx * tez_zz_z_1[j];

                tex_xyyy_x_0[j] = pa_x[j] * tex_yyy_x_0[j] - pc_x[j] * tex_yyy_x_1[j] + 0.5 * fl1_fx * tex_yyy_0_0[j] -
                                  0.5 * fl1_fx * tex_yyy_0_1[j] + ta_yyy_x_1[j];

                tey_xyyy_x_0[j] = pa_x[j] * tey_yyy_x_0[j] - pc_x[j] * tey_yyy_x_1[j] + 0.5 * fl1_fx * tey_yyy_0_0[j] - 0.5 * fl1_fx * tey_yyy_0_1[j];

                tez_xyyy_x_0[j] = pa_x[j] * tez_yyy_x_0[j] - pc_x[j] * tez_yyy_x_1[j] + 0.5 * fl1_fx * tez_yyy_0_0[j] - 0.5 * fl1_fx * tez_yyy_0_1[j];

                tex_xyyy_y_0[j] = pa_x[j] * tex_yyy_y_0[j] - pc_x[j] * tex_yyy_y_1[j] + ta_yyy_y_1[j];

                tey_xyyy_y_0[j] = pa_x[j] * tey_yyy_y_0[j] - pc_x[j] * tey_yyy_y_1[j];

                tez_xyyy_y_0[j] = pa_x[j] * tez_yyy_y_0[j] - pc_x[j] * tez_yyy_y_1[j];

                tex_xyyy_z_0[j] = pa_x[j] * tex_yyy_z_0[j] - pc_x[j] * tex_yyy_z_1[j] + ta_yyy_z_1[j];

                tey_xyyy_z_0[j] = pa_x[j] * tey_yyy_z_0[j] - pc_x[j] * tey_yyy_z_1[j];

                tez_xyyy_z_0[j] = pa_x[j] * tez_yyy_z_0[j] - pc_x[j] * tez_yyy_z_1[j];

                tex_xyyz_x_0[j] = pa_x[j] * tex_yyz_x_0[j] - pc_x[j] * tex_yyz_x_1[j] + 0.5 * fl1_fx * tex_yyz_0_0[j] -
                                  0.5 * fl1_fx * tex_yyz_0_1[j] + ta_yyz_x_1[j];

                tey_xyyz_x_0[j] = pa_x[j] * tey_yyz_x_0[j] - pc_x[j] * tey_yyz_x_1[j] + 0.5 * fl1_fx * tey_yyz_0_0[j] - 0.5 * fl1_fx * tey_yyz_0_1[j];

                tez_xyyz_x_0[j] = pa_x[j] * tez_yyz_x_0[j] - pc_x[j] * tez_yyz_x_1[j] + 0.5 * fl1_fx * tez_yyz_0_0[j] - 0.5 * fl1_fx * tez_yyz_0_1[j];

                tex_xyyz_y_0[j] = pa_x[j] * tex_yyz_y_0[j] - pc_x[j] * tex_yyz_y_1[j] + ta_yyz_y_1[j];

                tey_xyyz_y_0[j] = pa_x[j] * tey_yyz_y_0[j] - pc_x[j] * tey_yyz_y_1[j];

                tez_xyyz_y_0[j] = pa_x[j] * tez_yyz_y_0[j] - pc_x[j] * tez_yyz_y_1[j];

                tex_xyyz_z_0[j] = pa_x[j] * tex_yyz_z_0[j] - pc_x[j] * tex_yyz_z_1[j] + ta_yyz_z_1[j];

                tey_xyyz_z_0[j] = pa_x[j] * tey_yyz_z_0[j] - pc_x[j] * tey_yyz_z_1[j];

                tez_xyyz_z_0[j] = pa_x[j] * tez_yyz_z_0[j] - pc_x[j] * tez_yyz_z_1[j];

                tex_xyzz_x_0[j] = pa_x[j] * tex_yzz_x_0[j] - pc_x[j] * tex_yzz_x_1[j] + 0.5 * fl1_fx * tex_yzz_0_0[j] -
                                  0.5 * fl1_fx * tex_yzz_0_1[j] + ta_yzz_x_1[j];

                tey_xyzz_x_0[j] = pa_x[j] * tey_yzz_x_0[j] - pc_x[j] * tey_yzz_x_1[j] + 0.5 * fl1_fx * tey_yzz_0_0[j] - 0.5 * fl1_fx * tey_yzz_0_1[j];

                tez_xyzz_x_0[j] = pa_x[j] * tez_yzz_x_0[j] - pc_x[j] * tez_yzz_x_1[j] + 0.5 * fl1_fx * tez_yzz_0_0[j] - 0.5 * fl1_fx * tez_yzz_0_1[j];

                tex_xyzz_y_0[j] = pa_x[j] * tex_yzz_y_0[j] - pc_x[j] * tex_yzz_y_1[j] + ta_yzz_y_1[j];

                tey_xyzz_y_0[j] = pa_x[j] * tey_yzz_y_0[j] - pc_x[j] * tey_yzz_y_1[j];

                tez_xyzz_y_0[j] = pa_x[j] * tez_yzz_y_0[j] - pc_x[j] * tez_yzz_y_1[j];

                tex_xyzz_z_0[j] = pa_x[j] * tex_yzz_z_0[j] - pc_x[j] * tex_yzz_z_1[j] + ta_yzz_z_1[j];

                tey_xyzz_z_0[j] = pa_x[j] * tey_yzz_z_0[j] - pc_x[j] * tey_yzz_z_1[j];

                tez_xyzz_z_0[j] = pa_x[j] * tez_yzz_z_0[j] - pc_x[j] * tez_yzz_z_1[j];

                tex_xzzz_x_0[j] = pa_x[j] * tex_zzz_x_0[j] - pc_x[j] * tex_zzz_x_1[j] + 0.5 * fl1_fx * tex_zzz_0_0[j] -
                                  0.5 * fl1_fx * tex_zzz_0_1[j] + ta_zzz_x_1[j];

                tey_xzzz_x_0[j] = pa_x[j] * tey_zzz_x_0[j] - pc_x[j] * tey_zzz_x_1[j] + 0.5 * fl1_fx * tey_zzz_0_0[j] - 0.5 * fl1_fx * tey_zzz_0_1[j];

                tez_xzzz_x_0[j] = pa_x[j] * tez_zzz_x_0[j] - pc_x[j] * tez_zzz_x_1[j] + 0.5 * fl1_fx * tez_zzz_0_0[j] - 0.5 * fl1_fx * tez_zzz_0_1[j];

                tex_xzzz_y_0[j] = pa_x[j] * tex_zzz_y_0[j] - pc_x[j] * tex_zzz_y_1[j] + ta_zzz_y_1[j];

                tey_xzzz_y_0[j] = pa_x[j] * tey_zzz_y_0[j] - pc_x[j] * tey_zzz_y_1[j];

                tez_xzzz_y_0[j] = pa_x[j] * tez_zzz_y_0[j] - pc_x[j] * tez_zzz_y_1[j];

                tex_xzzz_z_0[j] = pa_x[j] * tex_zzz_z_0[j] - pc_x[j] * tex_zzz_z_1[j] + ta_zzz_z_1[j];

                tey_xzzz_z_0[j] = pa_x[j] * tey_zzz_z_0[j] - pc_x[j] * tey_zzz_z_1[j];

                tez_xzzz_z_0[j] = pa_x[j] * tez_zzz_z_0[j] - pc_x[j] * tez_zzz_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGP_90_135(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
                              const CMemBlock2D<double>& pcDistances,
                              const CGtoBlock&           braGtoBlock,
                              const CGtoBlock&           ketGtoBlock,
                              const int32_t              iContrGto)
{
    // Batch of Integrals (90,135)

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_3_1_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 18);

            auto tey_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 18);

            auto tez_yyy_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 18);

            auto tex_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 19);

            auto tey_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 19);

            auto tez_yyy_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 19);

            auto tex_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 20);

            auto tey_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 20);

            auto tez_yyy_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 20);

            auto tex_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 21);

            auto tey_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 21);

            auto tez_yyz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 21);

            auto tex_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 22);

            auto tey_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 22);

            auto tez_yyz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 22);

            auto tex_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 23);

            auto tey_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 23);

            auto tez_yyz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 23);

            auto tex_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 24);

            auto tey_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 24);

            auto tez_yzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 24);

            auto tex_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 25);

            auto tey_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 25);

            auto tez_yzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 25);

            auto tex_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 26);

            auto tey_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 26);

            auto tez_yzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 26);

            auto tex_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 27);

            auto tey_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 27);

            auto tez_zzz_x_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 27);

            auto tex_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 28);

            auto tey_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 28);

            auto tez_zzz_y_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 28);

            auto tex_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * idx + 29);

            auto tey_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 30 * bdim + 30 * idx + 29);

            auto tez_zzz_z_0 = primBuffer.data(pidx_e_3_1_m0 + 60 * bdim + 30 * idx + 29);

            auto tex_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 18);

            auto tey_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 18);

            auto tez_yyy_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 18);

            auto tex_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 19);

            auto tey_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 19);

            auto tez_yyy_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 19);

            auto tex_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 20);

            auto tey_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 20);

            auto tez_yyy_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 20);

            auto tex_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 21);

            auto tey_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 21);

            auto tez_yyz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 21);

            auto tex_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 22);

            auto tey_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 22);

            auto tez_yyz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 22);

            auto tex_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 23);

            auto tey_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 23);

            auto tez_yyz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 23);

            auto tex_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 24);

            auto tey_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 24);

            auto tez_yzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 24);

            auto tex_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 25);

            auto tey_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 25);

            auto tez_yzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 25);

            auto tex_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 26);

            auto tey_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 26);

            auto tez_yzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 26);

            auto tex_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 27);

            auto tey_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 27);

            auto tez_zzz_x_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 27);

            auto tex_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 28);

            auto tey_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 28);

            auto tez_zzz_y_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 28);

            auto tex_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * idx + 29);

            auto tey_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 30 * bdim + 30 * idx + 29);

            auto tez_zzz_z_1 = primBuffer.data(pidx_e_3_1_m1 + 60 * bdim + 30 * idx + 29);

            auto tex_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 9);

            auto tey_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 10);

            auto tey_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 11);

            auto tey_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 12);

            auto tey_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 13);

            auto tey_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 14);

            auto tey_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 14);

            auto tex_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 15);

            auto tey_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 16);

            auto tey_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * idx + 17);

            auto tey_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_0 = primBuffer.data(pidx_e_2_1_m0 + 36 * bdim + 18 * idx + 17);

            auto tex_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 9);

            auto tey_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 9);

            auto tez_yy_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 9);

            auto tex_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 10);

            auto tey_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 10);

            auto tez_yy_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 10);

            auto tex_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 11);

            auto tey_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 11);

            auto tez_yy_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 11);

            auto tex_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 12);

            auto tey_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 12);

            auto tez_yz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 12);

            auto tex_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 13);

            auto tey_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 13);

            auto tez_yz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 13);

            auto tex_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 14);

            auto tey_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 14);

            auto tez_yz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 14);

            auto tex_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 15);

            auto tey_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 15);

            auto tez_zz_x_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 15);

            auto tex_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 16);

            auto tey_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 16);

            auto tez_zz_y_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 16);

            auto tex_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * idx + 17);

            auto tey_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 18 * bdim + 18 * idx + 17);

            auto tez_zz_z_1 = primBuffer.data(pidx_e_2_1_m1 + 36 * bdim + 18 * idx + 17);

            auto tex_yyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 6);

            auto tey_yyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 6);

            auto tez_yyy_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 6);

            auto tex_yyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 7);

            auto tey_yyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 7);

            auto tez_yyz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 7);

            auto tex_yzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 8);

            auto tey_yzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 8);

            auto tez_yzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 8);

            auto tex_zzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * idx + 9);

            auto tey_zzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 10 * bdim + 10 * idx + 9);

            auto tez_zzz_0_0 = primBuffer.data(pidx_e_3_0_m0 + 20 * bdim + 10 * idx + 9);

            auto tex_yyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 6);

            auto tey_yyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 6);

            auto tez_yyy_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 6);

            auto tex_yyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 7);

            auto tey_yyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 7);

            auto tez_yyz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 7);

            auto tex_yzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 8);

            auto tey_yzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 8);

            auto tez_yzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 8);

            auto tex_zzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * idx + 9);

            auto tey_zzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 10 * bdim + 10 * idx + 9);

            auto tez_zzz_0_1 = primBuffer.data(pidx_e_3_0_m1 + 20 * bdim + 10 * idx + 9);

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

            auto tex_yyyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 30);

            auto tey_yyyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 30);

            auto tez_yyyy_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 30);

            auto tex_yyyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 31);

            auto tey_yyyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 31);

            auto tez_yyyy_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 31);

            auto tex_yyyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 32);

            auto tey_yyyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 32);

            auto tez_yyyy_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 32);

            auto tex_yyyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 33);

            auto tey_yyyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 33);

            auto tez_yyyz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 33);

            auto tex_yyyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 34);

            auto tey_yyyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 34);

            auto tez_yyyz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 34);

            auto tex_yyyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 35);

            auto tey_yyyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 35);

            auto tez_yyyz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 35);

            auto tex_yyzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 36);

            auto tey_yyzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 36);

            auto tez_yyzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 36);

            auto tex_yyzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 37);

            auto tey_yyzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 37);

            auto tez_yyzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 37);

            auto tex_yyzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 38);

            auto tey_yyzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 38);

            auto tez_yyzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 38);

            auto tex_yzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 39);

            auto tey_yzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 39);

            auto tez_yzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 39);

            auto tex_yzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 40);

            auto tey_yzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 40);

            auto tez_yzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 40);

            auto tex_yzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 41);

            auto tey_yzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 41);

            auto tez_yzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 41);

            auto tex_zzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 42);

            auto tey_zzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 42);

            auto tez_zzzz_x_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 42);

            auto tex_zzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 43);

            auto tey_zzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 43);

            auto tez_zzzz_y_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 43);

            auto tex_zzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * idx + 44);

            auto tey_zzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 45 * bdim + 45 * idx + 44);

            auto tez_zzzz_z_0 = primBuffer.data(pidx_e_4_1_m0 + 90 * bdim + 45 * idx + 44);

            // Batch of Integrals (90,135)

            #pragma omp simd aligned(fx, pa_y, pa_z, pc_y, pc_z, ta_yyy_x_1, ta_yyy_y_1, ta_yyy_z_1, ta_yyz_x_1, \
                                         ta_yyz_y_1, ta_yyz_z_1, ta_yzz_x_1, ta_yzz_y_1, ta_yzz_z_1, ta_zzz_x_1, ta_zzz_y_1, \
                                         ta_zzz_z_1, tex_yy_x_0, tex_yy_x_1, tex_yy_y_0, tex_yy_y_1, tex_yy_z_0, tex_yy_z_1, \
                                         tex_yyy_0_0, tex_yyy_0_1, tex_yyy_x_0, tex_yyy_x_1, tex_yyy_y_0, tex_yyy_y_1, \
                                         tex_yyy_z_0, tex_yyy_z_1, tex_yyyy_x_0, tex_yyyy_y_0, tex_yyyy_z_0, tex_yyyz_x_0, \
                                         tex_yyyz_y_0, tex_yyyz_z_0, tex_yyz_0_0, tex_yyz_0_1, tex_yyz_x_0, tex_yyz_x_1, \
                                         tex_yyz_y_0, tex_yyz_y_1, tex_yyz_z_0, tex_yyz_z_1, tex_yyzz_x_0, tex_yyzz_y_0, \
                                         tex_yyzz_z_0, tex_yz_x_0, tex_yz_x_1, tex_yz_y_0, tex_yz_y_1, tex_yz_z_0, tex_yz_z_1, \
                                         tex_yzz_0_0, tex_yzz_0_1, tex_yzz_x_0, tex_yzz_x_1, tex_yzz_y_0, tex_yzz_y_1, \
                                         tex_yzz_z_0, tex_yzz_z_1, tex_yzzz_x_0, tex_yzzz_y_0, tex_yzzz_z_0, tex_zz_x_0, \
                                         tex_zz_x_1, tex_zz_y_0, tex_zz_y_1, tex_zz_z_0, tex_zz_z_1, tex_zzz_0_0, \
                                         tex_zzz_0_1, tex_zzz_x_0, tex_zzz_x_1, tex_zzz_y_0, tex_zzz_y_1, tex_zzz_z_0, \
                                         tex_zzz_z_1, tex_zzzz_x_0, tex_zzzz_y_0, tex_zzzz_z_0, tey_yy_x_0, tey_yy_x_1, \
                                         tey_yy_y_0, tey_yy_y_1, tey_yy_z_0, tey_yy_z_1, tey_yyy_0_0, tey_yyy_0_1, \
                                         tey_yyy_x_0, tey_yyy_x_1, tey_yyy_y_0, tey_yyy_y_1, tey_yyy_z_0, tey_yyy_z_1, \
                                         tey_yyyy_x_0, tey_yyyy_y_0, tey_yyyy_z_0, tey_yyyz_x_0, tey_yyyz_y_0, tey_yyyz_z_0, \
                                         tey_yyz_0_0, tey_yyz_0_1, tey_yyz_x_0, tey_yyz_x_1, tey_yyz_y_0, tey_yyz_y_1, \
                                         tey_yyz_z_0, tey_yyz_z_1, tey_yyzz_x_0, tey_yyzz_y_0, tey_yyzz_z_0, tey_yz_x_0, \
                                         tey_yz_x_1, tey_yz_y_0, tey_yz_y_1, tey_yz_z_0, tey_yz_z_1, tey_yzz_0_0, \
                                         tey_yzz_0_1, tey_yzz_x_0, tey_yzz_x_1, tey_yzz_y_0, tey_yzz_y_1, tey_yzz_z_0, \
                                         tey_yzz_z_1, tey_yzzz_x_0, tey_yzzz_y_0, tey_yzzz_z_0, tey_zz_x_0, tey_zz_x_1, \
                                         tey_zz_y_0, tey_zz_y_1, tey_zz_z_0, tey_zz_z_1, tey_zzz_0_0, tey_zzz_0_1, \
                                         tey_zzz_x_0, tey_zzz_x_1, tey_zzz_y_0, tey_zzz_y_1, tey_zzz_z_0, tey_zzz_z_1, \
                                         tey_zzzz_x_0, tey_zzzz_y_0, tey_zzzz_z_0, tez_yy_x_0, tez_yy_x_1, tez_yy_y_0, \
                                         tez_yy_y_1, tez_yy_z_0, tez_yy_z_1, tez_yyy_0_0, tez_yyy_0_1, tez_yyy_x_0, \
                                         tez_yyy_x_1, tez_yyy_y_0, tez_yyy_y_1, tez_yyy_z_0, tez_yyy_z_1, tez_yyyy_x_0, \
                                         tez_yyyy_y_0, tez_yyyy_z_0, tez_yyyz_x_0, tez_yyyz_y_0, tez_yyyz_z_0, tez_yyz_0_0, \
                                         tez_yyz_0_1, tez_yyz_x_0, tez_yyz_x_1, tez_yyz_y_0, tez_yyz_y_1, tez_yyz_z_0, \
                                         tez_yyz_z_1, tez_yyzz_x_0, tez_yyzz_y_0, tez_yyzz_z_0, tez_yz_x_0, tez_yz_x_1, \
                                         tez_yz_y_0, tez_yz_y_1, tez_yz_z_0, tez_yz_z_1, tez_yzz_0_0, tez_yzz_0_1, \
                                         tez_yzz_x_0, tez_yzz_x_1, tez_yzz_y_0, tez_yzz_y_1, tez_yzz_z_0, tez_yzz_z_1, \
                                         tez_yzzz_x_0, tez_yzzz_y_0, tez_yzzz_z_0, tez_zz_x_0, tez_zz_x_1, tez_zz_y_0, \
                                         tez_zz_y_1, tez_zz_z_0, tez_zz_z_1, tez_zzz_0_0, tez_zzz_0_1, tez_zzz_x_0, \
                                         tez_zzz_x_1, tez_zzz_y_0, tez_zzz_y_1, tez_zzz_z_0, tez_zzz_z_1, tez_zzzz_x_0, \
                                         tez_zzzz_y_0, tez_zzzz_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_yyyy_x_0[j] = pa_y[j] * tex_yyy_x_0[j] - pc_y[j] * tex_yyy_x_1[j] + 1.5 * fl1_fx * tex_yy_x_0[j] - 1.5 * fl1_fx * tex_yy_x_1[j];

                tey_yyyy_x_0[j] =
                    pa_y[j] * tey_yyy_x_0[j] - pc_y[j] * tey_yyy_x_1[j] + 1.5 * fl1_fx * tey_yy_x_0[j] - 1.5 * fl1_fx * tey_yy_x_1[j] + ta_yyy_x_1[j];

                tez_yyyy_x_0[j] = pa_y[j] * tez_yyy_x_0[j] - pc_y[j] * tez_yyy_x_1[j] + 1.5 * fl1_fx * tez_yy_x_0[j] - 1.5 * fl1_fx * tez_yy_x_1[j];

                tex_yyyy_y_0[j] = pa_y[j] * tex_yyy_y_0[j] - pc_y[j] * tex_yyy_y_1[j] + 1.5 * fl1_fx * tex_yy_y_0[j] - 1.5 * fl1_fx * tex_yy_y_1[j] +
                                  0.5 * fl1_fx * tex_yyy_0_0[j] - 0.5 * fl1_fx * tex_yyy_0_1[j];

                tey_yyyy_y_0[j] = pa_y[j] * tey_yyy_y_0[j] - pc_y[j] * tey_yyy_y_1[j] + 1.5 * fl1_fx * tey_yy_y_0[j] - 1.5 * fl1_fx * tey_yy_y_1[j] +
                                  0.5 * fl1_fx * tey_yyy_0_0[j] - 0.5 * fl1_fx * tey_yyy_0_1[j] + ta_yyy_y_1[j];

                tez_yyyy_y_0[j] = pa_y[j] * tez_yyy_y_0[j] - pc_y[j] * tez_yyy_y_1[j] + 1.5 * fl1_fx * tez_yy_y_0[j] - 1.5 * fl1_fx * tez_yy_y_1[j] +
                                  0.5 * fl1_fx * tez_yyy_0_0[j] - 0.5 * fl1_fx * tez_yyy_0_1[j];

                tex_yyyy_z_0[j] = pa_y[j] * tex_yyy_z_0[j] - pc_y[j] * tex_yyy_z_1[j] + 1.5 * fl1_fx * tex_yy_z_0[j] - 1.5 * fl1_fx * tex_yy_z_1[j];

                tey_yyyy_z_0[j] =
                    pa_y[j] * tey_yyy_z_0[j] - pc_y[j] * tey_yyy_z_1[j] + 1.5 * fl1_fx * tey_yy_z_0[j] - 1.5 * fl1_fx * tey_yy_z_1[j] + ta_yyy_z_1[j];

                tez_yyyy_z_0[j] = pa_y[j] * tez_yyy_z_0[j] - pc_y[j] * tez_yyy_z_1[j] + 1.5 * fl1_fx * tez_yy_z_0[j] - 1.5 * fl1_fx * tez_yy_z_1[j];

                tex_yyyz_x_0[j] = pa_y[j] * tex_yyz_x_0[j] - pc_y[j] * tex_yyz_x_1[j] + fl1_fx * tex_yz_x_0[j] - fl1_fx * tex_yz_x_1[j];

                tey_yyyz_x_0[j] =
                    pa_y[j] * tey_yyz_x_0[j] - pc_y[j] * tey_yyz_x_1[j] + fl1_fx * tey_yz_x_0[j] - fl1_fx * tey_yz_x_1[j] + ta_yyz_x_1[j];

                tez_yyyz_x_0[j] = pa_y[j] * tez_yyz_x_0[j] - pc_y[j] * tez_yyz_x_1[j] + fl1_fx * tez_yz_x_0[j] - fl1_fx * tez_yz_x_1[j];

                tex_yyyz_y_0[j] = pa_y[j] * tex_yyz_y_0[j] - pc_y[j] * tex_yyz_y_1[j] + fl1_fx * tex_yz_y_0[j] - fl1_fx * tex_yz_y_1[j] +
                                  0.5 * fl1_fx * tex_yyz_0_0[j] - 0.5 * fl1_fx * tex_yyz_0_1[j];

                tey_yyyz_y_0[j] = pa_y[j] * tey_yyz_y_0[j] - pc_y[j] * tey_yyz_y_1[j] + fl1_fx * tey_yz_y_0[j] - fl1_fx * tey_yz_y_1[j] +
                                  0.5 * fl1_fx * tey_yyz_0_0[j] - 0.5 * fl1_fx * tey_yyz_0_1[j] + ta_yyz_y_1[j];

                tez_yyyz_y_0[j] = pa_y[j] * tez_yyz_y_0[j] - pc_y[j] * tez_yyz_y_1[j] + fl1_fx * tez_yz_y_0[j] - fl1_fx * tez_yz_y_1[j] +
                                  0.5 * fl1_fx * tez_yyz_0_0[j] - 0.5 * fl1_fx * tez_yyz_0_1[j];

                tex_yyyz_z_0[j] = pa_y[j] * tex_yyz_z_0[j] - pc_y[j] * tex_yyz_z_1[j] + fl1_fx * tex_yz_z_0[j] - fl1_fx * tex_yz_z_1[j];

                tey_yyyz_z_0[j] =
                    pa_y[j] * tey_yyz_z_0[j] - pc_y[j] * tey_yyz_z_1[j] + fl1_fx * tey_yz_z_0[j] - fl1_fx * tey_yz_z_1[j] + ta_yyz_z_1[j];

                tez_yyyz_z_0[j] = pa_y[j] * tez_yyz_z_0[j] - pc_y[j] * tez_yyz_z_1[j] + fl1_fx * tez_yz_z_0[j] - fl1_fx * tez_yz_z_1[j];

                tex_yyzz_x_0[j] = pa_y[j] * tex_yzz_x_0[j] - pc_y[j] * tex_yzz_x_1[j] + 0.5 * fl1_fx * tex_zz_x_0[j] - 0.5 * fl1_fx * tex_zz_x_1[j];

                tey_yyzz_x_0[j] =
                    pa_y[j] * tey_yzz_x_0[j] - pc_y[j] * tey_yzz_x_1[j] + 0.5 * fl1_fx * tey_zz_x_0[j] - 0.5 * fl1_fx * tey_zz_x_1[j] + ta_yzz_x_1[j];

                tez_yyzz_x_0[j] = pa_y[j] * tez_yzz_x_0[j] - pc_y[j] * tez_yzz_x_1[j] + 0.5 * fl1_fx * tez_zz_x_0[j] - 0.5 * fl1_fx * tez_zz_x_1[j];

                tex_yyzz_y_0[j] = pa_y[j] * tex_yzz_y_0[j] - pc_y[j] * tex_yzz_y_1[j] + 0.5 * fl1_fx * tex_zz_y_0[j] - 0.5 * fl1_fx * tex_zz_y_1[j] +
                                  0.5 * fl1_fx * tex_yzz_0_0[j] - 0.5 * fl1_fx * tex_yzz_0_1[j];

                tey_yyzz_y_0[j] = pa_y[j] * tey_yzz_y_0[j] - pc_y[j] * tey_yzz_y_1[j] + 0.5 * fl1_fx * tey_zz_y_0[j] - 0.5 * fl1_fx * tey_zz_y_1[j] +
                                  0.5 * fl1_fx * tey_yzz_0_0[j] - 0.5 * fl1_fx * tey_yzz_0_1[j] + ta_yzz_y_1[j];

                tez_yyzz_y_0[j] = pa_y[j] * tez_yzz_y_0[j] - pc_y[j] * tez_yzz_y_1[j] + 0.5 * fl1_fx * tez_zz_y_0[j] - 0.5 * fl1_fx * tez_zz_y_1[j] +
                                  0.5 * fl1_fx * tez_yzz_0_0[j] - 0.5 * fl1_fx * tez_yzz_0_1[j];

                tex_yyzz_z_0[j] = pa_y[j] * tex_yzz_z_0[j] - pc_y[j] * tex_yzz_z_1[j] + 0.5 * fl1_fx * tex_zz_z_0[j] - 0.5 * fl1_fx * tex_zz_z_1[j];

                tey_yyzz_z_0[j] =
                    pa_y[j] * tey_yzz_z_0[j] - pc_y[j] * tey_yzz_z_1[j] + 0.5 * fl1_fx * tey_zz_z_0[j] - 0.5 * fl1_fx * tey_zz_z_1[j] + ta_yzz_z_1[j];

                tez_yyzz_z_0[j] = pa_y[j] * tez_yzz_z_0[j] - pc_y[j] * tez_yzz_z_1[j] + 0.5 * fl1_fx * tez_zz_z_0[j] - 0.5 * fl1_fx * tez_zz_z_1[j];

                tex_yzzz_x_0[j] = pa_y[j] * tex_zzz_x_0[j] - pc_y[j] * tex_zzz_x_1[j];

                tey_yzzz_x_0[j] = pa_y[j] * tey_zzz_x_0[j] - pc_y[j] * tey_zzz_x_1[j] + ta_zzz_x_1[j];

                tez_yzzz_x_0[j] = pa_y[j] * tez_zzz_x_0[j] - pc_y[j] * tez_zzz_x_1[j];

                tex_yzzz_y_0[j] = pa_y[j] * tex_zzz_y_0[j] - pc_y[j] * tex_zzz_y_1[j] + 0.5 * fl1_fx * tex_zzz_0_0[j] - 0.5 * fl1_fx * tex_zzz_0_1[j];

                tey_yzzz_y_0[j] = pa_y[j] * tey_zzz_y_0[j] - pc_y[j] * tey_zzz_y_1[j] + 0.5 * fl1_fx * tey_zzz_0_0[j] -
                                  0.5 * fl1_fx * tey_zzz_0_1[j] + ta_zzz_y_1[j];

                tez_yzzz_y_0[j] = pa_y[j] * tez_zzz_y_0[j] - pc_y[j] * tez_zzz_y_1[j] + 0.5 * fl1_fx * tez_zzz_0_0[j] - 0.5 * fl1_fx * tez_zzz_0_1[j];

                tex_yzzz_z_0[j] = pa_y[j] * tex_zzz_z_0[j] - pc_y[j] * tex_zzz_z_1[j];

                tey_yzzz_z_0[j] = pa_y[j] * tey_zzz_z_0[j] - pc_y[j] * tey_zzz_z_1[j] + ta_zzz_z_1[j];

                tez_yzzz_z_0[j] = pa_y[j] * tez_zzz_z_0[j] - pc_y[j] * tez_zzz_z_1[j];

                tex_zzzz_x_0[j] = pa_z[j] * tex_zzz_x_0[j] - pc_z[j] * tex_zzz_x_1[j] + 1.5 * fl1_fx * tex_zz_x_0[j] - 1.5 * fl1_fx * tex_zz_x_1[j];

                tey_zzzz_x_0[j] = pa_z[j] * tey_zzz_x_0[j] - pc_z[j] * tey_zzz_x_1[j] + 1.5 * fl1_fx * tey_zz_x_0[j] - 1.5 * fl1_fx * tey_zz_x_1[j];

                tez_zzzz_x_0[j] =
                    pa_z[j] * tez_zzz_x_0[j] - pc_z[j] * tez_zzz_x_1[j] + 1.5 * fl1_fx * tez_zz_x_0[j] - 1.5 * fl1_fx * tez_zz_x_1[j] + ta_zzz_x_1[j];

                tex_zzzz_y_0[j] = pa_z[j] * tex_zzz_y_0[j] - pc_z[j] * tex_zzz_y_1[j] + 1.5 * fl1_fx * tex_zz_y_0[j] - 1.5 * fl1_fx * tex_zz_y_1[j];

                tey_zzzz_y_0[j] = pa_z[j] * tey_zzz_y_0[j] - pc_z[j] * tey_zzz_y_1[j] + 1.5 * fl1_fx * tey_zz_y_0[j] - 1.5 * fl1_fx * tey_zz_y_1[j];

                tez_zzzz_y_0[j] =
                    pa_z[j] * tez_zzz_y_0[j] - pc_z[j] * tez_zzz_y_1[j] + 1.5 * fl1_fx * tez_zz_y_0[j] - 1.5 * fl1_fx * tez_zz_y_1[j] + ta_zzz_y_1[j];

                tex_zzzz_z_0[j] = pa_z[j] * tex_zzz_z_0[j] - pc_z[j] * tex_zzz_z_1[j] + 1.5 * fl1_fx * tex_zz_z_0[j] - 1.5 * fl1_fx * tex_zz_z_1[j] +
                                  0.5 * fl1_fx * tex_zzz_0_0[j] - 0.5 * fl1_fx * tex_zzz_0_1[j];

                tey_zzzz_z_0[j] = pa_z[j] * tey_zzz_z_0[j] - pc_z[j] * tey_zzz_z_1[j] + 1.5 * fl1_fx * tey_zz_z_0[j] - 1.5 * fl1_fx * tey_zz_z_1[j] +
                                  0.5 * fl1_fx * tey_zzz_0_0[j] - 0.5 * fl1_fx * tey_zzz_0_1[j];

                tez_zzzz_z_0[j] = pa_z[j] * tez_zzz_z_0[j] - pc_z[j] * tez_zzz_z_1[j] + 1.5 * fl1_fx * tez_zz_z_0[j] - 1.5 * fl1_fx * tez_zz_z_1[j] +
                                  0.5 * fl1_fx * tez_zzz_0_0[j] - 0.5 * fl1_fx * tez_zzz_0_1[j] + ta_zzz_z_1[j];
            }

            idx++;
        }
    }
}

}  // namespace efieldrecfunc
