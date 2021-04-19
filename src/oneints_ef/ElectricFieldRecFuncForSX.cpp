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

#include "ElectricFieldRecFuncForSX.hpp"

namespace efieldrecfunc {  // efieldrecfunc namespace

void
compElectricFieldForSS(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx == -1) continue;

        // set up indexes of auxilary integral

        auto sidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up Obara-Saika prefactors

            auto fg = osFactors.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto fex = primBuffer.data(pidx + idx);

            auto fey = primBuffer.data(pidx + bdim + idx);

            auto fez = primBuffer.data(pidx + 2 * bdim + idx);

            auto fa = primBuffer.data(sidx + idx);

            #pragma omp simd aligned(fg, pc_x, pc_y, pc_z, fex, fey, fez, fa: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fx = 2.0 * fg[j] * fa[j];

                fex[j] = fx * pc_x[j];

                fey[j] = fx * pc_y[j];

                fez[j] = fx * pc_z[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForSP(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& pbDistances,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_0_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tex_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + idx);

            auto tey_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + bdim + idx);

            auto tez_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + 2 * bdim + idx);

            auto tex_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + idx);

            auto tey_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + bdim + idx);

            auto tez_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + 2 * bdim + idx);

            auto ta_0_0_1 = primBuffer.data(pidx_a_0_0_m1 + idx);

            // set up pointers to integrals

            auto tex_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx);

            auto tey_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx);

            auto tez_0_x_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx);

            auto tex_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 1);

            auto tey_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_0_y_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * idx + 2);

            auto tey_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_0_z_0 = primBuffer.data(pidx_e_0_1_m0 + 6 * bdim + 3 * idx + 2);

            #pragma omp simd aligned(pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_0_1, tex_0_0_0, tex_0_0_1, tex_0_x_0, \
                                         tex_0_y_0, tex_0_z_0, tey_0_0_0, tey_0_0_1, tey_0_x_0, tey_0_y_0, tey_0_z_0, \
                                         tez_0_0_0, tez_0_0_1, tez_0_x_0, tez_0_y_0, tez_0_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                tex_0_x_0[j] = pb_x[j] * tex_0_0_0[j] - pc_x[j] * tex_0_0_1[j] + ta_0_0_1[j];

                tey_0_x_0[j] = pb_x[j] * tey_0_0_0[j] - pc_x[j] * tey_0_0_1[j];

                tez_0_x_0[j] = pb_x[j] * tez_0_0_0[j] - pc_x[j] * tez_0_0_1[j];

                tex_0_y_0[j] = pb_y[j] * tex_0_0_0[j] - pc_y[j] * tex_0_0_1[j];

                tey_0_y_0[j] = pb_y[j] * tey_0_0_0[j] - pc_y[j] * tey_0_0_1[j] + ta_0_0_1[j];

                tez_0_y_0[j] = pb_y[j] * tez_0_0_0[j] - pc_y[j] * tez_0_0_1[j];

                tex_0_z_0[j] = pb_z[j] * tex_0_0_0[j] - pc_z[j] * tex_0_0_1[j];

                tey_0_z_0[j] = pb_z[j] * tey_0_0_0[j] - pc_z[j] * tey_0_0_1[j];

                tez_0_z_0[j] = pb_z[j] * tez_0_0_0[j] - pc_z[j] * tez_0_0_1[j] + ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForPS(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_1_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_a_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        // loop over contracted GTO on bra side

        int32_t idx = 0;

        for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
        {
            // set up pointers to tensors product of distances R(PA) = P - A

            auto pa_x = paDistances.data(3 * idx);

            auto pa_y = paDistances.data(3 * idx + 1);

            auto pa_z = paDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

            auto pc_y = pcDistances.data(3 * idx + 1);

            auto pc_z = pcDistances.data(3 * idx + 2);

            // set up pointers to auxilary integrals

            auto tex_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + idx);

            auto tey_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + bdim + idx);

            auto tez_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + 2 * bdim + idx);

            auto tex_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + idx);

            auto tey_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + bdim + idx);

            auto tez_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + 2 * bdim + idx);

            auto ta_0_0_1 = primBuffer.data(pidx_a_0_0_m1 + idx);

            // set up pointers to integrals

            auto tex_x_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx);

            auto tey_x_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx);

            auto tez_x_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx);

            auto tex_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx + 1);

            auto tey_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx + 1);

            auto tez_y_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx + 1);

            auto tex_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * idx + 2);

            auto tey_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 3 * bdim + 3 * idx + 2);

            auto tez_z_0_0 = primBuffer.data(pidx_e_1_0_m0 + 6 * bdim + 3 * idx + 2);

            #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_0_1, tex_0_0_0, tex_0_0_1, tex_x_0_0, \
                                         tex_y_0_0, tex_z_0_0, tey_0_0_0, tey_0_0_1, tey_x_0_0, tey_y_0_0, tey_z_0_0, \
                                         tez_0_0_0, tez_0_0_1, tez_x_0_0, tez_y_0_0, tez_z_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                tex_x_0_0[j] = pa_x[j] * tex_0_0_0[j] - pc_x[j] * tex_0_0_1[j] + ta_0_0_1[j];

                tey_x_0_0[j] = pa_x[j] * tey_0_0_0[j] - pc_x[j] * tey_0_0_1[j];

                tez_x_0_0[j] = pa_x[j] * tez_0_0_0[j] - pc_x[j] * tez_0_0_1[j];

                tex_y_0_0[j] = pa_y[j] * tex_0_0_0[j] - pc_y[j] * tex_0_0_1[j];

                tey_y_0_0[j] = pa_y[j] * tey_0_0_0[j] - pc_y[j] * tey_0_0_1[j] + ta_0_0_1[j];

                tez_y_0_0[j] = pa_y[j] * tez_0_0_0[j] - pc_y[j] * tez_0_0_1[j];

                tex_z_0_0[j] = pa_z[j] * tex_0_0_0[j] - pc_z[j] * tex_0_0_1[j];

                tey_z_0_0[j] = pa_z[j] * tey_0_0_0[j] - pc_z[j] * tey_0_0_1[j];

                tez_z_0_0[j] = pa_z[j] * tez_0_0_0[j] - pc_z[j] * tez_0_0_1[j] + ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForSD(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& pbDistances,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_0_2_m0 == -1) continue;

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_x_1, ta_0_y_1, ta_0_z_1, tex_0_0_0, \
                                         tex_0_0_1, tex_0_x_0, tex_0_x_1, tex_0_xx_0, tex_0_xy_0, tex_0_xz_0, tex_0_y_0, \
                                         tex_0_y_1, tex_0_yy_0, tex_0_yz_0, tex_0_z_0, tex_0_z_1, tex_0_zz_0, tey_0_0_0, \
                                         tey_0_0_1, tey_0_x_0, tey_0_x_1, tey_0_xx_0, tey_0_xy_0, tey_0_xz_0, tey_0_y_0, \
                                         tey_0_y_1, tey_0_yy_0, tey_0_yz_0, tey_0_z_0, tey_0_z_1, tey_0_zz_0, tez_0_0_0, \
                                         tez_0_0_1, tez_0_x_0, tez_0_x_1, tez_0_xx_0, tez_0_xy_0, tez_0_xz_0, tez_0_y_0, \
                                         tez_0_y_1, tez_0_yy_0, tez_0_yz_0, tez_0_z_0, tez_0_z_1, tez_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_0_xx_0[j] =
                    pb_x[j] * tex_0_x_0[j] - pc_x[j] * tex_0_x_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j] + ta_0_x_1[j];

                tey_0_xx_0[j] = pb_x[j] * tey_0_x_0[j] - pc_x[j] * tey_0_x_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j];

                tez_0_xx_0[j] = pb_x[j] * tez_0_x_0[j] - pc_x[j] * tez_0_x_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j];

                tex_0_xy_0[j] = pb_x[j] * tex_0_y_0[j] - pc_x[j] * tex_0_y_1[j] + ta_0_y_1[j];

                tey_0_xy_0[j] = pb_x[j] * tey_0_y_0[j] - pc_x[j] * tey_0_y_1[j];

                tez_0_xy_0[j] = pb_x[j] * tez_0_y_0[j] - pc_x[j] * tez_0_y_1[j];

                tex_0_xz_0[j] = pb_x[j] * tex_0_z_0[j] - pc_x[j] * tex_0_z_1[j] + ta_0_z_1[j];

                tey_0_xz_0[j] = pb_x[j] * tey_0_z_0[j] - pc_x[j] * tey_0_z_1[j];

                tez_0_xz_0[j] = pb_x[j] * tez_0_z_0[j] - pc_x[j] * tez_0_z_1[j];

                tex_0_yy_0[j] = pb_y[j] * tex_0_y_0[j] - pc_y[j] * tex_0_y_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j];

                tey_0_yy_0[j] =
                    pb_y[j] * tey_0_y_0[j] - pc_y[j] * tey_0_y_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j] + ta_0_y_1[j];

                tez_0_yy_0[j] = pb_y[j] * tez_0_y_0[j] - pc_y[j] * tez_0_y_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j];

                tex_0_yz_0[j] = pb_y[j] * tex_0_z_0[j] - pc_y[j] * tex_0_z_1[j];

                tey_0_yz_0[j] = pb_y[j] * tey_0_z_0[j] - pc_y[j] * tey_0_z_1[j] + ta_0_z_1[j];

                tez_0_yz_0[j] = pb_y[j] * tez_0_z_0[j] - pc_y[j] * tez_0_z_1[j];

                tex_0_zz_0[j] = pb_z[j] * tex_0_z_0[j] - pc_z[j] * tex_0_z_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j];

                tey_0_zz_0[j] = pb_z[j] * tey_0_z_0[j] - pc_z[j] * tey_0_z_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j];

                tez_0_zz_0[j] =
                    pb_z[j] * tez_0_z_0[j] - pc_z[j] * tez_0_z_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j] + ta_0_z_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForDS(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_2_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_0_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + idx);

            auto tey_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + bdim + idx);

            auto tez_0_0_0 = primBuffer.data(pidx_e_0_0_m0 + 2 * bdim + idx);

            auto tex_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + idx);

            auto tey_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + bdim + idx);

            auto tez_0_0_1 = primBuffer.data(pidx_e_0_0_m1 + 2 * bdim + idx);

            auto ta_x_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx);

            auto ta_y_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 1);

            auto ta_z_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 2);

            // set up pointers to integrals

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

            auto tex_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 5);

            auto tey_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 5);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_x_0_1, ta_y_0_1, ta_z_0_1, tex_0_0_0, \
                                         tex_0_0_1, tex_x_0_0, tex_x_0_1, tex_xx_0_0, tex_xy_0_0, tex_xz_0_0, tex_y_0_0, \
                                         tex_y_0_1, tex_yy_0_0, tex_yz_0_0, tex_z_0_0, tex_z_0_1, tex_zz_0_0, tey_0_0_0, \
                                         tey_0_0_1, tey_x_0_0, tey_x_0_1, tey_xx_0_0, tey_xy_0_0, tey_xz_0_0, tey_y_0_0, \
                                         tey_y_0_1, tey_yy_0_0, tey_yz_0_0, tey_z_0_0, tey_z_0_1, tey_zz_0_0, tez_0_0_0, \
                                         tez_0_0_1, tez_x_0_0, tez_x_0_1, tez_xx_0_0, tez_xy_0_0, tez_xz_0_0, tez_y_0_0, \
                                         tez_y_0_1, tez_yy_0_0, tez_yz_0_0, tez_z_0_0, tez_z_0_1, tez_zz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xx_0_0[j] =
                    pa_x[j] * tex_x_0_0[j] - pc_x[j] * tex_x_0_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j] + ta_x_0_1[j];

                tey_xx_0_0[j] = pa_x[j] * tey_x_0_0[j] - pc_x[j] * tey_x_0_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j];

                tez_xx_0_0[j] = pa_x[j] * tez_x_0_0[j] - pc_x[j] * tez_x_0_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j];

                tex_xy_0_0[j] = pa_x[j] * tex_y_0_0[j] - pc_x[j] * tex_y_0_1[j] + ta_y_0_1[j];

                tey_xy_0_0[j] = pa_x[j] * tey_y_0_0[j] - pc_x[j] * tey_y_0_1[j];

                tez_xy_0_0[j] = pa_x[j] * tez_y_0_0[j] - pc_x[j] * tez_y_0_1[j];

                tex_xz_0_0[j] = pa_x[j] * tex_z_0_0[j] - pc_x[j] * tex_z_0_1[j] + ta_z_0_1[j];

                tey_xz_0_0[j] = pa_x[j] * tey_z_0_0[j] - pc_x[j] * tey_z_0_1[j];

                tez_xz_0_0[j] = pa_x[j] * tez_z_0_0[j] - pc_x[j] * tez_z_0_1[j];

                tex_yy_0_0[j] = pa_y[j] * tex_y_0_0[j] - pc_y[j] * tex_y_0_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j];

                tey_yy_0_0[j] =
                    pa_y[j] * tey_y_0_0[j] - pc_y[j] * tey_y_0_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j] + ta_y_0_1[j];

                tez_yy_0_0[j] = pa_y[j] * tez_y_0_0[j] - pc_y[j] * tez_y_0_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j];

                tex_yz_0_0[j] = pa_y[j] * tex_z_0_0[j] - pc_y[j] * tex_z_0_1[j];

                tey_yz_0_0[j] = pa_y[j] * tey_z_0_0[j] - pc_y[j] * tey_z_0_1[j] + ta_z_0_1[j];

                tez_yz_0_0[j] = pa_y[j] * tez_z_0_0[j] - pc_y[j] * tez_z_0_1[j];

                tex_zz_0_0[j] = pa_z[j] * tex_z_0_0[j] - pc_z[j] * tex_z_0_1[j] + 0.5 * fl1_fx * tex_0_0_0[j] - 0.5 * fl1_fx * tex_0_0_1[j];

                tey_zz_0_0[j] = pa_z[j] * tey_z_0_0[j] - pc_z[j] * tey_z_0_1[j] + 0.5 * fl1_fx * tey_0_0_0[j] - 0.5 * fl1_fx * tey_0_0_1[j];

                tez_zz_0_0[j] =
                    pa_z[j] * tez_z_0_0[j] - pc_z[j] * tez_z_0_1[j] + 0.5 * fl1_fx * tez_0_0_0[j] - 0.5 * fl1_fx * tez_0_0_1[j] + ta_z_0_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForSF(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& pbDistances,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_0_3_m0 == -1) continue;

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_xx_1, ta_0_xy_1, ta_0_xz_1, \
                                         ta_0_yy_1, ta_0_yz_1, ta_0_zz_1, tex_0_x_0, tex_0_x_1, tex_0_xx_0, tex_0_xx_1, \
                                         tex_0_xxx_0, tex_0_xxy_0, tex_0_xxz_0, tex_0_xy_0, tex_0_xy_1, tex_0_xyy_0, \
                                         tex_0_xyz_0, tex_0_xz_0, tex_0_xz_1, tex_0_xzz_0, tex_0_y_0, tex_0_y_1, tex_0_yy_0, \
                                         tex_0_yy_1, tex_0_yyy_0, tex_0_yyz_0, tex_0_yz_0, tex_0_yz_1, tex_0_yzz_0, \
                                         tex_0_z_0, tex_0_z_1, tex_0_zz_0, tex_0_zz_1, tex_0_zzz_0, tey_0_x_0, tey_0_x_1, \
                                         tey_0_xx_0, tey_0_xx_1, tey_0_xxx_0, tey_0_xxy_0, tey_0_xxz_0, tey_0_xy_0, \
                                         tey_0_xy_1, tey_0_xyy_0, tey_0_xyz_0, tey_0_xz_0, tey_0_xz_1, tey_0_xzz_0, \
                                         tey_0_y_0, tey_0_y_1, tey_0_yy_0, tey_0_yy_1, tey_0_yyy_0, tey_0_yyz_0, \
                                         tey_0_yz_0, tey_0_yz_1, tey_0_yzz_0, tey_0_z_0, tey_0_z_1, tey_0_zz_0, tey_0_zz_1, \
                                         tey_0_zzz_0, tez_0_x_0, tez_0_x_1, tez_0_xx_0, tez_0_xx_1, tez_0_xxx_0, tez_0_xxy_0, \
                                         tez_0_xxz_0, tez_0_xy_0, tez_0_xy_1, tez_0_xyy_0, tez_0_xyz_0, tez_0_xz_0, \
                                         tez_0_xz_1, tez_0_xzz_0, tez_0_y_0, tez_0_y_1, tez_0_yy_0, tez_0_yy_1, tez_0_yyy_0, \
                                         tez_0_yyz_0, tez_0_yz_0, tez_0_yz_1, tez_0_yzz_0, tez_0_z_0, tez_0_z_1, tez_0_zz_0, \
                                         tez_0_zz_1, tez_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_0_xxx_0[j] = pb_x[j] * tex_0_xx_0[j] - pc_x[j] * tex_0_xx_1[j] + fl1_fx * tex_0_x_0[j] - fl1_fx * tex_0_x_1[j] + ta_0_xx_1[j];

                tey_0_xxx_0[j] = pb_x[j] * tey_0_xx_0[j] - pc_x[j] * tey_0_xx_1[j] + fl1_fx * tey_0_x_0[j] - fl1_fx * tey_0_x_1[j];

                tez_0_xxx_0[j] = pb_x[j] * tez_0_xx_0[j] - pc_x[j] * tez_0_xx_1[j] + fl1_fx * tez_0_x_0[j] - fl1_fx * tez_0_x_1[j];

                tex_0_xxy_0[j] =
                    pb_x[j] * tex_0_xy_0[j] - pc_x[j] * tex_0_xy_1[j] + 0.5 * fl1_fx * tex_0_y_0[j] - 0.5 * fl1_fx * tex_0_y_1[j] + ta_0_xy_1[j];

                tey_0_xxy_0[j] = pb_x[j] * tey_0_xy_0[j] - pc_x[j] * tey_0_xy_1[j] + 0.5 * fl1_fx * tey_0_y_0[j] - 0.5 * fl1_fx * tey_0_y_1[j];

                tez_0_xxy_0[j] = pb_x[j] * tez_0_xy_0[j] - pc_x[j] * tez_0_xy_1[j] + 0.5 * fl1_fx * tez_0_y_0[j] - 0.5 * fl1_fx * tez_0_y_1[j];

                tex_0_xxz_0[j] =
                    pb_x[j] * tex_0_xz_0[j] - pc_x[j] * tex_0_xz_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j] + ta_0_xz_1[j];

                tey_0_xxz_0[j] = pb_x[j] * tey_0_xz_0[j] - pc_x[j] * tey_0_xz_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j];

                tez_0_xxz_0[j] = pb_x[j] * tez_0_xz_0[j] - pc_x[j] * tez_0_xz_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j];

                tex_0_xyy_0[j] = pb_x[j] * tex_0_yy_0[j] - pc_x[j] * tex_0_yy_1[j] + ta_0_yy_1[j];

                tey_0_xyy_0[j] = pb_x[j] * tey_0_yy_0[j] - pc_x[j] * tey_0_yy_1[j];

                tez_0_xyy_0[j] = pb_x[j] * tez_0_yy_0[j] - pc_x[j] * tez_0_yy_1[j];

                tex_0_xyz_0[j] = pb_x[j] * tex_0_yz_0[j] - pc_x[j] * tex_0_yz_1[j] + ta_0_yz_1[j];

                tey_0_xyz_0[j] = pb_x[j] * tey_0_yz_0[j] - pc_x[j] * tey_0_yz_1[j];

                tez_0_xyz_0[j] = pb_x[j] * tez_0_yz_0[j] - pc_x[j] * tez_0_yz_1[j];

                tex_0_xzz_0[j] = pb_x[j] * tex_0_zz_0[j] - pc_x[j] * tex_0_zz_1[j] + ta_0_zz_1[j];

                tey_0_xzz_0[j] = pb_x[j] * tey_0_zz_0[j] - pc_x[j] * tey_0_zz_1[j];

                tez_0_xzz_0[j] = pb_x[j] * tez_0_zz_0[j] - pc_x[j] * tez_0_zz_1[j];

                tex_0_yyy_0[j] = pb_y[j] * tex_0_yy_0[j] - pc_y[j] * tex_0_yy_1[j] + fl1_fx * tex_0_y_0[j] - fl1_fx * tex_0_y_1[j];

                tey_0_yyy_0[j] = pb_y[j] * tey_0_yy_0[j] - pc_y[j] * tey_0_yy_1[j] + fl1_fx * tey_0_y_0[j] - fl1_fx * tey_0_y_1[j] + ta_0_yy_1[j];

                tez_0_yyy_0[j] = pb_y[j] * tez_0_yy_0[j] - pc_y[j] * tez_0_yy_1[j] + fl1_fx * tez_0_y_0[j] - fl1_fx * tez_0_y_1[j];

                tex_0_yyz_0[j] = pb_y[j] * tex_0_yz_0[j] - pc_y[j] * tex_0_yz_1[j] + 0.5 * fl1_fx * tex_0_z_0[j] - 0.5 * fl1_fx * tex_0_z_1[j];

                tey_0_yyz_0[j] =
                    pb_y[j] * tey_0_yz_0[j] - pc_y[j] * tey_0_yz_1[j] + 0.5 * fl1_fx * tey_0_z_0[j] - 0.5 * fl1_fx * tey_0_z_1[j] + ta_0_yz_1[j];

                tez_0_yyz_0[j] = pb_y[j] * tez_0_yz_0[j] - pc_y[j] * tez_0_yz_1[j] + 0.5 * fl1_fx * tez_0_z_0[j] - 0.5 * fl1_fx * tez_0_z_1[j];

                tex_0_yzz_0[j] = pb_y[j] * tex_0_zz_0[j] - pc_y[j] * tex_0_zz_1[j];

                tey_0_yzz_0[j] = pb_y[j] * tey_0_zz_0[j] - pc_y[j] * tey_0_zz_1[j] + ta_0_zz_1[j];

                tez_0_yzz_0[j] = pb_y[j] * tez_0_zz_0[j] - pc_y[j] * tez_0_zz_1[j];

                tex_0_zzz_0[j] = pb_z[j] * tex_0_zz_0[j] - pc_z[j] * tex_0_zz_1[j] + fl1_fx * tex_0_z_0[j] - fl1_fx * tex_0_z_1[j];

                tey_0_zzz_0[j] = pb_z[j] * tey_0_zz_0[j] - pc_z[j] * tey_0_zz_1[j] + fl1_fx * tey_0_z_0[j] - fl1_fx * tey_0_z_1[j];

                tez_0_zzz_0[j] = pb_z[j] * tez_0_zz_0[j] - pc_z[j] * tez_0_zz_1[j] + fl1_fx * tez_0_z_0[j] - fl1_fx * tez_0_z_1[j] + ta_0_zz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForFS(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_3_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_1_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 5);

            auto tey_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 5);

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

            auto tex_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 5);

            auto tey_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 5);

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

            auto ta_xx_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx);

            auto ta_xy_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 1);

            auto ta_xz_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 2);

            auto ta_yy_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 3);

            auto ta_yz_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 4);

            auto ta_zz_0_1 = primBuffer.data(pidx_a_2_0_m1 + 6 * idx + 5);

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xx_0_1, ta_xy_0_1, ta_xz_0_1, \
                                         ta_yy_0_1, ta_yz_0_1, ta_zz_0_1, tex_x_0_0, tex_x_0_1, tex_xx_0_0, tex_xx_0_1, \
                                         tex_xxx_0_0, tex_xxy_0_0, tex_xxz_0_0, tex_xy_0_0, tex_xy_0_1, tex_xyy_0_0, \
                                         tex_xyz_0_0, tex_xz_0_0, tex_xz_0_1, tex_xzz_0_0, tex_y_0_0, tex_y_0_1, tex_yy_0_0, \
                                         tex_yy_0_1, tex_yyy_0_0, tex_yyz_0_0, tex_yz_0_0, tex_yz_0_1, tex_yzz_0_0, \
                                         tex_z_0_0, tex_z_0_1, tex_zz_0_0, tex_zz_0_1, tex_zzz_0_0, tey_x_0_0, tey_x_0_1, \
                                         tey_xx_0_0, tey_xx_0_1, tey_xxx_0_0, tey_xxy_0_0, tey_xxz_0_0, tey_xy_0_0, \
                                         tey_xy_0_1, tey_xyy_0_0, tey_xyz_0_0, tey_xz_0_0, tey_xz_0_1, tey_xzz_0_0, \
                                         tey_y_0_0, tey_y_0_1, tey_yy_0_0, tey_yy_0_1, tey_yyy_0_0, tey_yyz_0_0, \
                                         tey_yz_0_0, tey_yz_0_1, tey_yzz_0_0, tey_z_0_0, tey_z_0_1, tey_zz_0_0, tey_zz_0_1, \
                                         tey_zzz_0_0, tez_x_0_0, tez_x_0_1, tez_xx_0_0, tez_xx_0_1, tez_xxx_0_0, tez_xxy_0_0, \
                                         tez_xxz_0_0, tez_xy_0_0, tez_xy_0_1, tez_xyy_0_0, tez_xyz_0_0, tez_xz_0_0, \
                                         tez_xz_0_1, tez_xzz_0_0, tez_y_0_0, tez_y_0_1, tez_yy_0_0, tez_yy_0_1, tez_yyy_0_0, \
                                         tez_yyz_0_0, tez_yz_0_0, tez_yz_0_1, tez_yzz_0_0, tez_z_0_0, tez_z_0_1, tez_zz_0_0, \
                                         tez_zz_0_1, tez_zzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxx_0_0[j] = pa_x[j] * tex_xx_0_0[j] - pc_x[j] * tex_xx_0_1[j] + fl1_fx * tex_x_0_0[j] - fl1_fx * tex_x_0_1[j] + ta_xx_0_1[j];

                tey_xxx_0_0[j] = pa_x[j] * tey_xx_0_0[j] - pc_x[j] * tey_xx_0_1[j] + fl1_fx * tey_x_0_0[j] - fl1_fx * tey_x_0_1[j];

                tez_xxx_0_0[j] = pa_x[j] * tez_xx_0_0[j] - pc_x[j] * tez_xx_0_1[j] + fl1_fx * tez_x_0_0[j] - fl1_fx * tez_x_0_1[j];

                tex_xxy_0_0[j] =
                    pa_x[j] * tex_xy_0_0[j] - pc_x[j] * tex_xy_0_1[j] + 0.5 * fl1_fx * tex_y_0_0[j] - 0.5 * fl1_fx * tex_y_0_1[j] + ta_xy_0_1[j];

                tey_xxy_0_0[j] = pa_x[j] * tey_xy_0_0[j] - pc_x[j] * tey_xy_0_1[j] + 0.5 * fl1_fx * tey_y_0_0[j] - 0.5 * fl1_fx * tey_y_0_1[j];

                tez_xxy_0_0[j] = pa_x[j] * tez_xy_0_0[j] - pc_x[j] * tez_xy_0_1[j] + 0.5 * fl1_fx * tez_y_0_0[j] - 0.5 * fl1_fx * tez_y_0_1[j];

                tex_xxz_0_0[j] =
                    pa_x[j] * tex_xz_0_0[j] - pc_x[j] * tex_xz_0_1[j] + 0.5 * fl1_fx * tex_z_0_0[j] - 0.5 * fl1_fx * tex_z_0_1[j] + ta_xz_0_1[j];

                tey_xxz_0_0[j] = pa_x[j] * tey_xz_0_0[j] - pc_x[j] * tey_xz_0_1[j] + 0.5 * fl1_fx * tey_z_0_0[j] - 0.5 * fl1_fx * tey_z_0_1[j];

                tez_xxz_0_0[j] = pa_x[j] * tez_xz_0_0[j] - pc_x[j] * tez_xz_0_1[j] + 0.5 * fl1_fx * tez_z_0_0[j] - 0.5 * fl1_fx * tez_z_0_1[j];

                tex_xyy_0_0[j] = pa_x[j] * tex_yy_0_0[j] - pc_x[j] * tex_yy_0_1[j] + ta_yy_0_1[j];

                tey_xyy_0_0[j] = pa_x[j] * tey_yy_0_0[j] - pc_x[j] * tey_yy_0_1[j];

                tez_xyy_0_0[j] = pa_x[j] * tez_yy_0_0[j] - pc_x[j] * tez_yy_0_1[j];

                tex_xyz_0_0[j] = pa_x[j] * tex_yz_0_0[j] - pc_x[j] * tex_yz_0_1[j] + ta_yz_0_1[j];

                tey_xyz_0_0[j] = pa_x[j] * tey_yz_0_0[j] - pc_x[j] * tey_yz_0_1[j];

                tez_xyz_0_0[j] = pa_x[j] * tez_yz_0_0[j] - pc_x[j] * tez_yz_0_1[j];

                tex_xzz_0_0[j] = pa_x[j] * tex_zz_0_0[j] - pc_x[j] * tex_zz_0_1[j] + ta_zz_0_1[j];

                tey_xzz_0_0[j] = pa_x[j] * tey_zz_0_0[j] - pc_x[j] * tey_zz_0_1[j];

                tez_xzz_0_0[j] = pa_x[j] * tez_zz_0_0[j] - pc_x[j] * tez_zz_0_1[j];

                tex_yyy_0_0[j] = pa_y[j] * tex_yy_0_0[j] - pc_y[j] * tex_yy_0_1[j] + fl1_fx * tex_y_0_0[j] - fl1_fx * tex_y_0_1[j];

                tey_yyy_0_0[j] = pa_y[j] * tey_yy_0_0[j] - pc_y[j] * tey_yy_0_1[j] + fl1_fx * tey_y_0_0[j] - fl1_fx * tey_y_0_1[j] + ta_yy_0_1[j];

                tez_yyy_0_0[j] = pa_y[j] * tez_yy_0_0[j] - pc_y[j] * tez_yy_0_1[j] + fl1_fx * tez_y_0_0[j] - fl1_fx * tez_y_0_1[j];

                tex_yyz_0_0[j] = pa_y[j] * tex_yz_0_0[j] - pc_y[j] * tex_yz_0_1[j] + 0.5 * fl1_fx * tex_z_0_0[j] - 0.5 * fl1_fx * tex_z_0_1[j];

                tey_yyz_0_0[j] =
                    pa_y[j] * tey_yz_0_0[j] - pc_y[j] * tey_yz_0_1[j] + 0.5 * fl1_fx * tey_z_0_0[j] - 0.5 * fl1_fx * tey_z_0_1[j] + ta_yz_0_1[j];

                tez_yyz_0_0[j] = pa_y[j] * tez_yz_0_0[j] - pc_y[j] * tez_yz_0_1[j] + 0.5 * fl1_fx * tez_z_0_0[j] - 0.5 * fl1_fx * tez_z_0_1[j];

                tex_yzz_0_0[j] = pa_y[j] * tex_zz_0_0[j] - pc_y[j] * tex_zz_0_1[j];

                tey_yzz_0_0[j] = pa_y[j] * tey_zz_0_0[j] - pc_y[j] * tey_zz_0_1[j] + ta_zz_0_1[j];

                tez_yzz_0_0[j] = pa_y[j] * tez_zz_0_0[j] - pc_y[j] * tez_zz_0_1[j];

                tex_zzz_0_0[j] = pa_z[j] * tex_zz_0_0[j] - pc_z[j] * tex_zz_0_1[j] + fl1_fx * tex_z_0_0[j] - fl1_fx * tex_z_0_1[j];

                tey_zzz_0_0[j] = pa_z[j] * tey_zz_0_0[j] - pc_z[j] * tey_zz_0_1[j] + fl1_fx * tey_z_0_0[j] - fl1_fx * tey_z_0_1[j];

                tez_zzz_0_0[j] = pa_z[j] * tez_zz_0_0[j] - pc_z[j] * tez_zz_0_1[j] + fl1_fx * tez_z_0_0[j] - fl1_fx * tez_z_0_1[j] + ta_zz_0_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForSG(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& pbDistances,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_0_4_m0 == -1) continue;

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

            // set up pointers to tensors product of distances R(PC) = P - C

            auto pc_x = pcDistances.data(3 * idx);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_xxx_1, ta_0_xxy_1, ta_0_xxz_1, \
                                         ta_0_xyy_1, ta_0_xyz_1, ta_0_xzz_1, ta_0_yyy_1, ta_0_yyz_1, ta_0_yzz_1, ta_0_zzz_1, \
                                         tex_0_xx_0, tex_0_xx_1, tex_0_xxx_0, tex_0_xxx_1, tex_0_xxxx_0, tex_0_xxxy_0, \
                                         tex_0_xxxz_0, tex_0_xxy_0, tex_0_xxy_1, tex_0_xxyy_0, tex_0_xxyz_0, tex_0_xxz_0, \
                                         tex_0_xxz_1, tex_0_xxzz_0, tex_0_xy_0, tex_0_xy_1, tex_0_xyy_0, tex_0_xyy_1, \
                                         tex_0_xyyy_0, tex_0_xyyz_0, tex_0_xyz_0, tex_0_xyz_1, tex_0_xyzz_0, tex_0_xz_0, \
                                         tex_0_xz_1, tex_0_xzz_0, tex_0_xzz_1, tex_0_xzzz_0, tex_0_yy_0, tex_0_yy_1, \
                                         tex_0_yyy_0, tex_0_yyy_1, tex_0_yyyy_0, tex_0_yyyz_0, tex_0_yyz_0, tex_0_yyz_1, \
                                         tex_0_yyzz_0, tex_0_yz_0, tex_0_yz_1, tex_0_yzz_0, tex_0_yzz_1, tex_0_yzzz_0, \
                                         tex_0_zz_0, tex_0_zz_1, tex_0_zzz_0, tex_0_zzz_1, tex_0_zzzz_0, tey_0_xx_0, \
                                         tey_0_xx_1, tey_0_xxx_0, tey_0_xxx_1, tey_0_xxxx_0, tey_0_xxxy_0, tey_0_xxxz_0, \
                                         tey_0_xxy_0, tey_0_xxy_1, tey_0_xxyy_0, tey_0_xxyz_0, tey_0_xxz_0, tey_0_xxz_1, \
                                         tey_0_xxzz_0, tey_0_xy_0, tey_0_xy_1, tey_0_xyy_0, tey_0_xyy_1, tey_0_xyyy_0, \
                                         tey_0_xyyz_0, tey_0_xyz_0, tey_0_xyz_1, tey_0_xyzz_0, tey_0_xz_0, tey_0_xz_1, \
                                         tey_0_xzz_0, tey_0_xzz_1, tey_0_xzzz_0, tey_0_yy_0, tey_0_yy_1, tey_0_yyy_0, \
                                         tey_0_yyy_1, tey_0_yyyy_0, tey_0_yyyz_0, tey_0_yyz_0, tey_0_yyz_1, tey_0_yyzz_0, \
                                         tey_0_yz_0, tey_0_yz_1, tey_0_yzz_0, tey_0_yzz_1, tey_0_yzzz_0, tey_0_zz_0, \
                                         tey_0_zz_1, tey_0_zzz_0, tey_0_zzz_1, tey_0_zzzz_0, tez_0_xx_0, tez_0_xx_1, \
                                         tez_0_xxx_0, tez_0_xxx_1, tez_0_xxxx_0, tez_0_xxxy_0, tez_0_xxxz_0, tez_0_xxy_0, \
                                         tez_0_xxy_1, tez_0_xxyy_0, tez_0_xxyz_0, tez_0_xxz_0, tez_0_xxz_1, tez_0_xxzz_0, \
                                         tez_0_xy_0, tez_0_xy_1, tez_0_xyy_0, tez_0_xyy_1, tez_0_xyyy_0, tez_0_xyyz_0, \
                                         tez_0_xyz_0, tez_0_xyz_1, tez_0_xyzz_0, tez_0_xz_0, tez_0_xz_1, tez_0_xzz_0, \
                                         tez_0_xzz_1, tez_0_xzzz_0, tez_0_yy_0, tez_0_yy_1, tez_0_yyy_0, tez_0_yyy_1, \
                                         tez_0_yyyy_0, tez_0_yyyz_0, tez_0_yyz_0, tez_0_yyz_1, tez_0_yyzz_0, tez_0_yz_0, \
                                         tez_0_yz_1, tez_0_yzz_0, tez_0_yzz_1, tez_0_yzzz_0, tez_0_zz_0, tez_0_zz_1, \
                                         tez_0_zzz_0, tez_0_zzz_1, tez_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_0_xxxx_0[j] =
                    pb_x[j] * tex_0_xxx_0[j] - pc_x[j] * tex_0_xxx_1[j] + 1.5 * fl1_fx * tex_0_xx_0[j] - 1.5 * fl1_fx * tex_0_xx_1[j] + ta_0_xxx_1[j];

                tey_0_xxxx_0[j] = pb_x[j] * tey_0_xxx_0[j] - pc_x[j] * tey_0_xxx_1[j] + 1.5 * fl1_fx * tey_0_xx_0[j] - 1.5 * fl1_fx * tey_0_xx_1[j];

                tez_0_xxxx_0[j] = pb_x[j] * tez_0_xxx_0[j] - pc_x[j] * tez_0_xxx_1[j] + 1.5 * fl1_fx * tez_0_xx_0[j] - 1.5 * fl1_fx * tez_0_xx_1[j];

                tex_0_xxxy_0[j] =
                    pb_x[j] * tex_0_xxy_0[j] - pc_x[j] * tex_0_xxy_1[j] + fl1_fx * tex_0_xy_0[j] - fl1_fx * tex_0_xy_1[j] + ta_0_xxy_1[j];

                tey_0_xxxy_0[j] = pb_x[j] * tey_0_xxy_0[j] - pc_x[j] * tey_0_xxy_1[j] + fl1_fx * tey_0_xy_0[j] - fl1_fx * tey_0_xy_1[j];

                tez_0_xxxy_0[j] = pb_x[j] * tez_0_xxy_0[j] - pc_x[j] * tez_0_xxy_1[j] + fl1_fx * tez_0_xy_0[j] - fl1_fx * tez_0_xy_1[j];

                tex_0_xxxz_0[j] =
                    pb_x[j] * tex_0_xxz_0[j] - pc_x[j] * tex_0_xxz_1[j] + fl1_fx * tex_0_xz_0[j] - fl1_fx * tex_0_xz_1[j] + ta_0_xxz_1[j];

                tey_0_xxxz_0[j] = pb_x[j] * tey_0_xxz_0[j] - pc_x[j] * tey_0_xxz_1[j] + fl1_fx * tey_0_xz_0[j] - fl1_fx * tey_0_xz_1[j];

                tez_0_xxxz_0[j] = pb_x[j] * tez_0_xxz_0[j] - pc_x[j] * tez_0_xxz_1[j] + fl1_fx * tez_0_xz_0[j] - fl1_fx * tez_0_xz_1[j];

                tex_0_xxyy_0[j] =
                    pb_x[j] * tex_0_xyy_0[j] - pc_x[j] * tex_0_xyy_1[j] + 0.5 * fl1_fx * tex_0_yy_0[j] - 0.5 * fl1_fx * tex_0_yy_1[j] + ta_0_xyy_1[j];

                tey_0_xxyy_0[j] = pb_x[j] * tey_0_xyy_0[j] - pc_x[j] * tey_0_xyy_1[j] + 0.5 * fl1_fx * tey_0_yy_0[j] - 0.5 * fl1_fx * tey_0_yy_1[j];

                tez_0_xxyy_0[j] = pb_x[j] * tez_0_xyy_0[j] - pc_x[j] * tez_0_xyy_1[j] + 0.5 * fl1_fx * tez_0_yy_0[j] - 0.5 * fl1_fx * tez_0_yy_1[j];

                tex_0_xxyz_0[j] =
                    pb_x[j] * tex_0_xyz_0[j] - pc_x[j] * tex_0_xyz_1[j] + 0.5 * fl1_fx * tex_0_yz_0[j] - 0.5 * fl1_fx * tex_0_yz_1[j] + ta_0_xyz_1[j];

                tey_0_xxyz_0[j] = pb_x[j] * tey_0_xyz_0[j] - pc_x[j] * tey_0_xyz_1[j] + 0.5 * fl1_fx * tey_0_yz_0[j] - 0.5 * fl1_fx * tey_0_yz_1[j];

                tez_0_xxyz_0[j] = pb_x[j] * tez_0_xyz_0[j] - pc_x[j] * tez_0_xyz_1[j] + 0.5 * fl1_fx * tez_0_yz_0[j] - 0.5 * fl1_fx * tez_0_yz_1[j];

                tex_0_xxzz_0[j] =
                    pb_x[j] * tex_0_xzz_0[j] - pc_x[j] * tex_0_xzz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j] + ta_0_xzz_1[j];

                tey_0_xxzz_0[j] = pb_x[j] * tey_0_xzz_0[j] - pc_x[j] * tey_0_xzz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j];

                tez_0_xxzz_0[j] = pb_x[j] * tez_0_xzz_0[j] - pc_x[j] * tez_0_xzz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j];

                tex_0_xyyy_0[j] = pb_x[j] * tex_0_yyy_0[j] - pc_x[j] * tex_0_yyy_1[j] + ta_0_yyy_1[j];

                tey_0_xyyy_0[j] = pb_x[j] * tey_0_yyy_0[j] - pc_x[j] * tey_0_yyy_1[j];

                tez_0_xyyy_0[j] = pb_x[j] * tez_0_yyy_0[j] - pc_x[j] * tez_0_yyy_1[j];

                tex_0_xyyz_0[j] = pb_x[j] * tex_0_yyz_0[j] - pc_x[j] * tex_0_yyz_1[j] + ta_0_yyz_1[j];

                tey_0_xyyz_0[j] = pb_x[j] * tey_0_yyz_0[j] - pc_x[j] * tey_0_yyz_1[j];

                tez_0_xyyz_0[j] = pb_x[j] * tez_0_yyz_0[j] - pc_x[j] * tez_0_yyz_1[j];

                tex_0_xyzz_0[j] = pb_x[j] * tex_0_yzz_0[j] - pc_x[j] * tex_0_yzz_1[j] + ta_0_yzz_1[j];

                tey_0_xyzz_0[j] = pb_x[j] * tey_0_yzz_0[j] - pc_x[j] * tey_0_yzz_1[j];

                tez_0_xyzz_0[j] = pb_x[j] * tez_0_yzz_0[j] - pc_x[j] * tez_0_yzz_1[j];

                tex_0_xzzz_0[j] = pb_x[j] * tex_0_zzz_0[j] - pc_x[j] * tex_0_zzz_1[j] + ta_0_zzz_1[j];

                tey_0_xzzz_0[j] = pb_x[j] * tey_0_zzz_0[j] - pc_x[j] * tey_0_zzz_1[j];

                tez_0_xzzz_0[j] = pb_x[j] * tez_0_zzz_0[j] - pc_x[j] * tez_0_zzz_1[j];

                tex_0_yyyy_0[j] = pb_y[j] * tex_0_yyy_0[j] - pc_y[j] * tex_0_yyy_1[j] + 1.5 * fl1_fx * tex_0_yy_0[j] - 1.5 * fl1_fx * tex_0_yy_1[j];

                tey_0_yyyy_0[j] =
                    pb_y[j] * tey_0_yyy_0[j] - pc_y[j] * tey_0_yyy_1[j] + 1.5 * fl1_fx * tey_0_yy_0[j] - 1.5 * fl1_fx * tey_0_yy_1[j] + ta_0_yyy_1[j];

                tez_0_yyyy_0[j] = pb_y[j] * tez_0_yyy_0[j] - pc_y[j] * tez_0_yyy_1[j] + 1.5 * fl1_fx * tez_0_yy_0[j] - 1.5 * fl1_fx * tez_0_yy_1[j];

                tex_0_yyyz_0[j] = pb_y[j] * tex_0_yyz_0[j] - pc_y[j] * tex_0_yyz_1[j] + fl1_fx * tex_0_yz_0[j] - fl1_fx * tex_0_yz_1[j];

                tey_0_yyyz_0[j] =
                    pb_y[j] * tey_0_yyz_0[j] - pc_y[j] * tey_0_yyz_1[j] + fl1_fx * tey_0_yz_0[j] - fl1_fx * tey_0_yz_1[j] + ta_0_yyz_1[j];

                tez_0_yyyz_0[j] = pb_y[j] * tez_0_yyz_0[j] - pc_y[j] * tez_0_yyz_1[j] + fl1_fx * tez_0_yz_0[j] - fl1_fx * tez_0_yz_1[j];

                tex_0_yyzz_0[j] = pb_y[j] * tex_0_yzz_0[j] - pc_y[j] * tex_0_yzz_1[j] + 0.5 * fl1_fx * tex_0_zz_0[j] - 0.5 * fl1_fx * tex_0_zz_1[j];

                tey_0_yyzz_0[j] =
                    pb_y[j] * tey_0_yzz_0[j] - pc_y[j] * tey_0_yzz_1[j] + 0.5 * fl1_fx * tey_0_zz_0[j] - 0.5 * fl1_fx * tey_0_zz_1[j] + ta_0_yzz_1[j];

                tez_0_yyzz_0[j] = pb_y[j] * tez_0_yzz_0[j] - pc_y[j] * tez_0_yzz_1[j] + 0.5 * fl1_fx * tez_0_zz_0[j] - 0.5 * fl1_fx * tez_0_zz_1[j];

                tex_0_yzzz_0[j] = pb_y[j] * tex_0_zzz_0[j] - pc_y[j] * tex_0_zzz_1[j];

                tey_0_yzzz_0[j] = pb_y[j] * tey_0_zzz_0[j] - pc_y[j] * tey_0_zzz_1[j] + ta_0_zzz_1[j];

                tez_0_yzzz_0[j] = pb_y[j] * tez_0_zzz_0[j] - pc_y[j] * tez_0_zzz_1[j];

                tex_0_zzzz_0[j] = pb_z[j] * tex_0_zzz_0[j] - pc_z[j] * tex_0_zzz_1[j] + 1.5 * fl1_fx * tex_0_zz_0[j] - 1.5 * fl1_fx * tex_0_zz_1[j];

                tey_0_zzzz_0[j] = pb_z[j] * tey_0_zzz_0[j] - pc_z[j] * tey_0_zzz_1[j] + 1.5 * fl1_fx * tey_0_zz_0[j] - 1.5 * fl1_fx * tey_0_zz_1[j];

                tez_0_zzzz_0[j] =
                    pb_z[j] * tez_0_zzz_0[j] - pc_z[j] * tez_0_zzz_1[j] + 1.5 * fl1_fx * tez_0_zz_0[j] - 1.5 * fl1_fx * tez_0_zz_1[j] + ta_0_zzz_1[j];
            }

            idx++;
        }
    }
}

void
compElectricFieldForGS(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Electric Field"}, {4, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_e_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {4, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_e_4_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_e_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_3_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

        auto pidx_e_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_e_2_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto tex_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * idx + 5);

            auto tey_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_0 = primBuffer.data(pidx_e_2_0_m0 + 12 * bdim + 6 * idx + 5);

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

            auto tex_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * idx + 5);

            auto tey_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 6 * bdim + 6 * idx + 5);

            auto tez_zz_0_1 = primBuffer.data(pidx_e_2_0_m1 + 12 * bdim + 6 * idx + 5);

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

            auto tex_xxxx_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx);

            auto tey_xxxx_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx);

            auto tez_xxxx_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx);

            auto tex_xxxy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 1);

            auto tey_xxxy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 1);

            auto tez_xxxy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 1);

            auto tex_xxxz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 2);

            auto tey_xxxz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 2);

            auto tez_xxxz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 2);

            auto tex_xxyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 3);

            auto tey_xxyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 3);

            auto tez_xxyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 3);

            auto tex_xxyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 4);

            auto tey_xxyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 4);

            auto tez_xxyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 4);

            auto tex_xxzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 5);

            auto tey_xxzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 5);

            auto tez_xxzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 5);

            auto tex_xyyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 6);

            auto tey_xyyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 6);

            auto tez_xyyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 6);

            auto tex_xyyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 7);

            auto tey_xyyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 7);

            auto tez_xyyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 7);

            auto tex_xyzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 8);

            auto tey_xyzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 8);

            auto tez_xyzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 8);

            auto tex_xzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 9);

            auto tey_xzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 9);

            auto tez_xzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 9);

            auto tex_yyyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 10);

            auto tey_yyyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 10);

            auto tez_yyyy_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 10);

            auto tex_yyyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 11);

            auto tey_yyyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 11);

            auto tez_yyyz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 11);

            auto tex_yyzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 12);

            auto tey_yyzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 12);

            auto tez_yyzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 12);

            auto tex_yzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 13);

            auto tey_yzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 13);

            auto tez_yzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 13);

            auto tex_zzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * idx + 14);

            auto tey_zzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 15 * bdim + 15 * idx + 14);

            auto tez_zzzz_0_0 = primBuffer.data(pidx_e_4_0_m0 + 30 * bdim + 15 * idx + 14);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xxx_0_1, ta_xxy_0_1, ta_xxz_0_1, \
                                         ta_xyy_0_1, ta_xyz_0_1, ta_xzz_0_1, ta_yyy_0_1, ta_yyz_0_1, ta_yzz_0_1, ta_zzz_0_1, \
                                         tex_xx_0_0, tex_xx_0_1, tex_xxx_0_0, tex_xxx_0_1, tex_xxxx_0_0, tex_xxxy_0_0, \
                                         tex_xxxz_0_0, tex_xxy_0_0, tex_xxy_0_1, tex_xxyy_0_0, tex_xxyz_0_0, tex_xxz_0_0, \
                                         tex_xxz_0_1, tex_xxzz_0_0, tex_xy_0_0, tex_xy_0_1, tex_xyy_0_0, tex_xyy_0_1, \
                                         tex_xyyy_0_0, tex_xyyz_0_0, tex_xyz_0_0, tex_xyz_0_1, tex_xyzz_0_0, tex_xz_0_0, \
                                         tex_xz_0_1, tex_xzz_0_0, tex_xzz_0_1, tex_xzzz_0_0, tex_yy_0_0, tex_yy_0_1, \
                                         tex_yyy_0_0, tex_yyy_0_1, tex_yyyy_0_0, tex_yyyz_0_0, tex_yyz_0_0, tex_yyz_0_1, \
                                         tex_yyzz_0_0, tex_yz_0_0, tex_yz_0_1, tex_yzz_0_0, tex_yzz_0_1, tex_yzzz_0_0, \
                                         tex_zz_0_0, tex_zz_0_1, tex_zzz_0_0, tex_zzz_0_1, tex_zzzz_0_0, tey_xx_0_0, \
                                         tey_xx_0_1, tey_xxx_0_0, tey_xxx_0_1, tey_xxxx_0_0, tey_xxxy_0_0, tey_xxxz_0_0, \
                                         tey_xxy_0_0, tey_xxy_0_1, tey_xxyy_0_0, tey_xxyz_0_0, tey_xxz_0_0, tey_xxz_0_1, \
                                         tey_xxzz_0_0, tey_xy_0_0, tey_xy_0_1, tey_xyy_0_0, tey_xyy_0_1, tey_xyyy_0_0, \
                                         tey_xyyz_0_0, tey_xyz_0_0, tey_xyz_0_1, tey_xyzz_0_0, tey_xz_0_0, tey_xz_0_1, \
                                         tey_xzz_0_0, tey_xzz_0_1, tey_xzzz_0_0, tey_yy_0_0, tey_yy_0_1, tey_yyy_0_0, \
                                         tey_yyy_0_1, tey_yyyy_0_0, tey_yyyz_0_0, tey_yyz_0_0, tey_yyz_0_1, tey_yyzz_0_0, \
                                         tey_yz_0_0, tey_yz_0_1, tey_yzz_0_0, tey_yzz_0_1, tey_yzzz_0_0, tey_zz_0_0, \
                                         tey_zz_0_1, tey_zzz_0_0, tey_zzz_0_1, tey_zzzz_0_0, tez_xx_0_0, tez_xx_0_1, \
                                         tez_xxx_0_0, tez_xxx_0_1, tez_xxxx_0_0, tez_xxxy_0_0, tez_xxxz_0_0, tez_xxy_0_0, \
                                         tez_xxy_0_1, tez_xxyy_0_0, tez_xxyz_0_0, tez_xxz_0_0, tez_xxz_0_1, tez_xxzz_0_0, \
                                         tez_xy_0_0, tez_xy_0_1, tez_xyy_0_0, tez_xyy_0_1, tez_xyyy_0_0, tez_xyyz_0_0, \
                                         tez_xyz_0_0, tez_xyz_0_1, tez_xyzz_0_0, tez_xz_0_0, tez_xz_0_1, tez_xzz_0_0, \
                                         tez_xzz_0_1, tez_xzzz_0_0, tez_yy_0_0, tez_yy_0_1, tez_yyy_0_0, tez_yyy_0_1, \
                                         tez_yyyy_0_0, tez_yyyz_0_0, tez_yyz_0_0, tez_yyz_0_1, tez_yyzz_0_0, tez_yz_0_0, \
                                         tez_yz_0_1, tez_yzz_0_0, tez_yzz_0_1, tez_yzzz_0_0, tez_zz_0_0, tez_zz_0_1, \
                                         tez_zzz_0_0, tez_zzz_0_1, tez_zzzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                tex_xxxx_0_0[j] =
                    pa_x[j] * tex_xxx_0_0[j] - pc_x[j] * tex_xxx_0_1[j] + 1.5 * fl1_fx * tex_xx_0_0[j] - 1.5 * fl1_fx * tex_xx_0_1[j] + ta_xxx_0_1[j];

                tey_xxxx_0_0[j] = pa_x[j] * tey_xxx_0_0[j] - pc_x[j] * tey_xxx_0_1[j] + 1.5 * fl1_fx * tey_xx_0_0[j] - 1.5 * fl1_fx * tey_xx_0_1[j];

                tez_xxxx_0_0[j] = pa_x[j] * tez_xxx_0_0[j] - pc_x[j] * tez_xxx_0_1[j] + 1.5 * fl1_fx * tez_xx_0_0[j] - 1.5 * fl1_fx * tez_xx_0_1[j];

                tex_xxxy_0_0[j] =
                    pa_x[j] * tex_xxy_0_0[j] - pc_x[j] * tex_xxy_0_1[j] + fl1_fx * tex_xy_0_0[j] - fl1_fx * tex_xy_0_1[j] + ta_xxy_0_1[j];

                tey_xxxy_0_0[j] = pa_x[j] * tey_xxy_0_0[j] - pc_x[j] * tey_xxy_0_1[j] + fl1_fx * tey_xy_0_0[j] - fl1_fx * tey_xy_0_1[j];

                tez_xxxy_0_0[j] = pa_x[j] * tez_xxy_0_0[j] - pc_x[j] * tez_xxy_0_1[j] + fl1_fx * tez_xy_0_0[j] - fl1_fx * tez_xy_0_1[j];

                tex_xxxz_0_0[j] =
                    pa_x[j] * tex_xxz_0_0[j] - pc_x[j] * tex_xxz_0_1[j] + fl1_fx * tex_xz_0_0[j] - fl1_fx * tex_xz_0_1[j] + ta_xxz_0_1[j];

                tey_xxxz_0_0[j] = pa_x[j] * tey_xxz_0_0[j] - pc_x[j] * tey_xxz_0_1[j] + fl1_fx * tey_xz_0_0[j] - fl1_fx * tey_xz_0_1[j];

                tez_xxxz_0_0[j] = pa_x[j] * tez_xxz_0_0[j] - pc_x[j] * tez_xxz_0_1[j] + fl1_fx * tez_xz_0_0[j] - fl1_fx * tez_xz_0_1[j];

                tex_xxyy_0_0[j] =
                    pa_x[j] * tex_xyy_0_0[j] - pc_x[j] * tex_xyy_0_1[j] + 0.5 * fl1_fx * tex_yy_0_0[j] - 0.5 * fl1_fx * tex_yy_0_1[j] + ta_xyy_0_1[j];

                tey_xxyy_0_0[j] = pa_x[j] * tey_xyy_0_0[j] - pc_x[j] * tey_xyy_0_1[j] + 0.5 * fl1_fx * tey_yy_0_0[j] - 0.5 * fl1_fx * tey_yy_0_1[j];

                tez_xxyy_0_0[j] = pa_x[j] * tez_xyy_0_0[j] - pc_x[j] * tez_xyy_0_1[j] + 0.5 * fl1_fx * tez_yy_0_0[j] - 0.5 * fl1_fx * tez_yy_0_1[j];

                tex_xxyz_0_0[j] =
                    pa_x[j] * tex_xyz_0_0[j] - pc_x[j] * tex_xyz_0_1[j] + 0.5 * fl1_fx * tex_yz_0_0[j] - 0.5 * fl1_fx * tex_yz_0_1[j] + ta_xyz_0_1[j];

                tey_xxyz_0_0[j] = pa_x[j] * tey_xyz_0_0[j] - pc_x[j] * tey_xyz_0_1[j] + 0.5 * fl1_fx * tey_yz_0_0[j] - 0.5 * fl1_fx * tey_yz_0_1[j];

                tez_xxyz_0_0[j] = pa_x[j] * tez_xyz_0_0[j] - pc_x[j] * tez_xyz_0_1[j] + 0.5 * fl1_fx * tez_yz_0_0[j] - 0.5 * fl1_fx * tez_yz_0_1[j];

                tex_xxzz_0_0[j] =
                    pa_x[j] * tex_xzz_0_0[j] - pc_x[j] * tex_xzz_0_1[j] + 0.5 * fl1_fx * tex_zz_0_0[j] - 0.5 * fl1_fx * tex_zz_0_1[j] + ta_xzz_0_1[j];

                tey_xxzz_0_0[j] = pa_x[j] * tey_xzz_0_0[j] - pc_x[j] * tey_xzz_0_1[j] + 0.5 * fl1_fx * tey_zz_0_0[j] - 0.5 * fl1_fx * tey_zz_0_1[j];

                tez_xxzz_0_0[j] = pa_x[j] * tez_xzz_0_0[j] - pc_x[j] * tez_xzz_0_1[j] + 0.5 * fl1_fx * tez_zz_0_0[j] - 0.5 * fl1_fx * tez_zz_0_1[j];

                tex_xyyy_0_0[j] = pa_x[j] * tex_yyy_0_0[j] - pc_x[j] * tex_yyy_0_1[j] + ta_yyy_0_1[j];

                tey_xyyy_0_0[j] = pa_x[j] * tey_yyy_0_0[j] - pc_x[j] * tey_yyy_0_1[j];

                tez_xyyy_0_0[j] = pa_x[j] * tez_yyy_0_0[j] - pc_x[j] * tez_yyy_0_1[j];

                tex_xyyz_0_0[j] = pa_x[j] * tex_yyz_0_0[j] - pc_x[j] * tex_yyz_0_1[j] + ta_yyz_0_1[j];

                tey_xyyz_0_0[j] = pa_x[j] * tey_yyz_0_0[j] - pc_x[j] * tey_yyz_0_1[j];

                tez_xyyz_0_0[j] = pa_x[j] * tez_yyz_0_0[j] - pc_x[j] * tez_yyz_0_1[j];

                tex_xyzz_0_0[j] = pa_x[j] * tex_yzz_0_0[j] - pc_x[j] * tex_yzz_0_1[j] + ta_yzz_0_1[j];

                tey_xyzz_0_0[j] = pa_x[j] * tey_yzz_0_0[j] - pc_x[j] * tey_yzz_0_1[j];

                tez_xyzz_0_0[j] = pa_x[j] * tez_yzz_0_0[j] - pc_x[j] * tez_yzz_0_1[j];

                tex_xzzz_0_0[j] = pa_x[j] * tex_zzz_0_0[j] - pc_x[j] * tex_zzz_0_1[j] + ta_zzz_0_1[j];

                tey_xzzz_0_0[j] = pa_x[j] * tey_zzz_0_0[j] - pc_x[j] * tey_zzz_0_1[j];

                tez_xzzz_0_0[j] = pa_x[j] * tez_zzz_0_0[j] - pc_x[j] * tez_zzz_0_1[j];

                tex_yyyy_0_0[j] = pa_y[j] * tex_yyy_0_0[j] - pc_y[j] * tex_yyy_0_1[j] + 1.5 * fl1_fx * tex_yy_0_0[j] - 1.5 * fl1_fx * tex_yy_0_1[j];

                tey_yyyy_0_0[j] =
                    pa_y[j] * tey_yyy_0_0[j] - pc_y[j] * tey_yyy_0_1[j] + 1.5 * fl1_fx * tey_yy_0_0[j] - 1.5 * fl1_fx * tey_yy_0_1[j] + ta_yyy_0_1[j];

                tez_yyyy_0_0[j] = pa_y[j] * tez_yyy_0_0[j] - pc_y[j] * tez_yyy_0_1[j] + 1.5 * fl1_fx * tez_yy_0_0[j] - 1.5 * fl1_fx * tez_yy_0_1[j];

                tex_yyyz_0_0[j] = pa_y[j] * tex_yyz_0_0[j] - pc_y[j] * tex_yyz_0_1[j] + fl1_fx * tex_yz_0_0[j] - fl1_fx * tex_yz_0_1[j];

                tey_yyyz_0_0[j] =
                    pa_y[j] * tey_yyz_0_0[j] - pc_y[j] * tey_yyz_0_1[j] + fl1_fx * tey_yz_0_0[j] - fl1_fx * tey_yz_0_1[j] + ta_yyz_0_1[j];

                tez_yyyz_0_0[j] = pa_y[j] * tez_yyz_0_0[j] - pc_y[j] * tez_yyz_0_1[j] + fl1_fx * tez_yz_0_0[j] - fl1_fx * tez_yz_0_1[j];

                tex_yyzz_0_0[j] = pa_y[j] * tex_yzz_0_0[j] - pc_y[j] * tex_yzz_0_1[j] + 0.5 * fl1_fx * tex_zz_0_0[j] - 0.5 * fl1_fx * tex_zz_0_1[j];

                tey_yyzz_0_0[j] =
                    pa_y[j] * tey_yzz_0_0[j] - pc_y[j] * tey_yzz_0_1[j] + 0.5 * fl1_fx * tey_zz_0_0[j] - 0.5 * fl1_fx * tey_zz_0_1[j] + ta_yzz_0_1[j];

                tez_yyzz_0_0[j] = pa_y[j] * tez_yzz_0_0[j] - pc_y[j] * tez_yzz_0_1[j] + 0.5 * fl1_fx * tez_zz_0_0[j] - 0.5 * fl1_fx * tez_zz_0_1[j];

                tex_yzzz_0_0[j] = pa_y[j] * tex_zzz_0_0[j] - pc_y[j] * tex_zzz_0_1[j];

                tey_yzzz_0_0[j] = pa_y[j] * tey_zzz_0_0[j] - pc_y[j] * tey_zzz_0_1[j] + ta_zzz_0_1[j];

                tez_yzzz_0_0[j] = pa_y[j] * tez_zzz_0_0[j] - pc_y[j] * tez_zzz_0_1[j];

                tex_zzzz_0_0[j] = pa_z[j] * tex_zzz_0_0[j] - pc_z[j] * tex_zzz_0_1[j] + 1.5 * fl1_fx * tex_zz_0_0[j] - 1.5 * fl1_fx * tex_zz_0_1[j];

                tey_zzzz_0_0[j] = pa_z[j] * tey_zzz_0_0[j] - pc_z[j] * tey_zzz_0_1[j] + 1.5 * fl1_fx * tey_zz_0_0[j] - 1.5 * fl1_fx * tey_zz_0_1[j];

                tez_zzzz_0_0[j] =
                    pa_z[j] * tez_zzz_0_0[j] - pc_z[j] * tez_zzz_0_1[j] + 1.5 * fl1_fx * tez_zz_0_0[j] - 1.5 * fl1_fx * tez_zz_0_1[j] + ta_zzz_0_1[j];
            }

            idx++;
        }
    }
}

}  // namespace efieldrecfunc
