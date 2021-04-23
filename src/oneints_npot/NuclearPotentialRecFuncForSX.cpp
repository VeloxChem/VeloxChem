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

#include "NuclearPotentialRecFuncForSX.hpp"

#include "MathConst.hpp"

namespace npotrecfunc {  // npotrecfunc namespace

void
compNuclearPotentialForSS(CMemBlock2D<double>&       primBuffer,
                          const CRecursionMap&       recursionMap,
                          const CBoysFunction&       bfTable,
                          CMemBlock<double>&         bfArguments,
                          CMemBlock2D<double>&       bfValues,
                          const int32_t              bfOrder,
                          const CMemBlock2D<double>& osFactors,
                          const CMemBlock2D<double>& abDistances,
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

    // set up pointer to overlap integrals

    auto sidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up Obara-Saika prefactors

        auto fg = osFactors.data(3 * idx + 2);

        // set up pointers to ditances R(PC)

        auto pcx = pcDistances.data(3 * idx);

        auto pcy = pcDistances.data(3 * idx + 1);

        auto pcz = pcDistances.data(3 * idx + 2);

        // compute Boys function argument

        auto fargs = bfArguments.data();

        #pragma omp simd aligned(fargs, fg, pcx, pcy, pcz: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            fargs[j] = fg[j] * (pcx[j] * pcx[j] + pcy[j] * pcy[j] +

                                pcz[j] * pcz[j]);
        }

        // evaluate Boys function values

        bfTable.compute(bfValues, bfArguments, bfOrder);

        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(3 * idx);

        // fetch up pi value

        auto fpi = 1.0 / mathconst::getPiValue();

        // compute overlap scaling factor

        #pragma omp simd aligned(fx, fargs: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            fargs[j] = 2.0 * std::sqrt(fpi / fx[j]);
        }

        // distribute (s|A(0)|s) integrals

        auto s_0_0 = primBuffer.data(sidx + idx);

        for (int32_t j = 0; j <= bfOrder; j++)
        {
            auto pidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, j));

            if (pidx == -1) continue;

            auto t_0_0 = primBuffer.data(pidx + idx);

            auto bvals = bfValues.data(j);

            #pragma omp simd aligned(t_0_0, bvals, fargs, s_0_0: VLX_ALIGN)
            for (int32_t k = 0; k < nprim; k++)
            {
                t_0_0[k] = fargs[k] * s_0_0[k] * bvals[k];
            }
        }

        idx++;
    }
}

void
compNuclearPotentialForSP(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_0_1_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_0_1_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_0_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

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

            auto ta_0_0_0 = primBuffer.data(pidx_a_0_0_m0 + idx);

            auto ta_0_0_1 = primBuffer.data(pidx_a_0_0_m1 + idx);

            // set up pointers to integrals

            auto ta_0_x_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx);

            auto ta_0_y_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 1);

            auto ta_0_z_0 = primBuffer.data(pidx_a_0_1_m0 + 3 * idx + 2);

            #pragma omp simd aligned(pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_0_0, ta_0_0_1, ta_0_x_0, ta_0_y_0, \
                                         ta_0_z_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                ta_0_x_0[j] = pb_x[j] * ta_0_0_0[j] - pc_x[j] * ta_0_0_1[j];

                ta_0_y_0[j] = pb_y[j] * ta_0_0_0[j] - pc_y[j] * ta_0_0_1[j];

                ta_0_z_0[j] = pb_z[j] * ta_0_0_0[j] - pc_z[j] * ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForPS(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_1_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_1_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_0_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

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

            auto ta_0_0_0 = primBuffer.data(pidx_a_0_0_m0 + idx);

            auto ta_0_0_1 = primBuffer.data(pidx_a_0_0_m1 + idx);

            // set up pointers to integrals

            auto ta_x_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx);

            auto ta_y_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 1);

            auto ta_z_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 2);

            #pragma omp simd aligned(pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_0_0, ta_0_0_1, ta_x_0_0, ta_y_0_0, \
                                         ta_z_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                ta_x_0_0[j] = pa_x[j] * ta_0_0_0[j] - pc_x[j] * ta_0_0_1[j];

                ta_y_0_0[j] = pa_y[j] * ta_0_0_0[j] - pc_y[j] * ta_0_0_1[j];

                ta_z_0_0[j] = pa_z[j] * ta_0_0_0[j] - pc_z[j] * ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForSD(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_0_2_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_0_2_m0 == -1) continue;

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            auto ta_0_xx_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx);

            auto ta_0_xy_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 1);

            auto ta_0_xz_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 2);

            auto ta_0_yy_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 3);

            auto ta_0_yz_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 4);

            auto ta_0_zz_0 = primBuffer.data(pidx_a_0_2_m0 + 6 * idx + 5);

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_0_0, ta_0_0_1, ta_0_x_0, ta_0_x_1, \
                                         ta_0_xx_0, ta_0_xy_0, ta_0_xz_0, ta_0_y_0, ta_0_y_1, ta_0_yy_0, ta_0_yz_0, ta_0_z_0, \
                                         ta_0_z_1, ta_0_zz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_0_xx_0[j] = pb_x[j] * ta_0_x_0[j] - pc_x[j] * ta_0_x_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];

                ta_0_xy_0[j] = pb_x[j] * ta_0_y_0[j] - pc_x[j] * ta_0_y_1[j];

                ta_0_xz_0[j] = pb_x[j] * ta_0_z_0[j] - pc_x[j] * ta_0_z_1[j];

                ta_0_yy_0[j] = pb_y[j] * ta_0_y_0[j] - pc_y[j] * ta_0_y_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];

                ta_0_yz_0[j] = pb_y[j] * ta_0_z_0[j] - pc_y[j] * ta_0_z_1[j];

                ta_0_zz_0[j] = pb_z[j] * ta_0_z_0[j] - pc_z[j] * ta_0_z_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForDS(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_2_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_2_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_1_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_1_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_x_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx);

            auto ta_y_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 1);

            auto ta_z_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 2);

            auto ta_x_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx);

            auto ta_y_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 1);

            auto ta_z_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 2);

            auto ta_0_0_0 = primBuffer.data(pidx_a_0_0_m0 + idx);

            auto ta_0_0_1 = primBuffer.data(pidx_a_0_0_m1 + idx);

            // set up pointers to integrals

            auto ta_xx_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx);

            auto ta_xy_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 1);

            auto ta_xz_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 2);

            auto ta_yy_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 3);

            auto ta_yz_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 4);

            auto ta_zz_0_0 = primBuffer.data(pidx_a_2_0_m0 + 6 * idx + 5);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_0_0_0, ta_0_0_1, ta_x_0_0, ta_x_0_1, \
                                         ta_xx_0_0, ta_xy_0_0, ta_xz_0_0, ta_y_0_0, ta_y_0_1, ta_yy_0_0, ta_yz_0_0, ta_z_0_0, \
                                         ta_z_0_1, ta_zz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xx_0_0[j] = pa_x[j] * ta_x_0_0[j] - pc_x[j] * ta_x_0_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];

                ta_xy_0_0[j] = pa_x[j] * ta_y_0_0[j] - pc_x[j] * ta_y_0_1[j];

                ta_xz_0_0[j] = pa_x[j] * ta_z_0_0[j] - pc_x[j] * ta_z_0_1[j];

                ta_yy_0_0[j] = pa_y[j] * ta_y_0_0[j] - pc_y[j] * ta_y_0_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];

                ta_yz_0_0[j] = pa_y[j] * ta_z_0_0[j] - pc_y[j] * ta_z_0_1[j];

                ta_zz_0_0[j] = pa_z[j] * ta_z_0_0[j] - pc_z[j] * ta_z_0_1[j] + 0.5 * fl1_fx * ta_0_0_0[j] - 0.5 * fl1_fx * ta_0_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForSF(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_0_3_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_0_3_m0 == -1) continue;

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_x_0, ta_0_x_1, ta_0_xx_0, ta_0_xx_1, \
                                         ta_0_xxx_0, ta_0_xxy_0, ta_0_xxz_0, ta_0_xy_0, ta_0_xy_1, ta_0_xyy_0, ta_0_xyz_0, \
                                         ta_0_xz_0, ta_0_xz_1, ta_0_xzz_0, ta_0_y_0, ta_0_y_1, ta_0_yy_0, ta_0_yy_1, \
                                         ta_0_yyy_0, ta_0_yyz_0, ta_0_yz_0, ta_0_yz_1, ta_0_yzz_0, ta_0_z_0, ta_0_z_1, \
                                         ta_0_zz_0, ta_0_zz_1, ta_0_zzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_0_xxx_0[j] = pb_x[j] * ta_0_xx_0[j] - pc_x[j] * ta_0_xx_1[j] + fl1_fx * ta_0_x_0[j] - fl1_fx * ta_0_x_1[j];

                ta_0_xxy_0[j] = pb_x[j] * ta_0_xy_0[j] - pc_x[j] * ta_0_xy_1[j] + 0.5 * fl1_fx * ta_0_y_0[j] - 0.5 * fl1_fx * ta_0_y_1[j];

                ta_0_xxz_0[j] = pb_x[j] * ta_0_xz_0[j] - pc_x[j] * ta_0_xz_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j];

                ta_0_xyy_0[j] = pb_x[j] * ta_0_yy_0[j] - pc_x[j] * ta_0_yy_1[j];

                ta_0_xyz_0[j] = pb_x[j] * ta_0_yz_0[j] - pc_x[j] * ta_0_yz_1[j];

                ta_0_xzz_0[j] = pb_x[j] * ta_0_zz_0[j] - pc_x[j] * ta_0_zz_1[j];

                ta_0_yyy_0[j] = pb_y[j] * ta_0_yy_0[j] - pc_y[j] * ta_0_yy_1[j] + fl1_fx * ta_0_y_0[j] - fl1_fx * ta_0_y_1[j];

                ta_0_yyz_0[j] = pb_y[j] * ta_0_yz_0[j] - pc_y[j] * ta_0_yz_1[j] + 0.5 * fl1_fx * ta_0_z_0[j] - 0.5 * fl1_fx * ta_0_z_1[j];

                ta_0_yzz_0[j] = pb_y[j] * ta_0_zz_0[j] - pc_y[j] * ta_0_zz_1[j];

                ta_0_zzz_0[j] = pb_z[j] * ta_0_zz_0[j] - pc_z[j] * ta_0_zz_1[j] + fl1_fx * ta_0_z_0[j] - fl1_fx * ta_0_z_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForFS(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_3_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_3_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_2_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_2_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_x_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx);

            auto ta_y_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 1);

            auto ta_z_0_0 = primBuffer.data(pidx_a_1_0_m0 + 3 * idx + 2);

            auto ta_x_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx);

            auto ta_y_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 1);

            auto ta_z_0_1 = primBuffer.data(pidx_a_1_0_m1 + 3 * idx + 2);

            // set up pointers to integrals

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

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_x_0_0, ta_x_0_1, ta_xx_0_0, ta_xx_0_1, \
                                         ta_xxx_0_0, ta_xxy_0_0, ta_xxz_0_0, ta_xy_0_0, ta_xy_0_1, ta_xyy_0_0, ta_xyz_0_0, \
                                         ta_xz_0_0, ta_xz_0_1, ta_xzz_0_0, ta_y_0_0, ta_y_0_1, ta_yy_0_0, ta_yy_0_1, \
                                         ta_yyy_0_0, ta_yyz_0_0, ta_yz_0_0, ta_yz_0_1, ta_yzz_0_0, ta_z_0_0, ta_z_0_1, \
                                         ta_zz_0_0, ta_zz_0_1, ta_zzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxx_0_0[j] = pa_x[j] * ta_xx_0_0[j] - pc_x[j] * ta_xx_0_1[j] + fl1_fx * ta_x_0_0[j] - fl1_fx * ta_x_0_1[j];

                ta_xxy_0_0[j] = pa_x[j] * ta_xy_0_0[j] - pc_x[j] * ta_xy_0_1[j] + 0.5 * fl1_fx * ta_y_0_0[j] - 0.5 * fl1_fx * ta_y_0_1[j];

                ta_xxz_0_0[j] = pa_x[j] * ta_xz_0_0[j] - pc_x[j] * ta_xz_0_1[j] + 0.5 * fl1_fx * ta_z_0_0[j] - 0.5 * fl1_fx * ta_z_0_1[j];

                ta_xyy_0_0[j] = pa_x[j] * ta_yy_0_0[j] - pc_x[j] * ta_yy_0_1[j];

                ta_xyz_0_0[j] = pa_x[j] * ta_yz_0_0[j] - pc_x[j] * ta_yz_0_1[j];

                ta_xzz_0_0[j] = pa_x[j] * ta_zz_0_0[j] - pc_x[j] * ta_zz_0_1[j];

                ta_yyy_0_0[j] = pa_y[j] * ta_yy_0_0[j] - pc_y[j] * ta_yy_0_1[j] + fl1_fx * ta_y_0_0[j] - fl1_fx * ta_y_0_1[j];

                ta_yyz_0_0[j] = pa_y[j] * ta_yz_0_0[j] - pc_y[j] * ta_yz_0_1[j] + 0.5 * fl1_fx * ta_z_0_0[j] - 0.5 * fl1_fx * ta_z_0_1[j];

                ta_yzz_0_0[j] = pa_y[j] * ta_zz_0_0[j] - pc_y[j] * ta_zz_0_1[j];

                ta_zzz_0_0[j] = pa_z[j] * ta_zz_0_0[j] - pc_z[j] * ta_zz_0_1[j] + fl1_fx * ta_z_0_0[j] - fl1_fx * ta_z_0_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForSG(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up maximum order of integral

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_0_4_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_0_4_m0 == -1) continue;

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

            // set up pointers to tensors product of distances R(PB) = P - B

            auto pb_x = pbDistances.data(3 * idx);

            auto pb_y = pbDistances.data(3 * idx + 1);

            auto pb_z = pbDistances.data(3 * idx + 2);

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

            #pragma omp simd aligned(fx, pb_x, pb_y, pb_z, pc_x, pc_y, pc_z, ta_0_xx_0, ta_0_xx_1, ta_0_xxx_0, \
                                         ta_0_xxx_1, ta_0_xxxx_0, ta_0_xxxy_0, ta_0_xxxz_0, ta_0_xxy_0, ta_0_xxy_1, \
                                         ta_0_xxyy_0, ta_0_xxyz_0, ta_0_xxz_0, ta_0_xxz_1, ta_0_xxzz_0, ta_0_xy_0, ta_0_xy_1, \
                                         ta_0_xyy_0, ta_0_xyy_1, ta_0_xyyy_0, ta_0_xyyz_0, ta_0_xyz_0, ta_0_xyz_1, \
                                         ta_0_xyzz_0, ta_0_xz_0, ta_0_xz_1, ta_0_xzz_0, ta_0_xzz_1, ta_0_xzzz_0, ta_0_yy_0, \
                                         ta_0_yy_1, ta_0_yyy_0, ta_0_yyy_1, ta_0_yyyy_0, ta_0_yyyz_0, ta_0_yyz_0, \
                                         ta_0_yyz_1, ta_0_yyzz_0, ta_0_yz_0, ta_0_yz_1, ta_0_yzz_0, ta_0_yzz_1, ta_0_yzzz_0, \
                                         ta_0_zz_0, ta_0_zz_1, ta_0_zzz_0, ta_0_zzz_1, ta_0_zzzz_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_0_xxxx_0[j] = pb_x[j] * ta_0_xxx_0[j] - pc_x[j] * ta_0_xxx_1[j] + 1.5 * fl1_fx * ta_0_xx_0[j] - 1.5 * fl1_fx * ta_0_xx_1[j];

                ta_0_xxxy_0[j] = pb_x[j] * ta_0_xxy_0[j] - pc_x[j] * ta_0_xxy_1[j] + fl1_fx * ta_0_xy_0[j] - fl1_fx * ta_0_xy_1[j];

                ta_0_xxxz_0[j] = pb_x[j] * ta_0_xxz_0[j] - pc_x[j] * ta_0_xxz_1[j] + fl1_fx * ta_0_xz_0[j] - fl1_fx * ta_0_xz_1[j];

                ta_0_xxyy_0[j] = pb_x[j] * ta_0_xyy_0[j] - pc_x[j] * ta_0_xyy_1[j] + 0.5 * fl1_fx * ta_0_yy_0[j] - 0.5 * fl1_fx * ta_0_yy_1[j];

                ta_0_xxyz_0[j] = pb_x[j] * ta_0_xyz_0[j] - pc_x[j] * ta_0_xyz_1[j] + 0.5 * fl1_fx * ta_0_yz_0[j] - 0.5 * fl1_fx * ta_0_yz_1[j];

                ta_0_xxzz_0[j] = pb_x[j] * ta_0_xzz_0[j] - pc_x[j] * ta_0_xzz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j];

                ta_0_xyyy_0[j] = pb_x[j] * ta_0_yyy_0[j] - pc_x[j] * ta_0_yyy_1[j];

                ta_0_xyyz_0[j] = pb_x[j] * ta_0_yyz_0[j] - pc_x[j] * ta_0_yyz_1[j];

                ta_0_xyzz_0[j] = pb_x[j] * ta_0_yzz_0[j] - pc_x[j] * ta_0_yzz_1[j];

                ta_0_xzzz_0[j] = pb_x[j] * ta_0_zzz_0[j] - pc_x[j] * ta_0_zzz_1[j];

                ta_0_yyyy_0[j] = pb_y[j] * ta_0_yyy_0[j] - pc_y[j] * ta_0_yyy_1[j] + 1.5 * fl1_fx * ta_0_yy_0[j] - 1.5 * fl1_fx * ta_0_yy_1[j];

                ta_0_yyyz_0[j] = pb_y[j] * ta_0_yyz_0[j] - pc_y[j] * ta_0_yyz_1[j] + fl1_fx * ta_0_yz_0[j] - fl1_fx * ta_0_yz_1[j];

                ta_0_yyzz_0[j] = pb_y[j] * ta_0_yzz_0[j] - pc_y[j] * ta_0_yzz_1[j] + 0.5 * fl1_fx * ta_0_zz_0[j] - 0.5 * fl1_fx * ta_0_zz_1[j];

                ta_0_yzzz_0[j] = pb_y[j] * ta_0_zzz_0[j] - pc_y[j] * ta_0_zzz_1[j];

                ta_0_zzzz_0[j] = pb_z[j] * ta_0_zzz_0[j] - pc_z[j] * ta_0_zzz_1[j] + 1.5 * fl1_fx * ta_0_zz_0[j] - 1.5 * fl1_fx * ta_0_zz_1[j];
            }

            idx++;
        }
    }
}

void
compNuclearPotentialForGS(CMemBlock2D<double>&       primBuffer,
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

    auto mord = recursionMap.getMaxOrder({"Nuclear Potential"}, {4, -1, -1, -1}, {0, -1, -1, -1}, 1, 1);

    for (int32_t iord = 0; iord <= mord; iord++)
    {
        // set up index of integral

        auto pidx_a_4_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {4, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        // check if integral is needed in recursion expansion

        if (pidx_a_4_0_m0 == -1) continue;

        // set up indexes of auxilary integral

        auto pidx_a_3_0_m0 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord));

        auto pidx_a_3_0_m1 =
            recursionMap.getIndexOfTerm(CRecursionTerm({"Nuclear Potential"}, 0, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, iord + 1));

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

            auto ta_xxxx_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx);

            auto ta_xxxy_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 1);

            auto ta_xxxz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 2);

            auto ta_xxyy_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 3);

            auto ta_xxyz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 4);

            auto ta_xxzz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 5);

            auto ta_xyyy_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 6);

            auto ta_xyyz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 7);

            auto ta_xyzz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 8);

            auto ta_xzzz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 9);

            auto ta_yyyy_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 10);

            auto ta_yyyz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 11);

            auto ta_yyzz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 12);

            auto ta_yzzz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 13);

            auto ta_zzzz_0_0 = primBuffer.data(pidx_a_4_0_m0 + 15 * idx + 14);

            #pragma omp simd aligned(fx, pa_x, pa_y, pa_z, pc_x, pc_y, pc_z, ta_xx_0_0, ta_xx_0_1, ta_xxx_0_0, \
                                         ta_xxx_0_1, ta_xxxx_0_0, ta_xxxy_0_0, ta_xxxz_0_0, ta_xxy_0_0, ta_xxy_0_1, \
                                         ta_xxyy_0_0, ta_xxyz_0_0, ta_xxz_0_0, ta_xxz_0_1, ta_xxzz_0_0, ta_xy_0_0, ta_xy_0_1, \
                                         ta_xyy_0_0, ta_xyy_0_1, ta_xyyy_0_0, ta_xyyz_0_0, ta_xyz_0_0, ta_xyz_0_1, \
                                         ta_xyzz_0_0, ta_xz_0_0, ta_xz_0_1, ta_xzz_0_0, ta_xzz_0_1, ta_xzzz_0_0, ta_yy_0_0, \
                                         ta_yy_0_1, ta_yyy_0_0, ta_yyy_0_1, ta_yyyy_0_0, ta_yyyz_0_0, ta_yyz_0_0, \
                                         ta_yyz_0_1, ta_yyzz_0_0, ta_yz_0_0, ta_yz_0_1, ta_yzz_0_0, ta_yzz_0_1, ta_yzzz_0_0, \
                                         ta_zz_0_0, ta_zz_0_1, ta_zzz_0_0, ta_zzz_0_1, ta_zzzz_0_0: VLX_ALIGN)
            for (int32_t j = 0; j < nprim; j++)
            {
                double fl1_fx = fx[j];

                ta_xxxx_0_0[j] = pa_x[j] * ta_xxx_0_0[j] - pc_x[j] * ta_xxx_0_1[j] + 1.5 * fl1_fx * ta_xx_0_0[j] - 1.5 * fl1_fx * ta_xx_0_1[j];

                ta_xxxy_0_0[j] = pa_x[j] * ta_xxy_0_0[j] - pc_x[j] * ta_xxy_0_1[j] + fl1_fx * ta_xy_0_0[j] - fl1_fx * ta_xy_0_1[j];

                ta_xxxz_0_0[j] = pa_x[j] * ta_xxz_0_0[j] - pc_x[j] * ta_xxz_0_1[j] + fl1_fx * ta_xz_0_0[j] - fl1_fx * ta_xz_0_1[j];

                ta_xxyy_0_0[j] = pa_x[j] * ta_xyy_0_0[j] - pc_x[j] * ta_xyy_0_1[j] + 0.5 * fl1_fx * ta_yy_0_0[j] - 0.5 * fl1_fx * ta_yy_0_1[j];

                ta_xxyz_0_0[j] = pa_x[j] * ta_xyz_0_0[j] - pc_x[j] * ta_xyz_0_1[j] + 0.5 * fl1_fx * ta_yz_0_0[j] - 0.5 * fl1_fx * ta_yz_0_1[j];

                ta_xxzz_0_0[j] = pa_x[j] * ta_xzz_0_0[j] - pc_x[j] * ta_xzz_0_1[j] + 0.5 * fl1_fx * ta_zz_0_0[j] - 0.5 * fl1_fx * ta_zz_0_1[j];

                ta_xyyy_0_0[j] = pa_x[j] * ta_yyy_0_0[j] - pc_x[j] * ta_yyy_0_1[j];

                ta_xyyz_0_0[j] = pa_x[j] * ta_yyz_0_0[j] - pc_x[j] * ta_yyz_0_1[j];

                ta_xyzz_0_0[j] = pa_x[j] * ta_yzz_0_0[j] - pc_x[j] * ta_yzz_0_1[j];

                ta_xzzz_0_0[j] = pa_x[j] * ta_zzz_0_0[j] - pc_x[j] * ta_zzz_0_1[j];

                ta_yyyy_0_0[j] = pa_y[j] * ta_yyy_0_0[j] - pc_y[j] * ta_yyy_0_1[j] + 1.5 * fl1_fx * ta_yy_0_0[j] - 1.5 * fl1_fx * ta_yy_0_1[j];

                ta_yyyz_0_0[j] = pa_y[j] * ta_yyz_0_0[j] - pc_y[j] * ta_yyz_0_1[j] + fl1_fx * ta_yz_0_0[j] - fl1_fx * ta_yz_0_1[j];

                ta_yyzz_0_0[j] = pa_y[j] * ta_yzz_0_0[j] - pc_y[j] * ta_yzz_0_1[j] + 0.5 * fl1_fx * ta_zz_0_0[j] - 0.5 * fl1_fx * ta_zz_0_1[j];

                ta_yzzz_0_0[j] = pa_y[j] * ta_zzz_0_0[j] - pc_y[j] * ta_zzz_0_1[j];

                ta_zzzz_0_0[j] = pa_z[j] * ta_zzz_0_0[j] - pc_z[j] * ta_zzz_0_1[j] + 1.5 * fl1_fx * ta_zz_0_0[j] - 1.5 * fl1_fx * ta_zz_0_1[j];
            }

            idx++;
        }
    }
}

}  // namespace npotrecfunc
