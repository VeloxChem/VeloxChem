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

#include "AngularMomentumRecFuncForSX.hpp"

namespace amomrecfunc {  // amomrecfunc namespace

void
compAngularMomentumForSS(CMemBlock2D<double>&       primBuffer,
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

    // set up pointer to overlap integrals

    auto sidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // set up pointer to kinetic energy integrals

    auto tidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to ditances R(PA)

        auto pax = paDistances.data(3 * idx);

        auto pay = paDistances.data(3 * idx + 1);

        auto paz = paDistances.data(3 * idx + 2);

        // set up pointers to ditances R(PC)

        auto pcx = pcDistances.data(3 * idx);

        auto pcy = pcDistances.data(3 * idx + 1);

        auto pcz = pcDistances.data(3 * idx + 2);

        // set up primitive integrals data

        auto fovl = primBuffer.data(sidx + idx);

        auto fmomx = primBuffer.data(tidx + idx);

        auto fmomy = primBuffer.data(tidx + bdim + idx);

        auto fmomz = primBuffer.data(tidx + 2 * bdim + idx);

        #pragma omp simd aligned(fovl, fmomx, fmomy, fmomz, fga, pax, pay, paz, pcx, pcy, pcz: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fx = 2.0 * fga[j] * fovl[j];

            fmomx[j] = fx * (pcz[j] * pay[j] - pcy[j] * paz[j]);

            fmomy[j] = fx * (pcx[j] * paz[j] - pcz[j] * pax[j]);

            fmomz[j] = fx * (pcy[j] * pax[j] - pcx[j] * pay[j]);
        }

        idx++;
    }
}

void
compAngularMomentumForSP(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& pbDistances,
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

    // set up index of integral

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_0_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + idx);

        auto tly_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + bdim + idx);

        auto tlz_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + 2 * bdim + idx);

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto tdx_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + idx);

        auto tdy_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + bdim + idx);

        auto tdz_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + 2 * bdim + idx);

        // set up pointers to integrals

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tdx_0_0_0, tdy_0_0_0, tdz_0_0_0, tlx_0_0_0, \
                                     tlx_0_x_0, tlx_0_y_0, tlx_0_z_0, tly_0_0_0, tly_0_x_0, tly_0_y_0, tly_0_z_0, \
                                     tlz_0_0_0, tlz_0_x_0, tlz_0_y_0, tlz_0_z_0, tpx_0_0_0, tpy_0_0_0, tpz_0_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tlx_0_x_0[j] = pb_x[j] * tlx_0_0_0[j];

            tly_0_x_0[j] = pb_x[j] * tly_0_0_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j] - fl1_fx * fl1_fga * tdz_0_0_0[j];

            tlz_0_x_0[j] = pb_x[j] * tlz_0_0_0[j] - 0.5 * fl1_fx * tpy_0_0_0[j] + fl1_fx * fl1_fga * tdy_0_0_0[j];

            tlx_0_y_0[j] = pb_y[j] * tlx_0_0_0[j] - 0.5 * fl1_fx * tpz_0_0_0[j] + fl1_fx * fl1_fga * tdz_0_0_0[j];

            tly_0_y_0[j] = pb_y[j] * tly_0_0_0[j];

            tlz_0_y_0[j] = pb_y[j] * tlz_0_0_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j] - fl1_fx * fl1_fga * tdx_0_0_0[j];

            tlx_0_z_0[j] = pb_z[j] * tlx_0_0_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j] - fl1_fx * fl1_fga * tdy_0_0_0[j];

            tly_0_z_0[j] = pb_z[j] * tly_0_0_0[j] - 0.5 * fl1_fx * tpx_0_0_0[j] + fl1_fx * fl1_fga * tdx_0_0_0[j];

            tlz_0_z_0[j] = pb_z[j] * tlz_0_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForPS(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_l_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_1_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + idx);

        auto tly_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + bdim + idx);

        auto tlz_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + 2 * bdim + idx);

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto tdx_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + idx);

        auto tdy_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + bdim + idx);

        auto tdz_0_0_0 = primBuffer.data(pidx_d_0_0_m0 + 2 * bdim + idx);

        // set up pointers to integrals

        auto tlx_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx);

        auto tly_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx);

        auto tlz_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx);

        auto tlx_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 1);

        auto tly_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 2);

        auto tly_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 2);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tdx_0_0_0, tdy_0_0_0, tdz_0_0_0, tlx_0_0_0, \
                                     tlx_x_0_0, tlx_y_0_0, tlx_z_0_0, tly_0_0_0, tly_x_0_0, tly_y_0_0, tly_z_0_0, \
                                     tlz_0_0_0, tlz_x_0_0, tlz_y_0_0, tlz_z_0_0, tpx_0_0_0, tpy_0_0_0, tpz_0_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_x_0_0[j] = pa_x[j] * tlx_0_0_0[j];

            tly_x_0_0[j] = pa_x[j] * tly_0_0_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j] + fl1_fx * fl1_fgb * tdz_0_0_0[j];

            tlz_x_0_0[j] = pa_x[j] * tlz_0_0_0[j] - 0.5 * fl1_fx * tpy_0_0_0[j] - fl1_fx * fl1_fgb * tdy_0_0_0[j];

            tlx_y_0_0[j] = pa_y[j] * tlx_0_0_0[j] - 0.5 * fl1_fx * tpz_0_0_0[j] - fl1_fx * fl1_fgb * tdz_0_0_0[j];

            tly_y_0_0[j] = pa_y[j] * tly_0_0_0[j];

            tlz_y_0_0[j] = pa_y[j] * tlz_0_0_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j] + fl1_fx * fl1_fgb * tdx_0_0_0[j];

            tlx_z_0_0[j] = pa_z[j] * tlx_0_0_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j] + fl1_fx * fl1_fgb * tdy_0_0_0[j];

            tly_z_0_0[j] = pa_z[j] * tly_0_0_0[j] - 0.5 * fl1_fx * tpx_0_0_0[j] - fl1_fx * fl1_fgb * tdx_0_0_0[j];

            tlz_z_0_0[j] = pa_z[j] * tlz_0_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForSD(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& pbDistances,
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

    // set up index of integral

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_0_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tlx_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + idx);

        auto tly_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + bdim + idx);

        auto tlz_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + 2 * bdim + idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tdy_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx);

        auto tdz_0_x_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx);

        auto tdx_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 1);

        auto tdy_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tdz_0_y_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tdx_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * idx + 2);

        auto tdy_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tdz_0_z_0 = primBuffer.data(pidx_d_0_1_m0 + 6 * bdim + 3 * idx + 2);

        // set up pointers to integrals

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tdx_0_y_0, tdx_0_z_0, tdy_0_x_0, tdy_0_y_0, \
                                     tdy_0_z_0, tdz_0_x_0, tdz_0_y_0, tdz_0_z_0, tlx_0_0_0, tlx_0_x_0, tlx_0_xx_0, \
                                     tlx_0_xy_0, tlx_0_xz_0, tlx_0_y_0, tlx_0_yy_0, tlx_0_yz_0, tlx_0_z_0, tlx_0_zz_0, \
                                     tly_0_0_0, tly_0_x_0, tly_0_xx_0, tly_0_xy_0, tly_0_xz_0, tly_0_y_0, tly_0_yy_0, \
                                     tly_0_yz_0, tly_0_z_0, tly_0_zz_0, tlz_0_0_0, tlz_0_x_0, tlz_0_xx_0, tlz_0_xy_0, \
                                     tlz_0_xz_0, tlz_0_y_0, tlz_0_yy_0, tlz_0_yz_0, tlz_0_z_0, tlz_0_zz_0, tpx_0_y_0, \
                                     tpx_0_z_0, tpy_0_x_0, tpy_0_y_0, tpy_0_z_0, tpz_0_x_0, tpz_0_y_0, tpz_0_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tlx_0_xx_0[j] = pb_x[j] * tlx_0_x_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j];

            tly_0_xx_0[j] = pb_x[j] * tly_0_x_0[j] + 0.5 * fl1_fx * tly_0_0_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j] - fl1_fx * fl1_fga * tdz_0_x_0[j];

            tlz_0_xx_0[j] = pb_x[j] * tlz_0_x_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j] - 0.5 * fl1_fx * tpy_0_x_0[j] + fl1_fx * fl1_fga * tdy_0_x_0[j];

            tlx_0_xy_0[j] = pb_x[j] * tlx_0_y_0[j];

            tly_0_xy_0[j] = pb_x[j] * tly_0_y_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j] - fl1_fx * fl1_fga * tdz_0_y_0[j];

            tlz_0_xy_0[j] = pb_x[j] * tlz_0_y_0[j] - 0.5 * fl1_fx * tpy_0_y_0[j] + fl1_fx * fl1_fga * tdy_0_y_0[j];

            tlx_0_xz_0[j] = pb_x[j] * tlx_0_z_0[j];

            tly_0_xz_0[j] = pb_x[j] * tly_0_z_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j] - fl1_fx * fl1_fga * tdz_0_z_0[j];

            tlz_0_xz_0[j] = pb_x[j] * tlz_0_z_0[j] - 0.5 * fl1_fx * tpy_0_z_0[j] + fl1_fx * fl1_fga * tdy_0_z_0[j];

            tlx_0_yy_0[j] = pb_y[j] * tlx_0_y_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j] - 0.5 * fl1_fx * tpz_0_y_0[j] + fl1_fx * fl1_fga * tdz_0_y_0[j];

            tly_0_yy_0[j] = pb_y[j] * tly_0_y_0[j] + 0.5 * fl1_fx * tly_0_0_0[j];

            tlz_0_yy_0[j] = pb_y[j] * tlz_0_y_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j] - fl1_fx * fl1_fga * tdx_0_y_0[j];

            tlx_0_yz_0[j] = pb_y[j] * tlx_0_z_0[j] - 0.5 * fl1_fx * tpz_0_z_0[j] + fl1_fx * fl1_fga * tdz_0_z_0[j];

            tly_0_yz_0[j] = pb_y[j] * tly_0_z_0[j];

            tlz_0_yz_0[j] = pb_y[j] * tlz_0_z_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j] - fl1_fx * fl1_fga * tdx_0_z_0[j];

            tlx_0_zz_0[j] = pb_z[j] * tlx_0_z_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j] - fl1_fx * fl1_fga * tdy_0_z_0[j];

            tly_0_zz_0[j] = pb_z[j] * tly_0_z_0[j] + 0.5 * fl1_fx * tly_0_0_0[j] - 0.5 * fl1_fx * tpx_0_z_0[j] + fl1_fx * fl1_fga * tdx_0_z_0[j];

            tlz_0_zz_0[j] = pb_z[j] * tlz_0_z_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForDS(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_l_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_2_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx);

        auto tly_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx);

        auto tlz_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx);

        auto tlx_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 1);

        auto tly_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 2);

        auto tly_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto tlx_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + idx);

        auto tly_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + bdim + idx);

        auto tlz_0_0_0 = primBuffer.data(pidx_l_0_0_m0 + 2 * bdim + idx);

        auto tpy_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx);

        auto tpz_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx);

        auto tpx_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 1);

        auto tpy_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 2);

        auto tpy_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto tdy_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx);

        auto tdz_x_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx);

        auto tdx_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 1);

        auto tdy_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tdz_y_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tdx_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * idx + 2);

        auto tdy_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tdz_z_0_0 = primBuffer.data(pidx_d_1_0_m0 + 6 * bdim + 3 * idx + 2);

        // set up pointers to integrals

        auto tlx_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx);

        auto tly_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx);

        auto tlz_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx);

        auto tlx_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 1);

        auto tly_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 2);

        auto tly_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 3);

        auto tly_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 4);

        auto tly_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 5);

        auto tly_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 5);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tdx_y_0_0, tdx_z_0_0, tdy_x_0_0, tdy_y_0_0, \
                                     tdy_z_0_0, tdz_x_0_0, tdz_y_0_0, tdz_z_0_0, tlx_0_0_0, tlx_x_0_0, tlx_xx_0_0, \
                                     tlx_xy_0_0, tlx_xz_0_0, tlx_y_0_0, tlx_yy_0_0, tlx_yz_0_0, tlx_z_0_0, tlx_zz_0_0, \
                                     tly_0_0_0, tly_x_0_0, tly_xx_0_0, tly_xy_0_0, tly_xz_0_0, tly_y_0_0, tly_yy_0_0, \
                                     tly_yz_0_0, tly_z_0_0, tly_zz_0_0, tlz_0_0_0, tlz_x_0_0, tlz_xx_0_0, tlz_xy_0_0, \
                                     tlz_xz_0_0, tlz_y_0_0, tlz_yy_0_0, tlz_yz_0_0, tlz_z_0_0, tlz_zz_0_0, tpx_y_0_0, \
                                     tpx_z_0_0, tpy_x_0_0, tpy_y_0_0, tpy_z_0_0, tpz_x_0_0, tpz_y_0_0, tpz_z_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xx_0_0[j] = pa_x[j] * tlx_x_0_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j];

            tly_xx_0_0[j] = pa_x[j] * tly_x_0_0[j] + 0.5 * fl1_fx * tly_0_0_0[j] + 0.5 * fl1_fx * tpz_x_0_0[j] + fl1_fx * fl1_fgb * tdz_x_0_0[j];

            tlz_xx_0_0[j] = pa_x[j] * tlz_x_0_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j] - 0.5 * fl1_fx * tpy_x_0_0[j] - fl1_fx * fl1_fgb * tdy_x_0_0[j];

            tlx_xy_0_0[j] = pa_x[j] * tlx_y_0_0[j];

            tly_xy_0_0[j] = pa_x[j] * tly_y_0_0[j] + 0.5 * fl1_fx * tpz_y_0_0[j] + fl1_fx * fl1_fgb * tdz_y_0_0[j];

            tlz_xy_0_0[j] = pa_x[j] * tlz_y_0_0[j] - 0.5 * fl1_fx * tpy_y_0_0[j] - fl1_fx * fl1_fgb * tdy_y_0_0[j];

            tlx_xz_0_0[j] = pa_x[j] * tlx_z_0_0[j];

            tly_xz_0_0[j] = pa_x[j] * tly_z_0_0[j] + 0.5 * fl1_fx * tpz_z_0_0[j] + fl1_fx * fl1_fgb * tdz_z_0_0[j];

            tlz_xz_0_0[j] = pa_x[j] * tlz_z_0_0[j] - 0.5 * fl1_fx * tpy_z_0_0[j] - fl1_fx * fl1_fgb * tdy_z_0_0[j];

            tlx_yy_0_0[j] = pa_y[j] * tlx_y_0_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j] - 0.5 * fl1_fx * tpz_y_0_0[j] - fl1_fx * fl1_fgb * tdz_y_0_0[j];

            tly_yy_0_0[j] = pa_y[j] * tly_y_0_0[j] + 0.5 * fl1_fx * tly_0_0_0[j];

            tlz_yy_0_0[j] = pa_y[j] * tlz_y_0_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j] + 0.5 * fl1_fx * tpx_y_0_0[j] + fl1_fx * fl1_fgb * tdx_y_0_0[j];

            tlx_yz_0_0[j] = pa_y[j] * tlx_z_0_0[j] - 0.5 * fl1_fx * tpz_z_0_0[j] - fl1_fx * fl1_fgb * tdz_z_0_0[j];

            tly_yz_0_0[j] = pa_y[j] * tly_z_0_0[j];

            tlz_yz_0_0[j] = pa_y[j] * tlz_z_0_0[j] + 0.5 * fl1_fx * tpx_z_0_0[j] + fl1_fx * fl1_fgb * tdx_z_0_0[j];

            tlx_zz_0_0[j] = pa_z[j] * tlx_z_0_0[j] + 0.5 * fl1_fx * tlx_0_0_0[j] + 0.5 * fl1_fx * tpy_z_0_0[j] + fl1_fx * fl1_fgb * tdy_z_0_0[j];

            tly_zz_0_0[j] = pa_z[j] * tly_z_0_0[j] + 0.5 * fl1_fx * tly_0_0_0[j] - 0.5 * fl1_fx * tpx_z_0_0[j] - fl1_fx * fl1_fgb * tdx_z_0_0[j];

            tlz_zz_0_0[j] = pa_z[j] * tlz_z_0_0[j] + 0.5 * fl1_fx * tlz_0_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForSF(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& pbDistances,
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

    // set up index of integral

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_0_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tlx_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx);

        auto tly_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx);

        auto tlz_0_x_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx);

        auto tlx_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 1);

        auto tly_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_0_y_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * idx + 2);

        auto tly_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_0_z_0 = primBuffer.data(pidx_l_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpy_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 3);

        auto tpy_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_0_yy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 4);

        auto tpy_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_0_yz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 5);

        auto tpy_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_0_zz_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tdy_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx);

        auto tdz_0_xx_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx);

        auto tdy_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tdz_0_xy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tdy_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tdz_0_xz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tdx_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 3);

        auto tdy_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tdz_0_yy_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tdx_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 4);

        auto tdy_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tdz_0_yz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tdx_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * idx + 5);

        auto tdy_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tdz_0_zz_0 = primBuffer.data(pidx_d_0_2_m0 + 12 * bdim + 6 * idx + 5);

        // set up pointers to integrals

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tdx_0_yy_0, tdx_0_yz_0, tdx_0_zz_0, tdy_0_xx_0, \
                                     tdy_0_xy_0, tdy_0_xz_0, tdy_0_yy_0, tdy_0_yz_0, tdy_0_zz_0, tdz_0_xx_0, tdz_0_xy_0, \
                                     tdz_0_xz_0, tdz_0_yy_0, tdz_0_yz_0, tdz_0_zz_0, tlx_0_x_0, tlx_0_xx_0, tlx_0_xxx_0, \
                                     tlx_0_xxy_0, tlx_0_xxz_0, tlx_0_xy_0, tlx_0_xyy_0, tlx_0_xyz_0, tlx_0_xz_0, \
                                     tlx_0_xzz_0, tlx_0_y_0, tlx_0_yy_0, tlx_0_yyy_0, tlx_0_yyz_0, tlx_0_yz_0, \
                                     tlx_0_yzz_0, tlx_0_z_0, tlx_0_zz_0, tlx_0_zzz_0, tly_0_x_0, tly_0_xx_0, tly_0_xxx_0, \
                                     tly_0_xxy_0, tly_0_xxz_0, tly_0_xy_0, tly_0_xyy_0, tly_0_xyz_0, tly_0_xz_0, \
                                     tly_0_xzz_0, tly_0_y_0, tly_0_yy_0, tly_0_yyy_0, tly_0_yyz_0, tly_0_yz_0, \
                                     tly_0_yzz_0, tly_0_z_0, tly_0_zz_0, tly_0_zzz_0, tlz_0_x_0, tlz_0_xx_0, tlz_0_xxx_0, \
                                     tlz_0_xxy_0, tlz_0_xxz_0, tlz_0_xy_0, tlz_0_xyy_0, tlz_0_xyz_0, tlz_0_xz_0, \
                                     tlz_0_xzz_0, tlz_0_y_0, tlz_0_yy_0, tlz_0_yyy_0, tlz_0_yyz_0, tlz_0_yz_0, \
                                     tlz_0_yzz_0, tlz_0_z_0, tlz_0_zz_0, tlz_0_zzz_0, tpx_0_yy_0, tpx_0_yz_0, tpx_0_zz_0, \
                                     tpy_0_xx_0, tpy_0_xy_0, tpy_0_xz_0, tpy_0_yy_0, tpy_0_yz_0, tpy_0_zz_0, tpz_0_xx_0, \
                                     tpz_0_xy_0, tpz_0_xz_0, tpz_0_yy_0, tpz_0_yz_0, tpz_0_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tlx_0_xxx_0[j] = pb_x[j] * tlx_0_xx_0[j] + fl1_fx * tlx_0_x_0[j];

            tly_0_xxx_0[j] = pb_x[j] * tly_0_xx_0[j] + fl1_fx * tly_0_x_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j] - fl1_fx * fl1_fga * tdz_0_xx_0[j];

            tlz_0_xxx_0[j] = pb_x[j] * tlz_0_xx_0[j] + fl1_fx * tlz_0_x_0[j] - 0.5 * fl1_fx * tpy_0_xx_0[j] + fl1_fx * fl1_fga * tdy_0_xx_0[j];

            tlx_0_xxy_0[j] = pb_x[j] * tlx_0_xy_0[j] + 0.5 * fl1_fx * tlx_0_y_0[j];

            tly_0_xxy_0[j] = pb_x[j] * tly_0_xy_0[j] + 0.5 * fl1_fx * tly_0_y_0[j] + 0.5 * fl1_fx * tpz_0_xy_0[j] - fl1_fx * fl1_fga * tdz_0_xy_0[j];

            tlz_0_xxy_0[j] = pb_x[j] * tlz_0_xy_0[j] + 0.5 * fl1_fx * tlz_0_y_0[j] - 0.5 * fl1_fx * tpy_0_xy_0[j] + fl1_fx * fl1_fga * tdy_0_xy_0[j];

            tlx_0_xxz_0[j] = pb_x[j] * tlx_0_xz_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j];

            tly_0_xxz_0[j] = pb_x[j] * tly_0_xz_0[j] + 0.5 * fl1_fx * tly_0_z_0[j] + 0.5 * fl1_fx * tpz_0_xz_0[j] - fl1_fx * fl1_fga * tdz_0_xz_0[j];

            tlz_0_xxz_0[j] = pb_x[j] * tlz_0_xz_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] - 0.5 * fl1_fx * tpy_0_xz_0[j] + fl1_fx * fl1_fga * tdy_0_xz_0[j];

            tlx_0_xyy_0[j] = pb_x[j] * tlx_0_yy_0[j];

            tly_0_xyy_0[j] = pb_x[j] * tly_0_yy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j] - fl1_fx * fl1_fga * tdz_0_yy_0[j];

            tlz_0_xyy_0[j] = pb_x[j] * tlz_0_yy_0[j] - 0.5 * fl1_fx * tpy_0_yy_0[j] + fl1_fx * fl1_fga * tdy_0_yy_0[j];

            tlx_0_xyz_0[j] = pb_x[j] * tlx_0_yz_0[j];

            tly_0_xyz_0[j] = pb_x[j] * tly_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j] - fl1_fx * fl1_fga * tdz_0_yz_0[j];

            tlz_0_xyz_0[j] = pb_x[j] * tlz_0_yz_0[j] - 0.5 * fl1_fx * tpy_0_yz_0[j] + fl1_fx * fl1_fga * tdy_0_yz_0[j];

            tlx_0_xzz_0[j] = pb_x[j] * tlx_0_zz_0[j];

            tly_0_xzz_0[j] = pb_x[j] * tly_0_zz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j] - fl1_fx * fl1_fga * tdz_0_zz_0[j];

            tlz_0_xzz_0[j] = pb_x[j] * tlz_0_zz_0[j] - 0.5 * fl1_fx * tpy_0_zz_0[j] + fl1_fx * fl1_fga * tdy_0_zz_0[j];

            tlx_0_yyy_0[j] = pb_y[j] * tlx_0_yy_0[j] + fl1_fx * tlx_0_y_0[j] - 0.5 * fl1_fx * tpz_0_yy_0[j] + fl1_fx * fl1_fga * tdz_0_yy_0[j];

            tly_0_yyy_0[j] = pb_y[j] * tly_0_yy_0[j] + fl1_fx * tly_0_y_0[j];

            tlz_0_yyy_0[j] = pb_y[j] * tlz_0_yy_0[j] + fl1_fx * tlz_0_y_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j] - fl1_fx * fl1_fga * tdx_0_yy_0[j];

            tlx_0_yyz_0[j] = pb_y[j] * tlx_0_yz_0[j] + 0.5 * fl1_fx * tlx_0_z_0[j] - 0.5 * fl1_fx * tpz_0_yz_0[j] + fl1_fx * fl1_fga * tdz_0_yz_0[j];

            tly_0_yyz_0[j] = pb_y[j] * tly_0_yz_0[j] + 0.5 * fl1_fx * tly_0_z_0[j];

            tlz_0_yyz_0[j] = pb_y[j] * tlz_0_yz_0[j] + 0.5 * fl1_fx * tlz_0_z_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] - fl1_fx * fl1_fga * tdx_0_yz_0[j];

            tlx_0_yzz_0[j] = pb_y[j] * tlx_0_zz_0[j] - 0.5 * fl1_fx * tpz_0_zz_0[j] + fl1_fx * fl1_fga * tdz_0_zz_0[j];

            tly_0_yzz_0[j] = pb_y[j] * tly_0_zz_0[j];

            tlz_0_yzz_0[j] = pb_y[j] * tlz_0_zz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j] - fl1_fx * fl1_fga * tdx_0_zz_0[j];

            tlx_0_zzz_0[j] = pb_z[j] * tlx_0_zz_0[j] + fl1_fx * tlx_0_z_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j] - fl1_fx * fl1_fga * tdy_0_zz_0[j];

            tly_0_zzz_0[j] = pb_z[j] * tly_0_zz_0[j] + fl1_fx * tly_0_z_0[j] - 0.5 * fl1_fx * tpx_0_zz_0[j] + fl1_fx * fl1_fga * tdx_0_zz_0[j];

            tlz_0_zzz_0[j] = pb_z[j] * tlz_0_zz_0[j] + fl1_fx * tlz_0_z_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForFS(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_l_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_3_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx);

        auto tly_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx);

        auto tlz_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx);

        auto tlx_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 1);

        auto tly_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 2);

        auto tly_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 3);

        auto tly_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 4);

        auto tly_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 5);

        auto tly_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 5);

        auto tlx_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx);

        auto tly_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx);

        auto tlz_x_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx);

        auto tlx_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 1);

        auto tly_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tlz_y_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tlx_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * idx + 2);

        auto tly_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tlz_z_0_0 = primBuffer.data(pidx_l_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto tpy_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx);

        auto tpz_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx);

        auto tpy_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tpy_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 3);

        auto tpy_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 4);

        auto tpy_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 5);

        auto tpy_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 5);

        auto tdy_xx_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx);

        auto tdz_xx_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx);

        auto tdy_xy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tdz_xy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tdy_xz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tdz_xz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tdx_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 3);

        auto tdy_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tdz_yy_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tdx_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 4);

        auto tdy_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tdz_yz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tdx_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * idx + 5);

        auto tdy_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tdz_zz_0_0 = primBuffer.data(pidx_d_2_0_m0 + 12 * bdim + 6 * idx + 5);

        // set up pointers to integrals

        auto tlx_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx);

        auto tly_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx);

        auto tlz_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx);

        auto tlx_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 1);

        auto tly_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 2);

        auto tly_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 3);

        auto tly_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 4);

        auto tly_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 5);

        auto tly_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 6);

        auto tly_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 7);

        auto tly_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 8);

        auto tly_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 9);

        auto tly_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 9);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tdx_yy_0_0, tdx_yz_0_0, tdx_zz_0_0, tdy_xx_0_0, \
                                     tdy_xy_0_0, tdy_xz_0_0, tdy_yy_0_0, tdy_yz_0_0, tdy_zz_0_0, tdz_xx_0_0, tdz_xy_0_0, \
                                     tdz_xz_0_0, tdz_yy_0_0, tdz_yz_0_0, tdz_zz_0_0, tlx_x_0_0, tlx_xx_0_0, tlx_xxx_0_0, \
                                     tlx_xxy_0_0, tlx_xxz_0_0, tlx_xy_0_0, tlx_xyy_0_0, tlx_xyz_0_0, tlx_xz_0_0, \
                                     tlx_xzz_0_0, tlx_y_0_0, tlx_yy_0_0, tlx_yyy_0_0, tlx_yyz_0_0, tlx_yz_0_0, \
                                     tlx_yzz_0_0, tlx_z_0_0, tlx_zz_0_0, tlx_zzz_0_0, tly_x_0_0, tly_xx_0_0, tly_xxx_0_0, \
                                     tly_xxy_0_0, tly_xxz_0_0, tly_xy_0_0, tly_xyy_0_0, tly_xyz_0_0, tly_xz_0_0, \
                                     tly_xzz_0_0, tly_y_0_0, tly_yy_0_0, tly_yyy_0_0, tly_yyz_0_0, tly_yz_0_0, \
                                     tly_yzz_0_0, tly_z_0_0, tly_zz_0_0, tly_zzz_0_0, tlz_x_0_0, tlz_xx_0_0, tlz_xxx_0_0, \
                                     tlz_xxy_0_0, tlz_xxz_0_0, tlz_xy_0_0, tlz_xyy_0_0, tlz_xyz_0_0, tlz_xz_0_0, \
                                     tlz_xzz_0_0, tlz_y_0_0, tlz_yy_0_0, tlz_yyy_0_0, tlz_yyz_0_0, tlz_yz_0_0, \
                                     tlz_yzz_0_0, tlz_z_0_0, tlz_zz_0_0, tlz_zzz_0_0, tpx_yy_0_0, tpx_yz_0_0, tpx_zz_0_0, \
                                     tpy_xx_0_0, tpy_xy_0_0, tpy_xz_0_0, tpy_yy_0_0, tpy_yz_0_0, tpy_zz_0_0, tpz_xx_0_0, \
                                     tpz_xy_0_0, tpz_xz_0_0, tpz_yy_0_0, tpz_yz_0_0, tpz_zz_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxx_0_0[j] = pa_x[j] * tlx_xx_0_0[j] + fl1_fx * tlx_x_0_0[j];

            tly_xxx_0_0[j] = pa_x[j] * tly_xx_0_0[j] + fl1_fx * tly_x_0_0[j] + 0.5 * fl1_fx * tpz_xx_0_0[j] + fl1_fx * fl1_fgb * tdz_xx_0_0[j];

            tlz_xxx_0_0[j] = pa_x[j] * tlz_xx_0_0[j] + fl1_fx * tlz_x_0_0[j] - 0.5 * fl1_fx * tpy_xx_0_0[j] - fl1_fx * fl1_fgb * tdy_xx_0_0[j];

            tlx_xxy_0_0[j] = pa_x[j] * tlx_xy_0_0[j] + 0.5 * fl1_fx * tlx_y_0_0[j];

            tly_xxy_0_0[j] = pa_x[j] * tly_xy_0_0[j] + 0.5 * fl1_fx * tly_y_0_0[j] + 0.5 * fl1_fx * tpz_xy_0_0[j] + fl1_fx * fl1_fgb * tdz_xy_0_0[j];

            tlz_xxy_0_0[j] = pa_x[j] * tlz_xy_0_0[j] + 0.5 * fl1_fx * tlz_y_0_0[j] - 0.5 * fl1_fx * tpy_xy_0_0[j] - fl1_fx * fl1_fgb * tdy_xy_0_0[j];

            tlx_xxz_0_0[j] = pa_x[j] * tlx_xz_0_0[j] + 0.5 * fl1_fx * tlx_z_0_0[j];

            tly_xxz_0_0[j] = pa_x[j] * tly_xz_0_0[j] + 0.5 * fl1_fx * tly_z_0_0[j] + 0.5 * fl1_fx * tpz_xz_0_0[j] + fl1_fx * fl1_fgb * tdz_xz_0_0[j];

            tlz_xxz_0_0[j] = pa_x[j] * tlz_xz_0_0[j] + 0.5 * fl1_fx * tlz_z_0_0[j] - 0.5 * fl1_fx * tpy_xz_0_0[j] - fl1_fx * fl1_fgb * tdy_xz_0_0[j];

            tlx_xyy_0_0[j] = pa_x[j] * tlx_yy_0_0[j];

            tly_xyy_0_0[j] = pa_x[j] * tly_yy_0_0[j] + 0.5 * fl1_fx * tpz_yy_0_0[j] + fl1_fx * fl1_fgb * tdz_yy_0_0[j];

            tlz_xyy_0_0[j] = pa_x[j] * tlz_yy_0_0[j] - 0.5 * fl1_fx * tpy_yy_0_0[j] - fl1_fx * fl1_fgb * tdy_yy_0_0[j];

            tlx_xyz_0_0[j] = pa_x[j] * tlx_yz_0_0[j];

            tly_xyz_0_0[j] = pa_x[j] * tly_yz_0_0[j] + 0.5 * fl1_fx * tpz_yz_0_0[j] + fl1_fx * fl1_fgb * tdz_yz_0_0[j];

            tlz_xyz_0_0[j] = pa_x[j] * tlz_yz_0_0[j] - 0.5 * fl1_fx * tpy_yz_0_0[j] - fl1_fx * fl1_fgb * tdy_yz_0_0[j];

            tlx_xzz_0_0[j] = pa_x[j] * tlx_zz_0_0[j];

            tly_xzz_0_0[j] = pa_x[j] * tly_zz_0_0[j] + 0.5 * fl1_fx * tpz_zz_0_0[j] + fl1_fx * fl1_fgb * tdz_zz_0_0[j];

            tlz_xzz_0_0[j] = pa_x[j] * tlz_zz_0_0[j] - 0.5 * fl1_fx * tpy_zz_0_0[j] - fl1_fx * fl1_fgb * tdy_zz_0_0[j];

            tlx_yyy_0_0[j] = pa_y[j] * tlx_yy_0_0[j] + fl1_fx * tlx_y_0_0[j] - 0.5 * fl1_fx * tpz_yy_0_0[j] - fl1_fx * fl1_fgb * tdz_yy_0_0[j];

            tly_yyy_0_0[j] = pa_y[j] * tly_yy_0_0[j] + fl1_fx * tly_y_0_0[j];

            tlz_yyy_0_0[j] = pa_y[j] * tlz_yy_0_0[j] + fl1_fx * tlz_y_0_0[j] + 0.5 * fl1_fx * tpx_yy_0_0[j] + fl1_fx * fl1_fgb * tdx_yy_0_0[j];

            tlx_yyz_0_0[j] = pa_y[j] * tlx_yz_0_0[j] + 0.5 * fl1_fx * tlx_z_0_0[j] - 0.5 * fl1_fx * tpz_yz_0_0[j] - fl1_fx * fl1_fgb * tdz_yz_0_0[j];

            tly_yyz_0_0[j] = pa_y[j] * tly_yz_0_0[j] + 0.5 * fl1_fx * tly_z_0_0[j];

            tlz_yyz_0_0[j] = pa_y[j] * tlz_yz_0_0[j] + 0.5 * fl1_fx * tlz_z_0_0[j] + 0.5 * fl1_fx * tpx_yz_0_0[j] + fl1_fx * fl1_fgb * tdx_yz_0_0[j];

            tlx_yzz_0_0[j] = pa_y[j] * tlx_zz_0_0[j] - 0.5 * fl1_fx * tpz_zz_0_0[j] - fl1_fx * fl1_fgb * tdz_zz_0_0[j];

            tly_yzz_0_0[j] = pa_y[j] * tly_zz_0_0[j];

            tlz_yzz_0_0[j] = pa_y[j] * tlz_zz_0_0[j] + 0.5 * fl1_fx * tpx_zz_0_0[j] + fl1_fx * fl1_fgb * tdx_zz_0_0[j];

            tlx_zzz_0_0[j] = pa_z[j] * tlx_zz_0_0[j] + fl1_fx * tlx_z_0_0[j] + 0.5 * fl1_fx * tpy_zz_0_0[j] + fl1_fx * fl1_fgb * tdy_zz_0_0[j];

            tly_zzz_0_0[j] = pa_z[j] * tly_zz_0_0[j] + fl1_fx * tly_z_0_0[j] - 0.5 * fl1_fx * tpx_zz_0_0[j] - fl1_fx * fl1_fgb * tdx_zz_0_0[j];

            tlz_zzz_0_0[j] = pa_z[j] * tlz_zz_0_0[j] + fl1_fx * tlz_z_0_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForSG(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& pbDistances,
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

    // set up index of integral

    auto pidx_l_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_0_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fga = osFactors.data(4 * idx + 2);

        // set up pointers to tensors product of distances R(PB) = P - B

        auto pb_x = pbDistances.data(3 * idx);

        auto pb_y = pbDistances.data(3 * idx + 1);

        auto pb_z = pbDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx);

        auto tly_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx);

        auto tlz_0_xxx_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx);

        auto tlx_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 1);

        auto tly_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_0_xxy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 2);

        auto tly_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_0_xxz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 3);

        auto tly_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_0_xyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 4);

        auto tly_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_0_xyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 5);

        auto tly_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_0_xzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 6);

        auto tly_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_0_yyy_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 7);

        auto tly_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_0_yyz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 8);

        auto tly_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_0_yzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * idx + 9);

        auto tly_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_0_zzz_0 = primBuffer.data(pidx_l_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tlx_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx);

        auto tly_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx);

        auto tlz_0_xx_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx);

        auto tlx_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 1);

        auto tly_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_0_xy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 2);

        auto tly_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_0_xz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 3);

        auto tly_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_0_yy_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 4);

        auto tly_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_0_yz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * idx + 5);

        auto tly_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_0_zz_0 = primBuffer.data(pidx_l_0_2_m0 + 12 * bdim + 6 * idx + 5);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpy_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tpz_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tpx_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 6);

        auto tpy_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tpz_0_yyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tpx_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 7);

        auto tpy_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tpz_0_yyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tpx_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 8);

        auto tpy_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tpz_0_yzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tpx_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 9);

        auto tpy_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tpz_0_zzz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 9);

        auto tdy_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx);

        auto tdz_0_xxx_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx);

        auto tdy_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tdz_0_xxy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tdy_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tdz_0_xxz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tdy_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tdz_0_xyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tdy_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tdz_0_xyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tdy_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 5);

        auto tdz_0_xzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 5);

        auto tdx_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 6);

        auto tdy_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 6);

        auto tdz_0_yyy_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 6);

        auto tdx_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 7);

        auto tdy_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 7);

        auto tdz_0_yyz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 7);

        auto tdx_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 8);

        auto tdy_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 8);

        auto tdz_0_yzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 8);

        auto tdx_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * idx + 9);

        auto tdy_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 10 * bdim + 10 * idx + 9);

        auto tdz_0_zzz_0 = primBuffer.data(pidx_d_0_3_m0 + 20 * bdim + 10 * idx + 9);

        // set up pointers to integrals

        auto tlx_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx);

        auto tly_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx);

        auto tlz_0_xxxx_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx);

        auto tlx_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 1);

        auto tly_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tlz_0_xxxy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tlx_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 2);

        auto tly_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tlz_0_xxxz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tlx_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 3);

        auto tly_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tlz_0_xxyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tlx_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 4);

        auto tly_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tlz_0_xxyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tlx_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 5);

        auto tly_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tlz_0_xxzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tlx_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 6);

        auto tly_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tlz_0_xyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tlx_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 7);

        auto tly_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tlz_0_xyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tlx_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 8);

        auto tly_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tlz_0_xyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tlx_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 9);

        auto tly_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tlz_0_xzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tlx_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 10);

        auto tly_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tlz_0_yyyy_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tlx_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 11);

        auto tly_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tlz_0_yyyz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tlx_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 12);

        auto tly_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tlz_0_yyzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tlx_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 13);

        auto tly_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tlz_0_yzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tlx_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * idx + 14);

        auto tly_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tlz_0_zzzz_0 = primBuffer.data(pidx_l_0_4_m0 + 30 * bdim + 15 * idx + 14);

        #pragma omp simd aligned(fga, fx, pb_x, pb_y, pb_z, tdx_0_yyy_0, tdx_0_yyz_0, tdx_0_yzz_0, \
                                     tdx_0_zzz_0, tdy_0_xxx_0, tdy_0_xxy_0, tdy_0_xxz_0, tdy_0_xyy_0, tdy_0_xyz_0, \
                                     tdy_0_xzz_0, tdy_0_yyy_0, tdy_0_yyz_0, tdy_0_yzz_0, tdy_0_zzz_0, tdz_0_xxx_0, \
                                     tdz_0_xxy_0, tdz_0_xxz_0, tdz_0_xyy_0, tdz_0_xyz_0, tdz_0_xzz_0, tdz_0_yyy_0, \
                                     tdz_0_yyz_0, tdz_0_yzz_0, tdz_0_zzz_0, tlx_0_xx_0, tlx_0_xxx_0, tlx_0_xxxx_0, \
                                     tlx_0_xxxy_0, tlx_0_xxxz_0, tlx_0_xxy_0, tlx_0_xxyy_0, tlx_0_xxyz_0, tlx_0_xxz_0, \
                                     tlx_0_xxzz_0, tlx_0_xy_0, tlx_0_xyy_0, tlx_0_xyyy_0, tlx_0_xyyz_0, tlx_0_xyz_0, \
                                     tlx_0_xyzz_0, tlx_0_xz_0, tlx_0_xzz_0, tlx_0_xzzz_0, tlx_0_yy_0, tlx_0_yyy_0, \
                                     tlx_0_yyyy_0, tlx_0_yyyz_0, tlx_0_yyz_0, tlx_0_yyzz_0, tlx_0_yz_0, tlx_0_yzz_0, \
                                     tlx_0_yzzz_0, tlx_0_zz_0, tlx_0_zzz_0, tlx_0_zzzz_0, tly_0_xx_0, tly_0_xxx_0, \
                                     tly_0_xxxx_0, tly_0_xxxy_0, tly_0_xxxz_0, tly_0_xxy_0, tly_0_xxyy_0, tly_0_xxyz_0, \
                                     tly_0_xxz_0, tly_0_xxzz_0, tly_0_xy_0, tly_0_xyy_0, tly_0_xyyy_0, tly_0_xyyz_0, \
                                     tly_0_xyz_0, tly_0_xyzz_0, tly_0_xz_0, tly_0_xzz_0, tly_0_xzzz_0, tly_0_yy_0, \
                                     tly_0_yyy_0, tly_0_yyyy_0, tly_0_yyyz_0, tly_0_yyz_0, tly_0_yyzz_0, tly_0_yz_0, \
                                     tly_0_yzz_0, tly_0_yzzz_0, tly_0_zz_0, tly_0_zzz_0, tly_0_zzzz_0, tlz_0_xx_0, \
                                     tlz_0_xxx_0, tlz_0_xxxx_0, tlz_0_xxxy_0, tlz_0_xxxz_0, tlz_0_xxy_0, tlz_0_xxyy_0, \
                                     tlz_0_xxyz_0, tlz_0_xxz_0, tlz_0_xxzz_0, tlz_0_xy_0, tlz_0_xyy_0, tlz_0_xyyy_0, \
                                     tlz_0_xyyz_0, tlz_0_xyz_0, tlz_0_xyzz_0, tlz_0_xz_0, tlz_0_xzz_0, tlz_0_xzzz_0, \
                                     tlz_0_yy_0, tlz_0_yyy_0, tlz_0_yyyy_0, tlz_0_yyyz_0, tlz_0_yyz_0, tlz_0_yyzz_0, \
                                     tlz_0_yz_0, tlz_0_yzz_0, tlz_0_yzzz_0, tlz_0_zz_0, tlz_0_zzz_0, tlz_0_zzzz_0, \
                                     tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yzz_0, tpx_0_zzz_0, tpy_0_xxx_0, tpy_0_xxy_0, \
                                     tpy_0_xxz_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xzz_0, tpy_0_yyy_0, tpy_0_yyz_0, \
                                     tpy_0_yzz_0, tpy_0_zzz_0, tpz_0_xxx_0, tpz_0_xxy_0, tpz_0_xxz_0, tpz_0_xyy_0, \
                                     tpz_0_xyz_0, tpz_0_xzz_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yzz_0, tpz_0_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            tlx_0_xxxx_0[j] = pb_x[j] * tlx_0_xxx_0[j] + 1.5 * fl1_fx * tlx_0_xx_0[j];

            tly_0_xxxx_0[j] =
                pb_x[j] * tly_0_xxx_0[j] + 1.5 * fl1_fx * tly_0_xx_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j] - fl1_fx * fl1_fga * tdz_0_xxx_0[j];

            tlz_0_xxxx_0[j] =
                pb_x[j] * tlz_0_xxx_0[j] + 1.5 * fl1_fx * tlz_0_xx_0[j] - 0.5 * fl1_fx * tpy_0_xxx_0[j] + fl1_fx * fl1_fga * tdy_0_xxx_0[j];

            tlx_0_xxxy_0[j] = pb_x[j] * tlx_0_xxy_0[j] + fl1_fx * tlx_0_xy_0[j];

            tly_0_xxxy_0[j] = pb_x[j] * tly_0_xxy_0[j] + fl1_fx * tly_0_xy_0[j] + 0.5 * fl1_fx * tpz_0_xxy_0[j] - fl1_fx * fl1_fga * tdz_0_xxy_0[j];

            tlz_0_xxxy_0[j] = pb_x[j] * tlz_0_xxy_0[j] + fl1_fx * tlz_0_xy_0[j] - 0.5 * fl1_fx * tpy_0_xxy_0[j] + fl1_fx * fl1_fga * tdy_0_xxy_0[j];

            tlx_0_xxxz_0[j] = pb_x[j] * tlx_0_xxz_0[j] + fl1_fx * tlx_0_xz_0[j];

            tly_0_xxxz_0[j] = pb_x[j] * tly_0_xxz_0[j] + fl1_fx * tly_0_xz_0[j] + 0.5 * fl1_fx * tpz_0_xxz_0[j] - fl1_fx * fl1_fga * tdz_0_xxz_0[j];

            tlz_0_xxxz_0[j] = pb_x[j] * tlz_0_xxz_0[j] + fl1_fx * tlz_0_xz_0[j] - 0.5 * fl1_fx * tpy_0_xxz_0[j] + fl1_fx * fl1_fga * tdy_0_xxz_0[j];

            tlx_0_xxyy_0[j] = pb_x[j] * tlx_0_xyy_0[j] + 0.5 * fl1_fx * tlx_0_yy_0[j];

            tly_0_xxyy_0[j] =
                pb_x[j] * tly_0_xyy_0[j] + 0.5 * fl1_fx * tly_0_yy_0[j] + 0.5 * fl1_fx * tpz_0_xyy_0[j] - fl1_fx * fl1_fga * tdz_0_xyy_0[j];

            tlz_0_xxyy_0[j] =
                pb_x[j] * tlz_0_xyy_0[j] + 0.5 * fl1_fx * tlz_0_yy_0[j] - 0.5 * fl1_fx * tpy_0_xyy_0[j] + fl1_fx * fl1_fga * tdy_0_xyy_0[j];

            tlx_0_xxyz_0[j] = pb_x[j] * tlx_0_xyz_0[j] + 0.5 * fl1_fx * tlx_0_yz_0[j];

            tly_0_xxyz_0[j] =
                pb_x[j] * tly_0_xyz_0[j] + 0.5 * fl1_fx * tly_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_xyz_0[j] - fl1_fx * fl1_fga * tdz_0_xyz_0[j];

            tlz_0_xxyz_0[j] =
                pb_x[j] * tlz_0_xyz_0[j] + 0.5 * fl1_fx * tlz_0_yz_0[j] - 0.5 * fl1_fx * tpy_0_xyz_0[j] + fl1_fx * fl1_fga * tdy_0_xyz_0[j];

            tlx_0_xxzz_0[j] = pb_x[j] * tlx_0_xzz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j];

            tly_0_xxzz_0[j] =
                pb_x[j] * tly_0_xzz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j] + 0.5 * fl1_fx * tpz_0_xzz_0[j] - fl1_fx * fl1_fga * tdz_0_xzz_0[j];

            tlz_0_xxzz_0[j] =
                pb_x[j] * tlz_0_xzz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] - 0.5 * fl1_fx * tpy_0_xzz_0[j] + fl1_fx * fl1_fga * tdy_0_xzz_0[j];

            tlx_0_xyyy_0[j] = pb_x[j] * tlx_0_yyy_0[j];

            tly_0_xyyy_0[j] = pb_x[j] * tly_0_yyy_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j] - fl1_fx * fl1_fga * tdz_0_yyy_0[j];

            tlz_0_xyyy_0[j] = pb_x[j] * tlz_0_yyy_0[j] - 0.5 * fl1_fx * tpy_0_yyy_0[j] + fl1_fx * fl1_fga * tdy_0_yyy_0[j];

            tlx_0_xyyz_0[j] = pb_x[j] * tlx_0_yyz_0[j];

            tly_0_xyyz_0[j] = pb_x[j] * tly_0_yyz_0[j] + 0.5 * fl1_fx * tpz_0_yyz_0[j] - fl1_fx * fl1_fga * tdz_0_yyz_0[j];

            tlz_0_xyyz_0[j] = pb_x[j] * tlz_0_yyz_0[j] - 0.5 * fl1_fx * tpy_0_yyz_0[j] + fl1_fx * fl1_fga * tdy_0_yyz_0[j];

            tlx_0_xyzz_0[j] = pb_x[j] * tlx_0_yzz_0[j];

            tly_0_xyzz_0[j] = pb_x[j] * tly_0_yzz_0[j] + 0.5 * fl1_fx * tpz_0_yzz_0[j] - fl1_fx * fl1_fga * tdz_0_yzz_0[j];

            tlz_0_xyzz_0[j] = pb_x[j] * tlz_0_yzz_0[j] - 0.5 * fl1_fx * tpy_0_yzz_0[j] + fl1_fx * fl1_fga * tdy_0_yzz_0[j];

            tlx_0_xzzz_0[j] = pb_x[j] * tlx_0_zzz_0[j];

            tly_0_xzzz_0[j] = pb_x[j] * tly_0_zzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j] - fl1_fx * fl1_fga * tdz_0_zzz_0[j];

            tlz_0_xzzz_0[j] = pb_x[j] * tlz_0_zzz_0[j] - 0.5 * fl1_fx * tpy_0_zzz_0[j] + fl1_fx * fl1_fga * tdy_0_zzz_0[j];

            tlx_0_yyyy_0[j] =
                pb_y[j] * tlx_0_yyy_0[j] + 1.5 * fl1_fx * tlx_0_yy_0[j] - 0.5 * fl1_fx * tpz_0_yyy_0[j] + fl1_fx * fl1_fga * tdz_0_yyy_0[j];

            tly_0_yyyy_0[j] = pb_y[j] * tly_0_yyy_0[j] + 1.5 * fl1_fx * tly_0_yy_0[j];

            tlz_0_yyyy_0[j] =
                pb_y[j] * tlz_0_yyy_0[j] + 1.5 * fl1_fx * tlz_0_yy_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j] - fl1_fx * fl1_fga * tdx_0_yyy_0[j];

            tlx_0_yyyz_0[j] = pb_y[j] * tlx_0_yyz_0[j] + fl1_fx * tlx_0_yz_0[j] - 0.5 * fl1_fx * tpz_0_yyz_0[j] + fl1_fx * fl1_fga * tdz_0_yyz_0[j];

            tly_0_yyyz_0[j] = pb_y[j] * tly_0_yyz_0[j] + fl1_fx * tly_0_yz_0[j];

            tlz_0_yyyz_0[j] = pb_y[j] * tlz_0_yyz_0[j] + fl1_fx * tlz_0_yz_0[j] + 0.5 * fl1_fx * tpx_0_yyz_0[j] - fl1_fx * fl1_fga * tdx_0_yyz_0[j];

            tlx_0_yyzz_0[j] =
                pb_y[j] * tlx_0_yzz_0[j] + 0.5 * fl1_fx * tlx_0_zz_0[j] - 0.5 * fl1_fx * tpz_0_yzz_0[j] + fl1_fx * fl1_fga * tdz_0_yzz_0[j];

            tly_0_yyzz_0[j] = pb_y[j] * tly_0_yzz_0[j] + 0.5 * fl1_fx * tly_0_zz_0[j];

            tlz_0_yyzz_0[j] =
                pb_y[j] * tlz_0_yzz_0[j] + 0.5 * fl1_fx * tlz_0_zz_0[j] + 0.5 * fl1_fx * tpx_0_yzz_0[j] - fl1_fx * fl1_fga * tdx_0_yzz_0[j];

            tlx_0_yzzz_0[j] = pb_y[j] * tlx_0_zzz_0[j] - 0.5 * fl1_fx * tpz_0_zzz_0[j] + fl1_fx * fl1_fga * tdz_0_zzz_0[j];

            tly_0_yzzz_0[j] = pb_y[j] * tly_0_zzz_0[j];

            tlz_0_yzzz_0[j] = pb_y[j] * tlz_0_zzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j] - fl1_fx * fl1_fga * tdx_0_zzz_0[j];

            tlx_0_zzzz_0[j] =
                pb_z[j] * tlx_0_zzz_0[j] + 1.5 * fl1_fx * tlx_0_zz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j] - fl1_fx * fl1_fga * tdy_0_zzz_0[j];

            tly_0_zzzz_0[j] =
                pb_z[j] * tly_0_zzz_0[j] + 1.5 * fl1_fx * tly_0_zz_0[j] - 0.5 * fl1_fx * tpx_0_zzz_0[j] + fl1_fx * fl1_fga * tdx_0_zzz_0[j];

            tlz_0_zzzz_0[j] = pb_z[j] * tlz_0_zzz_0[j] + 1.5 * fl1_fx * tlz_0_zz_0[j];
        }

        idx++;
    }
}

void
compAngularMomentumForGS(CMemBlock2D<double>&       primBuffer,
                         const CRecursionMap&       recursionMap,
                         const CMemBlock2D<double>& osFactors,
                         const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_l_4_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {4, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_l_4_0_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_l_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_l_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_d_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tlx_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx);

        auto tly_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx);

        auto tlz_xxx_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx);

        auto tlx_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 1);

        auto tly_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tlz_xxy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tlx_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 2);

        auto tly_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tlz_xxz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tlx_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 3);

        auto tly_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tlz_xyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tlx_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 4);

        auto tly_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tlz_xyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tlx_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 5);

        auto tly_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tlz_xzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tlx_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 6);

        auto tly_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tlz_yyy_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tlx_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 7);

        auto tly_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tlz_yyz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tlx_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 8);

        auto tly_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tlz_yzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tlx_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * idx + 9);

        auto tly_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tlz_zzz_0_0 = primBuffer.data(pidx_l_3_0_m0 + 20 * bdim + 10 * idx + 9);

        auto tlx_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx);

        auto tly_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx);

        auto tlz_xx_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx);

        auto tlx_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 1);

        auto tly_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tlz_xy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tlx_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 2);

        auto tly_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tlz_xz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tlx_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 3);

        auto tly_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tlz_yy_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tlx_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 4);

        auto tly_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tlz_yz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tlx_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * idx + 5);

        auto tly_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tlz_zz_0_0 = primBuffer.data(pidx_l_2_0_m0 + 12 * bdim + 6 * idx + 5);

        auto tpy_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx);

        auto tpz_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx);

        auto tpy_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tpy_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tpy_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tpy_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tpy_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tpz_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tpx_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 6);

        auto tpy_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tpz_yyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tpx_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 7);

        auto tpy_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tpz_yyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tpx_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 8);

        auto tpy_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tpz_yzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tpx_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 9);

        auto tpy_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tpz_zzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 9);

        auto tdy_xxx_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx);

        auto tdz_xxx_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx);

        auto tdy_xxy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tdz_xxy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tdy_xxz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tdz_xxz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tdy_xyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tdz_xyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tdy_xyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tdz_xyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto tdy_xzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 5);

        auto tdz_xzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 5);

        auto tdx_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 6);

        auto tdy_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 6);

        auto tdz_yyy_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 6);

        auto tdx_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 7);

        auto tdy_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 7);

        auto tdz_yyz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 7);

        auto tdx_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 8);

        auto tdy_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 8);

        auto tdz_yzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 8);

        auto tdx_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * idx + 9);

        auto tdy_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 10 * bdim + 10 * idx + 9);

        auto tdz_zzz_0_0 = primBuffer.data(pidx_d_3_0_m0 + 20 * bdim + 10 * idx + 9);

        // set up pointers to integrals

        auto tlx_xxxx_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx);

        auto tly_xxxx_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx);

        auto tlz_xxxx_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx);

        auto tlx_xxxy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 1);

        auto tly_xxxy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 1);

        auto tlz_xxxy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 1);

        auto tlx_xxxz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 2);

        auto tly_xxxz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 2);

        auto tlz_xxxz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 2);

        auto tlx_xxyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 3);

        auto tly_xxyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 3);

        auto tlz_xxyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 3);

        auto tlx_xxyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 4);

        auto tly_xxyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 4);

        auto tlz_xxyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 4);

        auto tlx_xxzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 5);

        auto tly_xxzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 5);

        auto tlz_xxzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 5);

        auto tlx_xyyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 6);

        auto tly_xyyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 6);

        auto tlz_xyyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 6);

        auto tlx_xyyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 7);

        auto tly_xyyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 7);

        auto tlz_xyyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 7);

        auto tlx_xyzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 8);

        auto tly_xyzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 8);

        auto tlz_xyzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 8);

        auto tlx_xzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 9);

        auto tly_xzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 9);

        auto tlz_xzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 9);

        auto tlx_yyyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 10);

        auto tly_yyyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 10);

        auto tlz_yyyy_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 10);

        auto tlx_yyyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 11);

        auto tly_yyyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 11);

        auto tlz_yyyz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 11);

        auto tlx_yyzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 12);

        auto tly_yyzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 12);

        auto tlz_yyzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 12);

        auto tlx_yzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 13);

        auto tly_yzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 13);

        auto tlz_yzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 13);

        auto tlx_zzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * idx + 14);

        auto tly_zzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 15 * bdim + 15 * idx + 14);

        auto tlz_zzzz_0_0 = primBuffer.data(pidx_l_4_0_m0 + 30 * bdim + 15 * idx + 14);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tdx_yyy_0_0, tdx_yyz_0_0, tdx_yzz_0_0, \
                                     tdx_zzz_0_0, tdy_xxx_0_0, tdy_xxy_0_0, tdy_xxz_0_0, tdy_xyy_0_0, tdy_xyz_0_0, \
                                     tdy_xzz_0_0, tdy_yyy_0_0, tdy_yyz_0_0, tdy_yzz_0_0, tdy_zzz_0_0, tdz_xxx_0_0, \
                                     tdz_xxy_0_0, tdz_xxz_0_0, tdz_xyy_0_0, tdz_xyz_0_0, tdz_xzz_0_0, tdz_yyy_0_0, \
                                     tdz_yyz_0_0, tdz_yzz_0_0, tdz_zzz_0_0, tlx_xx_0_0, tlx_xxx_0_0, tlx_xxxx_0_0, \
                                     tlx_xxxy_0_0, tlx_xxxz_0_0, tlx_xxy_0_0, tlx_xxyy_0_0, tlx_xxyz_0_0, tlx_xxz_0_0, \
                                     tlx_xxzz_0_0, tlx_xy_0_0, tlx_xyy_0_0, tlx_xyyy_0_0, tlx_xyyz_0_0, tlx_xyz_0_0, \
                                     tlx_xyzz_0_0, tlx_xz_0_0, tlx_xzz_0_0, tlx_xzzz_0_0, tlx_yy_0_0, tlx_yyy_0_0, \
                                     tlx_yyyy_0_0, tlx_yyyz_0_0, tlx_yyz_0_0, tlx_yyzz_0_0, tlx_yz_0_0, tlx_yzz_0_0, \
                                     tlx_yzzz_0_0, tlx_zz_0_0, tlx_zzz_0_0, tlx_zzzz_0_0, tly_xx_0_0, tly_xxx_0_0, \
                                     tly_xxxx_0_0, tly_xxxy_0_0, tly_xxxz_0_0, tly_xxy_0_0, tly_xxyy_0_0, tly_xxyz_0_0, \
                                     tly_xxz_0_0, tly_xxzz_0_0, tly_xy_0_0, tly_xyy_0_0, tly_xyyy_0_0, tly_xyyz_0_0, \
                                     tly_xyz_0_0, tly_xyzz_0_0, tly_xz_0_0, tly_xzz_0_0, tly_xzzz_0_0, tly_yy_0_0, \
                                     tly_yyy_0_0, tly_yyyy_0_0, tly_yyyz_0_0, tly_yyz_0_0, tly_yyzz_0_0, tly_yz_0_0, \
                                     tly_yzz_0_0, tly_yzzz_0_0, tly_zz_0_0, tly_zzz_0_0, tly_zzzz_0_0, tlz_xx_0_0, \
                                     tlz_xxx_0_0, tlz_xxxx_0_0, tlz_xxxy_0_0, tlz_xxxz_0_0, tlz_xxy_0_0, tlz_xxyy_0_0, \
                                     tlz_xxyz_0_0, tlz_xxz_0_0, tlz_xxzz_0_0, tlz_xy_0_0, tlz_xyy_0_0, tlz_xyyy_0_0, \
                                     tlz_xyyz_0_0, tlz_xyz_0_0, tlz_xyzz_0_0, tlz_xz_0_0, tlz_xzz_0_0, tlz_xzzz_0_0, \
                                     tlz_yy_0_0, tlz_yyy_0_0, tlz_yyyy_0_0, tlz_yyyz_0_0, tlz_yyz_0_0, tlz_yyzz_0_0, \
                                     tlz_yz_0_0, tlz_yzz_0_0, tlz_yzzz_0_0, tlz_zz_0_0, tlz_zzz_0_0, tlz_zzzz_0_0, \
                                     tpx_yyy_0_0, tpx_yyz_0_0, tpx_yzz_0_0, tpx_zzz_0_0, tpy_xxx_0_0, tpy_xxy_0_0, \
                                     tpy_xxz_0_0, tpy_xyy_0_0, tpy_xyz_0_0, tpy_xzz_0_0, tpy_yyy_0_0, tpy_yyz_0_0, \
                                     tpy_yzz_0_0, tpy_zzz_0_0, tpz_xxx_0_0, tpz_xxy_0_0, tpz_xxz_0_0, tpz_xyy_0_0, \
                                     tpz_xyz_0_0, tpz_xzz_0_0, tpz_yyy_0_0, tpz_yyz_0_0, tpz_yzz_0_0, tpz_zzz_0_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tlx_xxxx_0_0[j] = pa_x[j] * tlx_xxx_0_0[j] + 1.5 * fl1_fx * tlx_xx_0_0[j];

            tly_xxxx_0_0[j] =
                pa_x[j] * tly_xxx_0_0[j] + 1.5 * fl1_fx * tly_xx_0_0[j] + 0.5 * fl1_fx * tpz_xxx_0_0[j] + fl1_fx * fl1_fgb * tdz_xxx_0_0[j];

            tlz_xxxx_0_0[j] =
                pa_x[j] * tlz_xxx_0_0[j] + 1.5 * fl1_fx * tlz_xx_0_0[j] - 0.5 * fl1_fx * tpy_xxx_0_0[j] - fl1_fx * fl1_fgb * tdy_xxx_0_0[j];

            tlx_xxxy_0_0[j] = pa_x[j] * tlx_xxy_0_0[j] + fl1_fx * tlx_xy_0_0[j];

            tly_xxxy_0_0[j] = pa_x[j] * tly_xxy_0_0[j] + fl1_fx * tly_xy_0_0[j] + 0.5 * fl1_fx * tpz_xxy_0_0[j] + fl1_fx * fl1_fgb * tdz_xxy_0_0[j];

            tlz_xxxy_0_0[j] = pa_x[j] * tlz_xxy_0_0[j] + fl1_fx * tlz_xy_0_0[j] - 0.5 * fl1_fx * tpy_xxy_0_0[j] - fl1_fx * fl1_fgb * tdy_xxy_0_0[j];

            tlx_xxxz_0_0[j] = pa_x[j] * tlx_xxz_0_0[j] + fl1_fx * tlx_xz_0_0[j];

            tly_xxxz_0_0[j] = pa_x[j] * tly_xxz_0_0[j] + fl1_fx * tly_xz_0_0[j] + 0.5 * fl1_fx * tpz_xxz_0_0[j] + fl1_fx * fl1_fgb * tdz_xxz_0_0[j];

            tlz_xxxz_0_0[j] = pa_x[j] * tlz_xxz_0_0[j] + fl1_fx * tlz_xz_0_0[j] - 0.5 * fl1_fx * tpy_xxz_0_0[j] - fl1_fx * fl1_fgb * tdy_xxz_0_0[j];

            tlx_xxyy_0_0[j] = pa_x[j] * tlx_xyy_0_0[j] + 0.5 * fl1_fx * tlx_yy_0_0[j];

            tly_xxyy_0_0[j] =
                pa_x[j] * tly_xyy_0_0[j] + 0.5 * fl1_fx * tly_yy_0_0[j] + 0.5 * fl1_fx * tpz_xyy_0_0[j] + fl1_fx * fl1_fgb * tdz_xyy_0_0[j];

            tlz_xxyy_0_0[j] =
                pa_x[j] * tlz_xyy_0_0[j] + 0.5 * fl1_fx * tlz_yy_0_0[j] - 0.5 * fl1_fx * tpy_xyy_0_0[j] - fl1_fx * fl1_fgb * tdy_xyy_0_0[j];

            tlx_xxyz_0_0[j] = pa_x[j] * tlx_xyz_0_0[j] + 0.5 * fl1_fx * tlx_yz_0_0[j];

            tly_xxyz_0_0[j] =
                pa_x[j] * tly_xyz_0_0[j] + 0.5 * fl1_fx * tly_yz_0_0[j] + 0.5 * fl1_fx * tpz_xyz_0_0[j] + fl1_fx * fl1_fgb * tdz_xyz_0_0[j];

            tlz_xxyz_0_0[j] =
                pa_x[j] * tlz_xyz_0_0[j] + 0.5 * fl1_fx * tlz_yz_0_0[j] - 0.5 * fl1_fx * tpy_xyz_0_0[j] - fl1_fx * fl1_fgb * tdy_xyz_0_0[j];

            tlx_xxzz_0_0[j] = pa_x[j] * tlx_xzz_0_0[j] + 0.5 * fl1_fx * tlx_zz_0_0[j];

            tly_xxzz_0_0[j] =
                pa_x[j] * tly_xzz_0_0[j] + 0.5 * fl1_fx * tly_zz_0_0[j] + 0.5 * fl1_fx * tpz_xzz_0_0[j] + fl1_fx * fl1_fgb * tdz_xzz_0_0[j];

            tlz_xxzz_0_0[j] =
                pa_x[j] * tlz_xzz_0_0[j] + 0.5 * fl1_fx * tlz_zz_0_0[j] - 0.5 * fl1_fx * tpy_xzz_0_0[j] - fl1_fx * fl1_fgb * tdy_xzz_0_0[j];

            tlx_xyyy_0_0[j] = pa_x[j] * tlx_yyy_0_0[j];

            tly_xyyy_0_0[j] = pa_x[j] * tly_yyy_0_0[j] + 0.5 * fl1_fx * tpz_yyy_0_0[j] + fl1_fx * fl1_fgb * tdz_yyy_0_0[j];

            tlz_xyyy_0_0[j] = pa_x[j] * tlz_yyy_0_0[j] - 0.5 * fl1_fx * tpy_yyy_0_0[j] - fl1_fx * fl1_fgb * tdy_yyy_0_0[j];

            tlx_xyyz_0_0[j] = pa_x[j] * tlx_yyz_0_0[j];

            tly_xyyz_0_0[j] = pa_x[j] * tly_yyz_0_0[j] + 0.5 * fl1_fx * tpz_yyz_0_0[j] + fl1_fx * fl1_fgb * tdz_yyz_0_0[j];

            tlz_xyyz_0_0[j] = pa_x[j] * tlz_yyz_0_0[j] - 0.5 * fl1_fx * tpy_yyz_0_0[j] - fl1_fx * fl1_fgb * tdy_yyz_0_0[j];

            tlx_xyzz_0_0[j] = pa_x[j] * tlx_yzz_0_0[j];

            tly_xyzz_0_0[j] = pa_x[j] * tly_yzz_0_0[j] + 0.5 * fl1_fx * tpz_yzz_0_0[j] + fl1_fx * fl1_fgb * tdz_yzz_0_0[j];

            tlz_xyzz_0_0[j] = pa_x[j] * tlz_yzz_0_0[j] - 0.5 * fl1_fx * tpy_yzz_0_0[j] - fl1_fx * fl1_fgb * tdy_yzz_0_0[j];

            tlx_xzzz_0_0[j] = pa_x[j] * tlx_zzz_0_0[j];

            tly_xzzz_0_0[j] = pa_x[j] * tly_zzz_0_0[j] + 0.5 * fl1_fx * tpz_zzz_0_0[j] + fl1_fx * fl1_fgb * tdz_zzz_0_0[j];

            tlz_xzzz_0_0[j] = pa_x[j] * tlz_zzz_0_0[j] - 0.5 * fl1_fx * tpy_zzz_0_0[j] - fl1_fx * fl1_fgb * tdy_zzz_0_0[j];

            tlx_yyyy_0_0[j] =
                pa_y[j] * tlx_yyy_0_0[j] + 1.5 * fl1_fx * tlx_yy_0_0[j] - 0.5 * fl1_fx * tpz_yyy_0_0[j] - fl1_fx * fl1_fgb * tdz_yyy_0_0[j];

            tly_yyyy_0_0[j] = pa_y[j] * tly_yyy_0_0[j] + 1.5 * fl1_fx * tly_yy_0_0[j];

            tlz_yyyy_0_0[j] =
                pa_y[j] * tlz_yyy_0_0[j] + 1.5 * fl1_fx * tlz_yy_0_0[j] + 0.5 * fl1_fx * tpx_yyy_0_0[j] + fl1_fx * fl1_fgb * tdx_yyy_0_0[j];

            tlx_yyyz_0_0[j] = pa_y[j] * tlx_yyz_0_0[j] + fl1_fx * tlx_yz_0_0[j] - 0.5 * fl1_fx * tpz_yyz_0_0[j] - fl1_fx * fl1_fgb * tdz_yyz_0_0[j];

            tly_yyyz_0_0[j] = pa_y[j] * tly_yyz_0_0[j] + fl1_fx * tly_yz_0_0[j];

            tlz_yyyz_0_0[j] = pa_y[j] * tlz_yyz_0_0[j] + fl1_fx * tlz_yz_0_0[j] + 0.5 * fl1_fx * tpx_yyz_0_0[j] + fl1_fx * fl1_fgb * tdx_yyz_0_0[j];

            tlx_yyzz_0_0[j] =
                pa_y[j] * tlx_yzz_0_0[j] + 0.5 * fl1_fx * tlx_zz_0_0[j] - 0.5 * fl1_fx * tpz_yzz_0_0[j] - fl1_fx * fl1_fgb * tdz_yzz_0_0[j];

            tly_yyzz_0_0[j] = pa_y[j] * tly_yzz_0_0[j] + 0.5 * fl1_fx * tly_zz_0_0[j];

            tlz_yyzz_0_0[j] =
                pa_y[j] * tlz_yzz_0_0[j] + 0.5 * fl1_fx * tlz_zz_0_0[j] + 0.5 * fl1_fx * tpx_yzz_0_0[j] + fl1_fx * fl1_fgb * tdx_yzz_0_0[j];

            tlx_yzzz_0_0[j] = pa_y[j] * tlx_zzz_0_0[j] - 0.5 * fl1_fx * tpz_zzz_0_0[j] - fl1_fx * fl1_fgb * tdz_zzz_0_0[j];

            tly_yzzz_0_0[j] = pa_y[j] * tly_zzz_0_0[j];

            tlz_yzzz_0_0[j] = pa_y[j] * tlz_zzz_0_0[j] + 0.5 * fl1_fx * tpx_zzz_0_0[j] + fl1_fx * fl1_fgb * tdx_zzz_0_0[j];

            tlx_zzzz_0_0[j] =
                pa_z[j] * tlx_zzz_0_0[j] + 1.5 * fl1_fx * tlx_zz_0_0[j] + 0.5 * fl1_fx * tpy_zzz_0_0[j] + fl1_fx * fl1_fgb * tdy_zzz_0_0[j];

            tly_zzzz_0_0[j] =
                pa_z[j] * tly_zzz_0_0[j] + 1.5 * fl1_fx * tly_zz_0_0[j] - 0.5 * fl1_fx * tpx_zzz_0_0[j] - fl1_fx * fl1_fgb * tdx_zzz_0_0[j];

            tlz_zzzz_0_0[j] = pa_z[j] * tlz_zzz_0_0[j] + 1.5 * fl1_fx * tlz_zz_0_0[j];
        }

        idx++;
    }
}

}  // namespace amomrecfunc
