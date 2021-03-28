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

#include "LinearMomentumRecFuncForPX.hpp"

namespace lmomrecfunc {  // lmomrecfunc namespace

void
compLinearMomentumForPP(CMemBlock2D<double>&       primBuffer,
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

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + idx);

        auto tpy_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + bdim + idx);

        auto tpz_0_0_0 = primBuffer.data(pidx_p_0_0_m0 + 2 * bdim + idx);

        auto ts_0_x_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx);

        auto ts_0_y_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 1);

        auto ts_0_z_0 = primBuffer.data(pidx_s_0_1_m0 + 3 * idx + 2);

        // set up pointers to integrals

        auto tpx_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx);

        auto tpy_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx);

        auto tpz_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx);

        auto tpx_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 1);

        auto tpy_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tpz_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tpx_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 2);

        auto tpy_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tpz_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tpx_0_0_0, tpx_0_x_0, tpx_0_y_0, tpx_0_z_0, \
                                     tpx_x_x_0, tpx_x_y_0, tpx_x_z_0, tpx_y_x_0, tpx_y_y_0, tpx_y_z_0, tpx_z_x_0, \
                                     tpx_z_y_0, tpx_z_z_0, tpy_0_0_0, tpy_0_x_0, tpy_0_y_0, tpy_0_z_0, tpy_x_x_0, \
                                     tpy_x_y_0, tpy_x_z_0, tpy_y_x_0, tpy_y_y_0, tpy_y_z_0, tpy_z_x_0, tpy_z_y_0, \
                                     tpy_z_z_0, tpz_0_0_0, tpz_0_x_0, tpz_0_y_0, tpz_0_z_0, tpz_x_x_0, tpz_x_y_0, \
                                     tpz_x_z_0, tpz_y_x_0, tpz_y_y_0, tpz_y_z_0, tpz_z_x_0, tpz_z_y_0, tpz_z_z_0, \
                                     ts_0_x_0, ts_0_y_0, ts_0_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_x_x_0[j] = pa_x[j] * tpx_0_x_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j] - fl1_fgb * fl1_fx * ts_0_x_0[j];

            tpy_x_x_0[j] = pa_x[j] * tpy_0_x_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j];

            tpz_x_x_0[j] = pa_x[j] * tpz_0_x_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j];

            tpx_x_y_0[j] = pa_x[j] * tpx_0_y_0[j] - fl1_fgb * fl1_fx * ts_0_y_0[j];

            tpy_x_y_0[j] = pa_x[j] * tpy_0_y_0[j];

            tpz_x_y_0[j] = pa_x[j] * tpz_0_y_0[j];

            tpx_x_z_0[j] = pa_x[j] * tpx_0_z_0[j] - fl1_fgb * fl1_fx * ts_0_z_0[j];

            tpy_x_z_0[j] = pa_x[j] * tpy_0_z_0[j];

            tpz_x_z_0[j] = pa_x[j] * tpz_0_z_0[j];

            tpx_y_x_0[j] = pa_y[j] * tpx_0_x_0[j];

            tpy_y_x_0[j] = pa_y[j] * tpy_0_x_0[j] - fl1_fgb * fl1_fx * ts_0_x_0[j];

            tpz_y_x_0[j] = pa_y[j] * tpz_0_x_0[j];

            tpx_y_y_0[j] = pa_y[j] * tpx_0_y_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j];

            tpy_y_y_0[j] = pa_y[j] * tpy_0_y_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j] - fl1_fgb * fl1_fx * ts_0_y_0[j];

            tpz_y_y_0[j] = pa_y[j] * tpz_0_y_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j];

            tpx_y_z_0[j] = pa_y[j] * tpx_0_z_0[j];

            tpy_y_z_0[j] = pa_y[j] * tpy_0_z_0[j] - fl1_fgb * fl1_fx * ts_0_z_0[j];

            tpz_y_z_0[j] = pa_y[j] * tpz_0_z_0[j];

            tpx_z_x_0[j] = pa_z[j] * tpx_0_x_0[j];

            tpy_z_x_0[j] = pa_z[j] * tpy_0_x_0[j];

            tpz_z_x_0[j] = pa_z[j] * tpz_0_x_0[j] - fl1_fgb * fl1_fx * ts_0_x_0[j];

            tpx_z_y_0[j] = pa_z[j] * tpx_0_y_0[j];

            tpy_z_y_0[j] = pa_z[j] * tpy_0_y_0[j];

            tpz_z_y_0[j] = pa_z[j] * tpz_0_y_0[j] - fl1_fgb * fl1_fx * ts_0_y_0[j];

            tpx_z_z_0[j] = pa_z[j] * tpx_0_z_0[j] + 0.5 * fl1_fx * tpx_0_0_0[j];

            tpy_z_z_0[j] = pa_z[j] * tpy_0_z_0[j] + 0.5 * fl1_fx * tpy_0_0_0[j];

            tpz_z_z_0[j] = pa_z[j] * tpz_0_z_0[j] + 0.5 * fl1_fx * tpz_0_0_0[j] - fl1_fgb * fl1_fx * ts_0_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPD(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForPD_0_27(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForPD_27_54(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForPD_0_27(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

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

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx);

        auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1);

        auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2);

        auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3);

        auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4);

        auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5);

        // set up pointers to integrals

        auto tpx_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx);

        auto tpy_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx);

        auto tpz_x_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx);

        auto tpx_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 1);

        auto tpy_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_x_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 2);

        auto tpy_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_x_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 3);

        auto tpy_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_x_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 4);

        auto tpy_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_x_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 5);

        auto tpy_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_x_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 6);

        auto tpy_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_y_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 7);

        auto tpy_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_y_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 8);

        auto tpy_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_y_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 8);

        // Batch of Integrals (0,27)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_0_x_0, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, \
                                     tpx_0_y_0, tpx_0_yy_0, tpx_0_yz_0, tpx_0_z_0, tpx_0_zz_0, tpx_x_xx_0, tpx_x_xy_0, \
                                     tpx_x_xz_0, tpx_x_yy_0, tpx_x_yz_0, tpx_x_zz_0, tpx_y_xx_0, tpx_y_xy_0, tpx_y_xz_0, \
                                     tpy_0_x_0, tpy_0_xx_0, tpy_0_xy_0, tpy_0_xz_0, tpy_0_y_0, tpy_0_yy_0, tpy_0_yz_0, \
                                     tpy_0_z_0, tpy_0_zz_0, tpy_x_xx_0, tpy_x_xy_0, tpy_x_xz_0, tpy_x_yy_0, tpy_x_yz_0, \
                                     tpy_x_zz_0, tpy_y_xx_0, tpy_y_xy_0, tpy_y_xz_0, tpz_0_x_0, tpz_0_xx_0, tpz_0_xy_0, \
                                     tpz_0_xz_0, tpz_0_y_0, tpz_0_yy_0, tpz_0_yz_0, tpz_0_z_0, tpz_0_zz_0, tpz_x_xx_0, \
                                     tpz_x_xy_0, tpz_x_xz_0, tpz_x_yy_0, tpz_x_yz_0, tpz_x_zz_0, tpz_y_xx_0, tpz_y_xy_0, \
                                     tpz_y_xz_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, ts_0_yy_0, ts_0_yz_0, ts_0_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_x_xx_0[j] = pa_x[j] * tpx_0_xx_0[j] + fl1_fx * tpx_0_x_0[j] - fl1_fgb * fl1_fx * ts_0_xx_0[j];

            tpy_x_xx_0[j] = pa_x[j] * tpy_0_xx_0[j] + fl1_fx * tpy_0_x_0[j];

            tpz_x_xx_0[j] = pa_x[j] * tpz_0_xx_0[j] + fl1_fx * tpz_0_x_0[j];

            tpx_x_xy_0[j] = pa_x[j] * tpx_0_xy_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j] - fl1_fgb * fl1_fx * ts_0_xy_0[j];

            tpy_x_xy_0[j] = pa_x[j] * tpy_0_xy_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j];

            tpz_x_xy_0[j] = pa_x[j] * tpz_0_xy_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j];

            tpx_x_xz_0[j] = pa_x[j] * tpx_0_xz_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j] - fl1_fgb * fl1_fx * ts_0_xz_0[j];

            tpy_x_xz_0[j] = pa_x[j] * tpy_0_xz_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j];

            tpz_x_xz_0[j] = pa_x[j] * tpz_0_xz_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j];

            tpx_x_yy_0[j] = pa_x[j] * tpx_0_yy_0[j] - fl1_fgb * fl1_fx * ts_0_yy_0[j];

            tpy_x_yy_0[j] = pa_x[j] * tpy_0_yy_0[j];

            tpz_x_yy_0[j] = pa_x[j] * tpz_0_yy_0[j];

            tpx_x_yz_0[j] = pa_x[j] * tpx_0_yz_0[j] - fl1_fgb * fl1_fx * ts_0_yz_0[j];

            tpy_x_yz_0[j] = pa_x[j] * tpy_0_yz_0[j];

            tpz_x_yz_0[j] = pa_x[j] * tpz_0_yz_0[j];

            tpx_x_zz_0[j] = pa_x[j] * tpx_0_zz_0[j] - fl1_fgb * fl1_fx * ts_0_zz_0[j];

            tpy_x_zz_0[j] = pa_x[j] * tpy_0_zz_0[j];

            tpz_x_zz_0[j] = pa_x[j] * tpz_0_zz_0[j];

            tpx_y_xx_0[j] = pa_y[j] * tpx_0_xx_0[j];

            tpy_y_xx_0[j] = pa_y[j] * tpy_0_xx_0[j] - fl1_fgb * fl1_fx * ts_0_xx_0[j];

            tpz_y_xx_0[j] = pa_y[j] * tpz_0_xx_0[j];

            tpx_y_xy_0[j] = pa_y[j] * tpx_0_xy_0[j] + 0.5 * fl1_fx * tpx_0_x_0[j];

            tpy_y_xy_0[j] = pa_y[j] * tpy_0_xy_0[j] + 0.5 * fl1_fx * tpy_0_x_0[j] - fl1_fgb * fl1_fx * ts_0_xy_0[j];

            tpz_y_xy_0[j] = pa_y[j] * tpz_0_xy_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j];

            tpx_y_xz_0[j] = pa_y[j] * tpx_0_xz_0[j];

            tpy_y_xz_0[j] = pa_y[j] * tpy_0_xz_0[j] - fl1_fgb * fl1_fx * ts_0_xz_0[j];

            tpz_y_xz_0[j] = pa_y[j] * tpz_0_xz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPD_27_54(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

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

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx);

        auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1);

        auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2);

        auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3);

        auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4);

        auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5);

        // set up pointers to integrals

        auto tpx_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 9);

        auto tpy_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_y_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 10);

        auto tpy_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_y_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 11);

        auto tpy_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_y_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 12);

        auto tpy_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_z_xx_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 13);

        auto tpy_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_z_xy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 14);

        auto tpy_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_z_xz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 15);

        auto tpy_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_z_yy_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 16);

        auto tpy_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_z_yz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * idx + 17);

        auto tpy_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_z_zz_0 = primBuffer.data(pidx_p_1_2_m0 + 36 * bdim + 18 * idx + 17);

        // Batch of Integrals (27,54)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_0_x_0, tpx_0_xx_0, tpx_0_xy_0, tpx_0_xz_0, \
                                     tpx_0_y_0, tpx_0_yy_0, tpx_0_yz_0, tpx_0_z_0, tpx_0_zz_0, tpx_y_yy_0, tpx_y_yz_0, \
                                     tpx_y_zz_0, tpx_z_xx_0, tpx_z_xy_0, tpx_z_xz_0, tpx_z_yy_0, tpx_z_yz_0, tpx_z_zz_0, \
                                     tpy_0_x_0, tpy_0_xx_0, tpy_0_xy_0, tpy_0_xz_0, tpy_0_y_0, tpy_0_yy_0, tpy_0_yz_0, \
                                     tpy_0_z_0, tpy_0_zz_0, tpy_y_yy_0, tpy_y_yz_0, tpy_y_zz_0, tpy_z_xx_0, tpy_z_xy_0, \
                                     tpy_z_xz_0, tpy_z_yy_0, tpy_z_yz_0, tpy_z_zz_0, tpz_0_x_0, tpz_0_xx_0, tpz_0_xy_0, \
                                     tpz_0_xz_0, tpz_0_y_0, tpz_0_yy_0, tpz_0_yz_0, tpz_0_z_0, tpz_0_zz_0, tpz_y_yy_0, \
                                     tpz_y_yz_0, tpz_y_zz_0, tpz_z_xx_0, tpz_z_xy_0, tpz_z_xz_0, tpz_z_yy_0, tpz_z_yz_0, \
                                     tpz_z_zz_0, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, ts_0_yy_0, ts_0_yz_0, ts_0_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_y_yy_0[j] = pa_y[j] * tpx_0_yy_0[j] + fl1_fx * tpx_0_y_0[j];

            tpy_y_yy_0[j] = pa_y[j] * tpy_0_yy_0[j] + fl1_fx * tpy_0_y_0[j] - fl1_fgb * fl1_fx * ts_0_yy_0[j];

            tpz_y_yy_0[j] = pa_y[j] * tpz_0_yy_0[j] + fl1_fx * tpz_0_y_0[j];

            tpx_y_yz_0[j] = pa_y[j] * tpx_0_yz_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j];

            tpy_y_yz_0[j] = pa_y[j] * tpy_0_yz_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j] - fl1_fgb * fl1_fx * ts_0_yz_0[j];

            tpz_y_yz_0[j] = pa_y[j] * tpz_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j];

            tpx_y_zz_0[j] = pa_y[j] * tpx_0_zz_0[j];

            tpy_y_zz_0[j] = pa_y[j] * tpy_0_zz_0[j] - fl1_fgb * fl1_fx * ts_0_zz_0[j];

            tpz_y_zz_0[j] = pa_y[j] * tpz_0_zz_0[j];

            tpx_z_xx_0[j] = pa_z[j] * tpx_0_xx_0[j];

            tpy_z_xx_0[j] = pa_z[j] * tpy_0_xx_0[j];

            tpz_z_xx_0[j] = pa_z[j] * tpz_0_xx_0[j] - fl1_fgb * fl1_fx * ts_0_xx_0[j];

            tpx_z_xy_0[j] = pa_z[j] * tpx_0_xy_0[j];

            tpy_z_xy_0[j] = pa_z[j] * tpy_0_xy_0[j];

            tpz_z_xy_0[j] = pa_z[j] * tpz_0_xy_0[j] - fl1_fgb * fl1_fx * ts_0_xy_0[j];

            tpx_z_xz_0[j] = pa_z[j] * tpx_0_xz_0[j] + 0.5 * fl1_fx * tpx_0_x_0[j];

            tpy_z_xz_0[j] = pa_z[j] * tpy_0_xz_0[j] + 0.5 * fl1_fx * tpy_0_x_0[j];

            tpz_z_xz_0[j] = pa_z[j] * tpz_0_xz_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j] - fl1_fgb * fl1_fx * ts_0_xz_0[j];

            tpx_z_yy_0[j] = pa_z[j] * tpx_0_yy_0[j];

            tpy_z_yy_0[j] = pa_z[j] * tpy_0_yy_0[j];

            tpz_z_yy_0[j] = pa_z[j] * tpz_0_yy_0[j] - fl1_fgb * fl1_fx * ts_0_yy_0[j];

            tpx_z_yz_0[j] = pa_z[j] * tpx_0_yz_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j];

            tpy_z_yz_0[j] = pa_z[j] * tpy_0_yz_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j];

            tpz_z_yz_0[j] = pa_z[j] * tpz_0_yz_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j] - fl1_fgb * fl1_fx * ts_0_yz_0[j];

            tpx_z_zz_0[j] = pa_z[j] * tpx_0_zz_0[j] + fl1_fx * tpx_0_z_0[j];

            tpy_z_zz_0[j] = pa_z[j] * tpy_0_zz_0[j] + fl1_fx * tpy_0_z_0[j];

            tpz_z_zz_0[j] = pa_z[j] * tpz_0_zz_0[j] + fl1_fx * tpz_0_z_0[j] - fl1_fgb * fl1_fx * ts_0_zz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDP(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForDP_0_27(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForDP_27_54(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForDP_0_27(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx);

        auto tpy_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx);

        auto tpz_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx);

        auto tpx_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 1);

        auto tpy_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tpz_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tpx_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 2);

        auto tpy_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tpz_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx);

        auto tpy_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx);

        auto tpz_x_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx);

        auto tpx_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 1);

        auto tpy_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 2);

        auto tpy_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto ts_x_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx);

        auto ts_x_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 1);

        auto ts_x_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 2);

        auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3);

        auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4);

        auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5);

        auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6);

        auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7);

        auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8);

        // set up pointers to integrals

        auto tpx_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx);

        auto tpy_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx);

        auto tpz_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx);

        auto tpx_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 1);

        auto tpy_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 2);

        auto tpy_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 3);

        auto tpy_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 4);

        auto tpy_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 5);

        auto tpy_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 6);

        auto tpy_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 7);

        auto tpy_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 8);

        auto tpy_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 8);

        // Batch of Integrals (0,27)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_0_x_0, tpx_0_y_0, tpx_0_z_0, tpx_x_0_0, tpx_x_x_0, \
                                     tpx_x_y_0, tpx_x_z_0, tpx_xx_x_0, tpx_xx_y_0, tpx_xx_z_0, tpx_xy_x_0, tpx_xy_y_0, \
                                     tpx_xy_z_0, tpx_xz_x_0, tpx_xz_y_0, tpx_xz_z_0, tpx_y_0_0, tpx_y_x_0, tpx_y_y_0, \
                                     tpx_y_z_0, tpx_z_0_0, tpx_z_x_0, tpx_z_y_0, tpx_z_z_0, tpy_0_x_0, tpy_0_y_0, \
                                     tpy_0_z_0, tpy_x_0_0, tpy_x_x_0, tpy_x_y_0, tpy_x_z_0, tpy_xx_x_0, tpy_xx_y_0, \
                                     tpy_xx_z_0, tpy_xy_x_0, tpy_xy_y_0, tpy_xy_z_0, tpy_xz_x_0, tpy_xz_y_0, tpy_xz_z_0, \
                                     tpy_y_0_0, tpy_y_x_0, tpy_y_y_0, tpy_y_z_0, tpy_z_0_0, tpy_z_x_0, tpy_z_y_0, \
                                     tpy_z_z_0, tpz_0_x_0, tpz_0_y_0, tpz_0_z_0, tpz_x_0_0, tpz_x_x_0, tpz_x_y_0, \
                                     tpz_x_z_0, tpz_xx_x_0, tpz_xx_y_0, tpz_xx_z_0, tpz_xy_x_0, tpz_xy_y_0, tpz_xy_z_0, \
                                     tpz_xz_x_0, tpz_xz_y_0, tpz_xz_z_0, tpz_y_0_0, tpz_y_x_0, tpz_y_y_0, tpz_y_z_0, \
                                     tpz_z_0_0, tpz_z_x_0, tpz_z_y_0, tpz_z_z_0, ts_x_x_0, ts_x_y_0, ts_x_z_0, ts_y_x_0, \
                                     ts_y_y_0, ts_y_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xx_x_0[j] = pa_x[j] * tpx_x_x_0[j] + 0.5 * fl1_fx * tpx_0_x_0[j] + 0.5 * fl1_fx * tpx_x_0_0[j] - fl1_fgb * fl1_fx * ts_x_x_0[j];

            tpy_xx_x_0[j] = pa_x[j] * tpy_x_x_0[j] + 0.5 * fl1_fx * tpy_0_x_0[j] + 0.5 * fl1_fx * tpy_x_0_0[j];

            tpz_xx_x_0[j] = pa_x[j] * tpz_x_x_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j] + 0.5 * fl1_fx * tpz_x_0_0[j];

            tpx_xx_y_0[j] = pa_x[j] * tpx_x_y_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j] - fl1_fgb * fl1_fx * ts_x_y_0[j];

            tpy_xx_y_0[j] = pa_x[j] * tpy_x_y_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j];

            tpz_xx_y_0[j] = pa_x[j] * tpz_x_y_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j];

            tpx_xx_z_0[j] = pa_x[j] * tpx_x_z_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j] - fl1_fgb * fl1_fx * ts_x_z_0[j];

            tpy_xx_z_0[j] = pa_x[j] * tpy_x_z_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j];

            tpz_xx_z_0[j] = pa_x[j] * tpz_x_z_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j];

            tpx_xy_x_0[j] = pa_x[j] * tpx_y_x_0[j] + 0.5 * fl1_fx * tpx_y_0_0[j] - fl1_fgb * fl1_fx * ts_y_x_0[j];

            tpy_xy_x_0[j] = pa_x[j] * tpy_y_x_0[j] + 0.5 * fl1_fx * tpy_y_0_0[j];

            tpz_xy_x_0[j] = pa_x[j] * tpz_y_x_0[j] + 0.5 * fl1_fx * tpz_y_0_0[j];

            tpx_xy_y_0[j] = pa_x[j] * tpx_y_y_0[j] - fl1_fgb * fl1_fx * ts_y_y_0[j];

            tpy_xy_y_0[j] = pa_x[j] * tpy_y_y_0[j];

            tpz_xy_y_0[j] = pa_x[j] * tpz_y_y_0[j];

            tpx_xy_z_0[j] = pa_x[j] * tpx_y_z_0[j] - fl1_fgb * fl1_fx * ts_y_z_0[j];

            tpy_xy_z_0[j] = pa_x[j] * tpy_y_z_0[j];

            tpz_xy_z_0[j] = pa_x[j] * tpz_y_z_0[j];

            tpx_xz_x_0[j] = pa_x[j] * tpx_z_x_0[j] + 0.5 * fl1_fx * tpx_z_0_0[j] - fl1_fgb * fl1_fx * ts_z_x_0[j];

            tpy_xz_x_0[j] = pa_x[j] * tpy_z_x_0[j] + 0.5 * fl1_fx * tpy_z_0_0[j];

            tpz_xz_x_0[j] = pa_x[j] * tpz_z_x_0[j] + 0.5 * fl1_fx * tpz_z_0_0[j];

            tpx_xz_y_0[j] = pa_x[j] * tpx_z_y_0[j] - fl1_fgb * fl1_fx * ts_z_y_0[j];

            tpy_xz_y_0[j] = pa_x[j] * tpy_z_y_0[j];

            tpz_xz_y_0[j] = pa_x[j] * tpz_z_y_0[j];

            tpx_xz_z_0[j] = pa_x[j] * tpx_z_z_0[j] - fl1_fgb * fl1_fx * ts_z_z_0[j];

            tpy_xz_z_0[j] = pa_x[j] * tpy_z_z_0[j];

            tpz_xz_z_0[j] = pa_x[j] * tpz_z_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForDP_27_54(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_2_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tpx_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx);

        auto tpy_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx);

        auto tpz_0_x_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx);

        auto tpx_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 1);

        auto tpy_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_0_y_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * idx + 2);

        auto tpy_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_0_z_0 = primBuffer.data(pidx_p_0_1_m0 + 6 * bdim + 3 * idx + 2);

        auto tpx_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 1);

        auto tpy_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 1);

        auto tpz_y_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 1);

        auto tpx_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * idx + 2);

        auto tpy_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 3 * bdim + 3 * idx + 2);

        auto tpz_z_0_0 = primBuffer.data(pidx_p_1_0_m0 + 6 * bdim + 3 * idx + 2);

        auto ts_y_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 3);

        auto ts_y_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 4);

        auto ts_y_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 5);

        auto ts_z_x_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 6);

        auto ts_z_y_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 7);

        auto ts_z_z_0 = primBuffer.data(pidx_s_1_1_m0 + 9 * idx + 8);

        // set up pointers to integrals

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 15);

        auto tpy_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 16);

        auto tpy_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 17);

        auto tpy_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 17);

        // Batch of Integrals (27,54)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_0_x_0, tpx_0_y_0, tpx_0_z_0, tpx_y_0_0, tpx_y_x_0, \
                                     tpx_y_y_0, tpx_y_z_0, tpx_yy_x_0, tpx_yy_y_0, tpx_yy_z_0, tpx_yz_x_0, tpx_yz_y_0, \
                                     tpx_yz_z_0, tpx_z_0_0, tpx_z_x_0, tpx_z_y_0, tpx_z_z_0, tpx_zz_x_0, tpx_zz_y_0, \
                                     tpx_zz_z_0, tpy_0_x_0, tpy_0_y_0, tpy_0_z_0, tpy_y_0_0, tpy_y_x_0, tpy_y_y_0, \
                                     tpy_y_z_0, tpy_yy_x_0, tpy_yy_y_0, tpy_yy_z_0, tpy_yz_x_0, tpy_yz_y_0, tpy_yz_z_0, \
                                     tpy_z_0_0, tpy_z_x_0, tpy_z_y_0, tpy_z_z_0, tpy_zz_x_0, tpy_zz_y_0, tpy_zz_z_0, \
                                     tpz_0_x_0, tpz_0_y_0, tpz_0_z_0, tpz_y_0_0, tpz_y_x_0, tpz_y_y_0, tpz_y_z_0, \
                                     tpz_yy_x_0, tpz_yy_y_0, tpz_yy_z_0, tpz_yz_x_0, tpz_yz_y_0, tpz_yz_z_0, tpz_z_0_0, \
                                     tpz_z_x_0, tpz_z_y_0, tpz_z_z_0, tpz_zz_x_0, tpz_zz_y_0, tpz_zz_z_0, ts_y_x_0, \
                                     ts_y_y_0, ts_y_z_0, ts_z_x_0, ts_z_y_0, ts_z_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yy_x_0[j] = pa_y[j] * tpx_y_x_0[j] + 0.5 * fl1_fx * tpx_0_x_0[j];

            tpy_yy_x_0[j] = pa_y[j] * tpy_y_x_0[j] + 0.5 * fl1_fx * tpy_0_x_0[j] - fl1_fgb * fl1_fx * ts_y_x_0[j];

            tpz_yy_x_0[j] = pa_y[j] * tpz_y_x_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j];

            tpx_yy_y_0[j] = pa_y[j] * tpx_y_y_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j] + 0.5 * fl1_fx * tpx_y_0_0[j];

            tpy_yy_y_0[j] = pa_y[j] * tpy_y_y_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j] + 0.5 * fl1_fx * tpy_y_0_0[j] - fl1_fgb * fl1_fx * ts_y_y_0[j];

            tpz_yy_y_0[j] = pa_y[j] * tpz_y_y_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j] + 0.5 * fl1_fx * tpz_y_0_0[j];

            tpx_yy_z_0[j] = pa_y[j] * tpx_y_z_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j];

            tpy_yy_z_0[j] = pa_y[j] * tpy_y_z_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j] - fl1_fgb * fl1_fx * ts_y_z_0[j];

            tpz_yy_z_0[j] = pa_y[j] * tpz_y_z_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j];

            tpx_yz_x_0[j] = pa_y[j] * tpx_z_x_0[j];

            tpy_yz_x_0[j] = pa_y[j] * tpy_z_x_0[j] - fl1_fgb * fl1_fx * ts_z_x_0[j];

            tpz_yz_x_0[j] = pa_y[j] * tpz_z_x_0[j];

            tpx_yz_y_0[j] = pa_y[j] * tpx_z_y_0[j] + 0.5 * fl1_fx * tpx_z_0_0[j];

            tpy_yz_y_0[j] = pa_y[j] * tpy_z_y_0[j] + 0.5 * fl1_fx * tpy_z_0_0[j] - fl1_fgb * fl1_fx * ts_z_y_0[j];

            tpz_yz_y_0[j] = pa_y[j] * tpz_z_y_0[j] + 0.5 * fl1_fx * tpz_z_0_0[j];

            tpx_yz_z_0[j] = pa_y[j] * tpx_z_z_0[j];

            tpy_yz_z_0[j] = pa_y[j] * tpy_z_z_0[j] - fl1_fgb * fl1_fx * ts_z_z_0[j];

            tpz_yz_z_0[j] = pa_y[j] * tpz_z_z_0[j];

            tpx_zz_x_0[j] = pa_z[j] * tpx_z_x_0[j] + 0.5 * fl1_fx * tpx_0_x_0[j];

            tpy_zz_x_0[j] = pa_z[j] * tpy_z_x_0[j] + 0.5 * fl1_fx * tpy_0_x_0[j];

            tpz_zz_x_0[j] = pa_z[j] * tpz_z_x_0[j] + 0.5 * fl1_fx * tpz_0_x_0[j] - fl1_fgb * fl1_fx * ts_z_x_0[j];

            tpx_zz_y_0[j] = pa_z[j] * tpx_z_y_0[j] + 0.5 * fl1_fx * tpx_0_y_0[j];

            tpy_zz_y_0[j] = pa_z[j] * tpy_z_y_0[j] + 0.5 * fl1_fx * tpy_0_y_0[j];

            tpz_zz_y_0[j] = pa_z[j] * tpz_z_y_0[j] + 0.5 * fl1_fx * tpz_0_y_0[j] - fl1_fgb * fl1_fx * ts_z_y_0[j];

            tpx_zz_z_0[j] = pa_z[j] * tpx_z_z_0[j] + 0.5 * fl1_fx * tpx_0_z_0[j] + 0.5 * fl1_fx * tpx_z_0_0[j];

            tpy_zz_z_0[j] = pa_z[j] * tpy_z_z_0[j] + 0.5 * fl1_fx * tpy_0_z_0[j] + 0.5 * fl1_fx * tpy_z_0_0[j];

            tpz_zz_z_0[j] = pa_z[j] * tpz_z_z_0[j] + 0.5 * fl1_fx * tpz_0_z_0[j] + 0.5 * fl1_fx * tpz_z_0_0[j] - fl1_fgb * fl1_fx * ts_z_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPF(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForPF_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForPF_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForPF_0_45(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 5);

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

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

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

        auto ts_0_xxx_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx);

        auto ts_0_xxy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 1);

        auto ts_0_xxz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 2);

        auto ts_0_xyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 3);

        auto ts_0_xyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 4);

        auto ts_0_xzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 5);

        auto ts_0_yyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 6);

        auto ts_0_yyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 7);

        auto ts_0_yzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 8);

        auto ts_0_zzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 9);

        // set up pointers to integrals

        auto tpx_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx);

        auto tpy_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx);

        auto tpz_x_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx);

        auto tpx_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 1);

        auto tpy_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_x_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 2);

        auto tpy_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_x_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 3);

        auto tpy_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_x_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 4);

        auto tpy_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_x_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 5);

        auto tpy_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_x_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 6);

        auto tpy_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_x_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 7);

        auto tpy_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_x_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 8);

        auto tpy_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_x_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 8);

        auto tpx_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 9);

        auto tpy_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_x_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 9);

        auto tpx_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 10);

        auto tpy_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_y_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 11);

        auto tpy_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_y_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 12);

        auto tpy_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_y_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 13);

        auto tpy_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_y_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 14);

        auto tpy_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_y_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, tpx_0_xx_0, tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, \
                                     tpx_0_xy_0, tpx_0_xyy_0, tpx_0_xyz_0, tpx_0_xz_0, tpx_0_xzz_0, tpx_0_yy_0, \
                                     tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yz_0, tpx_0_yzz_0, tpx_0_zz_0, tpx_0_zzz_0, \
                                     tpx_x_xxx_0, tpx_x_xxy_0, tpx_x_xxz_0, tpx_x_xyy_0, tpx_x_xyz_0, tpx_x_xzz_0, \
                                     tpx_x_yyy_0, tpx_x_yyz_0, tpx_x_yzz_0, tpx_x_zzz_0, tpx_y_xxx_0, tpx_y_xxy_0, \
                                     tpx_y_xxz_0, tpx_y_xyy_0, tpx_y_xyz_0, tpy_0_xx_0, tpy_0_xxx_0, tpy_0_xxy_0, \
                                     tpy_0_xxz_0, tpy_0_xy_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xz_0, tpy_0_xzz_0, \
                                     tpy_0_yy_0, tpy_0_yyy_0, tpy_0_yyz_0, tpy_0_yz_0, tpy_0_yzz_0, tpy_0_zz_0, \
                                     tpy_0_zzz_0, tpy_x_xxx_0, tpy_x_xxy_0, tpy_x_xxz_0, tpy_x_xyy_0, tpy_x_xyz_0, \
                                     tpy_x_xzz_0, tpy_x_yyy_0, tpy_x_yyz_0, tpy_x_yzz_0, tpy_x_zzz_0, tpy_y_xxx_0, \
                                     tpy_y_xxy_0, tpy_y_xxz_0, tpy_y_xyy_0, tpy_y_xyz_0, tpz_0_xx_0, tpz_0_xxx_0, \
                                     tpz_0_xxy_0, tpz_0_xxz_0, tpz_0_xy_0, tpz_0_xyy_0, tpz_0_xyz_0, tpz_0_xz_0, \
                                     tpz_0_xzz_0, tpz_0_yy_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yz_0, tpz_0_yzz_0, \
                                     tpz_0_zz_0, tpz_0_zzz_0, tpz_x_xxx_0, tpz_x_xxy_0, tpz_x_xxz_0, tpz_x_xyy_0, \
                                     tpz_x_xyz_0, tpz_x_xzz_0, tpz_x_yyy_0, tpz_x_yyz_0, tpz_x_yzz_0, tpz_x_zzz_0, \
                                     tpz_y_xxx_0, tpz_y_xxy_0, tpz_y_xxz_0, tpz_y_xyy_0, tpz_y_xyz_0, ts_0_xxx_0, \
                                     ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, \
                                     ts_0_yzz_0, ts_0_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_x_xxx_0[j] = pa_x[j] * tpx_0_xxx_0[j] + 1.5 * fl1_fx * tpx_0_xx_0[j] - fl1_fgb * fl1_fx * ts_0_xxx_0[j];

            tpy_x_xxx_0[j] = pa_x[j] * tpy_0_xxx_0[j] + 1.5 * fl1_fx * tpy_0_xx_0[j];

            tpz_x_xxx_0[j] = pa_x[j] * tpz_0_xxx_0[j] + 1.5 * fl1_fx * tpz_0_xx_0[j];

            tpx_x_xxy_0[j] = pa_x[j] * tpx_0_xxy_0[j] + fl1_fx * tpx_0_xy_0[j] - fl1_fgb * fl1_fx * ts_0_xxy_0[j];

            tpy_x_xxy_0[j] = pa_x[j] * tpy_0_xxy_0[j] + fl1_fx * tpy_0_xy_0[j];

            tpz_x_xxy_0[j] = pa_x[j] * tpz_0_xxy_0[j] + fl1_fx * tpz_0_xy_0[j];

            tpx_x_xxz_0[j] = pa_x[j] * tpx_0_xxz_0[j] + fl1_fx * tpx_0_xz_0[j] - fl1_fgb * fl1_fx * ts_0_xxz_0[j];

            tpy_x_xxz_0[j] = pa_x[j] * tpy_0_xxz_0[j] + fl1_fx * tpy_0_xz_0[j];

            tpz_x_xxz_0[j] = pa_x[j] * tpz_0_xxz_0[j] + fl1_fx * tpz_0_xz_0[j];

            tpx_x_xyy_0[j] = pa_x[j] * tpx_0_xyy_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j] - fl1_fgb * fl1_fx * ts_0_xyy_0[j];

            tpy_x_xyy_0[j] = pa_x[j] * tpy_0_xyy_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j];

            tpz_x_xyy_0[j] = pa_x[j] * tpz_0_xyy_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j];

            tpx_x_xyz_0[j] = pa_x[j] * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_0_yz_0[j] - fl1_fgb * fl1_fx * ts_0_xyz_0[j];

            tpy_x_xyz_0[j] = pa_x[j] * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_0_yz_0[j];

            tpz_x_xyz_0[j] = pa_x[j] * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_0_yz_0[j];

            tpx_x_xzz_0[j] = pa_x[j] * tpx_0_xzz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j] - fl1_fgb * fl1_fx * ts_0_xzz_0[j];

            tpy_x_xzz_0[j] = pa_x[j] * tpy_0_xzz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j];

            tpz_x_xzz_0[j] = pa_x[j] * tpz_0_xzz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j];

            tpx_x_yyy_0[j] = pa_x[j] * tpx_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_0_yyy_0[j];

            tpy_x_yyy_0[j] = pa_x[j] * tpy_0_yyy_0[j];

            tpz_x_yyy_0[j] = pa_x[j] * tpz_0_yyy_0[j];

            tpx_x_yyz_0[j] = pa_x[j] * tpx_0_yyz_0[j] - fl1_fgb * fl1_fx * ts_0_yyz_0[j];

            tpy_x_yyz_0[j] = pa_x[j] * tpy_0_yyz_0[j];

            tpz_x_yyz_0[j] = pa_x[j] * tpz_0_yyz_0[j];

            tpx_x_yzz_0[j] = pa_x[j] * tpx_0_yzz_0[j] - fl1_fgb * fl1_fx * ts_0_yzz_0[j];

            tpy_x_yzz_0[j] = pa_x[j] * tpy_0_yzz_0[j];

            tpz_x_yzz_0[j] = pa_x[j] * tpz_0_yzz_0[j];

            tpx_x_zzz_0[j] = pa_x[j] * tpx_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_0_zzz_0[j];

            tpy_x_zzz_0[j] = pa_x[j] * tpy_0_zzz_0[j];

            tpz_x_zzz_0[j] = pa_x[j] * tpz_0_zzz_0[j];

            tpx_y_xxx_0[j] = pa_y[j] * tpx_0_xxx_0[j];

            tpy_y_xxx_0[j] = pa_y[j] * tpy_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxx_0[j];

            tpz_y_xxx_0[j] = pa_y[j] * tpz_0_xxx_0[j];

            tpx_y_xxy_0[j] = pa_y[j] * tpx_0_xxy_0[j] + 0.5 * fl1_fx * tpx_0_xx_0[j];

            tpy_y_xxy_0[j] = pa_y[j] * tpy_0_xxy_0[j] + 0.5 * fl1_fx * tpy_0_xx_0[j] - fl1_fgb * fl1_fx * ts_0_xxy_0[j];

            tpz_y_xxy_0[j] = pa_y[j] * tpz_0_xxy_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j];

            tpx_y_xxz_0[j] = pa_y[j] * tpx_0_xxz_0[j];

            tpy_y_xxz_0[j] = pa_y[j] * tpy_0_xxz_0[j] - fl1_fgb * fl1_fx * ts_0_xxz_0[j];

            tpz_y_xxz_0[j] = pa_y[j] * tpz_0_xxz_0[j];

            tpx_y_xyy_0[j] = pa_y[j] * tpx_0_xyy_0[j] + fl1_fx * tpx_0_xy_0[j];

            tpy_y_xyy_0[j] = pa_y[j] * tpy_0_xyy_0[j] + fl1_fx * tpy_0_xy_0[j] - fl1_fgb * fl1_fx * ts_0_xyy_0[j];

            tpz_y_xyy_0[j] = pa_y[j] * tpz_0_xyy_0[j] + fl1_fx * tpz_0_xy_0[j];

            tpx_y_xyz_0[j] = pa_y[j] * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_0_xz_0[j];

            tpy_y_xyz_0[j] = pa_y[j] * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_0_xz_0[j] - fl1_fgb * fl1_fx * ts_0_xyz_0[j];

            tpz_y_xyz_0[j] = pa_y[j] * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_0_xz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPF_45_90(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 5);

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

        auto tpx_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx);

        auto tpy_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx);

        auto tpz_0_xx_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx);

        auto tpx_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 1);

        auto tpy_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_0_xy_0 = primBuffer.data(pidx_p_0_2_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_0_xz_0 = primBuffer.data(pidx_p_0_2_m0 + 6 * idx + 2);

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

        auto ts_0_xxx_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx);

        auto ts_0_xxy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 1);

        auto ts_0_xxz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 2);

        auto ts_0_xyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 3);

        auto ts_0_xyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 4);

        auto ts_0_xzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 5);

        auto ts_0_yyy_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 6);

        auto ts_0_yyz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 7);

        auto ts_0_yzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 8);

        auto ts_0_zzz_0 = primBuffer.data(pidx_s_0_3_m0 + 10 * idx + 9);

        // set up pointers to integrals

        auto tpx_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 15);

        auto tpy_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_y_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 16);

        auto tpy_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_y_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 17);

        auto tpy_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_y_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 18);

        auto tpy_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_y_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 19);

        auto tpy_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_y_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 20);

        auto tpy_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_z_xxx_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 21);

        auto tpy_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_z_xxy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 22);

        auto tpy_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_z_xxz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 23);

        auto tpy_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_z_xyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 24);

        auto tpy_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_z_xyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 25);

        auto tpy_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_z_xzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 26);

        auto tpy_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_z_yyy_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 27);

        auto tpy_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_z_yyz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 28);

        auto tpy_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_z_yzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * idx + 29);

        auto tpy_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_z_zzz_0 = primBuffer.data(pidx_p_1_3_m0 + 60 * bdim + 30 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_0_xx_0, tpx_0_xxx_0, tpx_0_xxy_0, tpx_0_xxz_0, \
                                     tpx_0_xy_0, tpx_0_xyy_0, tpx_0_xyz_0, tpx_0_xz_0, tpx_0_xzz_0, tpx_0_yy_0, \
                                     tpx_0_yyy_0, tpx_0_yyz_0, tpx_0_yz_0, tpx_0_yzz_0, tpx_0_zz_0, tpx_0_zzz_0, \
                                     tpx_y_xzz_0, tpx_y_yyy_0, tpx_y_yyz_0, tpx_y_yzz_0, tpx_y_zzz_0, tpx_z_xxx_0, \
                                     tpx_z_xxy_0, tpx_z_xxz_0, tpx_z_xyy_0, tpx_z_xyz_0, tpx_z_xzz_0, tpx_z_yyy_0, \
                                     tpx_z_yyz_0, tpx_z_yzz_0, tpx_z_zzz_0, tpy_0_xx_0, tpy_0_xxx_0, tpy_0_xxy_0, \
                                     tpy_0_xxz_0, tpy_0_xy_0, tpy_0_xyy_0, tpy_0_xyz_0, tpy_0_xz_0, tpy_0_xzz_0, \
                                     tpy_0_yy_0, tpy_0_yyy_0, tpy_0_yyz_0, tpy_0_yz_0, tpy_0_yzz_0, tpy_0_zz_0, \
                                     tpy_0_zzz_0, tpy_y_xzz_0, tpy_y_yyy_0, tpy_y_yyz_0, tpy_y_yzz_0, tpy_y_zzz_0, \
                                     tpy_z_xxx_0, tpy_z_xxy_0, tpy_z_xxz_0, tpy_z_xyy_0, tpy_z_xyz_0, tpy_z_xzz_0, \
                                     tpy_z_yyy_0, tpy_z_yyz_0, tpy_z_yzz_0, tpy_z_zzz_0, tpz_0_xx_0, tpz_0_xxx_0, \
                                     tpz_0_xxy_0, tpz_0_xxz_0, tpz_0_xy_0, tpz_0_xyy_0, tpz_0_xyz_0, tpz_0_xz_0, \
                                     tpz_0_xzz_0, tpz_0_yy_0, tpz_0_yyy_0, tpz_0_yyz_0, tpz_0_yz_0, tpz_0_yzz_0, \
                                     tpz_0_zz_0, tpz_0_zzz_0, tpz_y_xzz_0, tpz_y_yyy_0, tpz_y_yyz_0, tpz_y_yzz_0, \
                                     tpz_y_zzz_0, tpz_z_xxx_0, tpz_z_xxy_0, tpz_z_xxz_0, tpz_z_xyy_0, tpz_z_xyz_0, \
                                     tpz_z_xzz_0, tpz_z_yyy_0, tpz_z_yyz_0, tpz_z_yzz_0, tpz_z_zzz_0, ts_0_xxx_0, \
                                     ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, \
                                     ts_0_yzz_0, ts_0_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_y_xzz_0[j] = pa_y[j] * tpx_0_xzz_0[j];

            tpy_y_xzz_0[j] = pa_y[j] * tpy_0_xzz_0[j] - fl1_fgb * fl1_fx * ts_0_xzz_0[j];

            tpz_y_xzz_0[j] = pa_y[j] * tpz_0_xzz_0[j];

            tpx_y_yyy_0[j] = pa_y[j] * tpx_0_yyy_0[j] + 1.5 * fl1_fx * tpx_0_yy_0[j];

            tpy_y_yyy_0[j] = pa_y[j] * tpy_0_yyy_0[j] + 1.5 * fl1_fx * tpy_0_yy_0[j] - fl1_fgb * fl1_fx * ts_0_yyy_0[j];

            tpz_y_yyy_0[j] = pa_y[j] * tpz_0_yyy_0[j] + 1.5 * fl1_fx * tpz_0_yy_0[j];

            tpx_y_yyz_0[j] = pa_y[j] * tpx_0_yyz_0[j] + fl1_fx * tpx_0_yz_0[j];

            tpy_y_yyz_0[j] = pa_y[j] * tpy_0_yyz_0[j] + fl1_fx * tpy_0_yz_0[j] - fl1_fgb * fl1_fx * ts_0_yyz_0[j];

            tpz_y_yyz_0[j] = pa_y[j] * tpz_0_yyz_0[j] + fl1_fx * tpz_0_yz_0[j];

            tpx_y_yzz_0[j] = pa_y[j] * tpx_0_yzz_0[j] + 0.5 * fl1_fx * tpx_0_zz_0[j];

            tpy_y_yzz_0[j] = pa_y[j] * tpy_0_yzz_0[j] + 0.5 * fl1_fx * tpy_0_zz_0[j] - fl1_fgb * fl1_fx * ts_0_yzz_0[j];

            tpz_y_yzz_0[j] = pa_y[j] * tpz_0_yzz_0[j] + 0.5 * fl1_fx * tpz_0_zz_0[j];

            tpx_y_zzz_0[j] = pa_y[j] * tpx_0_zzz_0[j];

            tpy_y_zzz_0[j] = pa_y[j] * tpy_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_0_zzz_0[j];

            tpz_y_zzz_0[j] = pa_y[j] * tpz_0_zzz_0[j];

            tpx_z_xxx_0[j] = pa_z[j] * tpx_0_xxx_0[j];

            tpy_z_xxx_0[j] = pa_z[j] * tpy_0_xxx_0[j];

            tpz_z_xxx_0[j] = pa_z[j] * tpz_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxx_0[j];

            tpx_z_xxy_0[j] = pa_z[j] * tpx_0_xxy_0[j];

            tpy_z_xxy_0[j] = pa_z[j] * tpy_0_xxy_0[j];

            tpz_z_xxy_0[j] = pa_z[j] * tpz_0_xxy_0[j] - fl1_fgb * fl1_fx * ts_0_xxy_0[j];

            tpx_z_xxz_0[j] = pa_z[j] * tpx_0_xxz_0[j] + 0.5 * fl1_fx * tpx_0_xx_0[j];

            tpy_z_xxz_0[j] = pa_z[j] * tpy_0_xxz_0[j] + 0.5 * fl1_fx * tpy_0_xx_0[j];

            tpz_z_xxz_0[j] = pa_z[j] * tpz_0_xxz_0[j] + 0.5 * fl1_fx * tpz_0_xx_0[j] - fl1_fgb * fl1_fx * ts_0_xxz_0[j];

            tpx_z_xyy_0[j] = pa_z[j] * tpx_0_xyy_0[j];

            tpy_z_xyy_0[j] = pa_z[j] * tpy_0_xyy_0[j];

            tpz_z_xyy_0[j] = pa_z[j] * tpz_0_xyy_0[j] - fl1_fgb * fl1_fx * ts_0_xyy_0[j];

            tpx_z_xyz_0[j] = pa_z[j] * tpx_0_xyz_0[j] + 0.5 * fl1_fx * tpx_0_xy_0[j];

            tpy_z_xyz_0[j] = pa_z[j] * tpy_0_xyz_0[j] + 0.5 * fl1_fx * tpy_0_xy_0[j];

            tpz_z_xyz_0[j] = pa_z[j] * tpz_0_xyz_0[j] + 0.5 * fl1_fx * tpz_0_xy_0[j] - fl1_fgb * fl1_fx * ts_0_xyz_0[j];

            tpx_z_xzz_0[j] = pa_z[j] * tpx_0_xzz_0[j] + fl1_fx * tpx_0_xz_0[j];

            tpy_z_xzz_0[j] = pa_z[j] * tpy_0_xzz_0[j] + fl1_fx * tpy_0_xz_0[j];

            tpz_z_xzz_0[j] = pa_z[j] * tpz_0_xzz_0[j] + fl1_fx * tpz_0_xz_0[j] - fl1_fgb * fl1_fx * ts_0_xzz_0[j];

            tpx_z_yyy_0[j] = pa_z[j] * tpx_0_yyy_0[j];

            tpy_z_yyy_0[j] = pa_z[j] * tpy_0_yyy_0[j];

            tpz_z_yyy_0[j] = pa_z[j] * tpz_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_0_yyy_0[j];

            tpx_z_yyz_0[j] = pa_z[j] * tpx_0_yyz_0[j] + 0.5 * fl1_fx * tpx_0_yy_0[j];

            tpy_z_yyz_0[j] = pa_z[j] * tpy_0_yyz_0[j] + 0.5 * fl1_fx * tpy_0_yy_0[j];

            tpz_z_yyz_0[j] = pa_z[j] * tpz_0_yyz_0[j] + 0.5 * fl1_fx * tpz_0_yy_0[j] - fl1_fgb * fl1_fx * ts_0_yyz_0[j];

            tpx_z_yzz_0[j] = pa_z[j] * tpx_0_yzz_0[j] + fl1_fx * tpx_0_yz_0[j];

            tpy_z_yzz_0[j] = pa_z[j] * tpy_0_yzz_0[j] + fl1_fx * tpy_0_yz_0[j];

            tpz_z_yzz_0[j] = pa_z[j] * tpz_0_yzz_0[j] + fl1_fx * tpz_0_yz_0[j] - fl1_fgb * fl1_fx * ts_0_yzz_0[j];

            tpx_z_zzz_0[j] = pa_z[j] * tpx_0_zzz_0[j] + 1.5 * fl1_fx * tpx_0_zz_0[j];

            tpy_z_zzz_0[j] = pa_z[j] * tpy_0_zzz_0[j] + 1.5 * fl1_fx * tpy_0_zz_0[j];

            tpz_z_zzz_0[j] = pa_z[j] * tpz_0_zzz_0[j] + 1.5 * fl1_fx * tpz_0_zz_0[j] - fl1_fgb * fl1_fx * ts_0_zzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFP(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForFP_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForFP_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForFP_0_45(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx);

        auto tpy_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx);

        auto tpz_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx);

        auto tpx_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 1);

        auto tpy_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 2);

        auto tpy_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 3);

        auto tpy_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 4);

        auto tpy_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 5);

        auto tpy_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 6);

        auto tpy_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 7);

        auto tpy_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 8);

        auto tpy_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx);

        auto tpy_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx);

        auto tpz_x_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx);

        auto tpx_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 1);

        auto tpy_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 1);

        auto tpz_x_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 1);

        auto tpx_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 2);

        auto tpy_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 2);

        auto tpz_x_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 2);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tpx_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx);

        auto tpy_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx);

        auto tpz_xx_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx);

        auto tpx_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 1);

        auto tpy_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 1);

        auto tpz_xy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 1);

        auto tpx_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 2);

        auto tpy_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 2);

        auto tpz_xz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 2);

        auto tpx_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 3);

        auto tpy_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 4);

        auto tpy_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto ts_xx_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx);

        auto ts_xx_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 1);

        auto ts_xx_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 2);

        auto ts_xy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 3);

        auto ts_xy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 4);

        auto ts_xy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 5);

        auto ts_xz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 6);

        auto ts_xz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 7);

        auto ts_xz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 8);

        auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9);

        auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10);

        auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11);

        auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12);

        auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13);

        auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14);

        // set up pointers to integrals

        auto tpx_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx);

        auto tpy_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx);

        auto tpz_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx);

        auto tpx_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 1);

        auto tpy_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 2);

        auto tpy_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 3);

        auto tpy_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 4);

        auto tpy_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 5);

        auto tpy_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 6);

        auto tpy_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 7);

        auto tpy_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 8);

        auto tpy_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto tpx_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 9);

        auto tpy_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tpx_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 10);

        auto tpy_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 11);

        auto tpy_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 12);

        auto tpy_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 13);

        auto tpy_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 14);

        auto tpy_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_x_x_0, tpx_x_y_0, tpx_x_z_0, tpx_xx_0_0, tpx_xx_x_0, \
                                     tpx_xx_y_0, tpx_xx_z_0, tpx_xxx_x_0, tpx_xxx_y_0, tpx_xxx_z_0, tpx_xxy_x_0, \
                                     tpx_xxy_y_0, tpx_xxy_z_0, tpx_xxz_x_0, tpx_xxz_y_0, tpx_xxz_z_0, tpx_xy_0_0, \
                                     tpx_xy_x_0, tpx_xy_y_0, tpx_xy_z_0, tpx_xyy_x_0, tpx_xyy_y_0, tpx_xyy_z_0, \
                                     tpx_xyz_x_0, tpx_xyz_y_0, tpx_xyz_z_0, tpx_xz_0_0, tpx_xz_x_0, tpx_xz_y_0, \
                                     tpx_xz_z_0, tpx_y_x_0, tpx_y_y_0, tpx_y_z_0, tpx_yy_0_0, tpx_yy_x_0, tpx_yy_y_0, \
                                     tpx_yy_z_0, tpx_yz_0_0, tpx_yz_x_0, tpx_yz_y_0, tpx_yz_z_0, tpx_z_x_0, tpx_z_y_0, \
                                     tpx_z_z_0, tpy_x_x_0, tpy_x_y_0, tpy_x_z_0, tpy_xx_0_0, tpy_xx_x_0, tpy_xx_y_0, \
                                     tpy_xx_z_0, tpy_xxx_x_0, tpy_xxx_y_0, tpy_xxx_z_0, tpy_xxy_x_0, tpy_xxy_y_0, \
                                     tpy_xxy_z_0, tpy_xxz_x_0, tpy_xxz_y_0, tpy_xxz_z_0, tpy_xy_0_0, tpy_xy_x_0, \
                                     tpy_xy_y_0, tpy_xy_z_0, tpy_xyy_x_0, tpy_xyy_y_0, tpy_xyy_z_0, tpy_xyz_x_0, \
                                     tpy_xyz_y_0, tpy_xyz_z_0, tpy_xz_0_0, tpy_xz_x_0, tpy_xz_y_0, tpy_xz_z_0, tpy_y_x_0, \
                                     tpy_y_y_0, tpy_y_z_0, tpy_yy_0_0, tpy_yy_x_0, tpy_yy_y_0, tpy_yy_z_0, tpy_yz_0_0, \
                                     tpy_yz_x_0, tpy_yz_y_0, tpy_yz_z_0, tpy_z_x_0, tpy_z_y_0, tpy_z_z_0, tpz_x_x_0, \
                                     tpz_x_y_0, tpz_x_z_0, tpz_xx_0_0, tpz_xx_x_0, tpz_xx_y_0, tpz_xx_z_0, tpz_xxx_x_0, \
                                     tpz_xxx_y_0, tpz_xxx_z_0, tpz_xxy_x_0, tpz_xxy_y_0, tpz_xxy_z_0, tpz_xxz_x_0, \
                                     tpz_xxz_y_0, tpz_xxz_z_0, tpz_xy_0_0, tpz_xy_x_0, tpz_xy_y_0, tpz_xy_z_0, \
                                     tpz_xyy_x_0, tpz_xyy_y_0, tpz_xyy_z_0, tpz_xyz_x_0, tpz_xyz_y_0, tpz_xyz_z_0, \
                                     tpz_xz_0_0, tpz_xz_x_0, tpz_xz_y_0, tpz_xz_z_0, tpz_y_x_0, tpz_y_y_0, tpz_y_z_0, \
                                     tpz_yy_0_0, tpz_yy_x_0, tpz_yy_y_0, tpz_yy_z_0, tpz_yz_0_0, tpz_yz_x_0, tpz_yz_y_0, \
                                     tpz_yz_z_0, tpz_z_x_0, tpz_z_y_0, tpz_z_z_0, ts_xx_x_0, ts_xx_y_0, ts_xx_z_0, \
                                     ts_xy_x_0, ts_xy_y_0, ts_xy_z_0, ts_xz_x_0, ts_xz_y_0, ts_xz_z_0, ts_yy_x_0, \
                                     ts_yy_y_0, ts_yy_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxx_x_0[j] = pa_x[j] * tpx_xx_x_0[j] + fl1_fx * tpx_x_x_0[j] + 0.5 * fl1_fx * tpx_xx_0_0[j] - fl1_fgb * fl1_fx * ts_xx_x_0[j];

            tpy_xxx_x_0[j] = pa_x[j] * tpy_xx_x_0[j] + fl1_fx * tpy_x_x_0[j] + 0.5 * fl1_fx * tpy_xx_0_0[j];

            tpz_xxx_x_0[j] = pa_x[j] * tpz_xx_x_0[j] + fl1_fx * tpz_x_x_0[j] + 0.5 * fl1_fx * tpz_xx_0_0[j];

            tpx_xxx_y_0[j] = pa_x[j] * tpx_xx_y_0[j] + fl1_fx * tpx_x_y_0[j] - fl1_fgb * fl1_fx * ts_xx_y_0[j];

            tpy_xxx_y_0[j] = pa_x[j] * tpy_xx_y_0[j] + fl1_fx * tpy_x_y_0[j];

            tpz_xxx_y_0[j] = pa_x[j] * tpz_xx_y_0[j] + fl1_fx * tpz_x_y_0[j];

            tpx_xxx_z_0[j] = pa_x[j] * tpx_xx_z_0[j] + fl1_fx * tpx_x_z_0[j] - fl1_fgb * fl1_fx * ts_xx_z_0[j];

            tpy_xxx_z_0[j] = pa_x[j] * tpy_xx_z_0[j] + fl1_fx * tpy_x_z_0[j];

            tpz_xxx_z_0[j] = pa_x[j] * tpz_xx_z_0[j] + fl1_fx * tpz_x_z_0[j];

            tpx_xxy_x_0[j] = pa_x[j] * tpx_xy_x_0[j] + 0.5 * fl1_fx * tpx_y_x_0[j] + 0.5 * fl1_fx * tpx_xy_0_0[j] - fl1_fgb * fl1_fx * ts_xy_x_0[j];

            tpy_xxy_x_0[j] = pa_x[j] * tpy_xy_x_0[j] + 0.5 * fl1_fx * tpy_y_x_0[j] + 0.5 * fl1_fx * tpy_xy_0_0[j];

            tpz_xxy_x_0[j] = pa_x[j] * tpz_xy_x_0[j] + 0.5 * fl1_fx * tpz_y_x_0[j] + 0.5 * fl1_fx * tpz_xy_0_0[j];

            tpx_xxy_y_0[j] = pa_x[j] * tpx_xy_y_0[j] + 0.5 * fl1_fx * tpx_y_y_0[j] - fl1_fgb * fl1_fx * ts_xy_y_0[j];

            tpy_xxy_y_0[j] = pa_x[j] * tpy_xy_y_0[j] + 0.5 * fl1_fx * tpy_y_y_0[j];

            tpz_xxy_y_0[j] = pa_x[j] * tpz_xy_y_0[j] + 0.5 * fl1_fx * tpz_y_y_0[j];

            tpx_xxy_z_0[j] = pa_x[j] * tpx_xy_z_0[j] + 0.5 * fl1_fx * tpx_y_z_0[j] - fl1_fgb * fl1_fx * ts_xy_z_0[j];

            tpy_xxy_z_0[j] = pa_x[j] * tpy_xy_z_0[j] + 0.5 * fl1_fx * tpy_y_z_0[j];

            tpz_xxy_z_0[j] = pa_x[j] * tpz_xy_z_0[j] + 0.5 * fl1_fx * tpz_y_z_0[j];

            tpx_xxz_x_0[j] = pa_x[j] * tpx_xz_x_0[j] + 0.5 * fl1_fx * tpx_z_x_0[j] + 0.5 * fl1_fx * tpx_xz_0_0[j] - fl1_fgb * fl1_fx * ts_xz_x_0[j];

            tpy_xxz_x_0[j] = pa_x[j] * tpy_xz_x_0[j] + 0.5 * fl1_fx * tpy_z_x_0[j] + 0.5 * fl1_fx * tpy_xz_0_0[j];

            tpz_xxz_x_0[j] = pa_x[j] * tpz_xz_x_0[j] + 0.5 * fl1_fx * tpz_z_x_0[j] + 0.5 * fl1_fx * tpz_xz_0_0[j];

            tpx_xxz_y_0[j] = pa_x[j] * tpx_xz_y_0[j] + 0.5 * fl1_fx * tpx_z_y_0[j] - fl1_fgb * fl1_fx * ts_xz_y_0[j];

            tpy_xxz_y_0[j] = pa_x[j] * tpy_xz_y_0[j] + 0.5 * fl1_fx * tpy_z_y_0[j];

            tpz_xxz_y_0[j] = pa_x[j] * tpz_xz_y_0[j] + 0.5 * fl1_fx * tpz_z_y_0[j];

            tpx_xxz_z_0[j] = pa_x[j] * tpx_xz_z_0[j] + 0.5 * fl1_fx * tpx_z_z_0[j] - fl1_fgb * fl1_fx * ts_xz_z_0[j];

            tpy_xxz_z_0[j] = pa_x[j] * tpy_xz_z_0[j] + 0.5 * fl1_fx * tpy_z_z_0[j];

            tpz_xxz_z_0[j] = pa_x[j] * tpz_xz_z_0[j] + 0.5 * fl1_fx * tpz_z_z_0[j];

            tpx_xyy_x_0[j] = pa_x[j] * tpx_yy_x_0[j] + 0.5 * fl1_fx * tpx_yy_0_0[j] - fl1_fgb * fl1_fx * ts_yy_x_0[j];

            tpy_xyy_x_0[j] = pa_x[j] * tpy_yy_x_0[j] + 0.5 * fl1_fx * tpy_yy_0_0[j];

            tpz_xyy_x_0[j] = pa_x[j] * tpz_yy_x_0[j] + 0.5 * fl1_fx * tpz_yy_0_0[j];

            tpx_xyy_y_0[j] = pa_x[j] * tpx_yy_y_0[j] - fl1_fgb * fl1_fx * ts_yy_y_0[j];

            tpy_xyy_y_0[j] = pa_x[j] * tpy_yy_y_0[j];

            tpz_xyy_y_0[j] = pa_x[j] * tpz_yy_y_0[j];

            tpx_xyy_z_0[j] = pa_x[j] * tpx_yy_z_0[j] - fl1_fgb * fl1_fx * ts_yy_z_0[j];

            tpy_xyy_z_0[j] = pa_x[j] * tpy_yy_z_0[j];

            tpz_xyy_z_0[j] = pa_x[j] * tpz_yy_z_0[j];

            tpx_xyz_x_0[j] = pa_x[j] * tpx_yz_x_0[j] + 0.5 * fl1_fx * tpx_yz_0_0[j] - fl1_fgb * fl1_fx * ts_yz_x_0[j];

            tpy_xyz_x_0[j] = pa_x[j] * tpy_yz_x_0[j] + 0.5 * fl1_fx * tpy_yz_0_0[j];

            tpz_xyz_x_0[j] = pa_x[j] * tpz_yz_x_0[j] + 0.5 * fl1_fx * tpz_yz_0_0[j];

            tpx_xyz_y_0[j] = pa_x[j] * tpx_yz_y_0[j] - fl1_fgb * fl1_fx * ts_yz_y_0[j];

            tpy_xyz_y_0[j] = pa_x[j] * tpy_yz_y_0[j];

            tpz_xyz_y_0[j] = pa_x[j] * tpz_yz_y_0[j];

            tpx_xyz_z_0[j] = pa_x[j] * tpx_yz_z_0[j] - fl1_fgb * fl1_fx * ts_yz_z_0[j];

            tpy_xyz_z_0[j] = pa_x[j] * tpy_yz_z_0[j];

            tpz_xyz_z_0[j] = pa_x[j] * tpz_yz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForFP_45_90(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_3_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 15);

        auto tpy_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 16);

        auto tpy_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 17);

        auto tpy_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto tpx_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 3);

        auto tpy_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 3);

        auto tpz_y_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 3);

        auto tpx_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 4);

        auto tpy_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 4);

        auto tpz_y_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 4);

        auto tpx_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 5);

        auto tpy_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 5);

        auto tpz_y_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 5);

        auto tpx_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 6);

        auto tpy_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 6);

        auto tpz_z_x_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 6);

        auto tpx_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 7);

        auto tpy_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 7);

        auto tpz_z_y_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 7);

        auto tpx_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * idx + 8);

        auto tpy_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 9 * bdim + 9 * idx + 8);

        auto tpz_z_z_0 = primBuffer.data(pidx_p_1_1_m0 + 18 * bdim + 9 * idx + 8);

        auto tpx_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 3);

        auto tpy_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 3);

        auto tpz_yy_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 3);

        auto tpx_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 4);

        auto tpy_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 4);

        auto tpz_yz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 4);

        auto tpx_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * idx + 5);

        auto tpy_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 6 * bdim + 6 * idx + 5);

        auto tpz_zz_0_0 = primBuffer.data(pidx_p_2_0_m0 + 12 * bdim + 6 * idx + 5);

        auto ts_yy_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 9);

        auto ts_yy_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 10);

        auto ts_yy_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 11);

        auto ts_yz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 12);

        auto ts_yz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 13);

        auto ts_yz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 14);

        auto ts_zz_x_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 15);

        auto ts_zz_y_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 16);

        auto ts_zz_z_0 = primBuffer.data(pidx_s_2_1_m0 + 18 * idx + 17);

        // set up pointers to integrals

        auto tpx_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 15);

        auto tpy_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 16);

        auto tpy_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 17);

        auto tpy_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 18);

        auto tpy_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 19);

        auto tpy_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 20);

        auto tpy_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 21);

        auto tpy_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 22);

        auto tpy_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 23);

        auto tpy_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 24);

        auto tpy_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 25);

        auto tpy_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 26);

        auto tpy_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 27);

        auto tpy_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 28);

        auto tpy_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 29);

        auto tpy_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, pa_y, pa_z, tpx_xzz_x_0, tpx_xzz_y_0, tpx_xzz_z_0, tpx_y_x_0, \
                                     tpx_y_y_0, tpx_y_z_0, tpx_yy_0_0, tpx_yy_x_0, tpx_yy_y_0, tpx_yy_z_0, tpx_yyy_x_0, \
                                     tpx_yyy_y_0, tpx_yyy_z_0, tpx_yyz_x_0, tpx_yyz_y_0, tpx_yyz_z_0, tpx_yz_0_0, \
                                     tpx_yz_x_0, tpx_yz_y_0, tpx_yz_z_0, tpx_yzz_x_0, tpx_yzz_y_0, tpx_yzz_z_0, \
                                     tpx_z_x_0, tpx_z_y_0, tpx_z_z_0, tpx_zz_0_0, tpx_zz_x_0, tpx_zz_y_0, tpx_zz_z_0, \
                                     tpx_zzz_x_0, tpx_zzz_y_0, tpx_zzz_z_0, tpy_xzz_x_0, tpy_xzz_y_0, tpy_xzz_z_0, \
                                     tpy_y_x_0, tpy_y_y_0, tpy_y_z_0, tpy_yy_0_0, tpy_yy_x_0, tpy_yy_y_0, tpy_yy_z_0, \
                                     tpy_yyy_x_0, tpy_yyy_y_0, tpy_yyy_z_0, tpy_yyz_x_0, tpy_yyz_y_0, tpy_yyz_z_0, \
                                     tpy_yz_0_0, tpy_yz_x_0, tpy_yz_y_0, tpy_yz_z_0, tpy_yzz_x_0, tpy_yzz_y_0, \
                                     tpy_yzz_z_0, tpy_z_x_0, tpy_z_y_0, tpy_z_z_0, tpy_zz_0_0, tpy_zz_x_0, tpy_zz_y_0, \
                                     tpy_zz_z_0, tpy_zzz_x_0, tpy_zzz_y_0, tpy_zzz_z_0, tpz_xzz_x_0, tpz_xzz_y_0, \
                                     tpz_xzz_z_0, tpz_y_x_0, tpz_y_y_0, tpz_y_z_0, tpz_yy_0_0, tpz_yy_x_0, tpz_yy_y_0, \
                                     tpz_yy_z_0, tpz_yyy_x_0, tpz_yyy_y_0, tpz_yyy_z_0, tpz_yyz_x_0, tpz_yyz_y_0, \
                                     tpz_yyz_z_0, tpz_yz_0_0, tpz_yz_x_0, tpz_yz_y_0, tpz_yz_z_0, tpz_yzz_x_0, \
                                     tpz_yzz_y_0, tpz_yzz_z_0, tpz_z_x_0, tpz_z_y_0, tpz_z_z_0, tpz_zz_0_0, tpz_zz_x_0, \
                                     tpz_zz_y_0, tpz_zz_z_0, tpz_zzz_x_0, tpz_zzz_y_0, tpz_zzz_z_0, ts_yy_x_0, \
                                     ts_yy_y_0, ts_yy_z_0, ts_yz_x_0, ts_yz_y_0, ts_yz_z_0, ts_zz_x_0, ts_zz_y_0, \
                                     ts_zz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xzz_x_0[j] = pa_x[j] * tpx_zz_x_0[j] + 0.5 * fl1_fx * tpx_zz_0_0[j] - fl1_fgb * fl1_fx * ts_zz_x_0[j];

            tpy_xzz_x_0[j] = pa_x[j] * tpy_zz_x_0[j] + 0.5 * fl1_fx * tpy_zz_0_0[j];

            tpz_xzz_x_0[j] = pa_x[j] * tpz_zz_x_0[j] + 0.5 * fl1_fx * tpz_zz_0_0[j];

            tpx_xzz_y_0[j] = pa_x[j] * tpx_zz_y_0[j] - fl1_fgb * fl1_fx * ts_zz_y_0[j];

            tpy_xzz_y_0[j] = pa_x[j] * tpy_zz_y_0[j];

            tpz_xzz_y_0[j] = pa_x[j] * tpz_zz_y_0[j];

            tpx_xzz_z_0[j] = pa_x[j] * tpx_zz_z_0[j] - fl1_fgb * fl1_fx * ts_zz_z_0[j];

            tpy_xzz_z_0[j] = pa_x[j] * tpy_zz_z_0[j];

            tpz_xzz_z_0[j] = pa_x[j] * tpz_zz_z_0[j];

            tpx_yyy_x_0[j] = pa_y[j] * tpx_yy_x_0[j] + fl1_fx * tpx_y_x_0[j];

            tpy_yyy_x_0[j] = pa_y[j] * tpy_yy_x_0[j] + fl1_fx * tpy_y_x_0[j] - fl1_fgb * fl1_fx * ts_yy_x_0[j];

            tpz_yyy_x_0[j] = pa_y[j] * tpz_yy_x_0[j] + fl1_fx * tpz_y_x_0[j];

            tpx_yyy_y_0[j] = pa_y[j] * tpx_yy_y_0[j] + fl1_fx * tpx_y_y_0[j] + 0.5 * fl1_fx * tpx_yy_0_0[j];

            tpy_yyy_y_0[j] = pa_y[j] * tpy_yy_y_0[j] + fl1_fx * tpy_y_y_0[j] + 0.5 * fl1_fx * tpy_yy_0_0[j] - fl1_fgb * fl1_fx * ts_yy_y_0[j];

            tpz_yyy_y_0[j] = pa_y[j] * tpz_yy_y_0[j] + fl1_fx * tpz_y_y_0[j] + 0.5 * fl1_fx * tpz_yy_0_0[j];

            tpx_yyy_z_0[j] = pa_y[j] * tpx_yy_z_0[j] + fl1_fx * tpx_y_z_0[j];

            tpy_yyy_z_0[j] = pa_y[j] * tpy_yy_z_0[j] + fl1_fx * tpy_y_z_0[j] - fl1_fgb * fl1_fx * ts_yy_z_0[j];

            tpz_yyy_z_0[j] = pa_y[j] * tpz_yy_z_0[j] + fl1_fx * tpz_y_z_0[j];

            tpx_yyz_x_0[j] = pa_y[j] * tpx_yz_x_0[j] + 0.5 * fl1_fx * tpx_z_x_0[j];

            tpy_yyz_x_0[j] = pa_y[j] * tpy_yz_x_0[j] + 0.5 * fl1_fx * tpy_z_x_0[j] - fl1_fgb * fl1_fx * ts_yz_x_0[j];

            tpz_yyz_x_0[j] = pa_y[j] * tpz_yz_x_0[j] + 0.5 * fl1_fx * tpz_z_x_0[j];

            tpx_yyz_y_0[j] = pa_y[j] * tpx_yz_y_0[j] + 0.5 * fl1_fx * tpx_z_y_0[j] + 0.5 * fl1_fx * tpx_yz_0_0[j];

            tpy_yyz_y_0[j] = pa_y[j] * tpy_yz_y_0[j] + 0.5 * fl1_fx * tpy_z_y_0[j] + 0.5 * fl1_fx * tpy_yz_0_0[j] - fl1_fgb * fl1_fx * ts_yz_y_0[j];

            tpz_yyz_y_0[j] = pa_y[j] * tpz_yz_y_0[j] + 0.5 * fl1_fx * tpz_z_y_0[j] + 0.5 * fl1_fx * tpz_yz_0_0[j];

            tpx_yyz_z_0[j] = pa_y[j] * tpx_yz_z_0[j] + 0.5 * fl1_fx * tpx_z_z_0[j];

            tpy_yyz_z_0[j] = pa_y[j] * tpy_yz_z_0[j] + 0.5 * fl1_fx * tpy_z_z_0[j] - fl1_fgb * fl1_fx * ts_yz_z_0[j];

            tpz_yyz_z_0[j] = pa_y[j] * tpz_yz_z_0[j] + 0.5 * fl1_fx * tpz_z_z_0[j];

            tpx_yzz_x_0[j] = pa_y[j] * tpx_zz_x_0[j];

            tpy_yzz_x_0[j] = pa_y[j] * tpy_zz_x_0[j] - fl1_fgb * fl1_fx * ts_zz_x_0[j];

            tpz_yzz_x_0[j] = pa_y[j] * tpz_zz_x_0[j];

            tpx_yzz_y_0[j] = pa_y[j] * tpx_zz_y_0[j] + 0.5 * fl1_fx * tpx_zz_0_0[j];

            tpy_yzz_y_0[j] = pa_y[j] * tpy_zz_y_0[j] + 0.5 * fl1_fx * tpy_zz_0_0[j] - fl1_fgb * fl1_fx * ts_zz_y_0[j];

            tpz_yzz_y_0[j] = pa_y[j] * tpz_zz_y_0[j] + 0.5 * fl1_fx * tpz_zz_0_0[j];

            tpx_yzz_z_0[j] = pa_y[j] * tpx_zz_z_0[j];

            tpy_yzz_z_0[j] = pa_y[j] * tpy_zz_z_0[j] - fl1_fgb * fl1_fx * ts_zz_z_0[j];

            tpz_yzz_z_0[j] = pa_y[j] * tpz_zz_z_0[j];

            tpx_zzz_x_0[j] = pa_z[j] * tpx_zz_x_0[j] + fl1_fx * tpx_z_x_0[j];

            tpy_zzz_x_0[j] = pa_z[j] * tpy_zz_x_0[j] + fl1_fx * tpy_z_x_0[j];

            tpz_zzz_x_0[j] = pa_z[j] * tpz_zz_x_0[j] + fl1_fx * tpz_z_x_0[j] - fl1_fgb * fl1_fx * ts_zz_x_0[j];

            tpx_zzz_y_0[j] = pa_z[j] * tpx_zz_y_0[j] + fl1_fx * tpx_z_y_0[j];

            tpy_zzz_y_0[j] = pa_z[j] * tpy_zz_y_0[j] + fl1_fx * tpy_z_y_0[j];

            tpz_zzz_y_0[j] = pa_z[j] * tpz_zz_y_0[j] + fl1_fx * tpz_z_y_0[j] - fl1_fgb * fl1_fx * ts_zz_y_0[j];

            tpx_zzz_z_0[j] = pa_z[j] * tpx_zz_z_0[j] + fl1_fx * tpx_z_z_0[j] + 0.5 * fl1_fx * tpx_zz_0_0[j];

            tpy_zzz_z_0[j] = pa_z[j] * tpy_zz_z_0[j] + fl1_fx * tpy_z_z_0[j] + 0.5 * fl1_fx * tpy_zz_0_0[j];

            tpz_zzz_z_0[j] = pa_z[j] * tpz_zz_z_0[j] + fl1_fx * tpz_z_z_0[j] + 0.5 * fl1_fx * tpz_zz_0_0[j] - fl1_fgb * fl1_fx * ts_zz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPG(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForPG_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForPG_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForPG_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForPG_0_45(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 5);

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

        auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx);

        auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1);

        auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2);

        auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3);

        auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4);

        auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5);

        auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6);

        auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7);

        auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8);

        auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9);

        auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10);

        auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11);

        auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12);

        auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13);

        auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14);

        // set up pointers to integrals

        auto tpx_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx);

        auto tpy_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx);

        auto tpz_x_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx);

        auto tpx_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 1);

        auto tpy_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 1);

        auto tpz_x_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 1);

        auto tpx_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 2);

        auto tpy_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 2);

        auto tpz_x_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 2);

        auto tpx_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 3);

        auto tpy_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 3);

        auto tpz_x_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 3);

        auto tpx_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 4);

        auto tpy_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 4);

        auto tpz_x_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 4);

        auto tpx_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 5);

        auto tpy_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 5);

        auto tpz_x_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 5);

        auto tpx_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 6);

        auto tpy_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 6);

        auto tpz_x_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 6);

        auto tpx_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 7);

        auto tpy_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 7);

        auto tpz_x_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 7);

        auto tpx_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 8);

        auto tpy_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 8);

        auto tpz_x_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 8);

        auto tpx_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 9);

        auto tpy_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 9);

        auto tpz_x_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 9);

        auto tpx_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 10);

        auto tpy_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 10);

        auto tpz_x_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 10);

        auto tpx_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 11);

        auto tpy_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 11);

        auto tpz_x_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 11);

        auto tpx_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 12);

        auto tpy_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 12);

        auto tpz_x_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 12);

        auto tpx_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 13);

        auto tpy_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 13);

        auto tpz_x_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 13);

        auto tpx_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 14);

        auto tpy_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 14);

        auto tpz_x_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_0_xxx_0, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, \
                                     tpx_0_xxy_0, tpx_0_xxyy_0, tpx_0_xxyz_0, tpx_0_xxz_0, tpx_0_xxzz_0, tpx_0_xyy_0, \
                                     tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyz_0, tpx_0_xyzz_0, tpx_0_xzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyy_0, tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyz_0, tpx_0_yyzz_0, tpx_0_yzz_0, \
                                     tpx_0_yzzz_0, tpx_0_zzz_0, tpx_0_zzzz_0, tpx_x_xxxx_0, tpx_x_xxxy_0, tpx_x_xxxz_0, \
                                     tpx_x_xxyy_0, tpx_x_xxyz_0, tpx_x_xxzz_0, tpx_x_xyyy_0, tpx_x_xyyz_0, tpx_x_xyzz_0, \
                                     tpx_x_xzzz_0, tpx_x_yyyy_0, tpx_x_yyyz_0, tpx_x_yyzz_0, tpx_x_yzzz_0, tpx_x_zzzz_0, \
                                     tpy_0_xxx_0, tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxy_0, tpy_0_xxyy_0, \
                                     tpy_0_xxyz_0, tpy_0_xxz_0, tpy_0_xxzz_0, tpy_0_xyy_0, tpy_0_xyyy_0, tpy_0_xyyz_0, \
                                     tpy_0_xyz_0, tpy_0_xyzz_0, tpy_0_xzz_0, tpy_0_xzzz_0, tpy_0_yyy_0, tpy_0_yyyy_0, \
                                     tpy_0_yyyz_0, tpy_0_yyz_0, tpy_0_yyzz_0, tpy_0_yzz_0, tpy_0_yzzz_0, tpy_0_zzz_0, \
                                     tpy_0_zzzz_0, tpy_x_xxxx_0, tpy_x_xxxy_0, tpy_x_xxxz_0, tpy_x_xxyy_0, tpy_x_xxyz_0, \
                                     tpy_x_xxzz_0, tpy_x_xyyy_0, tpy_x_xyyz_0, tpy_x_xyzz_0, tpy_x_xzzz_0, tpy_x_yyyy_0, \
                                     tpy_x_yyyz_0, tpy_x_yyzz_0, tpy_x_yzzz_0, tpy_x_zzzz_0, tpz_0_xxx_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxy_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxz_0, \
                                     tpz_0_xxzz_0, tpz_0_xyy_0, tpz_0_xyyy_0, tpz_0_xyyz_0, tpz_0_xyz_0, tpz_0_xyzz_0, \
                                     tpz_0_xzz_0, tpz_0_xzzz_0, tpz_0_yyy_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyz_0, \
                                     tpz_0_yyzz_0, tpz_0_yzz_0, tpz_0_yzzz_0, tpz_0_zzz_0, tpz_0_zzzz_0, tpz_x_xxxx_0, \
                                     tpz_x_xxxy_0, tpz_x_xxxz_0, tpz_x_xxyy_0, tpz_x_xxyz_0, tpz_x_xxzz_0, tpz_x_xyyy_0, \
                                     tpz_x_xyyz_0, tpz_x_xyzz_0, tpz_x_xzzz_0, tpz_x_yyyy_0, tpz_x_yyyz_0, tpz_x_yyzz_0, \
                                     tpz_x_yzzz_0, tpz_x_zzzz_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_x_xxxx_0[j] = pa_x[j] * tpx_0_xxxx_0[j] + 2.0 * fl1_fx * tpx_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxxx_0[j];

            tpy_x_xxxx_0[j] = pa_x[j] * tpy_0_xxxx_0[j] + 2.0 * fl1_fx * tpy_0_xxx_0[j];

            tpz_x_xxxx_0[j] = pa_x[j] * tpz_0_xxxx_0[j] + 2.0 * fl1_fx * tpz_0_xxx_0[j];

            tpx_x_xxxy_0[j] = pa_x[j] * tpx_0_xxxy_0[j] + 1.5 * fl1_fx * tpx_0_xxy_0[j] - fl1_fgb * fl1_fx * ts_0_xxxy_0[j];

            tpy_x_xxxy_0[j] = pa_x[j] * tpy_0_xxxy_0[j] + 1.5 * fl1_fx * tpy_0_xxy_0[j];

            tpz_x_xxxy_0[j] = pa_x[j] * tpz_0_xxxy_0[j] + 1.5 * fl1_fx * tpz_0_xxy_0[j];

            tpx_x_xxxz_0[j] = pa_x[j] * tpx_0_xxxz_0[j] + 1.5 * fl1_fx * tpx_0_xxz_0[j] - fl1_fgb * fl1_fx * ts_0_xxxz_0[j];

            tpy_x_xxxz_0[j] = pa_x[j] * tpy_0_xxxz_0[j] + 1.5 * fl1_fx * tpy_0_xxz_0[j];

            tpz_x_xxxz_0[j] = pa_x[j] * tpz_0_xxxz_0[j] + 1.5 * fl1_fx * tpz_0_xxz_0[j];

            tpx_x_xxyy_0[j] = pa_x[j] * tpx_0_xxyy_0[j] + fl1_fx * tpx_0_xyy_0[j] - fl1_fgb * fl1_fx * ts_0_xxyy_0[j];

            tpy_x_xxyy_0[j] = pa_x[j] * tpy_0_xxyy_0[j] + fl1_fx * tpy_0_xyy_0[j];

            tpz_x_xxyy_0[j] = pa_x[j] * tpz_0_xxyy_0[j] + fl1_fx * tpz_0_xyy_0[j];

            tpx_x_xxyz_0[j] = pa_x[j] * tpx_0_xxyz_0[j] + fl1_fx * tpx_0_xyz_0[j] - fl1_fgb * fl1_fx * ts_0_xxyz_0[j];

            tpy_x_xxyz_0[j] = pa_x[j] * tpy_0_xxyz_0[j] + fl1_fx * tpy_0_xyz_0[j];

            tpz_x_xxyz_0[j] = pa_x[j] * tpz_0_xxyz_0[j] + fl1_fx * tpz_0_xyz_0[j];

            tpx_x_xxzz_0[j] = pa_x[j] * tpx_0_xxzz_0[j] + fl1_fx * tpx_0_xzz_0[j] - fl1_fgb * fl1_fx * ts_0_xxzz_0[j];

            tpy_x_xxzz_0[j] = pa_x[j] * tpy_0_xxzz_0[j] + fl1_fx * tpy_0_xzz_0[j];

            tpz_x_xxzz_0[j] = pa_x[j] * tpz_0_xxzz_0[j] + fl1_fx * tpz_0_xzz_0[j];

            tpx_x_xyyy_0[j] = pa_x[j] * tpx_0_xyyy_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_0_xyyy_0[j];

            tpy_x_xyyy_0[j] = pa_x[j] * tpy_0_xyyy_0[j] + 0.5 * fl1_fx * tpy_0_yyy_0[j];

            tpz_x_xyyy_0[j] = pa_x[j] * tpz_0_xyyy_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j];

            tpx_x_xyyz_0[j] = pa_x[j] * tpx_0_xyyz_0[j] + 0.5 * fl1_fx * tpx_0_yyz_0[j] - fl1_fgb * fl1_fx * ts_0_xyyz_0[j];

            tpy_x_xyyz_0[j] = pa_x[j] * tpy_0_xyyz_0[j] + 0.5 * fl1_fx * tpy_0_yyz_0[j];

            tpz_x_xyyz_0[j] = pa_x[j] * tpz_0_xyyz_0[j] + 0.5 * fl1_fx * tpz_0_yyz_0[j];

            tpx_x_xyzz_0[j] = pa_x[j] * tpx_0_xyzz_0[j] + 0.5 * fl1_fx * tpx_0_yzz_0[j] - fl1_fgb * fl1_fx * ts_0_xyzz_0[j];

            tpy_x_xyzz_0[j] = pa_x[j] * tpy_0_xyzz_0[j] + 0.5 * fl1_fx * tpy_0_yzz_0[j];

            tpz_x_xyzz_0[j] = pa_x[j] * tpz_0_xyzz_0[j] + 0.5 * fl1_fx * tpz_0_yzz_0[j];

            tpx_x_xzzz_0[j] = pa_x[j] * tpx_0_xzzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_0_xzzz_0[j];

            tpy_x_xzzz_0[j] = pa_x[j] * tpy_0_xzzz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j];

            tpz_x_xzzz_0[j] = pa_x[j] * tpz_0_xzzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j];

            tpx_x_yyyy_0[j] = pa_x[j] * tpx_0_yyyy_0[j] - fl1_fgb * fl1_fx * ts_0_yyyy_0[j];

            tpy_x_yyyy_0[j] = pa_x[j] * tpy_0_yyyy_0[j];

            tpz_x_yyyy_0[j] = pa_x[j] * tpz_0_yyyy_0[j];

            tpx_x_yyyz_0[j] = pa_x[j] * tpx_0_yyyz_0[j] - fl1_fgb * fl1_fx * ts_0_yyyz_0[j];

            tpy_x_yyyz_0[j] = pa_x[j] * tpy_0_yyyz_0[j];

            tpz_x_yyyz_0[j] = pa_x[j] * tpz_0_yyyz_0[j];

            tpx_x_yyzz_0[j] = pa_x[j] * tpx_0_yyzz_0[j] - fl1_fgb * fl1_fx * ts_0_yyzz_0[j];

            tpy_x_yyzz_0[j] = pa_x[j] * tpy_0_yyzz_0[j];

            tpz_x_yyzz_0[j] = pa_x[j] * tpz_0_yyzz_0[j];

            tpx_x_yzzz_0[j] = pa_x[j] * tpx_0_yzzz_0[j] - fl1_fgb * fl1_fx * ts_0_yzzz_0[j];

            tpy_x_yzzz_0[j] = pa_x[j] * tpy_0_yzzz_0[j];

            tpz_x_yzzz_0[j] = pa_x[j] * tpz_0_yzzz_0[j];

            tpx_x_zzzz_0[j] = pa_x[j] * tpx_0_zzzz_0[j] - fl1_fgb * fl1_fx * ts_0_zzzz_0[j];

            tpy_x_zzzz_0[j] = pa_x[j] * tpy_0_zzzz_0[j];

            tpz_x_zzzz_0[j] = pa_x[j] * tpz_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPG_45_90(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 5);

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

        auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx);

        auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1);

        auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2);

        auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3);

        auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4);

        auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5);

        auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6);

        auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7);

        auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8);

        auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9);

        auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10);

        auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11);

        auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12);

        auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13);

        auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14);

        // set up pointers to integrals

        auto tpx_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 15);

        auto tpy_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 15);

        auto tpz_y_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 15);

        auto tpx_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 16);

        auto tpy_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 16);

        auto tpz_y_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 16);

        auto tpx_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 17);

        auto tpy_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 17);

        auto tpz_y_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 17);

        auto tpx_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 18);

        auto tpy_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 18);

        auto tpz_y_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 18);

        auto tpx_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 19);

        auto tpy_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 19);

        auto tpz_y_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 19);

        auto tpx_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 20);

        auto tpy_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 20);

        auto tpz_y_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 20);

        auto tpx_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 21);

        auto tpy_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 21);

        auto tpz_y_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 21);

        auto tpx_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 22);

        auto tpy_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 22);

        auto tpz_y_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 22);

        auto tpx_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 23);

        auto tpy_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 23);

        auto tpz_y_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 23);

        auto tpx_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 24);

        auto tpy_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 24);

        auto tpz_y_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 24);

        auto tpx_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 25);

        auto tpy_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 25);

        auto tpz_y_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 25);

        auto tpx_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 26);

        auto tpy_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 26);

        auto tpz_y_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 26);

        auto tpx_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 27);

        auto tpy_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 27);

        auto tpz_y_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 27);

        auto tpx_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 28);

        auto tpy_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 28);

        auto tpz_y_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 28);

        auto tpx_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 29);

        auto tpy_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 29);

        auto tpz_y_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_y, tpx_0_xxx_0, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, \
                                     tpx_0_xxy_0, tpx_0_xxyy_0, tpx_0_xxyz_0, tpx_0_xxz_0, tpx_0_xxzz_0, tpx_0_xyy_0, \
                                     tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyz_0, tpx_0_xyzz_0, tpx_0_xzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyy_0, tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyz_0, tpx_0_yyzz_0, tpx_0_yzz_0, \
                                     tpx_0_yzzz_0, tpx_0_zzz_0, tpx_0_zzzz_0, tpx_y_xxxx_0, tpx_y_xxxy_0, tpx_y_xxxz_0, \
                                     tpx_y_xxyy_0, tpx_y_xxyz_0, tpx_y_xxzz_0, tpx_y_xyyy_0, tpx_y_xyyz_0, tpx_y_xyzz_0, \
                                     tpx_y_xzzz_0, tpx_y_yyyy_0, tpx_y_yyyz_0, tpx_y_yyzz_0, tpx_y_yzzz_0, tpx_y_zzzz_0, \
                                     tpy_0_xxx_0, tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxy_0, tpy_0_xxyy_0, \
                                     tpy_0_xxyz_0, tpy_0_xxz_0, tpy_0_xxzz_0, tpy_0_xyy_0, tpy_0_xyyy_0, tpy_0_xyyz_0, \
                                     tpy_0_xyz_0, tpy_0_xyzz_0, tpy_0_xzz_0, tpy_0_xzzz_0, tpy_0_yyy_0, tpy_0_yyyy_0, \
                                     tpy_0_yyyz_0, tpy_0_yyz_0, tpy_0_yyzz_0, tpy_0_yzz_0, tpy_0_yzzz_0, tpy_0_zzz_0, \
                                     tpy_0_zzzz_0, tpy_y_xxxx_0, tpy_y_xxxy_0, tpy_y_xxxz_0, tpy_y_xxyy_0, tpy_y_xxyz_0, \
                                     tpy_y_xxzz_0, tpy_y_xyyy_0, tpy_y_xyyz_0, tpy_y_xyzz_0, tpy_y_xzzz_0, tpy_y_yyyy_0, \
                                     tpy_y_yyyz_0, tpy_y_yyzz_0, tpy_y_yzzz_0, tpy_y_zzzz_0, tpz_0_xxx_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxy_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxz_0, \
                                     tpz_0_xxzz_0, tpz_0_xyy_0, tpz_0_xyyy_0, tpz_0_xyyz_0, tpz_0_xyz_0, tpz_0_xyzz_0, \
                                     tpz_0_xzz_0, tpz_0_xzzz_0, tpz_0_yyy_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyz_0, \
                                     tpz_0_yyzz_0, tpz_0_yzz_0, tpz_0_yzzz_0, tpz_0_zzz_0, tpz_0_zzzz_0, tpz_y_xxxx_0, \
                                     tpz_y_xxxy_0, tpz_y_xxxz_0, tpz_y_xxyy_0, tpz_y_xxyz_0, tpz_y_xxzz_0, tpz_y_xyyy_0, \
                                     tpz_y_xyyz_0, tpz_y_xyzz_0, tpz_y_xzzz_0, tpz_y_yyyy_0, tpz_y_yyyz_0, tpz_y_yyzz_0, \
                                     tpz_y_yzzz_0, tpz_y_zzzz_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_y_xxxx_0[j] = pa_y[j] * tpx_0_xxxx_0[j];

            tpy_y_xxxx_0[j] = pa_y[j] * tpy_0_xxxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxxx_0[j];

            tpz_y_xxxx_0[j] = pa_y[j] * tpz_0_xxxx_0[j];

            tpx_y_xxxy_0[j] = pa_y[j] * tpx_0_xxxy_0[j] + 0.5 * fl1_fx * tpx_0_xxx_0[j];

            tpy_y_xxxy_0[j] = pa_y[j] * tpy_0_xxxy_0[j] + 0.5 * fl1_fx * tpy_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxxy_0[j];

            tpz_y_xxxy_0[j] = pa_y[j] * tpz_0_xxxy_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j];

            tpx_y_xxxz_0[j] = pa_y[j] * tpx_0_xxxz_0[j];

            tpy_y_xxxz_0[j] = pa_y[j] * tpy_0_xxxz_0[j] - fl1_fgb * fl1_fx * ts_0_xxxz_0[j];

            tpz_y_xxxz_0[j] = pa_y[j] * tpz_0_xxxz_0[j];

            tpx_y_xxyy_0[j] = pa_y[j] * tpx_0_xxyy_0[j] + fl1_fx * tpx_0_xxy_0[j];

            tpy_y_xxyy_0[j] = pa_y[j] * tpy_0_xxyy_0[j] + fl1_fx * tpy_0_xxy_0[j] - fl1_fgb * fl1_fx * ts_0_xxyy_0[j];

            tpz_y_xxyy_0[j] = pa_y[j] * tpz_0_xxyy_0[j] + fl1_fx * tpz_0_xxy_0[j];

            tpx_y_xxyz_0[j] = pa_y[j] * tpx_0_xxyz_0[j] + 0.5 * fl1_fx * tpx_0_xxz_0[j];

            tpy_y_xxyz_0[j] = pa_y[j] * tpy_0_xxyz_0[j] + 0.5 * fl1_fx * tpy_0_xxz_0[j] - fl1_fgb * fl1_fx * ts_0_xxyz_0[j];

            tpz_y_xxyz_0[j] = pa_y[j] * tpz_0_xxyz_0[j] + 0.5 * fl1_fx * tpz_0_xxz_0[j];

            tpx_y_xxzz_0[j] = pa_y[j] * tpx_0_xxzz_0[j];

            tpy_y_xxzz_0[j] = pa_y[j] * tpy_0_xxzz_0[j] - fl1_fgb * fl1_fx * ts_0_xxzz_0[j];

            tpz_y_xxzz_0[j] = pa_y[j] * tpz_0_xxzz_0[j];

            tpx_y_xyyy_0[j] = pa_y[j] * tpx_0_xyyy_0[j] + 1.5 * fl1_fx * tpx_0_xyy_0[j];

            tpy_y_xyyy_0[j] = pa_y[j] * tpy_0_xyyy_0[j] + 1.5 * fl1_fx * tpy_0_xyy_0[j] - fl1_fgb * fl1_fx * ts_0_xyyy_0[j];

            tpz_y_xyyy_0[j] = pa_y[j] * tpz_0_xyyy_0[j] + 1.5 * fl1_fx * tpz_0_xyy_0[j];

            tpx_y_xyyz_0[j] = pa_y[j] * tpx_0_xyyz_0[j] + fl1_fx * tpx_0_xyz_0[j];

            tpy_y_xyyz_0[j] = pa_y[j] * tpy_0_xyyz_0[j] + fl1_fx * tpy_0_xyz_0[j] - fl1_fgb * fl1_fx * ts_0_xyyz_0[j];

            tpz_y_xyyz_0[j] = pa_y[j] * tpz_0_xyyz_0[j] + fl1_fx * tpz_0_xyz_0[j];

            tpx_y_xyzz_0[j] = pa_y[j] * tpx_0_xyzz_0[j] + 0.5 * fl1_fx * tpx_0_xzz_0[j];

            tpy_y_xyzz_0[j] = pa_y[j] * tpy_0_xyzz_0[j] + 0.5 * fl1_fx * tpy_0_xzz_0[j] - fl1_fgb * fl1_fx * ts_0_xyzz_0[j];

            tpz_y_xyzz_0[j] = pa_y[j] * tpz_0_xyzz_0[j] + 0.5 * fl1_fx * tpz_0_xzz_0[j];

            tpx_y_xzzz_0[j] = pa_y[j] * tpx_0_xzzz_0[j];

            tpy_y_xzzz_0[j] = pa_y[j] * tpy_0_xzzz_0[j] - fl1_fgb * fl1_fx * ts_0_xzzz_0[j];

            tpz_y_xzzz_0[j] = pa_y[j] * tpz_0_xzzz_0[j];

            tpx_y_yyyy_0[j] = pa_y[j] * tpx_0_yyyy_0[j] + 2.0 * fl1_fx * tpx_0_yyy_0[j];

            tpy_y_yyyy_0[j] = pa_y[j] * tpy_0_yyyy_0[j] + 2.0 * fl1_fx * tpy_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_0_yyyy_0[j];

            tpz_y_yyyy_0[j] = pa_y[j] * tpz_0_yyyy_0[j] + 2.0 * fl1_fx * tpz_0_yyy_0[j];

            tpx_y_yyyz_0[j] = pa_y[j] * tpx_0_yyyz_0[j] + 1.5 * fl1_fx * tpx_0_yyz_0[j];

            tpy_y_yyyz_0[j] = pa_y[j] * tpy_0_yyyz_0[j] + 1.5 * fl1_fx * tpy_0_yyz_0[j] - fl1_fgb * fl1_fx * ts_0_yyyz_0[j];

            tpz_y_yyyz_0[j] = pa_y[j] * tpz_0_yyyz_0[j] + 1.5 * fl1_fx * tpz_0_yyz_0[j];

            tpx_y_yyzz_0[j] = pa_y[j] * tpx_0_yyzz_0[j] + fl1_fx * tpx_0_yzz_0[j];

            tpy_y_yyzz_0[j] = pa_y[j] * tpy_0_yyzz_0[j] + fl1_fx * tpy_0_yzz_0[j] - fl1_fgb * fl1_fx * ts_0_yyzz_0[j];

            tpz_y_yyzz_0[j] = pa_y[j] * tpz_0_yyzz_0[j] + fl1_fx * tpz_0_yzz_0[j];

            tpx_y_yzzz_0[j] = pa_y[j] * tpx_0_yzzz_0[j] + 0.5 * fl1_fx * tpx_0_zzz_0[j];

            tpy_y_yzzz_0[j] = pa_y[j] * tpy_0_yzzz_0[j] + 0.5 * fl1_fx * tpy_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_0_yzzz_0[j];

            tpz_y_yzzz_0[j] = pa_y[j] * tpz_0_yzzz_0[j] + 0.5 * fl1_fx * tpz_0_zzz_0[j];

            tpx_y_zzzz_0[j] = pa_y[j] * tpx_0_zzzz_0[j];

            tpy_y_zzzz_0[j] = pa_y[j] * tpy_0_zzzz_0[j] - fl1_fgb * fl1_fx * ts_0_zzzz_0[j];

            tpz_y_zzzz_0[j] = pa_y[j] * tpz_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForPG_90_135(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_1_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // loop over contracted GTO on bra side

    int32_t idx = 0;

    for (int32_t i = spos[iContrGto]; i < epos[iContrGto]; i++)
    {
        // set up pointers to Obara-Saika factors

        auto fx = osFactors.data(4 * idx);

        auto fgb = osFactors.data(4 * idx + 3);

        // set up pointers to tensors product of distances R(PA) = P - A

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tpx_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx);

        auto tpy_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx);

        auto tpz_0_xxxx_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx);

        auto tpx_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 1);

        auto tpy_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 1);

        auto tpz_0_xxxy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 1);

        auto tpx_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 2);

        auto tpy_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 2);

        auto tpz_0_xxxz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 2);

        auto tpx_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 3);

        auto tpy_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 3);

        auto tpz_0_xxyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 3);

        auto tpx_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 4);

        auto tpy_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 4);

        auto tpz_0_xxyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 4);

        auto tpx_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 5);

        auto tpy_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 5);

        auto tpz_0_xxzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 5);

        auto tpx_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 6);

        auto tpy_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 6);

        auto tpz_0_xyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 6);

        auto tpx_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 7);

        auto tpy_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 7);

        auto tpz_0_xyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 7);

        auto tpx_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 8);

        auto tpy_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 8);

        auto tpz_0_xyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 8);

        auto tpx_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 9);

        auto tpy_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 9);

        auto tpz_0_xzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 9);

        auto tpx_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 10);

        auto tpy_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 10);

        auto tpz_0_yyyy_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 10);

        auto tpx_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 11);

        auto tpy_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 11);

        auto tpz_0_yyyz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 11);

        auto tpx_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 12);

        auto tpy_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 12);

        auto tpz_0_yyzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 12);

        auto tpx_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 13);

        auto tpy_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 13);

        auto tpz_0_yzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 13);

        auto tpx_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * idx + 14);

        auto tpy_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 15 * bdim + 15 * idx + 14);

        auto tpz_0_zzzz_0 = primBuffer.data(pidx_p_0_4_m0 + 30 * bdim + 15 * idx + 14);

        auto tpx_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx);

        auto tpy_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx);

        auto tpz_0_xxx_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx);

        auto tpx_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 1);

        auto tpy_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_0_xxy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 2);

        auto tpy_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_0_xxz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 3);

        auto tpy_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_0_xyy_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 4);

        auto tpy_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_0_xyz_0 = primBuffer.data(pidx_p_0_3_m0 + 20 * bdim + 10 * idx + 4);

        auto tpx_0_xzz_0 = primBuffer.data(pidx_p_0_3_m0 + 10 * idx + 5);

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

        auto ts_0_xxxx_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx);

        auto ts_0_xxxy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 1);

        auto ts_0_xxxz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 2);

        auto ts_0_xxyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 3);

        auto ts_0_xxyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 4);

        auto ts_0_xxzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 5);

        auto ts_0_xyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 6);

        auto ts_0_xyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 7);

        auto ts_0_xyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 8);

        auto ts_0_xzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 9);

        auto ts_0_yyyy_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 10);

        auto ts_0_yyyz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 11);

        auto ts_0_yyzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 12);

        auto ts_0_yzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 13);

        auto ts_0_zzzz_0 = primBuffer.data(pidx_s_0_4_m0 + 15 * idx + 14);

        // set up pointers to integrals

        auto tpx_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 30);

        auto tpy_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_z_xxxx_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 31);

        auto tpy_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_z_xxxy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 32);

        auto tpy_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_z_xxxz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 33);

        auto tpy_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_z_xxyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 34);

        auto tpy_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_z_xxyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 35);

        auto tpy_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_z_xxzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 36);

        auto tpy_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_z_xyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 37);

        auto tpy_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_z_xyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 38);

        auto tpy_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_z_xyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 39);

        auto tpy_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_z_xzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 40);

        auto tpy_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_z_yyyy_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 41);

        auto tpy_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_z_yyyz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 42);

        auto tpy_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_z_yyzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 43);

        auto tpy_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_z_yzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * idx + 44);

        auto tpy_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_z_zzzz_0 = primBuffer.data(pidx_p_1_4_m0 + 90 * bdim + 45 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_z, tpx_0_xxx_0, tpx_0_xxxx_0, tpx_0_xxxy_0, tpx_0_xxxz_0, \
                                     tpx_0_xxy_0, tpx_0_xxyy_0, tpx_0_xxyz_0, tpx_0_xxz_0, tpx_0_xxzz_0, tpx_0_xyy_0, \
                                     tpx_0_xyyy_0, tpx_0_xyyz_0, tpx_0_xyz_0, tpx_0_xyzz_0, tpx_0_xzz_0, tpx_0_xzzz_0, \
                                     tpx_0_yyy_0, tpx_0_yyyy_0, tpx_0_yyyz_0, tpx_0_yyz_0, tpx_0_yyzz_0, tpx_0_yzz_0, \
                                     tpx_0_yzzz_0, tpx_0_zzz_0, tpx_0_zzzz_0, tpx_z_xxxx_0, tpx_z_xxxy_0, tpx_z_xxxz_0, \
                                     tpx_z_xxyy_0, tpx_z_xxyz_0, tpx_z_xxzz_0, tpx_z_xyyy_0, tpx_z_xyyz_0, tpx_z_xyzz_0, \
                                     tpx_z_xzzz_0, tpx_z_yyyy_0, tpx_z_yyyz_0, tpx_z_yyzz_0, tpx_z_yzzz_0, tpx_z_zzzz_0, \
                                     tpy_0_xxx_0, tpy_0_xxxx_0, tpy_0_xxxy_0, tpy_0_xxxz_0, tpy_0_xxy_0, tpy_0_xxyy_0, \
                                     tpy_0_xxyz_0, tpy_0_xxz_0, tpy_0_xxzz_0, tpy_0_xyy_0, tpy_0_xyyy_0, tpy_0_xyyz_0, \
                                     tpy_0_xyz_0, tpy_0_xyzz_0, tpy_0_xzz_0, tpy_0_xzzz_0, tpy_0_yyy_0, tpy_0_yyyy_0, \
                                     tpy_0_yyyz_0, tpy_0_yyz_0, tpy_0_yyzz_0, tpy_0_yzz_0, tpy_0_yzzz_0, tpy_0_zzz_0, \
                                     tpy_0_zzzz_0, tpy_z_xxxx_0, tpy_z_xxxy_0, tpy_z_xxxz_0, tpy_z_xxyy_0, tpy_z_xxyz_0, \
                                     tpy_z_xxzz_0, tpy_z_xyyy_0, tpy_z_xyyz_0, tpy_z_xyzz_0, tpy_z_xzzz_0, tpy_z_yyyy_0, \
                                     tpy_z_yyyz_0, tpy_z_yyzz_0, tpy_z_yzzz_0, tpy_z_zzzz_0, tpz_0_xxx_0, tpz_0_xxxx_0, \
                                     tpz_0_xxxy_0, tpz_0_xxxz_0, tpz_0_xxy_0, tpz_0_xxyy_0, tpz_0_xxyz_0, tpz_0_xxz_0, \
                                     tpz_0_xxzz_0, tpz_0_xyy_0, tpz_0_xyyy_0, tpz_0_xyyz_0, tpz_0_xyz_0, tpz_0_xyzz_0, \
                                     tpz_0_xzz_0, tpz_0_xzzz_0, tpz_0_yyy_0, tpz_0_yyyy_0, tpz_0_yyyz_0, tpz_0_yyz_0, \
                                     tpz_0_yyzz_0, tpz_0_yzz_0, tpz_0_yzzz_0, tpz_0_zzz_0, tpz_0_zzzz_0, tpz_z_xxxx_0, \
                                     tpz_z_xxxy_0, tpz_z_xxxz_0, tpz_z_xxyy_0, tpz_z_xxyz_0, tpz_z_xxzz_0, tpz_z_xyyy_0, \
                                     tpz_z_xyyz_0, tpz_z_xyzz_0, tpz_z_xzzz_0, tpz_z_yyyy_0, tpz_z_yyyz_0, tpz_z_yyzz_0, \
                                     tpz_z_yzzz_0, tpz_z_zzzz_0, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_z_xxxx_0[j] = pa_z[j] * tpx_0_xxxx_0[j];

            tpy_z_xxxx_0[j] = pa_z[j] * tpy_0_xxxx_0[j];

            tpz_z_xxxx_0[j] = pa_z[j] * tpz_0_xxxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxxx_0[j];

            tpx_z_xxxy_0[j] = pa_z[j] * tpx_0_xxxy_0[j];

            tpy_z_xxxy_0[j] = pa_z[j] * tpy_0_xxxy_0[j];

            tpz_z_xxxy_0[j] = pa_z[j] * tpz_0_xxxy_0[j] - fl1_fgb * fl1_fx * ts_0_xxxy_0[j];

            tpx_z_xxxz_0[j] = pa_z[j] * tpx_0_xxxz_0[j] + 0.5 * fl1_fx * tpx_0_xxx_0[j];

            tpy_z_xxxz_0[j] = pa_z[j] * tpy_0_xxxz_0[j] + 0.5 * fl1_fx * tpy_0_xxx_0[j];

            tpz_z_xxxz_0[j] = pa_z[j] * tpz_0_xxxz_0[j] + 0.5 * fl1_fx * tpz_0_xxx_0[j] - fl1_fgb * fl1_fx * ts_0_xxxz_0[j];

            tpx_z_xxyy_0[j] = pa_z[j] * tpx_0_xxyy_0[j];

            tpy_z_xxyy_0[j] = pa_z[j] * tpy_0_xxyy_0[j];

            tpz_z_xxyy_0[j] = pa_z[j] * tpz_0_xxyy_0[j] - fl1_fgb * fl1_fx * ts_0_xxyy_0[j];

            tpx_z_xxyz_0[j] = pa_z[j] * tpx_0_xxyz_0[j] + 0.5 * fl1_fx * tpx_0_xxy_0[j];

            tpy_z_xxyz_0[j] = pa_z[j] * tpy_0_xxyz_0[j] + 0.5 * fl1_fx * tpy_0_xxy_0[j];

            tpz_z_xxyz_0[j] = pa_z[j] * tpz_0_xxyz_0[j] + 0.5 * fl1_fx * tpz_0_xxy_0[j] - fl1_fgb * fl1_fx * ts_0_xxyz_0[j];

            tpx_z_xxzz_0[j] = pa_z[j] * tpx_0_xxzz_0[j] + fl1_fx * tpx_0_xxz_0[j];

            tpy_z_xxzz_0[j] = pa_z[j] * tpy_0_xxzz_0[j] + fl1_fx * tpy_0_xxz_0[j];

            tpz_z_xxzz_0[j] = pa_z[j] * tpz_0_xxzz_0[j] + fl1_fx * tpz_0_xxz_0[j] - fl1_fgb * fl1_fx * ts_0_xxzz_0[j];

            tpx_z_xyyy_0[j] = pa_z[j] * tpx_0_xyyy_0[j];

            tpy_z_xyyy_0[j] = pa_z[j] * tpy_0_xyyy_0[j];

            tpz_z_xyyy_0[j] = pa_z[j] * tpz_0_xyyy_0[j] - fl1_fgb * fl1_fx * ts_0_xyyy_0[j];

            tpx_z_xyyz_0[j] = pa_z[j] * tpx_0_xyyz_0[j] + 0.5 * fl1_fx * tpx_0_xyy_0[j];

            tpy_z_xyyz_0[j] = pa_z[j] * tpy_0_xyyz_0[j] + 0.5 * fl1_fx * tpy_0_xyy_0[j];

            tpz_z_xyyz_0[j] = pa_z[j] * tpz_0_xyyz_0[j] + 0.5 * fl1_fx * tpz_0_xyy_0[j] - fl1_fgb * fl1_fx * ts_0_xyyz_0[j];

            tpx_z_xyzz_0[j] = pa_z[j] * tpx_0_xyzz_0[j] + fl1_fx * tpx_0_xyz_0[j];

            tpy_z_xyzz_0[j] = pa_z[j] * tpy_0_xyzz_0[j] + fl1_fx * tpy_0_xyz_0[j];

            tpz_z_xyzz_0[j] = pa_z[j] * tpz_0_xyzz_0[j] + fl1_fx * tpz_0_xyz_0[j] - fl1_fgb * fl1_fx * ts_0_xyzz_0[j];

            tpx_z_xzzz_0[j] = pa_z[j] * tpx_0_xzzz_0[j] + 1.5 * fl1_fx * tpx_0_xzz_0[j];

            tpy_z_xzzz_0[j] = pa_z[j] * tpy_0_xzzz_0[j] + 1.5 * fl1_fx * tpy_0_xzz_0[j];

            tpz_z_xzzz_0[j] = pa_z[j] * tpz_0_xzzz_0[j] + 1.5 * fl1_fx * tpz_0_xzz_0[j] - fl1_fgb * fl1_fx * ts_0_xzzz_0[j];

            tpx_z_yyyy_0[j] = pa_z[j] * tpx_0_yyyy_0[j];

            tpy_z_yyyy_0[j] = pa_z[j] * tpy_0_yyyy_0[j];

            tpz_z_yyyy_0[j] = pa_z[j] * tpz_0_yyyy_0[j] - fl1_fgb * fl1_fx * ts_0_yyyy_0[j];

            tpx_z_yyyz_0[j] = pa_z[j] * tpx_0_yyyz_0[j] + 0.5 * fl1_fx * tpx_0_yyy_0[j];

            tpy_z_yyyz_0[j] = pa_z[j] * tpy_0_yyyz_0[j] + 0.5 * fl1_fx * tpy_0_yyy_0[j];

            tpz_z_yyyz_0[j] = pa_z[j] * tpz_0_yyyz_0[j] + 0.5 * fl1_fx * tpz_0_yyy_0[j] - fl1_fgb * fl1_fx * ts_0_yyyz_0[j];

            tpx_z_yyzz_0[j] = pa_z[j] * tpx_0_yyzz_0[j] + fl1_fx * tpx_0_yyz_0[j];

            tpy_z_yyzz_0[j] = pa_z[j] * tpy_0_yyzz_0[j] + fl1_fx * tpy_0_yyz_0[j];

            tpz_z_yyzz_0[j] = pa_z[j] * tpz_0_yyzz_0[j] + fl1_fx * tpz_0_yyz_0[j] - fl1_fgb * fl1_fx * ts_0_yyzz_0[j];

            tpx_z_yzzz_0[j] = pa_z[j] * tpx_0_yzzz_0[j] + 1.5 * fl1_fx * tpx_0_yzz_0[j];

            tpy_z_yzzz_0[j] = pa_z[j] * tpy_0_yzzz_0[j] + 1.5 * fl1_fx * tpy_0_yzz_0[j];

            tpz_z_yzzz_0[j] = pa_z[j] * tpz_0_yzzz_0[j] + 1.5 * fl1_fx * tpz_0_yzz_0[j] - fl1_fgb * fl1_fx * ts_0_yzzz_0[j];

            tpx_z_zzzz_0[j] = pa_z[j] * tpx_0_zzzz_0[j] + 2.0 * fl1_fx * tpx_0_zzz_0[j];

            tpy_z_zzzz_0[j] = pa_z[j] * tpy_0_zzzz_0[j] + 2.0 * fl1_fx * tpy_0_zzz_0[j];

            tpz_z_zzzz_0[j] = pa_z[j] * tpz_0_zzzz_0[j] + 2.0 * fl1_fx * tpz_0_zzz_0[j] - fl1_fgb * fl1_fx * ts_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGP(CMemBlock2D<double>&       primBuffer,
                        const CRecursionMap&       recursionMap,
                        const CMemBlock2D<double>& osFactors,
                        const CMemBlock2D<double>& paDistances,
                        const CGtoBlock&           braGtoBlock,
                        const CGtoBlock&           ketGtoBlock,
                        const int32_t              iContrGto)
{
    lmomrecfunc::compLinearMomentumForGP_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGP_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    lmomrecfunc::compLinearMomentumForGP_90_135(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compLinearMomentumForGP_0_45(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx);

        auto tpy_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx);

        auto tpz_xxx_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx);

        auto tpx_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 1);

        auto tpy_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 1);

        auto tpz_xxx_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 1);

        auto tpx_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 2);

        auto tpy_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 2);

        auto tpz_xxx_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 2);

        auto tpx_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 3);

        auto tpy_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 3);

        auto tpz_xxy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 3);

        auto tpx_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 4);

        auto tpy_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 4);

        auto tpz_xxy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 4);

        auto tpx_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 5);

        auto tpy_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 5);

        auto tpz_xxy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 5);

        auto tpx_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 6);

        auto tpy_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 6);

        auto tpz_xxz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 6);

        auto tpx_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 7);

        auto tpy_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 7);

        auto tpz_xxz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 7);

        auto tpx_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 8);

        auto tpy_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 8);

        auto tpz_xxz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 8);

        auto tpx_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 9);

        auto tpy_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 9);

        auto tpz_xyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 9);

        auto tpx_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 10);

        auto tpy_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 10);

        auto tpz_xyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 10);

        auto tpx_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 11);

        auto tpy_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 11);

        auto tpz_xyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 11);

        auto tpx_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 12);

        auto tpy_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 12);

        auto tpz_xyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 12);

        auto tpx_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 13);

        auto tpy_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 13);

        auto tpz_xyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 13);

        auto tpx_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 14);

        auto tpy_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 14);

        auto tpz_xyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 14);

        auto tpx_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx);

        auto tpy_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx);

        auto tpz_xx_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx);

        auto tpx_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 1);

        auto tpy_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 1);

        auto tpz_xx_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 1);

        auto tpx_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 2);

        auto tpy_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 2);

        auto tpz_xx_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 2);

        auto tpx_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 3);

        auto tpy_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 3);

        auto tpz_xy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 3);

        auto tpx_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 4);

        auto tpy_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 4);

        auto tpz_xy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 4);

        auto tpx_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 5);

        auto tpy_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 5);

        auto tpz_xy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 5);

        auto tpx_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 6);

        auto tpy_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 6);

        auto tpz_xz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 6);

        auto tpx_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 7);

        auto tpy_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 7);

        auto tpz_xz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 7);

        auto tpx_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 8);

        auto tpy_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 8);

        auto tpz_xz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 8);

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx);

        auto tpy_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx);

        auto tpz_xxx_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx);

        auto tpx_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 1);

        auto tpy_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 1);

        auto tpz_xxy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 1);

        auto tpx_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 2);

        auto tpy_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 2);

        auto tpz_xxz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 2);

        auto tpx_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 3);

        auto tpy_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 3);

        auto tpz_xyy_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 3);

        auto tpx_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 4);

        auto tpy_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * bdim + 10 * idx + 4);

        auto tpz_xyz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 20 * bdim + 10 * idx + 4);

        auto ts_xxx_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx);

        auto ts_xxx_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 1);

        auto ts_xxx_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 2);

        auto ts_xxy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 3);

        auto ts_xxy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 4);

        auto ts_xxy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 5);

        auto ts_xxz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 6);

        auto ts_xxz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 7);

        auto ts_xxz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 8);

        auto ts_xyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 9);

        auto ts_xyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 10);

        auto ts_xyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 11);

        auto ts_xyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 12);

        auto ts_xyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 13);

        auto ts_xyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 14);

        // set up pointers to integrals

        auto tpx_xxxx_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx);

        auto tpy_xxxx_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx);

        auto tpz_xxxx_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx);

        auto tpx_xxxx_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 1);

        auto tpy_xxxx_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 1);

        auto tpz_xxxx_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 1);

        auto tpx_xxxx_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 2);

        auto tpy_xxxx_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 2);

        auto tpz_xxxx_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 2);

        auto tpx_xxxy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 3);

        auto tpy_xxxy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 3);

        auto tpz_xxxy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 3);

        auto tpx_xxxy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 4);

        auto tpy_xxxy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 4);

        auto tpz_xxxy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 4);

        auto tpx_xxxy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 5);

        auto tpy_xxxy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 5);

        auto tpz_xxxy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 5);

        auto tpx_xxxz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 6);

        auto tpy_xxxz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 6);

        auto tpz_xxxz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 6);

        auto tpx_xxxz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 7);

        auto tpy_xxxz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 7);

        auto tpz_xxxz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 7);

        auto tpx_xxxz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 8);

        auto tpy_xxxz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 8);

        auto tpz_xxxz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 8);

        auto tpx_xxyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 9);

        auto tpy_xxyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 9);

        auto tpz_xxyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 9);

        auto tpx_xxyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 10);

        auto tpy_xxyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 10);

        auto tpz_xxyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 10);

        auto tpx_xxyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 11);

        auto tpy_xxyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 11);

        auto tpz_xxyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 11);

        auto tpx_xxyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 12);

        auto tpy_xxyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 12);

        auto tpz_xxyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 12);

        auto tpx_xxyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 13);

        auto tpy_xxyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 13);

        auto tpz_xxyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 13);

        auto tpx_xxyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 14);

        auto tpy_xxyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 14);

        auto tpz_xxyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 14);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xx_x_0, tpx_xx_y_0, tpx_xx_z_0, tpx_xxx_0_0, \
                                     tpx_xxx_x_0, tpx_xxx_y_0, tpx_xxx_z_0, tpx_xxxx_x_0, tpx_xxxx_y_0, tpx_xxxx_z_0, \
                                     tpx_xxxy_x_0, tpx_xxxy_y_0, tpx_xxxy_z_0, tpx_xxxz_x_0, tpx_xxxz_y_0, tpx_xxxz_z_0, \
                                     tpx_xxy_0_0, tpx_xxy_x_0, tpx_xxy_y_0, tpx_xxy_z_0, tpx_xxyy_x_0, tpx_xxyy_y_0, \
                                     tpx_xxyy_z_0, tpx_xxyz_x_0, tpx_xxyz_y_0, tpx_xxyz_z_0, tpx_xxz_0_0, tpx_xxz_x_0, \
                                     tpx_xxz_y_0, tpx_xxz_z_0, tpx_xy_x_0, tpx_xy_y_0, tpx_xy_z_0, tpx_xyy_0_0, \
                                     tpx_xyy_x_0, tpx_xyy_y_0, tpx_xyy_z_0, tpx_xyz_0_0, tpx_xyz_x_0, tpx_xyz_y_0, \
                                     tpx_xyz_z_0, tpx_xz_x_0, tpx_xz_y_0, tpx_xz_z_0, tpx_yy_x_0, tpx_yy_y_0, tpx_yy_z_0, \
                                     tpx_yz_x_0, tpx_yz_y_0, tpx_yz_z_0, tpy_xx_x_0, tpy_xx_y_0, tpy_xx_z_0, \
                                     tpy_xxx_0_0, tpy_xxx_x_0, tpy_xxx_y_0, tpy_xxx_z_0, tpy_xxxx_x_0, tpy_xxxx_y_0, \
                                     tpy_xxxx_z_0, tpy_xxxy_x_0, tpy_xxxy_y_0, tpy_xxxy_z_0, tpy_xxxz_x_0, tpy_xxxz_y_0, \
                                     tpy_xxxz_z_0, tpy_xxy_0_0, tpy_xxy_x_0, tpy_xxy_y_0, tpy_xxy_z_0, tpy_xxyy_x_0, \
                                     tpy_xxyy_y_0, tpy_xxyy_z_0, tpy_xxyz_x_0, tpy_xxyz_y_0, tpy_xxyz_z_0, tpy_xxz_0_0, \
                                     tpy_xxz_x_0, tpy_xxz_y_0, tpy_xxz_z_0, tpy_xy_x_0, tpy_xy_y_0, tpy_xy_z_0, \
                                     tpy_xyy_0_0, tpy_xyy_x_0, tpy_xyy_y_0, tpy_xyy_z_0, tpy_xyz_0_0, tpy_xyz_x_0, \
                                     tpy_xyz_y_0, tpy_xyz_z_0, tpy_xz_x_0, tpy_xz_y_0, tpy_xz_z_0, tpy_yy_x_0, \
                                     tpy_yy_y_0, tpy_yy_z_0, tpy_yz_x_0, tpy_yz_y_0, tpy_yz_z_0, tpz_xx_x_0, tpz_xx_y_0, \
                                     tpz_xx_z_0, tpz_xxx_0_0, tpz_xxx_x_0, tpz_xxx_y_0, tpz_xxx_z_0, tpz_xxxx_x_0, \
                                     tpz_xxxx_y_0, tpz_xxxx_z_0, tpz_xxxy_x_0, tpz_xxxy_y_0, tpz_xxxy_z_0, tpz_xxxz_x_0, \
                                     tpz_xxxz_y_0, tpz_xxxz_z_0, tpz_xxy_0_0, tpz_xxy_x_0, tpz_xxy_y_0, tpz_xxy_z_0, \
                                     tpz_xxyy_x_0, tpz_xxyy_y_0, tpz_xxyy_z_0, tpz_xxyz_x_0, tpz_xxyz_y_0, tpz_xxyz_z_0, \
                                     tpz_xxz_0_0, tpz_xxz_x_0, tpz_xxz_y_0, tpz_xxz_z_0, tpz_xy_x_0, tpz_xy_y_0, \
                                     tpz_xy_z_0, tpz_xyy_0_0, tpz_xyy_x_0, tpz_xyy_y_0, tpz_xyy_z_0, tpz_xyz_0_0, \
                                     tpz_xyz_x_0, tpz_xyz_y_0, tpz_xyz_z_0, tpz_xz_x_0, tpz_xz_y_0, tpz_xz_z_0, \
                                     tpz_yy_x_0, tpz_yy_y_0, tpz_yy_z_0, tpz_yz_x_0, tpz_yz_y_0, tpz_yz_z_0, ts_xxx_x_0, \
                                     ts_xxx_y_0, ts_xxx_z_0, ts_xxy_x_0, ts_xxy_y_0, ts_xxy_z_0, ts_xxz_x_0, ts_xxz_y_0, \
                                     ts_xxz_z_0, ts_xyy_x_0, ts_xyy_y_0, ts_xyy_z_0, ts_xyz_x_0, ts_xyz_y_0, ts_xyz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxxx_x_0[j] =
                pa_x[j] * tpx_xxx_x_0[j] + 1.5 * fl1_fx * tpx_xx_x_0[j] + 0.5 * fl1_fx * tpx_xxx_0_0[j] - fl1_fgb * fl1_fx * ts_xxx_x_0[j];

            tpy_xxxx_x_0[j] = pa_x[j] * tpy_xxx_x_0[j] + 1.5 * fl1_fx * tpy_xx_x_0[j] + 0.5 * fl1_fx * tpy_xxx_0_0[j];

            tpz_xxxx_x_0[j] = pa_x[j] * tpz_xxx_x_0[j] + 1.5 * fl1_fx * tpz_xx_x_0[j] + 0.5 * fl1_fx * tpz_xxx_0_0[j];

            tpx_xxxx_y_0[j] = pa_x[j] * tpx_xxx_y_0[j] + 1.5 * fl1_fx * tpx_xx_y_0[j] - fl1_fgb * fl1_fx * ts_xxx_y_0[j];

            tpy_xxxx_y_0[j] = pa_x[j] * tpy_xxx_y_0[j] + 1.5 * fl1_fx * tpy_xx_y_0[j];

            tpz_xxxx_y_0[j] = pa_x[j] * tpz_xxx_y_0[j] + 1.5 * fl1_fx * tpz_xx_y_0[j];

            tpx_xxxx_z_0[j] = pa_x[j] * tpx_xxx_z_0[j] + 1.5 * fl1_fx * tpx_xx_z_0[j] - fl1_fgb * fl1_fx * ts_xxx_z_0[j];

            tpy_xxxx_z_0[j] = pa_x[j] * tpy_xxx_z_0[j] + 1.5 * fl1_fx * tpy_xx_z_0[j];

            tpz_xxxx_z_0[j] = pa_x[j] * tpz_xxx_z_0[j] + 1.5 * fl1_fx * tpz_xx_z_0[j];

            tpx_xxxy_x_0[j] = pa_x[j] * tpx_xxy_x_0[j] + fl1_fx * tpx_xy_x_0[j] + 0.5 * fl1_fx * tpx_xxy_0_0[j] - fl1_fgb * fl1_fx * ts_xxy_x_0[j];

            tpy_xxxy_x_0[j] = pa_x[j] * tpy_xxy_x_0[j] + fl1_fx * tpy_xy_x_0[j] + 0.5 * fl1_fx * tpy_xxy_0_0[j];

            tpz_xxxy_x_0[j] = pa_x[j] * tpz_xxy_x_0[j] + fl1_fx * tpz_xy_x_0[j] + 0.5 * fl1_fx * tpz_xxy_0_0[j];

            tpx_xxxy_y_0[j] = pa_x[j] * tpx_xxy_y_0[j] + fl1_fx * tpx_xy_y_0[j] - fl1_fgb * fl1_fx * ts_xxy_y_0[j];

            tpy_xxxy_y_0[j] = pa_x[j] * tpy_xxy_y_0[j] + fl1_fx * tpy_xy_y_0[j];

            tpz_xxxy_y_0[j] = pa_x[j] * tpz_xxy_y_0[j] + fl1_fx * tpz_xy_y_0[j];

            tpx_xxxy_z_0[j] = pa_x[j] * tpx_xxy_z_0[j] + fl1_fx * tpx_xy_z_0[j] - fl1_fgb * fl1_fx * ts_xxy_z_0[j];

            tpy_xxxy_z_0[j] = pa_x[j] * tpy_xxy_z_0[j] + fl1_fx * tpy_xy_z_0[j];

            tpz_xxxy_z_0[j] = pa_x[j] * tpz_xxy_z_0[j] + fl1_fx * tpz_xy_z_0[j];

            tpx_xxxz_x_0[j] = pa_x[j] * tpx_xxz_x_0[j] + fl1_fx * tpx_xz_x_0[j] + 0.5 * fl1_fx * tpx_xxz_0_0[j] - fl1_fgb * fl1_fx * ts_xxz_x_0[j];

            tpy_xxxz_x_0[j] = pa_x[j] * tpy_xxz_x_0[j] + fl1_fx * tpy_xz_x_0[j] + 0.5 * fl1_fx * tpy_xxz_0_0[j];

            tpz_xxxz_x_0[j] = pa_x[j] * tpz_xxz_x_0[j] + fl1_fx * tpz_xz_x_0[j] + 0.5 * fl1_fx * tpz_xxz_0_0[j];

            tpx_xxxz_y_0[j] = pa_x[j] * tpx_xxz_y_0[j] + fl1_fx * tpx_xz_y_0[j] - fl1_fgb * fl1_fx * ts_xxz_y_0[j];

            tpy_xxxz_y_0[j] = pa_x[j] * tpy_xxz_y_0[j] + fl1_fx * tpy_xz_y_0[j];

            tpz_xxxz_y_0[j] = pa_x[j] * tpz_xxz_y_0[j] + fl1_fx * tpz_xz_y_0[j];

            tpx_xxxz_z_0[j] = pa_x[j] * tpx_xxz_z_0[j] + fl1_fx * tpx_xz_z_0[j] - fl1_fgb * fl1_fx * ts_xxz_z_0[j];

            tpy_xxxz_z_0[j] = pa_x[j] * tpy_xxz_z_0[j] + fl1_fx * tpy_xz_z_0[j];

            tpz_xxxz_z_0[j] = pa_x[j] * tpz_xxz_z_0[j] + fl1_fx * tpz_xz_z_0[j];

            tpx_xxyy_x_0[j] =
                pa_x[j] * tpx_xyy_x_0[j] + 0.5 * fl1_fx * tpx_yy_x_0[j] + 0.5 * fl1_fx * tpx_xyy_0_0[j] - fl1_fgb * fl1_fx * ts_xyy_x_0[j];

            tpy_xxyy_x_0[j] = pa_x[j] * tpy_xyy_x_0[j] + 0.5 * fl1_fx * tpy_yy_x_0[j] + 0.5 * fl1_fx * tpy_xyy_0_0[j];

            tpz_xxyy_x_0[j] = pa_x[j] * tpz_xyy_x_0[j] + 0.5 * fl1_fx * tpz_yy_x_0[j] + 0.5 * fl1_fx * tpz_xyy_0_0[j];

            tpx_xxyy_y_0[j] = pa_x[j] * tpx_xyy_y_0[j] + 0.5 * fl1_fx * tpx_yy_y_0[j] - fl1_fgb * fl1_fx * ts_xyy_y_0[j];

            tpy_xxyy_y_0[j] = pa_x[j] * tpy_xyy_y_0[j] + 0.5 * fl1_fx * tpy_yy_y_0[j];

            tpz_xxyy_y_0[j] = pa_x[j] * tpz_xyy_y_0[j] + 0.5 * fl1_fx * tpz_yy_y_0[j];

            tpx_xxyy_z_0[j] = pa_x[j] * tpx_xyy_z_0[j] + 0.5 * fl1_fx * tpx_yy_z_0[j] - fl1_fgb * fl1_fx * ts_xyy_z_0[j];

            tpy_xxyy_z_0[j] = pa_x[j] * tpy_xyy_z_0[j] + 0.5 * fl1_fx * tpy_yy_z_0[j];

            tpz_xxyy_z_0[j] = pa_x[j] * tpz_xyy_z_0[j] + 0.5 * fl1_fx * tpz_yy_z_0[j];

            tpx_xxyz_x_0[j] =
                pa_x[j] * tpx_xyz_x_0[j] + 0.5 * fl1_fx * tpx_yz_x_0[j] + 0.5 * fl1_fx * tpx_xyz_0_0[j] - fl1_fgb * fl1_fx * ts_xyz_x_0[j];

            tpy_xxyz_x_0[j] = pa_x[j] * tpy_xyz_x_0[j] + 0.5 * fl1_fx * tpy_yz_x_0[j] + 0.5 * fl1_fx * tpy_xyz_0_0[j];

            tpz_xxyz_x_0[j] = pa_x[j] * tpz_xyz_x_0[j] + 0.5 * fl1_fx * tpz_yz_x_0[j] + 0.5 * fl1_fx * tpz_xyz_0_0[j];

            tpx_xxyz_y_0[j] = pa_x[j] * tpx_xyz_y_0[j] + 0.5 * fl1_fx * tpx_yz_y_0[j] - fl1_fgb * fl1_fx * ts_xyz_y_0[j];

            tpy_xxyz_y_0[j] = pa_x[j] * tpy_xyz_y_0[j] + 0.5 * fl1_fx * tpy_yz_y_0[j];

            tpz_xxyz_y_0[j] = pa_x[j] * tpz_xyz_y_0[j] + 0.5 * fl1_fx * tpz_yz_y_0[j];

            tpx_xxyz_z_0[j] = pa_x[j] * tpx_xyz_z_0[j] + 0.5 * fl1_fx * tpx_yz_z_0[j] - fl1_fgb * fl1_fx * ts_xyz_z_0[j];

            tpy_xxyz_z_0[j] = pa_x[j] * tpy_xyz_z_0[j] + 0.5 * fl1_fx * tpy_yz_z_0[j];

            tpz_xxyz_z_0[j] = pa_x[j] * tpz_xyz_z_0[j] + 0.5 * fl1_fx * tpz_yz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGP_45_90(CMemBlock2D<double>&       primBuffer,
                              const CRecursionMap&       recursionMap,
                              const CMemBlock2D<double>& osFactors,
                              const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 15);

        auto tpy_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 15);

        auto tpz_xzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 15);

        auto tpx_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 16);

        auto tpy_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 16);

        auto tpz_xzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 16);

        auto tpx_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 17);

        auto tpy_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 17);

        auto tpz_xzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 17);

        auto tpx_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 18);

        auto tpy_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 19);

        auto tpy_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 20);

        auto tpy_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 21);

        auto tpy_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 22);

        auto tpy_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 23);

        auto tpy_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 24);

        auto tpy_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 25);

        auto tpy_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 26);

        auto tpy_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 27);

        auto tpy_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 28);

        auto tpy_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 29);

        auto tpy_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 29);

        auto tpx_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 15);

        auto tpy_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 16);

        auto tpy_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 17);

        auto tpy_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 17);

        auto tpx_xzz_0_0 = primBuffer.data(pidx_p_3_0_m0 + 10 * idx + 5);

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

        auto ts_xzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 15);

        auto ts_xzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 16);

        auto ts_xzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 17);

        auto ts_yyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 18);

        auto ts_yyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 19);

        auto ts_yyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 20);

        auto ts_yyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 21);

        auto ts_yyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 22);

        auto ts_yyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 23);

        auto ts_yzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 24);

        auto ts_yzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 25);

        auto ts_yzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 26);

        auto ts_zzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 27);

        auto ts_zzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 28);

        auto ts_zzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 29);

        // set up pointers to integrals

        auto tpx_xxzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 15);

        auto tpy_xxzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 15);

        auto tpz_xxzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 15);

        auto tpx_xxzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 16);

        auto tpy_xxzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 16);

        auto tpz_xxzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 16);

        auto tpx_xxzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 17);

        auto tpy_xxzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 17);

        auto tpz_xxzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 17);

        auto tpx_xyyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 18);

        auto tpy_xyyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 18);

        auto tpz_xyyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 18);

        auto tpx_xyyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 19);

        auto tpy_xyyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 19);

        auto tpz_xyyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 19);

        auto tpx_xyyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 20);

        auto tpy_xyyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 20);

        auto tpz_xyyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 20);

        auto tpx_xyyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 21);

        auto tpy_xyyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 21);

        auto tpz_xyyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 21);

        auto tpx_xyyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 22);

        auto tpy_xyyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 22);

        auto tpz_xyyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 22);

        auto tpx_xyyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 23);

        auto tpy_xyyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 23);

        auto tpz_xyyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 23);

        auto tpx_xyzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 24);

        auto tpy_xyzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 24);

        auto tpz_xyzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 24);

        auto tpx_xyzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 25);

        auto tpy_xyzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 25);

        auto tpz_xyzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 25);

        auto tpx_xyzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 26);

        auto tpy_xyzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 26);

        auto tpz_xyzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 26);

        auto tpx_xzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 27);

        auto tpy_xzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 27);

        auto tpz_xzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 27);

        auto tpx_xzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 28);

        auto tpy_xzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 28);

        auto tpz_xzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 28);

        auto tpx_xzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 29);

        auto tpy_xzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 29);

        auto tpz_xzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 29);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fgb, fx, pa_x, tpx_xxzz_x_0, tpx_xxzz_y_0, tpx_xxzz_z_0, tpx_xyyy_x_0, \
                                     tpx_xyyy_y_0, tpx_xyyy_z_0, tpx_xyyz_x_0, tpx_xyyz_y_0, tpx_xyyz_z_0, tpx_xyzz_x_0, \
                                     tpx_xyzz_y_0, tpx_xyzz_z_0, tpx_xzz_0_0, tpx_xzz_x_0, tpx_xzz_y_0, tpx_xzz_z_0, \
                                     tpx_xzzz_x_0, tpx_xzzz_y_0, tpx_xzzz_z_0, tpx_yyy_0_0, tpx_yyy_x_0, tpx_yyy_y_0, \
                                     tpx_yyy_z_0, tpx_yyz_0_0, tpx_yyz_x_0, tpx_yyz_y_0, tpx_yyz_z_0, tpx_yzz_0_0, \
                                     tpx_yzz_x_0, tpx_yzz_y_0, tpx_yzz_z_0, tpx_zz_x_0, tpx_zz_y_0, tpx_zz_z_0, \
                                     tpx_zzz_0_0, tpx_zzz_x_0, tpx_zzz_y_0, tpx_zzz_z_0, tpy_xxzz_x_0, tpy_xxzz_y_0, \
                                     tpy_xxzz_z_0, tpy_xyyy_x_0, tpy_xyyy_y_0, tpy_xyyy_z_0, tpy_xyyz_x_0, tpy_xyyz_y_0, \
                                     tpy_xyyz_z_0, tpy_xyzz_x_0, tpy_xyzz_y_0, tpy_xyzz_z_0, tpy_xzz_0_0, tpy_xzz_x_0, \
                                     tpy_xzz_y_0, tpy_xzz_z_0, tpy_xzzz_x_0, tpy_xzzz_y_0, tpy_xzzz_z_0, tpy_yyy_0_0, \
                                     tpy_yyy_x_0, tpy_yyy_y_0, tpy_yyy_z_0, tpy_yyz_0_0, tpy_yyz_x_0, tpy_yyz_y_0, \
                                     tpy_yyz_z_0, tpy_yzz_0_0, tpy_yzz_x_0, tpy_yzz_y_0, tpy_yzz_z_0, tpy_zz_x_0, \
                                     tpy_zz_y_0, tpy_zz_z_0, tpy_zzz_0_0, tpy_zzz_x_0, tpy_zzz_y_0, tpy_zzz_z_0, \
                                     tpz_xxzz_x_0, tpz_xxzz_y_0, tpz_xxzz_z_0, tpz_xyyy_x_0, tpz_xyyy_y_0, tpz_xyyy_z_0, \
                                     tpz_xyyz_x_0, tpz_xyyz_y_0, tpz_xyyz_z_0, tpz_xyzz_x_0, tpz_xyzz_y_0, tpz_xyzz_z_0, \
                                     tpz_xzz_0_0, tpz_xzz_x_0, tpz_xzz_y_0, tpz_xzz_z_0, tpz_xzzz_x_0, tpz_xzzz_y_0, \
                                     tpz_xzzz_z_0, tpz_yyy_0_0, tpz_yyy_x_0, tpz_yyy_y_0, tpz_yyy_z_0, tpz_yyz_0_0, \
                                     tpz_yyz_x_0, tpz_yyz_y_0, tpz_yyz_z_0, tpz_yzz_0_0, tpz_yzz_x_0, tpz_yzz_y_0, \
                                     tpz_yzz_z_0, tpz_zz_x_0, tpz_zz_y_0, tpz_zz_z_0, tpz_zzz_0_0, tpz_zzz_x_0, \
                                     tpz_zzz_y_0, tpz_zzz_z_0, ts_xzz_x_0, ts_xzz_y_0, ts_xzz_z_0, ts_yyy_x_0, \
                                     ts_yyy_y_0, ts_yyy_z_0, ts_yyz_x_0, ts_yyz_y_0, ts_yyz_z_0, ts_yzz_x_0, ts_yzz_y_0, \
                                     ts_yzz_z_0, ts_zzz_x_0, ts_zzz_y_0, ts_zzz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_xxzz_x_0[j] =
                pa_x[j] * tpx_xzz_x_0[j] + 0.5 * fl1_fx * tpx_zz_x_0[j] + 0.5 * fl1_fx * tpx_xzz_0_0[j] - fl1_fgb * fl1_fx * ts_xzz_x_0[j];

            tpy_xxzz_x_0[j] = pa_x[j] * tpy_xzz_x_0[j] + 0.5 * fl1_fx * tpy_zz_x_0[j] + 0.5 * fl1_fx * tpy_xzz_0_0[j];

            tpz_xxzz_x_0[j] = pa_x[j] * tpz_xzz_x_0[j] + 0.5 * fl1_fx * tpz_zz_x_0[j] + 0.5 * fl1_fx * tpz_xzz_0_0[j];

            tpx_xxzz_y_0[j] = pa_x[j] * tpx_xzz_y_0[j] + 0.5 * fl1_fx * tpx_zz_y_0[j] - fl1_fgb * fl1_fx * ts_xzz_y_0[j];

            tpy_xxzz_y_0[j] = pa_x[j] * tpy_xzz_y_0[j] + 0.5 * fl1_fx * tpy_zz_y_0[j];

            tpz_xxzz_y_0[j] = pa_x[j] * tpz_xzz_y_0[j] + 0.5 * fl1_fx * tpz_zz_y_0[j];

            tpx_xxzz_z_0[j] = pa_x[j] * tpx_xzz_z_0[j] + 0.5 * fl1_fx * tpx_zz_z_0[j] - fl1_fgb * fl1_fx * ts_xzz_z_0[j];

            tpy_xxzz_z_0[j] = pa_x[j] * tpy_xzz_z_0[j] + 0.5 * fl1_fx * tpy_zz_z_0[j];

            tpz_xxzz_z_0[j] = pa_x[j] * tpz_xzz_z_0[j] + 0.5 * fl1_fx * tpz_zz_z_0[j];

            tpx_xyyy_x_0[j] = pa_x[j] * tpx_yyy_x_0[j] + 0.5 * fl1_fx * tpx_yyy_0_0[j] - fl1_fgb * fl1_fx * ts_yyy_x_0[j];

            tpy_xyyy_x_0[j] = pa_x[j] * tpy_yyy_x_0[j] + 0.5 * fl1_fx * tpy_yyy_0_0[j];

            tpz_xyyy_x_0[j] = pa_x[j] * tpz_yyy_x_0[j] + 0.5 * fl1_fx * tpz_yyy_0_0[j];

            tpx_xyyy_y_0[j] = pa_x[j] * tpx_yyy_y_0[j] - fl1_fgb * fl1_fx * ts_yyy_y_0[j];

            tpy_xyyy_y_0[j] = pa_x[j] * tpy_yyy_y_0[j];

            tpz_xyyy_y_0[j] = pa_x[j] * tpz_yyy_y_0[j];

            tpx_xyyy_z_0[j] = pa_x[j] * tpx_yyy_z_0[j] - fl1_fgb * fl1_fx * ts_yyy_z_0[j];

            tpy_xyyy_z_0[j] = pa_x[j] * tpy_yyy_z_0[j];

            tpz_xyyy_z_0[j] = pa_x[j] * tpz_yyy_z_0[j];

            tpx_xyyz_x_0[j] = pa_x[j] * tpx_yyz_x_0[j] + 0.5 * fl1_fx * tpx_yyz_0_0[j] - fl1_fgb * fl1_fx * ts_yyz_x_0[j];

            tpy_xyyz_x_0[j] = pa_x[j] * tpy_yyz_x_0[j] + 0.5 * fl1_fx * tpy_yyz_0_0[j];

            tpz_xyyz_x_0[j] = pa_x[j] * tpz_yyz_x_0[j] + 0.5 * fl1_fx * tpz_yyz_0_0[j];

            tpx_xyyz_y_0[j] = pa_x[j] * tpx_yyz_y_0[j] - fl1_fgb * fl1_fx * ts_yyz_y_0[j];

            tpy_xyyz_y_0[j] = pa_x[j] * tpy_yyz_y_0[j];

            tpz_xyyz_y_0[j] = pa_x[j] * tpz_yyz_y_0[j];

            tpx_xyyz_z_0[j] = pa_x[j] * tpx_yyz_z_0[j] - fl1_fgb * fl1_fx * ts_yyz_z_0[j];

            tpy_xyyz_z_0[j] = pa_x[j] * tpy_yyz_z_0[j];

            tpz_xyyz_z_0[j] = pa_x[j] * tpz_yyz_z_0[j];

            tpx_xyzz_x_0[j] = pa_x[j] * tpx_yzz_x_0[j] + 0.5 * fl1_fx * tpx_yzz_0_0[j] - fl1_fgb * fl1_fx * ts_yzz_x_0[j];

            tpy_xyzz_x_0[j] = pa_x[j] * tpy_yzz_x_0[j] + 0.5 * fl1_fx * tpy_yzz_0_0[j];

            tpz_xyzz_x_0[j] = pa_x[j] * tpz_yzz_x_0[j] + 0.5 * fl1_fx * tpz_yzz_0_0[j];

            tpx_xyzz_y_0[j] = pa_x[j] * tpx_yzz_y_0[j] - fl1_fgb * fl1_fx * ts_yzz_y_0[j];

            tpy_xyzz_y_0[j] = pa_x[j] * tpy_yzz_y_0[j];

            tpz_xyzz_y_0[j] = pa_x[j] * tpz_yzz_y_0[j];

            tpx_xyzz_z_0[j] = pa_x[j] * tpx_yzz_z_0[j] - fl1_fgb * fl1_fx * ts_yzz_z_0[j];

            tpy_xyzz_z_0[j] = pa_x[j] * tpy_yzz_z_0[j];

            tpz_xyzz_z_0[j] = pa_x[j] * tpz_yzz_z_0[j];

            tpx_xzzz_x_0[j] = pa_x[j] * tpx_zzz_x_0[j] + 0.5 * fl1_fx * tpx_zzz_0_0[j] - fl1_fgb * fl1_fx * ts_zzz_x_0[j];

            tpy_xzzz_x_0[j] = pa_x[j] * tpy_zzz_x_0[j] + 0.5 * fl1_fx * tpy_zzz_0_0[j];

            tpz_xzzz_x_0[j] = pa_x[j] * tpz_zzz_x_0[j] + 0.5 * fl1_fx * tpz_zzz_0_0[j];

            tpx_xzzz_y_0[j] = pa_x[j] * tpx_zzz_y_0[j] - fl1_fgb * fl1_fx * ts_zzz_y_0[j];

            tpy_xzzz_y_0[j] = pa_x[j] * tpy_zzz_y_0[j];

            tpz_xzzz_y_0[j] = pa_x[j] * tpz_zzz_y_0[j];

            tpx_xzzz_z_0[j] = pa_x[j] * tpx_zzz_z_0[j] - fl1_fgb * fl1_fx * ts_zzz_z_0[j];

            tpy_xzzz_z_0[j] = pa_x[j] * tpy_zzz_z_0[j];

            tpz_xzzz_z_0[j] = pa_x[j] * tpz_zzz_z_0[j];
        }

        idx++;
    }
}

void
compLinearMomentumForGP_90_135(CMemBlock2D<double>&       primBuffer,
                               const CRecursionMap&       recursionMap,
                               const CMemBlock2D<double>& osFactors,
                               const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_p_4_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {4, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_p_4_1_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_p_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_p_3_0_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Linear Momentum"}, 1, true, {3, -1, -1, -1}, {0, -1, -1, -1}, 1, 1, 0));

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

        auto tpx_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 18);

        auto tpy_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 18);

        auto tpz_yyy_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 18);

        auto tpx_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 19);

        auto tpy_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 19);

        auto tpz_yyy_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 19);

        auto tpx_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 20);

        auto tpy_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 20);

        auto tpz_yyy_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 20);

        auto tpx_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 21);

        auto tpy_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 21);

        auto tpz_yyz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 21);

        auto tpx_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 22);

        auto tpy_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 22);

        auto tpz_yyz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 22);

        auto tpx_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 23);

        auto tpy_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 23);

        auto tpz_yyz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 23);

        auto tpx_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 24);

        auto tpy_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 24);

        auto tpz_yzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 24);

        auto tpx_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 25);

        auto tpy_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 25);

        auto tpz_yzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 25);

        auto tpx_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 26);

        auto tpy_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 26);

        auto tpz_yzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 26);

        auto tpx_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 27);

        auto tpy_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 27);

        auto tpz_zzz_x_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 27);

        auto tpx_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 28);

        auto tpy_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 28);

        auto tpz_zzz_y_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 28);

        auto tpx_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * idx + 29);

        auto tpy_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 30 * bdim + 30 * idx + 29);

        auto tpz_zzz_z_0 = primBuffer.data(pidx_p_3_1_m0 + 60 * bdim + 30 * idx + 29);

        auto tpx_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 9);

        auto tpy_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 9);

        auto tpz_yy_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 9);

        auto tpx_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 10);

        auto tpy_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 10);

        auto tpz_yy_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 10);

        auto tpx_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 11);

        auto tpy_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 11);

        auto tpz_yy_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 11);

        auto tpx_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 12);

        auto tpy_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 12);

        auto tpz_yz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 12);

        auto tpx_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 13);

        auto tpy_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 13);

        auto tpz_yz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 13);

        auto tpx_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 14);

        auto tpy_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 14);

        auto tpz_yz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 14);

        auto tpx_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 15);

        auto tpy_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 15);

        auto tpz_zz_x_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 15);

        auto tpx_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 16);

        auto tpy_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 16);

        auto tpz_zz_y_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 16);

        auto tpx_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * idx + 17);

        auto tpy_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 18 * bdim + 18 * idx + 17);

        auto tpz_zz_z_0 = primBuffer.data(pidx_p_2_1_m0 + 36 * bdim + 18 * idx + 17);

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

        auto ts_yyy_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 18);

        auto ts_yyy_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 19);

        auto ts_yyy_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 20);

        auto ts_yyz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 21);

        auto ts_yyz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 22);

        auto ts_yyz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 23);

        auto ts_yzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 24);

        auto ts_yzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 25);

        auto ts_yzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 26);

        auto ts_zzz_x_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 27);

        auto ts_zzz_y_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 28);

        auto ts_zzz_z_0 = primBuffer.data(pidx_s_3_1_m0 + 30 * idx + 29);

        // set up pointers to integrals

        auto tpx_yyyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 30);

        auto tpy_yyyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 30);

        auto tpz_yyyy_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 30);

        auto tpx_yyyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 31);

        auto tpy_yyyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 31);

        auto tpz_yyyy_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 31);

        auto tpx_yyyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 32);

        auto tpy_yyyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 32);

        auto tpz_yyyy_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 32);

        auto tpx_yyyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 33);

        auto tpy_yyyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 33);

        auto tpz_yyyz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 33);

        auto tpx_yyyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 34);

        auto tpy_yyyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 34);

        auto tpz_yyyz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 34);

        auto tpx_yyyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 35);

        auto tpy_yyyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 35);

        auto tpz_yyyz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 35);

        auto tpx_yyzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 36);

        auto tpy_yyzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 36);

        auto tpz_yyzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 36);

        auto tpx_yyzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 37);

        auto tpy_yyzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 37);

        auto tpz_yyzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 37);

        auto tpx_yyzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 38);

        auto tpy_yyzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 38);

        auto tpz_yyzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 38);

        auto tpx_yzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 39);

        auto tpy_yzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 39);

        auto tpz_yzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 39);

        auto tpx_yzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 40);

        auto tpy_yzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 40);

        auto tpz_yzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 40);

        auto tpx_yzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 41);

        auto tpy_yzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 41);

        auto tpz_yzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 41);

        auto tpx_zzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 42);

        auto tpy_zzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 42);

        auto tpz_zzzz_x_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 42);

        auto tpx_zzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 43);

        auto tpy_zzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 43);

        auto tpz_zzzz_y_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 43);

        auto tpx_zzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * idx + 44);

        auto tpy_zzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 45 * bdim + 45 * idx + 44);

        auto tpz_zzzz_z_0 = primBuffer.data(pidx_p_4_1_m0 + 90 * bdim + 45 * idx + 44);

        // Batch of Integrals (90,135)

        #pragma omp simd aligned(fgb, fx, pa_y, pa_z, tpx_yy_x_0, tpx_yy_y_0, tpx_yy_z_0, tpx_yyy_0_0, \
                                     tpx_yyy_x_0, tpx_yyy_y_0, tpx_yyy_z_0, tpx_yyyy_x_0, tpx_yyyy_y_0, tpx_yyyy_z_0, \
                                     tpx_yyyz_x_0, tpx_yyyz_y_0, tpx_yyyz_z_0, tpx_yyz_0_0, tpx_yyz_x_0, tpx_yyz_y_0, \
                                     tpx_yyz_z_0, tpx_yyzz_x_0, tpx_yyzz_y_0, tpx_yyzz_z_0, tpx_yz_x_0, tpx_yz_y_0, \
                                     tpx_yz_z_0, tpx_yzz_0_0, tpx_yzz_x_0, tpx_yzz_y_0, tpx_yzz_z_0, tpx_yzzz_x_0, \
                                     tpx_yzzz_y_0, tpx_yzzz_z_0, tpx_zz_x_0, tpx_zz_y_0, tpx_zz_z_0, tpx_zzz_0_0, \
                                     tpx_zzz_x_0, tpx_zzz_y_0, tpx_zzz_z_0, tpx_zzzz_x_0, tpx_zzzz_y_0, tpx_zzzz_z_0, \
                                     tpy_yy_x_0, tpy_yy_y_0, tpy_yy_z_0, tpy_yyy_0_0, tpy_yyy_x_0, tpy_yyy_y_0, \
                                     tpy_yyy_z_0, tpy_yyyy_x_0, tpy_yyyy_y_0, tpy_yyyy_z_0, tpy_yyyz_x_0, tpy_yyyz_y_0, \
                                     tpy_yyyz_z_0, tpy_yyz_0_0, tpy_yyz_x_0, tpy_yyz_y_0, tpy_yyz_z_0, tpy_yyzz_x_0, \
                                     tpy_yyzz_y_0, tpy_yyzz_z_0, tpy_yz_x_0, tpy_yz_y_0, tpy_yz_z_0, tpy_yzz_0_0, \
                                     tpy_yzz_x_0, tpy_yzz_y_0, tpy_yzz_z_0, tpy_yzzz_x_0, tpy_yzzz_y_0, tpy_yzzz_z_0, \
                                     tpy_zz_x_0, tpy_zz_y_0, tpy_zz_z_0, tpy_zzz_0_0, tpy_zzz_x_0, tpy_zzz_y_0, \
                                     tpy_zzz_z_0, tpy_zzzz_x_0, tpy_zzzz_y_0, tpy_zzzz_z_0, tpz_yy_x_0, tpz_yy_y_0, \
                                     tpz_yy_z_0, tpz_yyy_0_0, tpz_yyy_x_0, tpz_yyy_y_0, tpz_yyy_z_0, tpz_yyyy_x_0, \
                                     tpz_yyyy_y_0, tpz_yyyy_z_0, tpz_yyyz_x_0, tpz_yyyz_y_0, tpz_yyyz_z_0, tpz_yyz_0_0, \
                                     tpz_yyz_x_0, tpz_yyz_y_0, tpz_yyz_z_0, tpz_yyzz_x_0, tpz_yyzz_y_0, tpz_yyzz_z_0, \
                                     tpz_yz_x_0, tpz_yz_y_0, tpz_yz_z_0, tpz_yzz_0_0, tpz_yzz_x_0, tpz_yzz_y_0, \
                                     tpz_yzz_z_0, tpz_yzzz_x_0, tpz_yzzz_y_0, tpz_yzzz_z_0, tpz_zz_x_0, tpz_zz_y_0, \
                                     tpz_zz_z_0, tpz_zzz_0_0, tpz_zzz_x_0, tpz_zzz_y_0, tpz_zzz_z_0, tpz_zzzz_x_0, \
                                     tpz_zzzz_y_0, tpz_zzzz_z_0, ts_yyy_x_0, ts_yyy_y_0, ts_yyy_z_0, ts_yyz_x_0, \
                                     ts_yyz_y_0, ts_yyz_z_0, ts_yzz_x_0, ts_yzz_y_0, ts_yzz_z_0, ts_zzz_x_0, ts_zzz_y_0, \
                                     ts_zzz_z_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fgb = fgb[j];

            double fl1_fx = fx[j];

            tpx_yyyy_x_0[j] = pa_y[j] * tpx_yyy_x_0[j] + 1.5 * fl1_fx * tpx_yy_x_0[j];

            tpy_yyyy_x_0[j] = pa_y[j] * tpy_yyy_x_0[j] + 1.5 * fl1_fx * tpy_yy_x_0[j] - fl1_fgb * fl1_fx * ts_yyy_x_0[j];

            tpz_yyyy_x_0[j] = pa_y[j] * tpz_yyy_x_0[j] + 1.5 * fl1_fx * tpz_yy_x_0[j];

            tpx_yyyy_y_0[j] = pa_y[j] * tpx_yyy_y_0[j] + 1.5 * fl1_fx * tpx_yy_y_0[j] + 0.5 * fl1_fx * tpx_yyy_0_0[j];

            tpy_yyyy_y_0[j] =
                pa_y[j] * tpy_yyy_y_0[j] + 1.5 * fl1_fx * tpy_yy_y_0[j] + 0.5 * fl1_fx * tpy_yyy_0_0[j] - fl1_fgb * fl1_fx * ts_yyy_y_0[j];

            tpz_yyyy_y_0[j] = pa_y[j] * tpz_yyy_y_0[j] + 1.5 * fl1_fx * tpz_yy_y_0[j] + 0.5 * fl1_fx * tpz_yyy_0_0[j];

            tpx_yyyy_z_0[j] = pa_y[j] * tpx_yyy_z_0[j] + 1.5 * fl1_fx * tpx_yy_z_0[j];

            tpy_yyyy_z_0[j] = pa_y[j] * tpy_yyy_z_0[j] + 1.5 * fl1_fx * tpy_yy_z_0[j] - fl1_fgb * fl1_fx * ts_yyy_z_0[j];

            tpz_yyyy_z_0[j] = pa_y[j] * tpz_yyy_z_0[j] + 1.5 * fl1_fx * tpz_yy_z_0[j];

            tpx_yyyz_x_0[j] = pa_y[j] * tpx_yyz_x_0[j] + fl1_fx * tpx_yz_x_0[j];

            tpy_yyyz_x_0[j] = pa_y[j] * tpy_yyz_x_0[j] + fl1_fx * tpy_yz_x_0[j] - fl1_fgb * fl1_fx * ts_yyz_x_0[j];

            tpz_yyyz_x_0[j] = pa_y[j] * tpz_yyz_x_0[j] + fl1_fx * tpz_yz_x_0[j];

            tpx_yyyz_y_0[j] = pa_y[j] * tpx_yyz_y_0[j] + fl1_fx * tpx_yz_y_0[j] + 0.5 * fl1_fx * tpx_yyz_0_0[j];

            tpy_yyyz_y_0[j] = pa_y[j] * tpy_yyz_y_0[j] + fl1_fx * tpy_yz_y_0[j] + 0.5 * fl1_fx * tpy_yyz_0_0[j] - fl1_fgb * fl1_fx * ts_yyz_y_0[j];

            tpz_yyyz_y_0[j] = pa_y[j] * tpz_yyz_y_0[j] + fl1_fx * tpz_yz_y_0[j] + 0.5 * fl1_fx * tpz_yyz_0_0[j];

            tpx_yyyz_z_0[j] = pa_y[j] * tpx_yyz_z_0[j] + fl1_fx * tpx_yz_z_0[j];

            tpy_yyyz_z_0[j] = pa_y[j] * tpy_yyz_z_0[j] + fl1_fx * tpy_yz_z_0[j] - fl1_fgb * fl1_fx * ts_yyz_z_0[j];

            tpz_yyyz_z_0[j] = pa_y[j] * tpz_yyz_z_0[j] + fl1_fx * tpz_yz_z_0[j];

            tpx_yyzz_x_0[j] = pa_y[j] * tpx_yzz_x_0[j] + 0.5 * fl1_fx * tpx_zz_x_0[j];

            tpy_yyzz_x_0[j] = pa_y[j] * tpy_yzz_x_0[j] + 0.5 * fl1_fx * tpy_zz_x_0[j] - fl1_fgb * fl1_fx * ts_yzz_x_0[j];

            tpz_yyzz_x_0[j] = pa_y[j] * tpz_yzz_x_0[j] + 0.5 * fl1_fx * tpz_zz_x_0[j];

            tpx_yyzz_y_0[j] = pa_y[j] * tpx_yzz_y_0[j] + 0.5 * fl1_fx * tpx_zz_y_0[j] + 0.5 * fl1_fx * tpx_yzz_0_0[j];

            tpy_yyzz_y_0[j] =
                pa_y[j] * tpy_yzz_y_0[j] + 0.5 * fl1_fx * tpy_zz_y_0[j] + 0.5 * fl1_fx * tpy_yzz_0_0[j] - fl1_fgb * fl1_fx * ts_yzz_y_0[j];

            tpz_yyzz_y_0[j] = pa_y[j] * tpz_yzz_y_0[j] + 0.5 * fl1_fx * tpz_zz_y_0[j] + 0.5 * fl1_fx * tpz_yzz_0_0[j];

            tpx_yyzz_z_0[j] = pa_y[j] * tpx_yzz_z_0[j] + 0.5 * fl1_fx * tpx_zz_z_0[j];

            tpy_yyzz_z_0[j] = pa_y[j] * tpy_yzz_z_0[j] + 0.5 * fl1_fx * tpy_zz_z_0[j] - fl1_fgb * fl1_fx * ts_yzz_z_0[j];

            tpz_yyzz_z_0[j] = pa_y[j] * tpz_yzz_z_0[j] + 0.5 * fl1_fx * tpz_zz_z_0[j];

            tpx_yzzz_x_0[j] = pa_y[j] * tpx_zzz_x_0[j];

            tpy_yzzz_x_0[j] = pa_y[j] * tpy_zzz_x_0[j] - fl1_fgb * fl1_fx * ts_zzz_x_0[j];

            tpz_yzzz_x_0[j] = pa_y[j] * tpz_zzz_x_0[j];

            tpx_yzzz_y_0[j] = pa_y[j] * tpx_zzz_y_0[j] + 0.5 * fl1_fx * tpx_zzz_0_0[j];

            tpy_yzzz_y_0[j] = pa_y[j] * tpy_zzz_y_0[j] + 0.5 * fl1_fx * tpy_zzz_0_0[j] - fl1_fgb * fl1_fx * ts_zzz_y_0[j];

            tpz_yzzz_y_0[j] = pa_y[j] * tpz_zzz_y_0[j] + 0.5 * fl1_fx * tpz_zzz_0_0[j];

            tpx_yzzz_z_0[j] = pa_y[j] * tpx_zzz_z_0[j];

            tpy_yzzz_z_0[j] = pa_y[j] * tpy_zzz_z_0[j] - fl1_fgb * fl1_fx * ts_zzz_z_0[j];

            tpz_yzzz_z_0[j] = pa_y[j] * tpz_zzz_z_0[j];

            tpx_zzzz_x_0[j] = pa_z[j] * tpx_zzz_x_0[j] + 1.5 * fl1_fx * tpx_zz_x_0[j];

            tpy_zzzz_x_0[j] = pa_z[j] * tpy_zzz_x_0[j] + 1.5 * fl1_fx * tpy_zz_x_0[j];

            tpz_zzzz_x_0[j] = pa_z[j] * tpz_zzz_x_0[j] + 1.5 * fl1_fx * tpz_zz_x_0[j] - fl1_fgb * fl1_fx * ts_zzz_x_0[j];

            tpx_zzzz_y_0[j] = pa_z[j] * tpx_zzz_y_0[j] + 1.5 * fl1_fx * tpx_zz_y_0[j];

            tpy_zzzz_y_0[j] = pa_z[j] * tpy_zzz_y_0[j] + 1.5 * fl1_fx * tpy_zz_y_0[j];

            tpz_zzzz_y_0[j] = pa_z[j] * tpz_zzz_y_0[j] + 1.5 * fl1_fx * tpz_zz_y_0[j] - fl1_fgb * fl1_fx * ts_zzz_y_0[j];

            tpx_zzzz_z_0[j] = pa_z[j] * tpx_zzz_z_0[j] + 1.5 * fl1_fx * tpx_zz_z_0[j] + 0.5 * fl1_fx * tpx_zzz_0_0[j];

            tpy_zzzz_z_0[j] = pa_z[j] * tpy_zzz_z_0[j] + 1.5 * fl1_fx * tpy_zz_z_0[j] + 0.5 * fl1_fx * tpy_zzz_0_0[j];

            tpz_zzzz_z_0[j] =
                pa_z[j] * tpz_zzz_z_0[j] + 1.5 * fl1_fx * tpz_zz_z_0[j] + 0.5 * fl1_fx * tpz_zzz_0_0[j] - fl1_fgb * fl1_fx * ts_zzz_z_0[j];
        }

        idx++;
    }
}

}  // namespace lmomrecfunc
