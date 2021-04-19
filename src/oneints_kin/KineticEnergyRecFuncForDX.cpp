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

#include "KineticEnergyRecFuncForDX.hpp"

namespace kinrecfunc {  // kinrecfunc namespace

void
compKineticEnergyForDD(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_2_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tt_x_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx);

        auto tt_x_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 1);

        auto tt_x_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 2);

        auto tt_x_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 3);

        auto tt_x_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 4);

        auto tt_x_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 5);

        auto tt_y_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 6);

        auto tt_y_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 7);

        auto tt_y_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 8);

        auto tt_y_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 9);

        auto tt_y_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 10);

        auto tt_y_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 11);

        auto tt_z_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 12);

        auto tt_z_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 13);

        auto tt_z_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 14);

        auto tt_z_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 15);

        auto tt_z_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 16);

        auto tt_z_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 17);

        auto tt_0_xx_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx);

        auto tt_0_xy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 1);

        auto tt_0_xz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 2);

        auto tt_0_yy_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 3);

        auto tt_0_yz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 4);

        auto tt_0_zz_0 = primBuffer.data(pidx_t_0_2_m0 + 6 * idx + 5);

        auto tt_x_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx);

        auto tt_x_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 1);

        auto tt_x_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 2);

        auto tt_y_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 3);

        auto tt_y_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 4);

        auto tt_y_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 5);

        auto tt_z_x_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 6);

        auto tt_z_y_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 7);

        auto tt_z_z_0 = primBuffer.data(pidx_t_1_1_m0 + 9 * idx + 8);

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

        auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30);

        auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31);

        auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32);

        auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33);

        auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34);

        auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35);

        auto ts_0_xx_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx);

        auto ts_0_xy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 1);

        auto ts_0_xz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 2);

        auto ts_0_yy_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 3);

        auto ts_0_yz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 4);

        auto ts_0_zz_0 = primBuffer.data(pidx_s_0_2_m0 + 6 * idx + 5);

        // set up pointers to integrals

        auto tt_xx_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx);

        auto tt_xx_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 1);

        auto tt_xx_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 2);

        auto tt_xx_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 3);

        auto tt_xx_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 4);

        auto tt_xx_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 5);

        auto tt_xy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 6);

        auto tt_xy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 7);

        auto tt_xy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 8);

        auto tt_xy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 9);

        auto tt_xy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 10);

        auto tt_xy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 11);

        auto tt_xz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 12);

        auto tt_xz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 13);

        auto tt_xz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 14);

        auto tt_xz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 15);

        auto tt_xz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 16);

        auto tt_xz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 17);

        auto tt_yy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 18);

        auto tt_yy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 19);

        auto tt_yy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 20);

        auto tt_yy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 21);

        auto tt_yy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 22);

        auto tt_yy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 23);

        auto tt_yz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 24);

        auto tt_yz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 25);

        auto tt_yz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 26);

        auto tt_yz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 27);

        auto tt_yz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 28);

        auto tt_yz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 29);

        auto tt_zz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 30);

        auto tt_zz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 31);

        auto tt_zz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 32);

        auto tt_zz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 33);

        auto tt_zz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 34);

        auto tt_zz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 35);

        #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_0_xx_0, ts_0_xy_0, ts_0_xz_0, ts_0_yy_0, \
                                     ts_0_yz_0, ts_0_zz_0, ts_xx_xx_0, ts_xx_xy_0, ts_xx_xz_0, ts_xx_yy_0, ts_xx_yz_0, \
                                     ts_xx_zz_0, ts_xy_xx_0, ts_xy_xy_0, ts_xy_xz_0, ts_xy_yy_0, ts_xy_yz_0, ts_xy_zz_0, \
                                     ts_xz_xx_0, ts_xz_xy_0, ts_xz_xz_0, ts_xz_yy_0, ts_xz_yz_0, ts_xz_zz_0, ts_yy_xx_0, \
                                     ts_yy_xy_0, ts_yy_xz_0, ts_yy_yy_0, ts_yy_yz_0, ts_yy_zz_0, ts_yz_xx_0, ts_yz_xy_0, \
                                     ts_yz_xz_0, ts_yz_yy_0, ts_yz_yz_0, ts_yz_zz_0, ts_zz_xx_0, ts_zz_xy_0, ts_zz_xz_0, \
                                     ts_zz_yy_0, ts_zz_yz_0, ts_zz_zz_0, tt_0_xx_0, tt_0_xy_0, tt_0_xz_0, tt_0_yy_0, \
                                     tt_0_yz_0, tt_0_zz_0, tt_x_x_0, tt_x_xx_0, tt_x_xy_0, tt_x_xz_0, tt_x_y_0, \
                                     tt_x_yy_0, tt_x_yz_0, tt_x_z_0, tt_x_zz_0, tt_xx_xx_0, tt_xx_xy_0, tt_xx_xz_0, \
                                     tt_xx_yy_0, tt_xx_yz_0, tt_xx_zz_0, tt_xy_xx_0, tt_xy_xy_0, tt_xy_xz_0, tt_xy_yy_0, \
                                     tt_xy_yz_0, tt_xy_zz_0, tt_xz_xx_0, tt_xz_xy_0, tt_xz_xz_0, tt_xz_yy_0, tt_xz_yz_0, \
                                     tt_xz_zz_0, tt_y_x_0, tt_y_xx_0, tt_y_xy_0, tt_y_xz_0, tt_y_y_0, tt_y_yy_0, \
                                     tt_y_yz_0, tt_y_z_0, tt_y_zz_0, tt_yy_xx_0, tt_yy_xy_0, tt_yy_xz_0, tt_yy_yy_0, \
                                     tt_yy_yz_0, tt_yy_zz_0, tt_yz_xx_0, tt_yz_xy_0, tt_yz_xz_0, tt_yz_yy_0, tt_yz_yz_0, \
                                     tt_yz_zz_0, tt_z_x_0, tt_z_xx_0, tt_z_xy_0, tt_z_xz_0, tt_z_y_0, tt_z_yy_0, \
                                     tt_z_yz_0, tt_z_z_0, tt_z_zz_0, tt_zz_xx_0, tt_zz_xy_0, tt_zz_xz_0, tt_zz_yy_0, \
                                     tt_zz_yz_0, tt_zz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xx_xx_0[j] = pa_x[j] * tt_x_xx_0[j] + 0.5 * fl1_fx * tt_0_xx_0[j] + fl1_fx * tt_x_x_0[j] + 2.0 * fl1_fz * ts_xx_xx_0[j] -
                            fl1_fz * fl1_fga * ts_0_xx_0[j];

            tt_xx_xy_0[j] = pa_x[j] * tt_x_xy_0[j] + 0.5 * fl1_fx * tt_0_xy_0[j] + 0.5 * fl1_fx * tt_x_y_0[j] + 2.0 * fl1_fz * ts_xx_xy_0[j] -
                            fl1_fz * fl1_fga * ts_0_xy_0[j];

            tt_xx_xz_0[j] = pa_x[j] * tt_x_xz_0[j] + 0.5 * fl1_fx * tt_0_xz_0[j] + 0.5 * fl1_fx * tt_x_z_0[j] + 2.0 * fl1_fz * ts_xx_xz_0[j] -
                            fl1_fz * fl1_fga * ts_0_xz_0[j];

            tt_xx_yy_0[j] = pa_x[j] * tt_x_yy_0[j] + 0.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_xx_yy_0[j] - fl1_fz * fl1_fga * ts_0_yy_0[j];

            tt_xx_yz_0[j] = pa_x[j] * tt_x_yz_0[j] + 0.5 * fl1_fx * tt_0_yz_0[j] + 2.0 * fl1_fz * ts_xx_yz_0[j] - fl1_fz * fl1_fga * ts_0_yz_0[j];

            tt_xx_zz_0[j] = pa_x[j] * tt_x_zz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_xx_zz_0[j] - fl1_fz * fl1_fga * ts_0_zz_0[j];

            tt_xy_xx_0[j] = pa_x[j] * tt_y_xx_0[j] + fl1_fx * tt_y_x_0[j] + 2.0 * fl1_fz * ts_xy_xx_0[j];

            tt_xy_xy_0[j] = pa_x[j] * tt_y_xy_0[j] + 0.5 * fl1_fx * tt_y_y_0[j] + 2.0 * fl1_fz * ts_xy_xy_0[j];

            tt_xy_xz_0[j] = pa_x[j] * tt_y_xz_0[j] + 0.5 * fl1_fx * tt_y_z_0[j] + 2.0 * fl1_fz * ts_xy_xz_0[j];

            tt_xy_yy_0[j] = pa_x[j] * tt_y_yy_0[j] + 2.0 * fl1_fz * ts_xy_yy_0[j];

            tt_xy_yz_0[j] = pa_x[j] * tt_y_yz_0[j] + 2.0 * fl1_fz * ts_xy_yz_0[j];

            tt_xy_zz_0[j] = pa_x[j] * tt_y_zz_0[j] + 2.0 * fl1_fz * ts_xy_zz_0[j];

            tt_xz_xx_0[j] = pa_x[j] * tt_z_xx_0[j] + fl1_fx * tt_z_x_0[j] + 2.0 * fl1_fz * ts_xz_xx_0[j];

            tt_xz_xy_0[j] = pa_x[j] * tt_z_xy_0[j] + 0.5 * fl1_fx * tt_z_y_0[j] + 2.0 * fl1_fz * ts_xz_xy_0[j];

            tt_xz_xz_0[j] = pa_x[j] * tt_z_xz_0[j] + 0.5 * fl1_fx * tt_z_z_0[j] + 2.0 * fl1_fz * ts_xz_xz_0[j];

            tt_xz_yy_0[j] = pa_x[j] * tt_z_yy_0[j] + 2.0 * fl1_fz * ts_xz_yy_0[j];

            tt_xz_yz_0[j] = pa_x[j] * tt_z_yz_0[j] + 2.0 * fl1_fz * ts_xz_yz_0[j];

            tt_xz_zz_0[j] = pa_x[j] * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_xz_zz_0[j];

            tt_yy_xx_0[j] = pa_y[j] * tt_y_xx_0[j] + 0.5 * fl1_fx * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_yy_xx_0[j] - fl1_fz * fl1_fga * ts_0_xx_0[j];

            tt_yy_xy_0[j] = pa_y[j] * tt_y_xy_0[j] + 0.5 * fl1_fx * tt_0_xy_0[j] + 0.5 * fl1_fx * tt_y_x_0[j] + 2.0 * fl1_fz * ts_yy_xy_0[j] -
                            fl1_fz * fl1_fga * ts_0_xy_0[j];

            tt_yy_xz_0[j] = pa_y[j] * tt_y_xz_0[j] + 0.5 * fl1_fx * tt_0_xz_0[j] + 2.0 * fl1_fz * ts_yy_xz_0[j] - fl1_fz * fl1_fga * ts_0_xz_0[j];

            tt_yy_yy_0[j] = pa_y[j] * tt_y_yy_0[j] + 0.5 * fl1_fx * tt_0_yy_0[j] + fl1_fx * tt_y_y_0[j] + 2.0 * fl1_fz * ts_yy_yy_0[j] -
                            fl1_fz * fl1_fga * ts_0_yy_0[j];

            tt_yy_yz_0[j] = pa_y[j] * tt_y_yz_0[j] + 0.5 * fl1_fx * tt_0_yz_0[j] + 0.5 * fl1_fx * tt_y_z_0[j] + 2.0 * fl1_fz * ts_yy_yz_0[j] -
                            fl1_fz * fl1_fga * ts_0_yz_0[j];

            tt_yy_zz_0[j] = pa_y[j] * tt_y_zz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + 2.0 * fl1_fz * ts_yy_zz_0[j] - fl1_fz * fl1_fga * ts_0_zz_0[j];

            tt_yz_xx_0[j] = pa_y[j] * tt_z_xx_0[j] + 2.0 * fl1_fz * ts_yz_xx_0[j];

            tt_yz_xy_0[j] = pa_y[j] * tt_z_xy_0[j] + 0.5 * fl1_fx * tt_z_x_0[j] + 2.0 * fl1_fz * ts_yz_xy_0[j];

            tt_yz_xz_0[j] = pa_y[j] * tt_z_xz_0[j] + 2.0 * fl1_fz * ts_yz_xz_0[j];

            tt_yz_yy_0[j] = pa_y[j] * tt_z_yy_0[j] + fl1_fx * tt_z_y_0[j] + 2.0 * fl1_fz * ts_yz_yy_0[j];

            tt_yz_yz_0[j] = pa_y[j] * tt_z_yz_0[j] + 0.5 * fl1_fx * tt_z_z_0[j] + 2.0 * fl1_fz * ts_yz_yz_0[j];

            tt_yz_zz_0[j] = pa_y[j] * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_yz_zz_0[j];

            tt_zz_xx_0[j] = pa_z[j] * tt_z_xx_0[j] + 0.5 * fl1_fx * tt_0_xx_0[j] + 2.0 * fl1_fz * ts_zz_xx_0[j] - fl1_fz * fl1_fga * ts_0_xx_0[j];

            tt_zz_xy_0[j] = pa_z[j] * tt_z_xy_0[j] + 0.5 * fl1_fx * tt_0_xy_0[j] + 2.0 * fl1_fz * ts_zz_xy_0[j] - fl1_fz * fl1_fga * ts_0_xy_0[j];

            tt_zz_xz_0[j] = pa_z[j] * tt_z_xz_0[j] + 0.5 * fl1_fx * tt_0_xz_0[j] + 0.5 * fl1_fx * tt_z_x_0[j] + 2.0 * fl1_fz * ts_zz_xz_0[j] -
                            fl1_fz * fl1_fga * ts_0_xz_0[j];

            tt_zz_yy_0[j] = pa_z[j] * tt_z_yy_0[j] + 0.5 * fl1_fx * tt_0_yy_0[j] + 2.0 * fl1_fz * ts_zz_yy_0[j] - fl1_fz * fl1_fga * ts_0_yy_0[j];

            tt_zz_yz_0[j] = pa_z[j] * tt_z_yz_0[j] + 0.5 * fl1_fx * tt_0_yz_0[j] + 0.5 * fl1_fx * tt_z_y_0[j] + 2.0 * fl1_fz * ts_zz_yz_0[j] -
                            fl1_fz * fl1_fga * ts_0_yz_0[j];

            tt_zz_zz_0[j] = pa_z[j] * tt_z_zz_0[j] + 0.5 * fl1_fx * tt_0_zz_0[j] + fl1_fx * tt_z_z_0[j] + 2.0 * fl1_fz * ts_zz_zz_0[j] -
                            fl1_fz * fl1_fga * ts_0_zz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForDF(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    kinrecfunc::compKineticEnergyForDF_0_30(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDF_30_60(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compKineticEnergyForDF_0_30(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_t_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_2_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tt_x_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx);

        auto tt_x_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 1);

        auto tt_x_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 2);

        auto tt_x_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 3);

        auto tt_x_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 4);

        auto tt_x_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 5);

        auto tt_x_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 6);

        auto tt_x_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 7);

        auto tt_x_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 8);

        auto tt_x_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 9);

        auto tt_y_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 10);

        auto tt_y_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 11);

        auto tt_y_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 12);

        auto tt_y_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 13);

        auto tt_y_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 14);

        auto tt_y_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 15);

        auto tt_y_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 16);

        auto tt_y_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 17);

        auto tt_y_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 18);

        auto tt_y_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 19);

        auto tt_z_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 20);

        auto tt_z_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 21);

        auto tt_z_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 22);

        auto tt_z_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 23);

        auto tt_z_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 24);

        auto tt_z_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 25);

        auto tt_z_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 26);

        auto tt_z_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 27);

        auto tt_z_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 28);

        auto tt_z_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 29);

        auto tt_0_xxx_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx);

        auto tt_0_xxy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 1);

        auto tt_0_xxz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 2);

        auto tt_0_xyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 3);

        auto tt_0_xyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 4);

        auto tt_0_xzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 5);

        auto tt_0_yyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 6);

        auto tt_0_yyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 7);

        auto tt_0_yzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 8);

        auto tt_0_zzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 9);

        auto tt_x_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx);

        auto tt_x_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 1);

        auto tt_x_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 2);

        auto tt_x_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 3);

        auto tt_x_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 4);

        auto tt_x_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 5);

        auto tt_y_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 6);

        auto tt_y_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 7);

        auto tt_y_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 8);

        auto tt_y_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 9);

        auto tt_y_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 10);

        auto tt_y_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 11);

        auto tt_z_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 12);

        auto tt_z_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 13);

        auto tt_z_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 14);

        auto tt_z_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 15);

        auto tt_z_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 16);

        auto tt_z_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 17);

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

        // Batch of Integrals (0,30)

        #pragma omp simd aligned(fga, fx, fz, pa_x, ts_0_xxx_0, ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, \
                                     ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, ts_0_yzz_0, ts_0_zzz_0, \
                                     ts_xx_xxx_0, ts_xx_xxy_0, ts_xx_xxz_0, ts_xx_xyy_0, ts_xx_xyz_0, ts_xx_xzz_0, \
                                     ts_xx_yyy_0, ts_xx_yyz_0, ts_xx_yzz_0, ts_xx_zzz_0, ts_xy_xxx_0, ts_xy_xxy_0, \
                                     ts_xy_xxz_0, ts_xy_xyy_0, ts_xy_xyz_0, ts_xy_xzz_0, ts_xy_yyy_0, ts_xy_yyz_0, \
                                     ts_xy_yzz_0, ts_xy_zzz_0, ts_xz_xxx_0, ts_xz_xxy_0, ts_xz_xxz_0, ts_xz_xyy_0, \
                                     ts_xz_xyz_0, ts_xz_xzz_0, ts_xz_yyy_0, ts_xz_yyz_0, ts_xz_yzz_0, ts_xz_zzz_0, \
                                     tt_0_xxx_0, tt_0_xxy_0, tt_0_xxz_0, tt_0_xyy_0, tt_0_xyz_0, tt_0_xzz_0, tt_0_yyy_0, \
                                     tt_0_yyz_0, tt_0_yzz_0, tt_0_zzz_0, tt_x_xx_0, tt_x_xxx_0, tt_x_xxy_0, tt_x_xxz_0, \
                                     tt_x_xy_0, tt_x_xyy_0, tt_x_xyz_0, tt_x_xz_0, tt_x_xzz_0, tt_x_yy_0, tt_x_yyy_0, \
                                     tt_x_yyz_0, tt_x_yz_0, tt_x_yzz_0, tt_x_zz_0, tt_x_zzz_0, tt_xx_xxx_0, tt_xx_xxy_0, \
                                     tt_xx_xxz_0, tt_xx_xyy_0, tt_xx_xyz_0, tt_xx_xzz_0, tt_xx_yyy_0, tt_xx_yyz_0, \
                                     tt_xx_yzz_0, tt_xx_zzz_0, tt_xy_xxx_0, tt_xy_xxy_0, tt_xy_xxz_0, tt_xy_xyy_0, \
                                     tt_xy_xyz_0, tt_xy_xzz_0, tt_xy_yyy_0, tt_xy_yyz_0, tt_xy_yzz_0, tt_xy_zzz_0, \
                                     tt_xz_xxx_0, tt_xz_xxy_0, tt_xz_xxz_0, tt_xz_xyy_0, tt_xz_xyz_0, tt_xz_xzz_0, \
                                     tt_xz_yyy_0, tt_xz_yyz_0, tt_xz_yzz_0, tt_xz_zzz_0, tt_y_xx_0, tt_y_xxx_0, \
                                     tt_y_xxy_0, tt_y_xxz_0, tt_y_xy_0, tt_y_xyy_0, tt_y_xyz_0, tt_y_xz_0, tt_y_xzz_0, \
                                     tt_y_yy_0, tt_y_yyy_0, tt_y_yyz_0, tt_y_yz_0, tt_y_yzz_0, tt_y_zz_0, tt_y_zzz_0, \
                                     tt_z_xx_0, tt_z_xxx_0, tt_z_xxy_0, tt_z_xxz_0, tt_z_xy_0, tt_z_xyy_0, tt_z_xyz_0, \
                                     tt_z_xz_0, tt_z_xzz_0, tt_z_yy_0, tt_z_yyy_0, tt_z_yyz_0, tt_z_yz_0, tt_z_yzz_0, \
                                     tt_z_zz_0, tt_z_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xx_xxx_0[j] = pa_x[j] * tt_x_xxx_0[j] + 0.5 * fl1_fx * tt_0_xxx_0[j] + 1.5 * fl1_fx * tt_x_xx_0[j] + 2.0 * fl1_fz * ts_xx_xxx_0[j] -
                             fl1_fz * fl1_fga * ts_0_xxx_0[j];

            tt_xx_xxy_0[j] = pa_x[j] * tt_x_xxy_0[j] + 0.5 * fl1_fx * tt_0_xxy_0[j] + fl1_fx * tt_x_xy_0[j] + 2.0 * fl1_fz * ts_xx_xxy_0[j] -
                             fl1_fz * fl1_fga * ts_0_xxy_0[j];

            tt_xx_xxz_0[j] = pa_x[j] * tt_x_xxz_0[j] + 0.5 * fl1_fx * tt_0_xxz_0[j] + fl1_fx * tt_x_xz_0[j] + 2.0 * fl1_fz * ts_xx_xxz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xxz_0[j];

            tt_xx_xyy_0[j] = pa_x[j] * tt_x_xyy_0[j] + 0.5 * fl1_fx * tt_0_xyy_0[j] + 0.5 * fl1_fx * tt_x_yy_0[j] + 2.0 * fl1_fz * ts_xx_xyy_0[j] -
                             fl1_fz * fl1_fga * ts_0_xyy_0[j];

            tt_xx_xyz_0[j] = pa_x[j] * tt_x_xyz_0[j] + 0.5 * fl1_fx * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_x_yz_0[j] + 2.0 * fl1_fz * ts_xx_xyz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xyz_0[j];

            tt_xx_xzz_0[j] = pa_x[j] * tt_x_xzz_0[j] + 0.5 * fl1_fx * tt_0_xzz_0[j] + 0.5 * fl1_fx * tt_x_zz_0[j] + 2.0 * fl1_fz * ts_xx_xzz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xzz_0[j];

            tt_xx_yyy_0[j] =
                pa_x[j] * tt_x_yyy_0[j] + 0.5 * fl1_fx * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_xx_yyy_0[j] - fl1_fz * fl1_fga * ts_0_yyy_0[j];

            tt_xx_yyz_0[j] =
                pa_x[j] * tt_x_yyz_0[j] + 0.5 * fl1_fx * tt_0_yyz_0[j] + 2.0 * fl1_fz * ts_xx_yyz_0[j] - fl1_fz * fl1_fga * ts_0_yyz_0[j];

            tt_xx_yzz_0[j] =
                pa_x[j] * tt_x_yzz_0[j] + 0.5 * fl1_fx * tt_0_yzz_0[j] + 2.0 * fl1_fz * ts_xx_yzz_0[j] - fl1_fz * fl1_fga * ts_0_yzz_0[j];

            tt_xx_zzz_0[j] =
                pa_x[j] * tt_x_zzz_0[j] + 0.5 * fl1_fx * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_xx_zzz_0[j] - fl1_fz * fl1_fga * ts_0_zzz_0[j];

            tt_xy_xxx_0[j] = pa_x[j] * tt_y_xxx_0[j] + 1.5 * fl1_fx * tt_y_xx_0[j] + 2.0 * fl1_fz * ts_xy_xxx_0[j];

            tt_xy_xxy_0[j] = pa_x[j] * tt_y_xxy_0[j] + fl1_fx * tt_y_xy_0[j] + 2.0 * fl1_fz * ts_xy_xxy_0[j];

            tt_xy_xxz_0[j] = pa_x[j] * tt_y_xxz_0[j] + fl1_fx * tt_y_xz_0[j] + 2.0 * fl1_fz * ts_xy_xxz_0[j];

            tt_xy_xyy_0[j] = pa_x[j] * tt_y_xyy_0[j] + 0.5 * fl1_fx * tt_y_yy_0[j] + 2.0 * fl1_fz * ts_xy_xyy_0[j];

            tt_xy_xyz_0[j] = pa_x[j] * tt_y_xyz_0[j] + 0.5 * fl1_fx * tt_y_yz_0[j] + 2.0 * fl1_fz * ts_xy_xyz_0[j];

            tt_xy_xzz_0[j] = pa_x[j] * tt_y_xzz_0[j] + 0.5 * fl1_fx * tt_y_zz_0[j] + 2.0 * fl1_fz * ts_xy_xzz_0[j];

            tt_xy_yyy_0[j] = pa_x[j] * tt_y_yyy_0[j] + 2.0 * fl1_fz * ts_xy_yyy_0[j];

            tt_xy_yyz_0[j] = pa_x[j] * tt_y_yyz_0[j] + 2.0 * fl1_fz * ts_xy_yyz_0[j];

            tt_xy_yzz_0[j] = pa_x[j] * tt_y_yzz_0[j] + 2.0 * fl1_fz * ts_xy_yzz_0[j];

            tt_xy_zzz_0[j] = pa_x[j] * tt_y_zzz_0[j] + 2.0 * fl1_fz * ts_xy_zzz_0[j];

            tt_xz_xxx_0[j] = pa_x[j] * tt_z_xxx_0[j] + 1.5 * fl1_fx * tt_z_xx_0[j] + 2.0 * fl1_fz * ts_xz_xxx_0[j];

            tt_xz_xxy_0[j] = pa_x[j] * tt_z_xxy_0[j] + fl1_fx * tt_z_xy_0[j] + 2.0 * fl1_fz * ts_xz_xxy_0[j];

            tt_xz_xxz_0[j] = pa_x[j] * tt_z_xxz_0[j] + fl1_fx * tt_z_xz_0[j] + 2.0 * fl1_fz * ts_xz_xxz_0[j];

            tt_xz_xyy_0[j] = pa_x[j] * tt_z_xyy_0[j] + 0.5 * fl1_fx * tt_z_yy_0[j] + 2.0 * fl1_fz * ts_xz_xyy_0[j];

            tt_xz_xyz_0[j] = pa_x[j] * tt_z_xyz_0[j] + 0.5 * fl1_fx * tt_z_yz_0[j] + 2.0 * fl1_fz * ts_xz_xyz_0[j];

            tt_xz_xzz_0[j] = pa_x[j] * tt_z_xzz_0[j] + 0.5 * fl1_fx * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_xz_xzz_0[j];

            tt_xz_yyy_0[j] = pa_x[j] * tt_z_yyy_0[j] + 2.0 * fl1_fz * ts_xz_yyy_0[j];

            tt_xz_yyz_0[j] = pa_x[j] * tt_z_yyz_0[j] + 2.0 * fl1_fz * ts_xz_yyz_0[j];

            tt_xz_yzz_0[j] = pa_x[j] * tt_z_yzz_0[j] + 2.0 * fl1_fz * ts_xz_yzz_0[j];

            tt_xz_zzz_0[j] = pa_x[j] * tt_z_zzz_0[j] + 2.0 * fl1_fz * ts_xz_zzz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForDF_30_60(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_t_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_2_3_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tt_y_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 10);

        auto tt_y_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 11);

        auto tt_y_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 12);

        auto tt_y_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 13);

        auto tt_y_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 14);

        auto tt_y_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 15);

        auto tt_y_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 16);

        auto tt_y_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 17);

        auto tt_y_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 18);

        auto tt_y_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 19);

        auto tt_z_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 20);

        auto tt_z_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 21);

        auto tt_z_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 22);

        auto tt_z_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 23);

        auto tt_z_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 24);

        auto tt_z_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 25);

        auto tt_z_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 26);

        auto tt_z_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 27);

        auto tt_z_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 28);

        auto tt_z_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 29);

        auto tt_0_xxx_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx);

        auto tt_0_xxy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 1);

        auto tt_0_xxz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 2);

        auto tt_0_xyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 3);

        auto tt_0_xyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 4);

        auto tt_0_xzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 5);

        auto tt_0_yyy_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 6);

        auto tt_0_yyz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 7);

        auto tt_0_yzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 8);

        auto tt_0_zzz_0 = primBuffer.data(pidx_t_0_3_m0 + 10 * idx + 9);

        auto tt_y_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 6);

        auto tt_y_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 7);

        auto tt_y_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 8);

        auto tt_y_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 9);

        auto tt_y_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 10);

        auto tt_y_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 11);

        auto tt_z_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 12);

        auto tt_z_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 13);

        auto tt_z_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 14);

        auto tt_z_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 15);

        auto tt_z_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 16);

        auto tt_z_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 17);

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

        // Batch of Integrals (30,60)

        #pragma omp simd aligned(fga, fx, fz, pa_y, pa_z, ts_0_xxx_0, ts_0_xxy_0, ts_0_xxz_0, ts_0_xyy_0, \
                                     ts_0_xyz_0, ts_0_xzz_0, ts_0_yyy_0, ts_0_yyz_0, ts_0_yzz_0, ts_0_zzz_0, \
                                     ts_yy_xxx_0, ts_yy_xxy_0, ts_yy_xxz_0, ts_yy_xyy_0, ts_yy_xyz_0, ts_yy_xzz_0, \
                                     ts_yy_yyy_0, ts_yy_yyz_0, ts_yy_yzz_0, ts_yy_zzz_0, ts_yz_xxx_0, ts_yz_xxy_0, \
                                     ts_yz_xxz_0, ts_yz_xyy_0, ts_yz_xyz_0, ts_yz_xzz_0, ts_yz_yyy_0, ts_yz_yyz_0, \
                                     ts_yz_yzz_0, ts_yz_zzz_0, ts_zz_xxx_0, ts_zz_xxy_0, ts_zz_xxz_0, ts_zz_xyy_0, \
                                     ts_zz_xyz_0, ts_zz_xzz_0, ts_zz_yyy_0, ts_zz_yyz_0, ts_zz_yzz_0, ts_zz_zzz_0, \
                                     tt_0_xxx_0, tt_0_xxy_0, tt_0_xxz_0, tt_0_xyy_0, tt_0_xyz_0, tt_0_xzz_0, tt_0_yyy_0, \
                                     tt_0_yyz_0, tt_0_yzz_0, tt_0_zzz_0, tt_y_xx_0, tt_y_xxx_0, tt_y_xxy_0, tt_y_xxz_0, \
                                     tt_y_xy_0, tt_y_xyy_0, tt_y_xyz_0, tt_y_xz_0, tt_y_xzz_0, tt_y_yy_0, tt_y_yyy_0, \
                                     tt_y_yyz_0, tt_y_yz_0, tt_y_yzz_0, tt_y_zz_0, tt_y_zzz_0, tt_yy_xxx_0, tt_yy_xxy_0, \
                                     tt_yy_xxz_0, tt_yy_xyy_0, tt_yy_xyz_0, tt_yy_xzz_0, tt_yy_yyy_0, tt_yy_yyz_0, \
                                     tt_yy_yzz_0, tt_yy_zzz_0, tt_yz_xxx_0, tt_yz_xxy_0, tt_yz_xxz_0, tt_yz_xyy_0, \
                                     tt_yz_xyz_0, tt_yz_xzz_0, tt_yz_yyy_0, tt_yz_yyz_0, tt_yz_yzz_0, tt_yz_zzz_0, \
                                     tt_z_xx_0, tt_z_xxx_0, tt_z_xxy_0, tt_z_xxz_0, tt_z_xy_0, tt_z_xyy_0, tt_z_xyz_0, \
                                     tt_z_xz_0, tt_z_xzz_0, tt_z_yy_0, tt_z_yyy_0, tt_z_yyz_0, tt_z_yz_0, tt_z_yzz_0, \
                                     tt_z_zz_0, tt_z_zzz_0, tt_zz_xxx_0, tt_zz_xxy_0, tt_zz_xxz_0, tt_zz_xyy_0, \
                                     tt_zz_xyz_0, tt_zz_xzz_0, tt_zz_yyy_0, tt_zz_yyz_0, tt_zz_yzz_0, tt_zz_zzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_yy_xxx_0[j] =
                pa_y[j] * tt_y_xxx_0[j] + 0.5 * fl1_fx * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_yy_xxx_0[j] - fl1_fz * fl1_fga * ts_0_xxx_0[j];

            tt_yy_xxy_0[j] = pa_y[j] * tt_y_xxy_0[j] + 0.5 * fl1_fx * tt_0_xxy_0[j] + 0.5 * fl1_fx * tt_y_xx_0[j] + 2.0 * fl1_fz * ts_yy_xxy_0[j] -
                             fl1_fz * fl1_fga * ts_0_xxy_0[j];

            tt_yy_xxz_0[j] =
                pa_y[j] * tt_y_xxz_0[j] + 0.5 * fl1_fx * tt_0_xxz_0[j] + 2.0 * fl1_fz * ts_yy_xxz_0[j] - fl1_fz * fl1_fga * ts_0_xxz_0[j];

            tt_yy_xyy_0[j] = pa_y[j] * tt_y_xyy_0[j] + 0.5 * fl1_fx * tt_0_xyy_0[j] + fl1_fx * tt_y_xy_0[j] + 2.0 * fl1_fz * ts_yy_xyy_0[j] -
                             fl1_fz * fl1_fga * ts_0_xyy_0[j];

            tt_yy_xyz_0[j] = pa_y[j] * tt_y_xyz_0[j] + 0.5 * fl1_fx * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_y_xz_0[j] + 2.0 * fl1_fz * ts_yy_xyz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xyz_0[j];

            tt_yy_xzz_0[j] =
                pa_y[j] * tt_y_xzz_0[j] + 0.5 * fl1_fx * tt_0_xzz_0[j] + 2.0 * fl1_fz * ts_yy_xzz_0[j] - fl1_fz * fl1_fga * ts_0_xzz_0[j];

            tt_yy_yyy_0[j] = pa_y[j] * tt_y_yyy_0[j] + 0.5 * fl1_fx * tt_0_yyy_0[j] + 1.5 * fl1_fx * tt_y_yy_0[j] + 2.0 * fl1_fz * ts_yy_yyy_0[j] -
                             fl1_fz * fl1_fga * ts_0_yyy_0[j];

            tt_yy_yyz_0[j] = pa_y[j] * tt_y_yyz_0[j] + 0.5 * fl1_fx * tt_0_yyz_0[j] + fl1_fx * tt_y_yz_0[j] + 2.0 * fl1_fz * ts_yy_yyz_0[j] -
                             fl1_fz * fl1_fga * ts_0_yyz_0[j];

            tt_yy_yzz_0[j] = pa_y[j] * tt_y_yzz_0[j] + 0.5 * fl1_fx * tt_0_yzz_0[j] + 0.5 * fl1_fx * tt_y_zz_0[j] + 2.0 * fl1_fz * ts_yy_yzz_0[j] -
                             fl1_fz * fl1_fga * ts_0_yzz_0[j];

            tt_yy_zzz_0[j] =
                pa_y[j] * tt_y_zzz_0[j] + 0.5 * fl1_fx * tt_0_zzz_0[j] + 2.0 * fl1_fz * ts_yy_zzz_0[j] - fl1_fz * fl1_fga * ts_0_zzz_0[j];

            tt_yz_xxx_0[j] = pa_y[j] * tt_z_xxx_0[j] + 2.0 * fl1_fz * ts_yz_xxx_0[j];

            tt_yz_xxy_0[j] = pa_y[j] * tt_z_xxy_0[j] + 0.5 * fl1_fx * tt_z_xx_0[j] + 2.0 * fl1_fz * ts_yz_xxy_0[j];

            tt_yz_xxz_0[j] = pa_y[j] * tt_z_xxz_0[j] + 2.0 * fl1_fz * ts_yz_xxz_0[j];

            tt_yz_xyy_0[j] = pa_y[j] * tt_z_xyy_0[j] + fl1_fx * tt_z_xy_0[j] + 2.0 * fl1_fz * ts_yz_xyy_0[j];

            tt_yz_xyz_0[j] = pa_y[j] * tt_z_xyz_0[j] + 0.5 * fl1_fx * tt_z_xz_0[j] + 2.0 * fl1_fz * ts_yz_xyz_0[j];

            tt_yz_xzz_0[j] = pa_y[j] * tt_z_xzz_0[j] + 2.0 * fl1_fz * ts_yz_xzz_0[j];

            tt_yz_yyy_0[j] = pa_y[j] * tt_z_yyy_0[j] + 1.5 * fl1_fx * tt_z_yy_0[j] + 2.0 * fl1_fz * ts_yz_yyy_0[j];

            tt_yz_yyz_0[j] = pa_y[j] * tt_z_yyz_0[j] + fl1_fx * tt_z_yz_0[j] + 2.0 * fl1_fz * ts_yz_yyz_0[j];

            tt_yz_yzz_0[j] = pa_y[j] * tt_z_yzz_0[j] + 0.5 * fl1_fx * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_yz_yzz_0[j];

            tt_yz_zzz_0[j] = pa_y[j] * tt_z_zzz_0[j] + 2.0 * fl1_fz * ts_yz_zzz_0[j];

            tt_zz_xxx_0[j] =
                pa_z[j] * tt_z_xxx_0[j] + 0.5 * fl1_fx * tt_0_xxx_0[j] + 2.0 * fl1_fz * ts_zz_xxx_0[j] - fl1_fz * fl1_fga * ts_0_xxx_0[j];

            tt_zz_xxy_0[j] =
                pa_z[j] * tt_z_xxy_0[j] + 0.5 * fl1_fx * tt_0_xxy_0[j] + 2.0 * fl1_fz * ts_zz_xxy_0[j] - fl1_fz * fl1_fga * ts_0_xxy_0[j];

            tt_zz_xxz_0[j] = pa_z[j] * tt_z_xxz_0[j] + 0.5 * fl1_fx * tt_0_xxz_0[j] + 0.5 * fl1_fx * tt_z_xx_0[j] + 2.0 * fl1_fz * ts_zz_xxz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xxz_0[j];

            tt_zz_xyy_0[j] =
                pa_z[j] * tt_z_xyy_0[j] + 0.5 * fl1_fx * tt_0_xyy_0[j] + 2.0 * fl1_fz * ts_zz_xyy_0[j] - fl1_fz * fl1_fga * ts_0_xyy_0[j];

            tt_zz_xyz_0[j] = pa_z[j] * tt_z_xyz_0[j] + 0.5 * fl1_fx * tt_0_xyz_0[j] + 0.5 * fl1_fx * tt_z_xy_0[j] + 2.0 * fl1_fz * ts_zz_xyz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xyz_0[j];

            tt_zz_xzz_0[j] = pa_z[j] * tt_z_xzz_0[j] + 0.5 * fl1_fx * tt_0_xzz_0[j] + fl1_fx * tt_z_xz_0[j] + 2.0 * fl1_fz * ts_zz_xzz_0[j] -
                             fl1_fz * fl1_fga * ts_0_xzz_0[j];

            tt_zz_yyy_0[j] =
                pa_z[j] * tt_z_yyy_0[j] + 0.5 * fl1_fx * tt_0_yyy_0[j] + 2.0 * fl1_fz * ts_zz_yyy_0[j] - fl1_fz * fl1_fga * ts_0_yyy_0[j];

            tt_zz_yyz_0[j] = pa_z[j] * tt_z_yyz_0[j] + 0.5 * fl1_fx * tt_0_yyz_0[j] + 0.5 * fl1_fx * tt_z_yy_0[j] + 2.0 * fl1_fz * ts_zz_yyz_0[j] -
                             fl1_fz * fl1_fga * ts_0_yyz_0[j];

            tt_zz_yzz_0[j] = pa_z[j] * tt_z_yzz_0[j] + 0.5 * fl1_fx * tt_0_yzz_0[j] + fl1_fx * tt_z_yz_0[j] + 2.0 * fl1_fz * ts_zz_yzz_0[j] -
                             fl1_fz * fl1_fga * ts_0_yzz_0[j];

            tt_zz_zzz_0[j] = pa_z[j] * tt_z_zzz_0[j] + 0.5 * fl1_fx * tt_0_zzz_0[j] + 1.5 * fl1_fx * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_zz_zzz_0[j] -
                             fl1_fz * fl1_fga * ts_0_zzz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForFD(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    kinrecfunc::compKineticEnergyForFD_0_30(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFD_30_60(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compKineticEnergyForFD_0_30(CMemBlock2D<double>&       primBuffer,
                            const CRecursionMap&       recursionMap,
                            const CMemBlock2D<double>& osFactors,
                            const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_t_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_3_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tt_xx_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx);

        auto tt_xx_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 1);

        auto tt_xx_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 2);

        auto tt_xx_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 3);

        auto tt_xx_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 4);

        auto tt_xx_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 5);

        auto tt_xy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 6);

        auto tt_xy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 7);

        auto tt_xy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 8);

        auto tt_xy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 9);

        auto tt_xy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 10);

        auto tt_xy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 11);

        auto tt_xz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 12);

        auto tt_xz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 13);

        auto tt_xz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 14);

        auto tt_xz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 15);

        auto tt_xz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 16);

        auto tt_xz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 17);

        auto tt_yy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 18);

        auto tt_yy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 19);

        auto tt_yy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 20);

        auto tt_yy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 21);

        auto tt_yy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 22);

        auto tt_yy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 23);

        auto tt_yz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 24);

        auto tt_yz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 25);

        auto tt_yz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 26);

        auto tt_yz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 27);

        auto tt_yz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 28);

        auto tt_yz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 29);

        auto tt_x_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx);

        auto tt_x_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 1);

        auto tt_x_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 2);

        auto tt_x_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 3);

        auto tt_x_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 4);

        auto tt_x_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 5);

        auto tt_y_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 6);

        auto tt_y_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 7);

        auto tt_y_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 8);

        auto tt_y_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 9);

        auto tt_y_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 10);

        auto tt_y_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 11);

        auto tt_z_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 12);

        auto tt_z_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 13);

        auto tt_z_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 14);

        auto tt_z_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 15);

        auto tt_z_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 16);

        auto tt_z_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 17);

        auto tt_xx_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx);

        auto tt_xx_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 1);

        auto tt_xx_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 2);

        auto tt_xy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 3);

        auto tt_xy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 4);

        auto tt_xy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 5);

        auto tt_xz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 6);

        auto tt_xz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 7);

        auto tt_xz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 8);

        auto tt_yy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 9);

        auto tt_yy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 10);

        auto tt_yy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 11);

        auto tt_yz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 12);

        auto tt_yz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 13);

        auto tt_yz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 14);

        auto ts_xxx_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx);

        auto ts_xxx_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 1);

        auto ts_xxx_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 2);

        auto ts_xxx_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 3);

        auto ts_xxx_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 4);

        auto ts_xxx_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 5);

        auto ts_xxy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 6);

        auto ts_xxy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 7);

        auto ts_xxy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 8);

        auto ts_xxy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 9);

        auto ts_xxy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 10);

        auto ts_xxy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 11);

        auto ts_xxz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 12);

        auto ts_xxz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 13);

        auto ts_xxz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 14);

        auto ts_xxz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 15);

        auto ts_xxz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 16);

        auto ts_xxz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 17);

        auto ts_xyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 18);

        auto ts_xyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 19);

        auto ts_xyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 20);

        auto ts_xyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 21);

        auto ts_xyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 22);

        auto ts_xyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 23);

        auto ts_xyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 24);

        auto ts_xyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 25);

        auto ts_xyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 26);

        auto ts_xyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 27);

        auto ts_xyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 28);

        auto ts_xyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 29);

        auto ts_x_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx);

        auto ts_x_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 1);

        auto ts_x_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 2);

        auto ts_x_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 3);

        auto ts_x_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 4);

        auto ts_x_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 5);

        auto ts_y_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 6);

        auto ts_y_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 7);

        auto ts_y_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 8);

        auto ts_y_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 9);

        auto ts_y_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 10);

        auto ts_y_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 11);

        auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12);

        auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13);

        auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14);

        auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15);

        auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16);

        auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17);

        // set up pointers to integrals

        auto tt_xxx_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx);

        auto tt_xxx_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 1);

        auto tt_xxx_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 2);

        auto tt_xxx_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 3);

        auto tt_xxx_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 4);

        auto tt_xxx_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 5);

        auto tt_xxy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 6);

        auto tt_xxy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 7);

        auto tt_xxy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 8);

        auto tt_xxy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 9);

        auto tt_xxy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 10);

        auto tt_xxy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 11);

        auto tt_xxz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 12);

        auto tt_xxz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 13);

        auto tt_xxz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 14);

        auto tt_xxz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 15);

        auto tt_xxz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 16);

        auto tt_xxz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 17);

        auto tt_xyy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 18);

        auto tt_xyy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 19);

        auto tt_xyy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 20);

        auto tt_xyy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 21);

        auto tt_xyy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 22);

        auto tt_xyy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 23);

        auto tt_xyz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 24);

        auto tt_xyz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 25);

        auto tt_xyz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 26);

        auto tt_xyz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 27);

        auto tt_xyz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 28);

        auto tt_xyz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 29);

        // Batch of Integrals (0,30)

        #pragma omp simd aligned(fga, fx, fz, pa_x, ts_x_xx_0, ts_x_xy_0, ts_x_xz_0, ts_x_yy_0, ts_x_yz_0, \
                                     ts_x_zz_0, ts_xxx_xx_0, ts_xxx_xy_0, ts_xxx_xz_0, ts_xxx_yy_0, ts_xxx_yz_0, \
                                     ts_xxx_zz_0, ts_xxy_xx_0, ts_xxy_xy_0, ts_xxy_xz_0, ts_xxy_yy_0, ts_xxy_yz_0, \
                                     ts_xxy_zz_0, ts_xxz_xx_0, ts_xxz_xy_0, ts_xxz_xz_0, ts_xxz_yy_0, ts_xxz_yz_0, \
                                     ts_xxz_zz_0, ts_xyy_xx_0, ts_xyy_xy_0, ts_xyy_xz_0, ts_xyy_yy_0, ts_xyy_yz_0, \
                                     ts_xyy_zz_0, ts_xyz_xx_0, ts_xyz_xy_0, ts_xyz_xz_0, ts_xyz_yy_0, ts_xyz_yz_0, \
                                     ts_xyz_zz_0, ts_y_xx_0, ts_y_xy_0, ts_y_xz_0, ts_y_yy_0, ts_y_yz_0, ts_y_zz_0, \
                                     ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, ts_z_zz_0, tt_x_xx_0, \
                                     tt_x_xy_0, tt_x_xz_0, tt_x_yy_0, tt_x_yz_0, tt_x_zz_0, tt_xx_x_0, tt_xx_xx_0, \
                                     tt_xx_xy_0, tt_xx_xz_0, tt_xx_y_0, tt_xx_yy_0, tt_xx_yz_0, tt_xx_z_0, tt_xx_zz_0, \
                                     tt_xxx_xx_0, tt_xxx_xy_0, tt_xxx_xz_0, tt_xxx_yy_0, tt_xxx_yz_0, tt_xxx_zz_0, \
                                     tt_xxy_xx_0, tt_xxy_xy_0, tt_xxy_xz_0, tt_xxy_yy_0, tt_xxy_yz_0, tt_xxy_zz_0, \
                                     tt_xxz_xx_0, tt_xxz_xy_0, tt_xxz_xz_0, tt_xxz_yy_0, tt_xxz_yz_0, tt_xxz_zz_0, \
                                     tt_xy_x_0, tt_xy_xx_0, tt_xy_xy_0, tt_xy_xz_0, tt_xy_y_0, tt_xy_yy_0, tt_xy_yz_0, \
                                     tt_xy_z_0, tt_xy_zz_0, tt_xyy_xx_0, tt_xyy_xy_0, tt_xyy_xz_0, tt_xyy_yy_0, \
                                     tt_xyy_yz_0, tt_xyy_zz_0, tt_xyz_xx_0, tt_xyz_xy_0, tt_xyz_xz_0, tt_xyz_yy_0, \
                                     tt_xyz_yz_0, tt_xyz_zz_0, tt_xz_x_0, tt_xz_xx_0, tt_xz_xy_0, tt_xz_xz_0, tt_xz_y_0, \
                                     tt_xz_yy_0, tt_xz_yz_0, tt_xz_z_0, tt_xz_zz_0, tt_y_xx_0, tt_y_xy_0, tt_y_xz_0, \
                                     tt_y_yy_0, tt_y_yz_0, tt_y_zz_0, tt_yy_x_0, tt_yy_xx_0, tt_yy_xy_0, tt_yy_xz_0, \
                                     tt_yy_y_0, tt_yy_yy_0, tt_yy_yz_0, tt_yy_z_0, tt_yy_zz_0, tt_yz_x_0, tt_yz_xx_0, \
                                     tt_yz_xy_0, tt_yz_xz_0, tt_yz_y_0, tt_yz_yy_0, tt_yz_yz_0, tt_yz_z_0, tt_yz_zz_0, \
                                     tt_z_xx_0, tt_z_xy_0, tt_z_xz_0, tt_z_yy_0, tt_z_yz_0, tt_z_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xxx_xx_0[j] = pa_x[j] * tt_xx_xx_0[j] + fl1_fx * tt_x_xx_0[j] + fl1_fx * tt_xx_x_0[j] + 2.0 * fl1_fz * ts_xxx_xx_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_x_xx_0[j];

            tt_xxx_xy_0[j] = pa_x[j] * tt_xx_xy_0[j] + fl1_fx * tt_x_xy_0[j] + 0.5 * fl1_fx * tt_xx_y_0[j] + 2.0 * fl1_fz * ts_xxx_xy_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_x_xy_0[j];

            tt_xxx_xz_0[j] = pa_x[j] * tt_xx_xz_0[j] + fl1_fx * tt_x_xz_0[j] + 0.5 * fl1_fx * tt_xx_z_0[j] + 2.0 * fl1_fz * ts_xxx_xz_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_x_xz_0[j];

            tt_xxx_yy_0[j] = pa_x[j] * tt_xx_yy_0[j] + fl1_fx * tt_x_yy_0[j] + 2.0 * fl1_fz * ts_xxx_yy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_yy_0[j];

            tt_xxx_yz_0[j] = pa_x[j] * tt_xx_yz_0[j] + fl1_fx * tt_x_yz_0[j] + 2.0 * fl1_fz * ts_xxx_yz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_yz_0[j];

            tt_xxx_zz_0[j] = pa_x[j] * tt_xx_zz_0[j] + fl1_fx * tt_x_zz_0[j] + 2.0 * fl1_fz * ts_xxx_zz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_x_zz_0[j];

            tt_xxy_xx_0[j] = pa_x[j] * tt_xy_xx_0[j] + 0.5 * fl1_fx * tt_y_xx_0[j] + fl1_fx * tt_xy_x_0[j] + 2.0 * fl1_fz * ts_xxy_xx_0[j] -
                             fl1_fz * fl1_fga * ts_y_xx_0[j];

            tt_xxy_xy_0[j] = pa_x[j] * tt_xy_xy_0[j] + 0.5 * fl1_fx * tt_y_xy_0[j] + 0.5 * fl1_fx * tt_xy_y_0[j] + 2.0 * fl1_fz * ts_xxy_xy_0[j] -
                             fl1_fz * fl1_fga * ts_y_xy_0[j];

            tt_xxy_xz_0[j] = pa_x[j] * tt_xy_xz_0[j] + 0.5 * fl1_fx * tt_y_xz_0[j] + 0.5 * fl1_fx * tt_xy_z_0[j] + 2.0 * fl1_fz * ts_xxy_xz_0[j] -
                             fl1_fz * fl1_fga * ts_y_xz_0[j];

            tt_xxy_yy_0[j] = pa_x[j] * tt_xy_yy_0[j] + 0.5 * fl1_fx * tt_y_yy_0[j] + 2.0 * fl1_fz * ts_xxy_yy_0[j] - fl1_fz * fl1_fga * ts_y_yy_0[j];

            tt_xxy_yz_0[j] = pa_x[j] * tt_xy_yz_0[j] + 0.5 * fl1_fx * tt_y_yz_0[j] + 2.0 * fl1_fz * ts_xxy_yz_0[j] - fl1_fz * fl1_fga * ts_y_yz_0[j];

            tt_xxy_zz_0[j] = pa_x[j] * tt_xy_zz_0[j] + 0.5 * fl1_fx * tt_y_zz_0[j] + 2.0 * fl1_fz * ts_xxy_zz_0[j] - fl1_fz * fl1_fga * ts_y_zz_0[j];

            tt_xxz_xx_0[j] = pa_x[j] * tt_xz_xx_0[j] + 0.5 * fl1_fx * tt_z_xx_0[j] + fl1_fx * tt_xz_x_0[j] + 2.0 * fl1_fz * ts_xxz_xx_0[j] -
                             fl1_fz * fl1_fga * ts_z_xx_0[j];

            tt_xxz_xy_0[j] = pa_x[j] * tt_xz_xy_0[j] + 0.5 * fl1_fx * tt_z_xy_0[j] + 0.5 * fl1_fx * tt_xz_y_0[j] + 2.0 * fl1_fz * ts_xxz_xy_0[j] -
                             fl1_fz * fl1_fga * ts_z_xy_0[j];

            tt_xxz_xz_0[j] = pa_x[j] * tt_xz_xz_0[j] + 0.5 * fl1_fx * tt_z_xz_0[j] + 0.5 * fl1_fx * tt_xz_z_0[j] + 2.0 * fl1_fz * ts_xxz_xz_0[j] -
                             fl1_fz * fl1_fga * ts_z_xz_0[j];

            tt_xxz_yy_0[j] = pa_x[j] * tt_xz_yy_0[j] + 0.5 * fl1_fx * tt_z_yy_0[j] + 2.0 * fl1_fz * ts_xxz_yy_0[j] - fl1_fz * fl1_fga * ts_z_yy_0[j];

            tt_xxz_yz_0[j] = pa_x[j] * tt_xz_yz_0[j] + 0.5 * fl1_fx * tt_z_yz_0[j] + 2.0 * fl1_fz * ts_xxz_yz_0[j] - fl1_fz * fl1_fga * ts_z_yz_0[j];

            tt_xxz_zz_0[j] = pa_x[j] * tt_xz_zz_0[j] + 0.5 * fl1_fx * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_xxz_zz_0[j] - fl1_fz * fl1_fga * ts_z_zz_0[j];

            tt_xyy_xx_0[j] = pa_x[j] * tt_yy_xx_0[j] + fl1_fx * tt_yy_x_0[j] + 2.0 * fl1_fz * ts_xyy_xx_0[j];

            tt_xyy_xy_0[j] = pa_x[j] * tt_yy_xy_0[j] + 0.5 * fl1_fx * tt_yy_y_0[j] + 2.0 * fl1_fz * ts_xyy_xy_0[j];

            tt_xyy_xz_0[j] = pa_x[j] * tt_yy_xz_0[j] + 0.5 * fl1_fx * tt_yy_z_0[j] + 2.0 * fl1_fz * ts_xyy_xz_0[j];

            tt_xyy_yy_0[j] = pa_x[j] * tt_yy_yy_0[j] + 2.0 * fl1_fz * ts_xyy_yy_0[j];

            tt_xyy_yz_0[j] = pa_x[j] * tt_yy_yz_0[j] + 2.0 * fl1_fz * ts_xyy_yz_0[j];

            tt_xyy_zz_0[j] = pa_x[j] * tt_yy_zz_0[j] + 2.0 * fl1_fz * ts_xyy_zz_0[j];

            tt_xyz_xx_0[j] = pa_x[j] * tt_yz_xx_0[j] + fl1_fx * tt_yz_x_0[j] + 2.0 * fl1_fz * ts_xyz_xx_0[j];

            tt_xyz_xy_0[j] = pa_x[j] * tt_yz_xy_0[j] + 0.5 * fl1_fx * tt_yz_y_0[j] + 2.0 * fl1_fz * ts_xyz_xy_0[j];

            tt_xyz_xz_0[j] = pa_x[j] * tt_yz_xz_0[j] + 0.5 * fl1_fx * tt_yz_z_0[j] + 2.0 * fl1_fz * ts_xyz_xz_0[j];

            tt_xyz_yy_0[j] = pa_x[j] * tt_yz_yy_0[j] + 2.0 * fl1_fz * ts_xyz_yy_0[j];

            tt_xyz_yz_0[j] = pa_x[j] * tt_yz_yz_0[j] + 2.0 * fl1_fz * ts_xyz_yz_0[j];

            tt_xyz_zz_0[j] = pa_x[j] * tt_yz_zz_0[j] + 2.0 * fl1_fz * ts_xyz_zz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForFD_30_60(CMemBlock2D<double>&       primBuffer,
                             const CRecursionMap&       recursionMap,
                             const CMemBlock2D<double>& osFactors,
                             const CMemBlock2D<double>& paDistances,
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

    // set up index of integral

    auto pidx_t_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_3_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_1_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {1, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tt_yy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 18);

        auto tt_yy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 19);

        auto tt_yy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 20);

        auto tt_yy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 21);

        auto tt_yy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 22);

        auto tt_yy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 23);

        auto tt_yz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 24);

        auto tt_yz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 25);

        auto tt_yz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 26);

        auto tt_yz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 27);

        auto tt_yz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 28);

        auto tt_yz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 29);

        auto tt_zz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 30);

        auto tt_zz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 31);

        auto tt_zz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 32);

        auto tt_zz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 33);

        auto tt_zz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 34);

        auto tt_zz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 35);

        auto tt_y_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 6);

        auto tt_y_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 7);

        auto tt_y_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 8);

        auto tt_y_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 9);

        auto tt_y_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 10);

        auto tt_y_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 11);

        auto tt_z_xx_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 12);

        auto tt_z_xy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 13);

        auto tt_z_xz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 14);

        auto tt_z_yy_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 15);

        auto tt_z_yz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 16);

        auto tt_z_zz_0 = primBuffer.data(pidx_t_1_2_m0 + 18 * idx + 17);

        auto tt_yy_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 9);

        auto tt_yy_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 10);

        auto tt_yy_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 11);

        auto tt_yz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 12);

        auto tt_yz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 13);

        auto tt_yz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 14);

        auto tt_zz_x_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 15);

        auto tt_zz_y_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 16);

        auto tt_zz_z_0 = primBuffer.data(pidx_t_2_1_m0 + 18 * idx + 17);

        auto ts_xzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 30);

        auto ts_xzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 31);

        auto ts_xzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 32);

        auto ts_xzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 33);

        auto ts_xzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 34);

        auto ts_xzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 35);

        auto ts_yyy_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 36);

        auto ts_yyy_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 37);

        auto ts_yyy_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 38);

        auto ts_yyy_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 39);

        auto ts_yyy_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 40);

        auto ts_yyy_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 41);

        auto ts_yyz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 42);

        auto ts_yyz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 43);

        auto ts_yyz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 44);

        auto ts_yyz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 45);

        auto ts_yyz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 46);

        auto ts_yyz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 47);

        auto ts_yzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 48);

        auto ts_yzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 49);

        auto ts_yzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 50);

        auto ts_yzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 51);

        auto ts_yzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 52);

        auto ts_yzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 53);

        auto ts_zzz_xx_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 54);

        auto ts_zzz_xy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 55);

        auto ts_zzz_xz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 56);

        auto ts_zzz_yy_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 57);

        auto ts_zzz_yz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 58);

        auto ts_zzz_zz_0 = primBuffer.data(pidx_s_3_2_m0 + 60 * idx + 59);

        auto ts_y_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 6);

        auto ts_y_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 7);

        auto ts_y_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 8);

        auto ts_y_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 9);

        auto ts_y_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 10);

        auto ts_y_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 11);

        auto ts_z_xx_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 12);

        auto ts_z_xy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 13);

        auto ts_z_xz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 14);

        auto ts_z_yy_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 15);

        auto ts_z_yz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 16);

        auto ts_z_zz_0 = primBuffer.data(pidx_s_1_2_m0 + 18 * idx + 17);

        // set up pointers to integrals

        auto tt_xzz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 30);

        auto tt_xzz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 31);

        auto tt_xzz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 32);

        auto tt_xzz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 33);

        auto tt_xzz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 34);

        auto tt_xzz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 35);

        auto tt_yyy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 36);

        auto tt_yyy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 37);

        auto tt_yyy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 38);

        auto tt_yyy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 39);

        auto tt_yyy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 40);

        auto tt_yyy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 41);

        auto tt_yyz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 42);

        auto tt_yyz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 43);

        auto tt_yyz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 44);

        auto tt_yyz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 45);

        auto tt_yyz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 46);

        auto tt_yyz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 47);

        auto tt_yzz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 48);

        auto tt_yzz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 49);

        auto tt_yzz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 50);

        auto tt_yzz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 51);

        auto tt_yzz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 52);

        auto tt_yzz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 53);

        auto tt_zzz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 54);

        auto tt_zzz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 55);

        auto tt_zzz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 56);

        auto tt_zzz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 57);

        auto tt_zzz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 58);

        auto tt_zzz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 59);

        // Batch of Integrals (30,60)

        #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_xzz_xx_0, ts_xzz_xy_0, ts_xzz_xz_0, \
                                     ts_xzz_yy_0, ts_xzz_yz_0, ts_xzz_zz_0, ts_y_xx_0, ts_y_xy_0, ts_y_xz_0, ts_y_yy_0, \
                                     ts_y_yz_0, ts_y_zz_0, ts_yyy_xx_0, ts_yyy_xy_0, ts_yyy_xz_0, ts_yyy_yy_0, \
                                     ts_yyy_yz_0, ts_yyy_zz_0, ts_yyz_xx_0, ts_yyz_xy_0, ts_yyz_xz_0, ts_yyz_yy_0, \
                                     ts_yyz_yz_0, ts_yyz_zz_0, ts_yzz_xx_0, ts_yzz_xy_0, ts_yzz_xz_0, ts_yzz_yy_0, \
                                     ts_yzz_yz_0, ts_yzz_zz_0, ts_z_xx_0, ts_z_xy_0, ts_z_xz_0, ts_z_yy_0, ts_z_yz_0, \
                                     ts_z_zz_0, ts_zzz_xx_0, ts_zzz_xy_0, ts_zzz_xz_0, ts_zzz_yy_0, ts_zzz_yz_0, \
                                     ts_zzz_zz_0, tt_xzz_xx_0, tt_xzz_xy_0, tt_xzz_xz_0, tt_xzz_yy_0, tt_xzz_yz_0, \
                                     tt_xzz_zz_0, tt_y_xx_0, tt_y_xy_0, tt_y_xz_0, tt_y_yy_0, tt_y_yz_0, tt_y_zz_0, \
                                     tt_yy_x_0, tt_yy_xx_0, tt_yy_xy_0, tt_yy_xz_0, tt_yy_y_0, tt_yy_yy_0, tt_yy_yz_0, \
                                     tt_yy_z_0, tt_yy_zz_0, tt_yyy_xx_0, tt_yyy_xy_0, tt_yyy_xz_0, tt_yyy_yy_0, \
                                     tt_yyy_yz_0, tt_yyy_zz_0, tt_yyz_xx_0, tt_yyz_xy_0, tt_yyz_xz_0, tt_yyz_yy_0, \
                                     tt_yyz_yz_0, tt_yyz_zz_0, tt_yz_x_0, tt_yz_xx_0, tt_yz_xy_0, tt_yz_xz_0, tt_yz_y_0, \
                                     tt_yz_yy_0, tt_yz_yz_0, tt_yz_z_0, tt_yz_zz_0, tt_yzz_xx_0, tt_yzz_xy_0, \
                                     tt_yzz_xz_0, tt_yzz_yy_0, tt_yzz_yz_0, tt_yzz_zz_0, tt_z_xx_0, tt_z_xy_0, tt_z_xz_0, \
                                     tt_z_yy_0, tt_z_yz_0, tt_z_zz_0, tt_zz_x_0, tt_zz_xx_0, tt_zz_xy_0, tt_zz_xz_0, \
                                     tt_zz_y_0, tt_zz_yy_0, tt_zz_yz_0, tt_zz_z_0, tt_zz_zz_0, tt_zzz_xx_0, \
                                     tt_zzz_xy_0, tt_zzz_xz_0, tt_zzz_yy_0, tt_zzz_yz_0, tt_zzz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xzz_xx_0[j] = pa_x[j] * tt_zz_xx_0[j] + fl1_fx * tt_zz_x_0[j] + 2.0 * fl1_fz * ts_xzz_xx_0[j];

            tt_xzz_xy_0[j] = pa_x[j] * tt_zz_xy_0[j] + 0.5 * fl1_fx * tt_zz_y_0[j] + 2.0 * fl1_fz * ts_xzz_xy_0[j];

            tt_xzz_xz_0[j] = pa_x[j] * tt_zz_xz_0[j] + 0.5 * fl1_fx * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_xzz_xz_0[j];

            tt_xzz_yy_0[j] = pa_x[j] * tt_zz_yy_0[j] + 2.0 * fl1_fz * ts_xzz_yy_0[j];

            tt_xzz_yz_0[j] = pa_x[j] * tt_zz_yz_0[j] + 2.0 * fl1_fz * ts_xzz_yz_0[j];

            tt_xzz_zz_0[j] = pa_x[j] * tt_zz_zz_0[j] + 2.0 * fl1_fz * ts_xzz_zz_0[j];

            tt_yyy_xx_0[j] = pa_y[j] * tt_yy_xx_0[j] + fl1_fx * tt_y_xx_0[j] + 2.0 * fl1_fz * ts_yyy_xx_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_xx_0[j];

            tt_yyy_xy_0[j] = pa_y[j] * tt_yy_xy_0[j] + fl1_fx * tt_y_xy_0[j] + 0.5 * fl1_fx * tt_yy_x_0[j] + 2.0 * fl1_fz * ts_yyy_xy_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_y_xy_0[j];

            tt_yyy_xz_0[j] = pa_y[j] * tt_yy_xz_0[j] + fl1_fx * tt_y_xz_0[j] + 2.0 * fl1_fz * ts_yyy_xz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_xz_0[j];

            tt_yyy_yy_0[j] = pa_y[j] * tt_yy_yy_0[j] + fl1_fx * tt_y_yy_0[j] + fl1_fx * tt_yy_y_0[j] + 2.0 * fl1_fz * ts_yyy_yy_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_y_yy_0[j];

            tt_yyy_yz_0[j] = pa_y[j] * tt_yy_yz_0[j] + fl1_fx * tt_y_yz_0[j] + 0.5 * fl1_fx * tt_yy_z_0[j] + 2.0 * fl1_fz * ts_yyy_yz_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_y_yz_0[j];

            tt_yyy_zz_0[j] = pa_y[j] * tt_yy_zz_0[j] + fl1_fx * tt_y_zz_0[j] + 2.0 * fl1_fz * ts_yyy_zz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_y_zz_0[j];

            tt_yyz_xx_0[j] = pa_y[j] * tt_yz_xx_0[j] + 0.5 * fl1_fx * tt_z_xx_0[j] + 2.0 * fl1_fz * ts_yyz_xx_0[j] - fl1_fz * fl1_fga * ts_z_xx_0[j];

            tt_yyz_xy_0[j] = pa_y[j] * tt_yz_xy_0[j] + 0.5 * fl1_fx * tt_z_xy_0[j] + 0.5 * fl1_fx * tt_yz_x_0[j] + 2.0 * fl1_fz * ts_yyz_xy_0[j] -
                             fl1_fz * fl1_fga * ts_z_xy_0[j];

            tt_yyz_xz_0[j] = pa_y[j] * tt_yz_xz_0[j] + 0.5 * fl1_fx * tt_z_xz_0[j] + 2.0 * fl1_fz * ts_yyz_xz_0[j] - fl1_fz * fl1_fga * ts_z_xz_0[j];

            tt_yyz_yy_0[j] = pa_y[j] * tt_yz_yy_0[j] + 0.5 * fl1_fx * tt_z_yy_0[j] + fl1_fx * tt_yz_y_0[j] + 2.0 * fl1_fz * ts_yyz_yy_0[j] -
                             fl1_fz * fl1_fga * ts_z_yy_0[j];

            tt_yyz_yz_0[j] = pa_y[j] * tt_yz_yz_0[j] + 0.5 * fl1_fx * tt_z_yz_0[j] + 0.5 * fl1_fx * tt_yz_z_0[j] + 2.0 * fl1_fz * ts_yyz_yz_0[j] -
                             fl1_fz * fl1_fga * ts_z_yz_0[j];

            tt_yyz_zz_0[j] = pa_y[j] * tt_yz_zz_0[j] + 0.5 * fl1_fx * tt_z_zz_0[j] + 2.0 * fl1_fz * ts_yyz_zz_0[j] - fl1_fz * fl1_fga * ts_z_zz_0[j];

            tt_yzz_xx_0[j] = pa_y[j] * tt_zz_xx_0[j] + 2.0 * fl1_fz * ts_yzz_xx_0[j];

            tt_yzz_xy_0[j] = pa_y[j] * tt_zz_xy_0[j] + 0.5 * fl1_fx * tt_zz_x_0[j] + 2.0 * fl1_fz * ts_yzz_xy_0[j];

            tt_yzz_xz_0[j] = pa_y[j] * tt_zz_xz_0[j] + 2.0 * fl1_fz * ts_yzz_xz_0[j];

            tt_yzz_yy_0[j] = pa_y[j] * tt_zz_yy_0[j] + fl1_fx * tt_zz_y_0[j] + 2.0 * fl1_fz * ts_yzz_yy_0[j];

            tt_yzz_yz_0[j] = pa_y[j] * tt_zz_yz_0[j] + 0.5 * fl1_fx * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_yzz_yz_0[j];

            tt_yzz_zz_0[j] = pa_y[j] * tt_zz_zz_0[j] + 2.0 * fl1_fz * ts_yzz_zz_0[j];

            tt_zzz_xx_0[j] = pa_z[j] * tt_zz_xx_0[j] + fl1_fx * tt_z_xx_0[j] + 2.0 * fl1_fz * ts_zzz_xx_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_xx_0[j];

            tt_zzz_xy_0[j] = pa_z[j] * tt_zz_xy_0[j] + fl1_fx * tt_z_xy_0[j] + 2.0 * fl1_fz * ts_zzz_xy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_xy_0[j];

            tt_zzz_xz_0[j] = pa_z[j] * tt_zz_xz_0[j] + fl1_fx * tt_z_xz_0[j] + 0.5 * fl1_fx * tt_zz_x_0[j] + 2.0 * fl1_fz * ts_zzz_xz_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_z_xz_0[j];

            tt_zzz_yy_0[j] = pa_z[j] * tt_zz_yy_0[j] + fl1_fx * tt_z_yy_0[j] + 2.0 * fl1_fz * ts_zzz_yy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_z_yy_0[j];

            tt_zzz_yz_0[j] = pa_z[j] * tt_zz_yz_0[j] + fl1_fx * tt_z_yz_0[j] + 0.5 * fl1_fx * tt_zz_y_0[j] + 2.0 * fl1_fz * ts_zzz_yz_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_z_yz_0[j];

            tt_zzz_zz_0[j] = pa_z[j] * tt_zz_zz_0[j] + fl1_fx * tt_z_zz_0[j] + fl1_fx * tt_zz_z_0[j] + 2.0 * fl1_fz * ts_zzz_zz_0[j] -
                             2.0 * fl1_fz * fl1_fga * ts_z_zz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForDG(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    kinrecfunc::compKineticEnergyForDG_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDG_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compKineticEnergyForDG_0_45(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

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

        auto tt_0_xxxx_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx);

        auto tt_0_xxxy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 1);

        auto tt_0_xxxz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 2);

        auto tt_0_xxyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 3);

        auto tt_0_xxyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 4);

        auto tt_0_xxzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 5);

        auto tt_0_xyyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 6);

        auto tt_0_xyyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 7);

        auto tt_0_xyzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 8);

        auto tt_0_xzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 9);

        auto tt_0_yyyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 10);

        auto tt_0_yyyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 11);

        auto tt_0_yyzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 12);

        auto tt_0_yzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 13);

        auto tt_0_zzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 14);

        auto tt_x_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx);

        auto tt_x_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 1);

        auto tt_x_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 2);

        auto tt_x_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 3);

        auto tt_x_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 4);

        auto tt_x_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 5);

        auto tt_x_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 6);

        auto tt_x_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 7);

        auto tt_x_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 8);

        auto tt_x_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 9);

        auto tt_y_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 10);

        auto tt_y_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 11);

        auto tt_y_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 12);

        auto tt_y_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 13);

        auto tt_y_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 14);

        auto tt_y_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 15);

        auto tt_y_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 16);

        auto tt_y_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 17);

        auto tt_y_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 18);

        auto tt_y_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 19);

        auto tt_z_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 20);

        auto tt_z_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 21);

        auto tt_z_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 22);

        auto tt_z_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 23);

        auto tt_z_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 24);

        auto tt_z_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 25);

        auto tt_z_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 26);

        auto tt_z_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 27);

        auto tt_z_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 28);

        auto tt_z_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 29);

        auto ts_xx_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx);

        auto ts_xx_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 1);

        auto ts_xx_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 2);

        auto ts_xx_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 3);

        auto ts_xx_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 4);

        auto ts_xx_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 5);

        auto ts_xx_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 6);

        auto ts_xx_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 7);

        auto ts_xx_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 8);

        auto ts_xx_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 9);

        auto ts_xx_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 10);

        auto ts_xx_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 11);

        auto ts_xx_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 12);

        auto ts_xx_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 13);

        auto ts_xx_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 14);

        auto ts_xy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 15);

        auto ts_xy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 16);

        auto ts_xy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 17);

        auto ts_xy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 18);

        auto ts_xy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 19);

        auto ts_xy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 20);

        auto ts_xy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 21);

        auto ts_xy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 22);

        auto ts_xy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 23);

        auto ts_xy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 24);

        auto ts_xy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 25);

        auto ts_xy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 26);

        auto ts_xy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 27);

        auto ts_xy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 28);

        auto ts_xy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 29);

        auto ts_xz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 30);

        auto ts_xz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 31);

        auto ts_xz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 32);

        auto ts_xz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 33);

        auto ts_xz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 34);

        auto ts_xz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 35);

        auto ts_xz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 36);

        auto ts_xz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 37);

        auto ts_xz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 38);

        auto ts_xz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 39);

        auto ts_xz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 40);

        auto ts_xz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 41);

        auto ts_xz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 42);

        auto ts_xz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 43);

        auto ts_xz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 44);

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

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fga, fx, fz, pa_x, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0, ts_xx_xxxx_0, \
                                     ts_xx_xxxy_0, ts_xx_xxxz_0, ts_xx_xxyy_0, ts_xx_xxyz_0, ts_xx_xxzz_0, ts_xx_xyyy_0, \
                                     ts_xx_xyyz_0, ts_xx_xyzz_0, ts_xx_xzzz_0, ts_xx_yyyy_0, ts_xx_yyyz_0, ts_xx_yyzz_0, \
                                     ts_xx_yzzz_0, ts_xx_zzzz_0, ts_xy_xxxx_0, ts_xy_xxxy_0, ts_xy_xxxz_0, ts_xy_xxyy_0, \
                                     ts_xy_xxyz_0, ts_xy_xxzz_0, ts_xy_xyyy_0, ts_xy_xyyz_0, ts_xy_xyzz_0, ts_xy_xzzz_0, \
                                     ts_xy_yyyy_0, ts_xy_yyyz_0, ts_xy_yyzz_0, ts_xy_yzzz_0, ts_xy_zzzz_0, ts_xz_xxxx_0, \
                                     ts_xz_xxxy_0, ts_xz_xxxz_0, ts_xz_xxyy_0, ts_xz_xxyz_0, ts_xz_xxzz_0, ts_xz_xyyy_0, \
                                     ts_xz_xyyz_0, ts_xz_xyzz_0, ts_xz_xzzz_0, ts_xz_yyyy_0, ts_xz_yyyz_0, ts_xz_yyzz_0, \
                                     ts_xz_yzzz_0, ts_xz_zzzz_0, tt_0_xxxx_0, tt_0_xxxy_0, tt_0_xxxz_0, tt_0_xxyy_0, \
                                     tt_0_xxyz_0, tt_0_xxzz_0, tt_0_xyyy_0, tt_0_xyyz_0, tt_0_xyzz_0, tt_0_xzzz_0, \
                                     tt_0_yyyy_0, tt_0_yyyz_0, tt_0_yyzz_0, tt_0_yzzz_0, tt_0_zzzz_0, tt_x_xxx_0, \
                                     tt_x_xxxx_0, tt_x_xxxy_0, tt_x_xxxz_0, tt_x_xxy_0, tt_x_xxyy_0, tt_x_xxyz_0, \
                                     tt_x_xxz_0, tt_x_xxzz_0, tt_x_xyy_0, tt_x_xyyy_0, tt_x_xyyz_0, tt_x_xyz_0, \
                                     tt_x_xyzz_0, tt_x_xzz_0, tt_x_xzzz_0, tt_x_yyy_0, tt_x_yyyy_0, tt_x_yyyz_0, \
                                     tt_x_yyz_0, tt_x_yyzz_0, tt_x_yzz_0, tt_x_yzzz_0, tt_x_zzz_0, tt_x_zzzz_0, \
                                     tt_xx_xxxx_0, tt_xx_xxxy_0, tt_xx_xxxz_0, tt_xx_xxyy_0, tt_xx_xxyz_0, tt_xx_xxzz_0, \
                                     tt_xx_xyyy_0, tt_xx_xyyz_0, tt_xx_xyzz_0, tt_xx_xzzz_0, tt_xx_yyyy_0, tt_xx_yyyz_0, \
                                     tt_xx_yyzz_0, tt_xx_yzzz_0, tt_xx_zzzz_0, tt_xy_xxxx_0, tt_xy_xxxy_0, tt_xy_xxxz_0, \
                                     tt_xy_xxyy_0, tt_xy_xxyz_0, tt_xy_xxzz_0, tt_xy_xyyy_0, tt_xy_xyyz_0, tt_xy_xyzz_0, \
                                     tt_xy_xzzz_0, tt_xy_yyyy_0, tt_xy_yyyz_0, tt_xy_yyzz_0, tt_xy_yzzz_0, tt_xy_zzzz_0, \
                                     tt_xz_xxxx_0, tt_xz_xxxy_0, tt_xz_xxxz_0, tt_xz_xxyy_0, tt_xz_xxyz_0, tt_xz_xxzz_0, \
                                     tt_xz_xyyy_0, tt_xz_xyyz_0, tt_xz_xyzz_0, tt_xz_xzzz_0, tt_xz_yyyy_0, tt_xz_yyyz_0, \
                                     tt_xz_yyzz_0, tt_xz_yzzz_0, tt_xz_zzzz_0, tt_y_xxx_0, tt_y_xxxx_0, tt_y_xxxy_0, \
                                     tt_y_xxxz_0, tt_y_xxy_0, tt_y_xxyy_0, tt_y_xxyz_0, tt_y_xxz_0, tt_y_xxzz_0, \
                                     tt_y_xyy_0, tt_y_xyyy_0, tt_y_xyyz_0, tt_y_xyz_0, tt_y_xyzz_0, tt_y_xzz_0, \
                                     tt_y_xzzz_0, tt_y_yyy_0, tt_y_yyyy_0, tt_y_yyyz_0, tt_y_yyz_0, tt_y_yyzz_0, \
                                     tt_y_yzz_0, tt_y_yzzz_0, tt_y_zzz_0, tt_y_zzzz_0, tt_z_xxx_0, tt_z_xxxx_0, \
                                     tt_z_xxxy_0, tt_z_xxxz_0, tt_z_xxy_0, tt_z_xxyy_0, tt_z_xxyz_0, tt_z_xxz_0, \
                                     tt_z_xxzz_0, tt_z_xyy_0, tt_z_xyyy_0, tt_z_xyyz_0, tt_z_xyz_0, tt_z_xyzz_0, \
                                     tt_z_xzz_0, tt_z_xzzz_0, tt_z_yyy_0, tt_z_yyyy_0, tt_z_yyyz_0, tt_z_yyz_0, \
                                     tt_z_yyzz_0, tt_z_yzz_0, tt_z_yzzz_0, tt_z_zzz_0, tt_z_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xx_xxxx_0[j] = pa_x[j] * tt_x_xxxx_0[j] + 0.5 * fl1_fx * tt_0_xxxx_0[j] + 2.0 * fl1_fx * tt_x_xxx_0[j] +
                              2.0 * fl1_fz * ts_xx_xxxx_0[j] - fl1_fz * fl1_fga * ts_0_xxxx_0[j];

            tt_xx_xxxy_0[j] = pa_x[j] * tt_x_xxxy_0[j] + 0.5 * fl1_fx * tt_0_xxxy_0[j] + 1.5 * fl1_fx * tt_x_xxy_0[j] +
                              2.0 * fl1_fz * ts_xx_xxxy_0[j] - fl1_fz * fl1_fga * ts_0_xxxy_0[j];

            tt_xx_xxxz_0[j] = pa_x[j] * tt_x_xxxz_0[j] + 0.5 * fl1_fx * tt_0_xxxz_0[j] + 1.5 * fl1_fx * tt_x_xxz_0[j] +
                              2.0 * fl1_fz * ts_xx_xxxz_0[j] - fl1_fz * fl1_fga * ts_0_xxxz_0[j];

            tt_xx_xxyy_0[j] = pa_x[j] * tt_x_xxyy_0[j] + 0.5 * fl1_fx * tt_0_xxyy_0[j] + fl1_fx * tt_x_xyy_0[j] + 2.0 * fl1_fz * ts_xx_xxyy_0[j] -
                              fl1_fz * fl1_fga * ts_0_xxyy_0[j];

            tt_xx_xxyz_0[j] = pa_x[j] * tt_x_xxyz_0[j] + 0.5 * fl1_fx * tt_0_xxyz_0[j] + fl1_fx * tt_x_xyz_0[j] + 2.0 * fl1_fz * ts_xx_xxyz_0[j] -
                              fl1_fz * fl1_fga * ts_0_xxyz_0[j];

            tt_xx_xxzz_0[j] = pa_x[j] * tt_x_xxzz_0[j] + 0.5 * fl1_fx * tt_0_xxzz_0[j] + fl1_fx * tt_x_xzz_0[j] + 2.0 * fl1_fz * ts_xx_xxzz_0[j] -
                              fl1_fz * fl1_fga * ts_0_xxzz_0[j];

            tt_xx_xyyy_0[j] = pa_x[j] * tt_x_xyyy_0[j] + 0.5 * fl1_fx * tt_0_xyyy_0[j] + 0.5 * fl1_fx * tt_x_yyy_0[j] +
                              2.0 * fl1_fz * ts_xx_xyyy_0[j] - fl1_fz * fl1_fga * ts_0_xyyy_0[j];

            tt_xx_xyyz_0[j] = pa_x[j] * tt_x_xyyz_0[j] + 0.5 * fl1_fx * tt_0_xyyz_0[j] + 0.5 * fl1_fx * tt_x_yyz_0[j] +
                              2.0 * fl1_fz * ts_xx_xyyz_0[j] - fl1_fz * fl1_fga * ts_0_xyyz_0[j];

            tt_xx_xyzz_0[j] = pa_x[j] * tt_x_xyzz_0[j] + 0.5 * fl1_fx * tt_0_xyzz_0[j] + 0.5 * fl1_fx * tt_x_yzz_0[j] +
                              2.0 * fl1_fz * ts_xx_xyzz_0[j] - fl1_fz * fl1_fga * ts_0_xyzz_0[j];

            tt_xx_xzzz_0[j] = pa_x[j] * tt_x_xzzz_0[j] + 0.5 * fl1_fx * tt_0_xzzz_0[j] + 0.5 * fl1_fx * tt_x_zzz_0[j] +
                              2.0 * fl1_fz * ts_xx_xzzz_0[j] - fl1_fz * fl1_fga * ts_0_xzzz_0[j];

            tt_xx_yyyy_0[j] =
                pa_x[j] * tt_x_yyyy_0[j] + 0.5 * fl1_fx * tt_0_yyyy_0[j] + 2.0 * fl1_fz * ts_xx_yyyy_0[j] - fl1_fz * fl1_fga * ts_0_yyyy_0[j];

            tt_xx_yyyz_0[j] =
                pa_x[j] * tt_x_yyyz_0[j] + 0.5 * fl1_fx * tt_0_yyyz_0[j] + 2.0 * fl1_fz * ts_xx_yyyz_0[j] - fl1_fz * fl1_fga * ts_0_yyyz_0[j];

            tt_xx_yyzz_0[j] =
                pa_x[j] * tt_x_yyzz_0[j] + 0.5 * fl1_fx * tt_0_yyzz_0[j] + 2.0 * fl1_fz * ts_xx_yyzz_0[j] - fl1_fz * fl1_fga * ts_0_yyzz_0[j];

            tt_xx_yzzz_0[j] =
                pa_x[j] * tt_x_yzzz_0[j] + 0.5 * fl1_fx * tt_0_yzzz_0[j] + 2.0 * fl1_fz * ts_xx_yzzz_0[j] - fl1_fz * fl1_fga * ts_0_yzzz_0[j];

            tt_xx_zzzz_0[j] =
                pa_x[j] * tt_x_zzzz_0[j] + 0.5 * fl1_fx * tt_0_zzzz_0[j] + 2.0 * fl1_fz * ts_xx_zzzz_0[j] - fl1_fz * fl1_fga * ts_0_zzzz_0[j];

            tt_xy_xxxx_0[j] = pa_x[j] * tt_y_xxxx_0[j] + 2.0 * fl1_fx * tt_y_xxx_0[j] + 2.0 * fl1_fz * ts_xy_xxxx_0[j];

            tt_xy_xxxy_0[j] = pa_x[j] * tt_y_xxxy_0[j] + 1.5 * fl1_fx * tt_y_xxy_0[j] + 2.0 * fl1_fz * ts_xy_xxxy_0[j];

            tt_xy_xxxz_0[j] = pa_x[j] * tt_y_xxxz_0[j] + 1.5 * fl1_fx * tt_y_xxz_0[j] + 2.0 * fl1_fz * ts_xy_xxxz_0[j];

            tt_xy_xxyy_0[j] = pa_x[j] * tt_y_xxyy_0[j] + fl1_fx * tt_y_xyy_0[j] + 2.0 * fl1_fz * ts_xy_xxyy_0[j];

            tt_xy_xxyz_0[j] = pa_x[j] * tt_y_xxyz_0[j] + fl1_fx * tt_y_xyz_0[j] + 2.0 * fl1_fz * ts_xy_xxyz_0[j];

            tt_xy_xxzz_0[j] = pa_x[j] * tt_y_xxzz_0[j] + fl1_fx * tt_y_xzz_0[j] + 2.0 * fl1_fz * ts_xy_xxzz_0[j];

            tt_xy_xyyy_0[j] = pa_x[j] * tt_y_xyyy_0[j] + 0.5 * fl1_fx * tt_y_yyy_0[j] + 2.0 * fl1_fz * ts_xy_xyyy_0[j];

            tt_xy_xyyz_0[j] = pa_x[j] * tt_y_xyyz_0[j] + 0.5 * fl1_fx * tt_y_yyz_0[j] + 2.0 * fl1_fz * ts_xy_xyyz_0[j];

            tt_xy_xyzz_0[j] = pa_x[j] * tt_y_xyzz_0[j] + 0.5 * fl1_fx * tt_y_yzz_0[j] + 2.0 * fl1_fz * ts_xy_xyzz_0[j];

            tt_xy_xzzz_0[j] = pa_x[j] * tt_y_xzzz_0[j] + 0.5 * fl1_fx * tt_y_zzz_0[j] + 2.0 * fl1_fz * ts_xy_xzzz_0[j];

            tt_xy_yyyy_0[j] = pa_x[j] * tt_y_yyyy_0[j] + 2.0 * fl1_fz * ts_xy_yyyy_0[j];

            tt_xy_yyyz_0[j] = pa_x[j] * tt_y_yyyz_0[j] + 2.0 * fl1_fz * ts_xy_yyyz_0[j];

            tt_xy_yyzz_0[j] = pa_x[j] * tt_y_yyzz_0[j] + 2.0 * fl1_fz * ts_xy_yyzz_0[j];

            tt_xy_yzzz_0[j] = pa_x[j] * tt_y_yzzz_0[j] + 2.0 * fl1_fz * ts_xy_yzzz_0[j];

            tt_xy_zzzz_0[j] = pa_x[j] * tt_y_zzzz_0[j] + 2.0 * fl1_fz * ts_xy_zzzz_0[j];

            tt_xz_xxxx_0[j] = pa_x[j] * tt_z_xxxx_0[j] + 2.0 * fl1_fx * tt_z_xxx_0[j] + 2.0 * fl1_fz * ts_xz_xxxx_0[j];

            tt_xz_xxxy_0[j] = pa_x[j] * tt_z_xxxy_0[j] + 1.5 * fl1_fx * tt_z_xxy_0[j] + 2.0 * fl1_fz * ts_xz_xxxy_0[j];

            tt_xz_xxxz_0[j] = pa_x[j] * tt_z_xxxz_0[j] + 1.5 * fl1_fx * tt_z_xxz_0[j] + 2.0 * fl1_fz * ts_xz_xxxz_0[j];

            tt_xz_xxyy_0[j] = pa_x[j] * tt_z_xxyy_0[j] + fl1_fx * tt_z_xyy_0[j] + 2.0 * fl1_fz * ts_xz_xxyy_0[j];

            tt_xz_xxyz_0[j] = pa_x[j] * tt_z_xxyz_0[j] + fl1_fx * tt_z_xyz_0[j] + 2.0 * fl1_fz * ts_xz_xxyz_0[j];

            tt_xz_xxzz_0[j] = pa_x[j] * tt_z_xxzz_0[j] + fl1_fx * tt_z_xzz_0[j] + 2.0 * fl1_fz * ts_xz_xxzz_0[j];

            tt_xz_xyyy_0[j] = pa_x[j] * tt_z_xyyy_0[j] + 0.5 * fl1_fx * tt_z_yyy_0[j] + 2.0 * fl1_fz * ts_xz_xyyy_0[j];

            tt_xz_xyyz_0[j] = pa_x[j] * tt_z_xyyz_0[j] + 0.5 * fl1_fx * tt_z_yyz_0[j] + 2.0 * fl1_fz * ts_xz_xyyz_0[j];

            tt_xz_xyzz_0[j] = pa_x[j] * tt_z_xyzz_0[j] + 0.5 * fl1_fx * tt_z_yzz_0[j] + 2.0 * fl1_fz * ts_xz_xyzz_0[j];

            tt_xz_xzzz_0[j] = pa_x[j] * tt_z_xzzz_0[j] + 0.5 * fl1_fx * tt_z_zzz_0[j] + 2.0 * fl1_fz * ts_xz_xzzz_0[j];

            tt_xz_yyyy_0[j] = pa_x[j] * tt_z_yyyy_0[j] + 2.0 * fl1_fz * ts_xz_yyyy_0[j];

            tt_xz_yyyz_0[j] = pa_x[j] * tt_z_yyyz_0[j] + 2.0 * fl1_fz * ts_xz_yyyz_0[j];

            tt_xz_yyzz_0[j] = pa_x[j] * tt_z_yyzz_0[j] + 2.0 * fl1_fz * ts_xz_yyzz_0[j];

            tt_xz_yzzz_0[j] = pa_x[j] * tt_z_yzzz_0[j] + 2.0 * fl1_fz * ts_xz_yzzz_0[j];

            tt_xz_zzzz_0[j] = pa_x[j] * tt_z_zzzz_0[j] + 2.0 * fl1_fz * ts_xz_zzzz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForDG_45_90(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_2_4_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_1_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_0_4_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {0, -1, -1, -1}, {4, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_1_3_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {1, -1, -1, -1}, {3, -1, -1, -1}, 1, 1, 0));

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

        auto tt_0_xxxx_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx);

        auto tt_0_xxxy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 1);

        auto tt_0_xxxz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 2);

        auto tt_0_xxyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 3);

        auto tt_0_xxyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 4);

        auto tt_0_xxzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 5);

        auto tt_0_xyyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 6);

        auto tt_0_xyyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 7);

        auto tt_0_xyzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 8);

        auto tt_0_xzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 9);

        auto tt_0_yyyy_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 10);

        auto tt_0_yyyz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 11);

        auto tt_0_yyzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 12);

        auto tt_0_yzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 13);

        auto tt_0_zzzz_0 = primBuffer.data(pidx_t_0_4_m0 + 15 * idx + 14);

        auto tt_y_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 10);

        auto tt_y_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 11);

        auto tt_y_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 12);

        auto tt_y_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 13);

        auto tt_y_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 14);

        auto tt_y_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 15);

        auto tt_y_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 16);

        auto tt_y_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 17);

        auto tt_y_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 18);

        auto tt_y_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 19);

        auto tt_z_xxx_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 20);

        auto tt_z_xxy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 21);

        auto tt_z_xxz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 22);

        auto tt_z_xyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 23);

        auto tt_z_xyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 24);

        auto tt_z_xzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 25);

        auto tt_z_yyy_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 26);

        auto tt_z_yyz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 27);

        auto tt_z_yzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 28);

        auto tt_z_zzz_0 = primBuffer.data(pidx_t_1_3_m0 + 30 * idx + 29);

        auto ts_yy_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 45);

        auto ts_yy_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 46);

        auto ts_yy_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 47);

        auto ts_yy_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 48);

        auto ts_yy_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 49);

        auto ts_yy_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 50);

        auto ts_yy_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 51);

        auto ts_yy_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 52);

        auto ts_yy_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 53);

        auto ts_yy_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 54);

        auto ts_yy_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 55);

        auto ts_yy_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 56);

        auto ts_yy_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 57);

        auto ts_yy_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 58);

        auto ts_yy_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 59);

        auto ts_yz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 60);

        auto ts_yz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 61);

        auto ts_yz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 62);

        auto ts_yz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 63);

        auto ts_yz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 64);

        auto ts_yz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 65);

        auto ts_yz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 66);

        auto ts_yz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 67);

        auto ts_yz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 68);

        auto ts_yz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 69);

        auto ts_yz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 70);

        auto ts_yz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 71);

        auto ts_yz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 72);

        auto ts_yz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 73);

        auto ts_yz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 74);

        auto ts_zz_xxxx_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 75);

        auto ts_zz_xxxy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 76);

        auto ts_zz_xxxz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 77);

        auto ts_zz_xxyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 78);

        auto ts_zz_xxyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 79);

        auto ts_zz_xxzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 80);

        auto ts_zz_xyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 81);

        auto ts_zz_xyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 82);

        auto ts_zz_xyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 83);

        auto ts_zz_xzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 84);

        auto ts_zz_yyyy_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 85);

        auto ts_zz_yyyz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 86);

        auto ts_zz_yyzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 87);

        auto ts_zz_yzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 88);

        auto ts_zz_zzzz_0 = primBuffer.data(pidx_s_2_4_m0 + 90 * idx + 89);

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

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fga, fx, fz, pa_y, pa_z, ts_0_xxxx_0, ts_0_xxxy_0, ts_0_xxxz_0, ts_0_xxyy_0, \
                                     ts_0_xxyz_0, ts_0_xxzz_0, ts_0_xyyy_0, ts_0_xyyz_0, ts_0_xyzz_0, ts_0_xzzz_0, \
                                     ts_0_yyyy_0, ts_0_yyyz_0, ts_0_yyzz_0, ts_0_yzzz_0, ts_0_zzzz_0, ts_yy_xxxx_0, \
                                     ts_yy_xxxy_0, ts_yy_xxxz_0, ts_yy_xxyy_0, ts_yy_xxyz_0, ts_yy_xxzz_0, ts_yy_xyyy_0, \
                                     ts_yy_xyyz_0, ts_yy_xyzz_0, ts_yy_xzzz_0, ts_yy_yyyy_0, ts_yy_yyyz_0, ts_yy_yyzz_0, \
                                     ts_yy_yzzz_0, ts_yy_zzzz_0, ts_yz_xxxx_0, ts_yz_xxxy_0, ts_yz_xxxz_0, ts_yz_xxyy_0, \
                                     ts_yz_xxyz_0, ts_yz_xxzz_0, ts_yz_xyyy_0, ts_yz_xyyz_0, ts_yz_xyzz_0, ts_yz_xzzz_0, \
                                     ts_yz_yyyy_0, ts_yz_yyyz_0, ts_yz_yyzz_0, ts_yz_yzzz_0, ts_yz_zzzz_0, ts_zz_xxxx_0, \
                                     ts_zz_xxxy_0, ts_zz_xxxz_0, ts_zz_xxyy_0, ts_zz_xxyz_0, ts_zz_xxzz_0, ts_zz_xyyy_0, \
                                     ts_zz_xyyz_0, ts_zz_xyzz_0, ts_zz_xzzz_0, ts_zz_yyyy_0, ts_zz_yyyz_0, ts_zz_yyzz_0, \
                                     ts_zz_yzzz_0, ts_zz_zzzz_0, tt_0_xxxx_0, tt_0_xxxy_0, tt_0_xxxz_0, tt_0_xxyy_0, \
                                     tt_0_xxyz_0, tt_0_xxzz_0, tt_0_xyyy_0, tt_0_xyyz_0, tt_0_xyzz_0, tt_0_xzzz_0, \
                                     tt_0_yyyy_0, tt_0_yyyz_0, tt_0_yyzz_0, tt_0_yzzz_0, tt_0_zzzz_0, tt_y_xxx_0, \
                                     tt_y_xxxx_0, tt_y_xxxy_0, tt_y_xxxz_0, tt_y_xxy_0, tt_y_xxyy_0, tt_y_xxyz_0, \
                                     tt_y_xxz_0, tt_y_xxzz_0, tt_y_xyy_0, tt_y_xyyy_0, tt_y_xyyz_0, tt_y_xyz_0, \
                                     tt_y_xyzz_0, tt_y_xzz_0, tt_y_xzzz_0, tt_y_yyy_0, tt_y_yyyy_0, tt_y_yyyz_0, \
                                     tt_y_yyz_0, tt_y_yyzz_0, tt_y_yzz_0, tt_y_yzzz_0, tt_y_zzz_0, tt_y_zzzz_0, \
                                     tt_yy_xxxx_0, tt_yy_xxxy_0, tt_yy_xxxz_0, tt_yy_xxyy_0, tt_yy_xxyz_0, tt_yy_xxzz_0, \
                                     tt_yy_xyyy_0, tt_yy_xyyz_0, tt_yy_xyzz_0, tt_yy_xzzz_0, tt_yy_yyyy_0, tt_yy_yyyz_0, \
                                     tt_yy_yyzz_0, tt_yy_yzzz_0, tt_yy_zzzz_0, tt_yz_xxxx_0, tt_yz_xxxy_0, tt_yz_xxxz_0, \
                                     tt_yz_xxyy_0, tt_yz_xxyz_0, tt_yz_xxzz_0, tt_yz_xyyy_0, tt_yz_xyyz_0, tt_yz_xyzz_0, \
                                     tt_yz_xzzz_0, tt_yz_yyyy_0, tt_yz_yyyz_0, tt_yz_yyzz_0, tt_yz_yzzz_0, tt_yz_zzzz_0, \
                                     tt_z_xxx_0, tt_z_xxxx_0, tt_z_xxxy_0, tt_z_xxxz_0, tt_z_xxy_0, tt_z_xxyy_0, \
                                     tt_z_xxyz_0, tt_z_xxz_0, tt_z_xxzz_0, tt_z_xyy_0, tt_z_xyyy_0, tt_z_xyyz_0, \
                                     tt_z_xyz_0, tt_z_xyzz_0, tt_z_xzz_0, tt_z_xzzz_0, tt_z_yyy_0, tt_z_yyyy_0, \
                                     tt_z_yyyz_0, tt_z_yyz_0, tt_z_yyzz_0, tt_z_yzz_0, tt_z_yzzz_0, tt_z_zzz_0, \
                                     tt_z_zzzz_0, tt_zz_xxxx_0, tt_zz_xxxy_0, tt_zz_xxxz_0, tt_zz_xxyy_0, tt_zz_xxyz_0, \
                                     tt_zz_xxzz_0, tt_zz_xyyy_0, tt_zz_xyyz_0, tt_zz_xyzz_0, tt_zz_xzzz_0, tt_zz_yyyy_0, \
                                     tt_zz_yyyz_0, tt_zz_yyzz_0, tt_zz_yzzz_0, tt_zz_zzzz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_yy_xxxx_0[j] =
                pa_y[j] * tt_y_xxxx_0[j] + 0.5 * fl1_fx * tt_0_xxxx_0[j] + 2.0 * fl1_fz * ts_yy_xxxx_0[j] - fl1_fz * fl1_fga * ts_0_xxxx_0[j];

            tt_yy_xxxy_0[j] = pa_y[j] * tt_y_xxxy_0[j] + 0.5 * fl1_fx * tt_0_xxxy_0[j] + 0.5 * fl1_fx * tt_y_xxx_0[j] +
                              2.0 * fl1_fz * ts_yy_xxxy_0[j] - fl1_fz * fl1_fga * ts_0_xxxy_0[j];

            tt_yy_xxxz_0[j] =
                pa_y[j] * tt_y_xxxz_0[j] + 0.5 * fl1_fx * tt_0_xxxz_0[j] + 2.0 * fl1_fz * ts_yy_xxxz_0[j] - fl1_fz * fl1_fga * ts_0_xxxz_0[j];

            tt_yy_xxyy_0[j] = pa_y[j] * tt_y_xxyy_0[j] + 0.5 * fl1_fx * tt_0_xxyy_0[j] + fl1_fx * tt_y_xxy_0[j] + 2.0 * fl1_fz * ts_yy_xxyy_0[j] -
                              fl1_fz * fl1_fga * ts_0_xxyy_0[j];

            tt_yy_xxyz_0[j] = pa_y[j] * tt_y_xxyz_0[j] + 0.5 * fl1_fx * tt_0_xxyz_0[j] + 0.5 * fl1_fx * tt_y_xxz_0[j] +
                              2.0 * fl1_fz * ts_yy_xxyz_0[j] - fl1_fz * fl1_fga * ts_0_xxyz_0[j];

            tt_yy_xxzz_0[j] =
                pa_y[j] * tt_y_xxzz_0[j] + 0.5 * fl1_fx * tt_0_xxzz_0[j] + 2.0 * fl1_fz * ts_yy_xxzz_0[j] - fl1_fz * fl1_fga * ts_0_xxzz_0[j];

            tt_yy_xyyy_0[j] = pa_y[j] * tt_y_xyyy_0[j] + 0.5 * fl1_fx * tt_0_xyyy_0[j] + 1.5 * fl1_fx * tt_y_xyy_0[j] +
                              2.0 * fl1_fz * ts_yy_xyyy_0[j] - fl1_fz * fl1_fga * ts_0_xyyy_0[j];

            tt_yy_xyyz_0[j] = pa_y[j] * tt_y_xyyz_0[j] + 0.5 * fl1_fx * tt_0_xyyz_0[j] + fl1_fx * tt_y_xyz_0[j] + 2.0 * fl1_fz * ts_yy_xyyz_0[j] -
                              fl1_fz * fl1_fga * ts_0_xyyz_0[j];

            tt_yy_xyzz_0[j] = pa_y[j] * tt_y_xyzz_0[j] + 0.5 * fl1_fx * tt_0_xyzz_0[j] + 0.5 * fl1_fx * tt_y_xzz_0[j] +
                              2.0 * fl1_fz * ts_yy_xyzz_0[j] - fl1_fz * fl1_fga * ts_0_xyzz_0[j];

            tt_yy_xzzz_0[j] =
                pa_y[j] * tt_y_xzzz_0[j] + 0.5 * fl1_fx * tt_0_xzzz_0[j] + 2.0 * fl1_fz * ts_yy_xzzz_0[j] - fl1_fz * fl1_fga * ts_0_xzzz_0[j];

            tt_yy_yyyy_0[j] = pa_y[j] * tt_y_yyyy_0[j] + 0.5 * fl1_fx * tt_0_yyyy_0[j] + 2.0 * fl1_fx * tt_y_yyy_0[j] +
                              2.0 * fl1_fz * ts_yy_yyyy_0[j] - fl1_fz * fl1_fga * ts_0_yyyy_0[j];

            tt_yy_yyyz_0[j] = pa_y[j] * tt_y_yyyz_0[j] + 0.5 * fl1_fx * tt_0_yyyz_0[j] + 1.5 * fl1_fx * tt_y_yyz_0[j] +
                              2.0 * fl1_fz * ts_yy_yyyz_0[j] - fl1_fz * fl1_fga * ts_0_yyyz_0[j];

            tt_yy_yyzz_0[j] = pa_y[j] * tt_y_yyzz_0[j] + 0.5 * fl1_fx * tt_0_yyzz_0[j] + fl1_fx * tt_y_yzz_0[j] + 2.0 * fl1_fz * ts_yy_yyzz_0[j] -
                              fl1_fz * fl1_fga * ts_0_yyzz_0[j];

            tt_yy_yzzz_0[j] = pa_y[j] * tt_y_yzzz_0[j] + 0.5 * fl1_fx * tt_0_yzzz_0[j] + 0.5 * fl1_fx * tt_y_zzz_0[j] +
                              2.0 * fl1_fz * ts_yy_yzzz_0[j] - fl1_fz * fl1_fga * ts_0_yzzz_0[j];

            tt_yy_zzzz_0[j] =
                pa_y[j] * tt_y_zzzz_0[j] + 0.5 * fl1_fx * tt_0_zzzz_0[j] + 2.0 * fl1_fz * ts_yy_zzzz_0[j] - fl1_fz * fl1_fga * ts_0_zzzz_0[j];

            tt_yz_xxxx_0[j] = pa_y[j] * tt_z_xxxx_0[j] + 2.0 * fl1_fz * ts_yz_xxxx_0[j];

            tt_yz_xxxy_0[j] = pa_y[j] * tt_z_xxxy_0[j] + 0.5 * fl1_fx * tt_z_xxx_0[j] + 2.0 * fl1_fz * ts_yz_xxxy_0[j];

            tt_yz_xxxz_0[j] = pa_y[j] * tt_z_xxxz_0[j] + 2.0 * fl1_fz * ts_yz_xxxz_0[j];

            tt_yz_xxyy_0[j] = pa_y[j] * tt_z_xxyy_0[j] + fl1_fx * tt_z_xxy_0[j] + 2.0 * fl1_fz * ts_yz_xxyy_0[j];

            tt_yz_xxyz_0[j] = pa_y[j] * tt_z_xxyz_0[j] + 0.5 * fl1_fx * tt_z_xxz_0[j] + 2.0 * fl1_fz * ts_yz_xxyz_0[j];

            tt_yz_xxzz_0[j] = pa_y[j] * tt_z_xxzz_0[j] + 2.0 * fl1_fz * ts_yz_xxzz_0[j];

            tt_yz_xyyy_0[j] = pa_y[j] * tt_z_xyyy_0[j] + 1.5 * fl1_fx * tt_z_xyy_0[j] + 2.0 * fl1_fz * ts_yz_xyyy_0[j];

            tt_yz_xyyz_0[j] = pa_y[j] * tt_z_xyyz_0[j] + fl1_fx * tt_z_xyz_0[j] + 2.0 * fl1_fz * ts_yz_xyyz_0[j];

            tt_yz_xyzz_0[j] = pa_y[j] * tt_z_xyzz_0[j] + 0.5 * fl1_fx * tt_z_xzz_0[j] + 2.0 * fl1_fz * ts_yz_xyzz_0[j];

            tt_yz_xzzz_0[j] = pa_y[j] * tt_z_xzzz_0[j] + 2.0 * fl1_fz * ts_yz_xzzz_0[j];

            tt_yz_yyyy_0[j] = pa_y[j] * tt_z_yyyy_0[j] + 2.0 * fl1_fx * tt_z_yyy_0[j] + 2.0 * fl1_fz * ts_yz_yyyy_0[j];

            tt_yz_yyyz_0[j] = pa_y[j] * tt_z_yyyz_0[j] + 1.5 * fl1_fx * tt_z_yyz_0[j] + 2.0 * fl1_fz * ts_yz_yyyz_0[j];

            tt_yz_yyzz_0[j] = pa_y[j] * tt_z_yyzz_0[j] + fl1_fx * tt_z_yzz_0[j] + 2.0 * fl1_fz * ts_yz_yyzz_0[j];

            tt_yz_yzzz_0[j] = pa_y[j] * tt_z_yzzz_0[j] + 0.5 * fl1_fx * tt_z_zzz_0[j] + 2.0 * fl1_fz * ts_yz_yzzz_0[j];

            tt_yz_zzzz_0[j] = pa_y[j] * tt_z_zzzz_0[j] + 2.0 * fl1_fz * ts_yz_zzzz_0[j];

            tt_zz_xxxx_0[j] =
                pa_z[j] * tt_z_xxxx_0[j] + 0.5 * fl1_fx * tt_0_xxxx_0[j] + 2.0 * fl1_fz * ts_zz_xxxx_0[j] - fl1_fz * fl1_fga * ts_0_xxxx_0[j];

            tt_zz_xxxy_0[j] =
                pa_z[j] * tt_z_xxxy_0[j] + 0.5 * fl1_fx * tt_0_xxxy_0[j] + 2.0 * fl1_fz * ts_zz_xxxy_0[j] - fl1_fz * fl1_fga * ts_0_xxxy_0[j];

            tt_zz_xxxz_0[j] = pa_z[j] * tt_z_xxxz_0[j] + 0.5 * fl1_fx * tt_0_xxxz_0[j] + 0.5 * fl1_fx * tt_z_xxx_0[j] +
                              2.0 * fl1_fz * ts_zz_xxxz_0[j] - fl1_fz * fl1_fga * ts_0_xxxz_0[j];

            tt_zz_xxyy_0[j] =
                pa_z[j] * tt_z_xxyy_0[j] + 0.5 * fl1_fx * tt_0_xxyy_0[j] + 2.0 * fl1_fz * ts_zz_xxyy_0[j] - fl1_fz * fl1_fga * ts_0_xxyy_0[j];

            tt_zz_xxyz_0[j] = pa_z[j] * tt_z_xxyz_0[j] + 0.5 * fl1_fx * tt_0_xxyz_0[j] + 0.5 * fl1_fx * tt_z_xxy_0[j] +
                              2.0 * fl1_fz * ts_zz_xxyz_0[j] - fl1_fz * fl1_fga * ts_0_xxyz_0[j];

            tt_zz_xxzz_0[j] = pa_z[j] * tt_z_xxzz_0[j] + 0.5 * fl1_fx * tt_0_xxzz_0[j] + fl1_fx * tt_z_xxz_0[j] + 2.0 * fl1_fz * ts_zz_xxzz_0[j] -
                              fl1_fz * fl1_fga * ts_0_xxzz_0[j];

            tt_zz_xyyy_0[j] =
                pa_z[j] * tt_z_xyyy_0[j] + 0.5 * fl1_fx * tt_0_xyyy_0[j] + 2.0 * fl1_fz * ts_zz_xyyy_0[j] - fl1_fz * fl1_fga * ts_0_xyyy_0[j];

            tt_zz_xyyz_0[j] = pa_z[j] * tt_z_xyyz_0[j] + 0.5 * fl1_fx * tt_0_xyyz_0[j] + 0.5 * fl1_fx * tt_z_xyy_0[j] +
                              2.0 * fl1_fz * ts_zz_xyyz_0[j] - fl1_fz * fl1_fga * ts_0_xyyz_0[j];

            tt_zz_xyzz_0[j] = pa_z[j] * tt_z_xyzz_0[j] + 0.5 * fl1_fx * tt_0_xyzz_0[j] + fl1_fx * tt_z_xyz_0[j] + 2.0 * fl1_fz * ts_zz_xyzz_0[j] -
                              fl1_fz * fl1_fga * ts_0_xyzz_0[j];

            tt_zz_xzzz_0[j] = pa_z[j] * tt_z_xzzz_0[j] + 0.5 * fl1_fx * tt_0_xzzz_0[j] + 1.5 * fl1_fx * tt_z_xzz_0[j] +
                              2.0 * fl1_fz * ts_zz_xzzz_0[j] - fl1_fz * fl1_fga * ts_0_xzzz_0[j];

            tt_zz_yyyy_0[j] =
                pa_z[j] * tt_z_yyyy_0[j] + 0.5 * fl1_fx * tt_0_yyyy_0[j] + 2.0 * fl1_fz * ts_zz_yyyy_0[j] - fl1_fz * fl1_fga * ts_0_yyyy_0[j];

            tt_zz_yyyz_0[j] = pa_z[j] * tt_z_yyyz_0[j] + 0.5 * fl1_fx * tt_0_yyyz_0[j] + 0.5 * fl1_fx * tt_z_yyy_0[j] +
                              2.0 * fl1_fz * ts_zz_yyyz_0[j] - fl1_fz * fl1_fga * ts_0_yyyz_0[j];

            tt_zz_yyzz_0[j] = pa_z[j] * tt_z_yyzz_0[j] + 0.5 * fl1_fx * tt_0_yyzz_0[j] + fl1_fx * tt_z_yyz_0[j] + 2.0 * fl1_fz * ts_zz_yyzz_0[j] -
                              fl1_fz * fl1_fga * ts_0_yyzz_0[j];

            tt_zz_yzzz_0[j] = pa_z[j] * tt_z_yzzz_0[j] + 0.5 * fl1_fx * tt_0_yzzz_0[j] + 1.5 * fl1_fx * tt_z_yzz_0[j] +
                              2.0 * fl1_fz * ts_zz_yzzz_0[j] - fl1_fz * fl1_fga * ts_0_yzzz_0[j];

            tt_zz_zzzz_0[j] = pa_z[j] * tt_z_zzzz_0[j] + 0.5 * fl1_fx * tt_0_zzzz_0[j] + 2.0 * fl1_fx * tt_z_zzz_0[j] +
                              2.0 * fl1_fz * ts_zz_zzzz_0[j] - fl1_fz * fl1_fga * ts_0_zzzz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForGD(CMemBlock2D<double>&       primBuffer,
                       const CRecursionMap&       recursionMap,
                       const CMemBlock2D<double>& osFactors,
                       const CMemBlock2D<double>& paDistances,
                       const CGtoBlock&           braGtoBlock,
                       const CGtoBlock&           ketGtoBlock,
                       const int32_t              iContrGto)
{
    kinrecfunc::compKineticEnergyForGD_0_45(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForGD_45_90(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

void
compKineticEnergyForGD_0_45(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto tt_xxx_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx);

        auto tt_xxx_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 1);

        auto tt_xxx_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 2);

        auto tt_xxx_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 3);

        auto tt_xxx_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 4);

        auto tt_xxx_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 5);

        auto tt_xxy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 6);

        auto tt_xxy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 7);

        auto tt_xxy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 8);

        auto tt_xxy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 9);

        auto tt_xxy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 10);

        auto tt_xxy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 11);

        auto tt_xxz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 12);

        auto tt_xxz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 13);

        auto tt_xxz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 14);

        auto tt_xxz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 15);

        auto tt_xxz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 16);

        auto tt_xxz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 17);

        auto tt_xyy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 18);

        auto tt_xyy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 19);

        auto tt_xyy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 20);

        auto tt_xyy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 21);

        auto tt_xyy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 22);

        auto tt_xyy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 23);

        auto tt_xyz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 24);

        auto tt_xyz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 25);

        auto tt_xyz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 26);

        auto tt_xyz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 27);

        auto tt_xyz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 28);

        auto tt_xyz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 29);

        auto tt_xzz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 30);

        auto tt_xzz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 31);

        auto tt_xzz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 32);

        auto tt_xzz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 33);

        auto tt_xzz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 34);

        auto tt_xzz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 35);

        auto tt_yyy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 36);

        auto tt_yyy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 37);

        auto tt_yyy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 38);

        auto tt_yyy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 39);

        auto tt_yyy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 40);

        auto tt_yyy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 41);

        auto tt_yyz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 42);

        auto tt_yyz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 43);

        auto tt_yyz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 44);

        auto tt_xx_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx);

        auto tt_xx_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 1);

        auto tt_xx_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 2);

        auto tt_xx_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 3);

        auto tt_xx_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 4);

        auto tt_xx_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 5);

        auto tt_xy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 6);

        auto tt_xy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 7);

        auto tt_xy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 8);

        auto tt_xy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 9);

        auto tt_xy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 10);

        auto tt_xy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 11);

        auto tt_xz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 12);

        auto tt_xz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 13);

        auto tt_xz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 14);

        auto tt_xz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 15);

        auto tt_xz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 16);

        auto tt_xz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 17);

        auto tt_yy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 18);

        auto tt_yy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 19);

        auto tt_yy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 20);

        auto tt_yy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 21);

        auto tt_yy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 22);

        auto tt_yy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 23);

        auto tt_yz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 24);

        auto tt_yz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 25);

        auto tt_yz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 26);

        auto tt_yz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 27);

        auto tt_yz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 28);

        auto tt_yz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 29);

        auto tt_zz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 30);

        auto tt_zz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 31);

        auto tt_zz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 32);

        auto tt_zz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 33);

        auto tt_zz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 34);

        auto tt_zz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 35);

        auto tt_xxx_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx);

        auto tt_xxx_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 1);

        auto tt_xxx_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 2);

        auto tt_xxy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 3);

        auto tt_xxy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 4);

        auto tt_xxy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 5);

        auto tt_xxz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 6);

        auto tt_xxz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 7);

        auto tt_xxz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 8);

        auto tt_xyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 9);

        auto tt_xyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 10);

        auto tt_xyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 11);

        auto tt_xyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 12);

        auto tt_xyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 13);

        auto tt_xyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 14);

        auto tt_xzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 15);

        auto tt_xzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 16);

        auto tt_xzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 17);

        auto tt_yyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 18);

        auto tt_yyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 19);

        auto tt_yyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 20);

        auto tt_yyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 21);

        auto tt_yyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 22);

        auto tt_yyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 23);

        auto ts_xxxx_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx);

        auto ts_xxxx_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 1);

        auto ts_xxxx_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 2);

        auto ts_xxxx_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 3);

        auto ts_xxxx_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 4);

        auto ts_xxxx_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 5);

        auto ts_xxxy_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 6);

        auto ts_xxxy_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 7);

        auto ts_xxxy_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 8);

        auto ts_xxxy_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 9);

        auto ts_xxxy_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 10);

        auto ts_xxxy_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 11);

        auto ts_xxxz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 12);

        auto ts_xxxz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 13);

        auto ts_xxxz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 14);

        auto ts_xxxz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 15);

        auto ts_xxxz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 16);

        auto ts_xxxz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 17);

        auto ts_xxyy_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 18);

        auto ts_xxyy_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 19);

        auto ts_xxyy_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 20);

        auto ts_xxyy_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 21);

        auto ts_xxyy_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 22);

        auto ts_xxyy_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 23);

        auto ts_xxyz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 24);

        auto ts_xxyz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 25);

        auto ts_xxyz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 26);

        auto ts_xxyz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 27);

        auto ts_xxyz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 28);

        auto ts_xxyz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 29);

        auto ts_xxzz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 30);

        auto ts_xxzz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 31);

        auto ts_xxzz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 32);

        auto ts_xxzz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 33);

        auto ts_xxzz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 34);

        auto ts_xxzz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 35);

        auto ts_xyyy_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 36);

        auto ts_xyyy_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 37);

        auto ts_xyyy_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 38);

        auto ts_xyyy_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 39);

        auto ts_xyyy_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 40);

        auto ts_xyyy_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 41);

        auto ts_xyyz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 42);

        auto ts_xyyz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 43);

        auto ts_xyyz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 44);

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

        auto ts_zz_xx_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 30);

        auto ts_zz_xy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 31);

        auto ts_zz_xz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 32);

        auto ts_zz_yy_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 33);

        auto ts_zz_yz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 34);

        auto ts_zz_zz_0 = primBuffer.data(pidx_s_2_2_m0 + 36 * idx + 35);

        // set up pointers to integrals

        auto tt_xxxx_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx);

        auto tt_xxxx_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 1);

        auto tt_xxxx_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 2);

        auto tt_xxxx_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 3);

        auto tt_xxxx_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 4);

        auto tt_xxxx_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 5);

        auto tt_xxxy_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 6);

        auto tt_xxxy_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 7);

        auto tt_xxxy_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 8);

        auto tt_xxxy_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 9);

        auto tt_xxxy_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 10);

        auto tt_xxxy_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 11);

        auto tt_xxxz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 12);

        auto tt_xxxz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 13);

        auto tt_xxxz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 14);

        auto tt_xxxz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 15);

        auto tt_xxxz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 16);

        auto tt_xxxz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 17);

        auto tt_xxyy_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 18);

        auto tt_xxyy_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 19);

        auto tt_xxyy_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 20);

        auto tt_xxyy_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 21);

        auto tt_xxyy_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 22);

        auto tt_xxyy_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 23);

        auto tt_xxyz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 24);

        auto tt_xxyz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 25);

        auto tt_xxyz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 26);

        auto tt_xxyz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 27);

        auto tt_xxyz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 28);

        auto tt_xxyz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 29);

        auto tt_xxzz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 30);

        auto tt_xxzz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 31);

        auto tt_xxzz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 32);

        auto tt_xxzz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 33);

        auto tt_xxzz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 34);

        auto tt_xxzz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 35);

        auto tt_xyyy_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 36);

        auto tt_xyyy_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 37);

        auto tt_xyyy_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 38);

        auto tt_xyyy_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 39);

        auto tt_xyyy_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 40);

        auto tt_xyyy_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 41);

        auto tt_xyyz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 42);

        auto tt_xyyz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 43);

        auto tt_xyyz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 44);

        // Batch of Integrals (0,45)

        #pragma omp simd aligned(fga, fx, fz, pa_x, ts_xx_xx_0, ts_xx_xy_0, ts_xx_xz_0, ts_xx_yy_0, \
                                     ts_xx_yz_0, ts_xx_zz_0, ts_xxxx_xx_0, ts_xxxx_xy_0, ts_xxxx_xz_0, ts_xxxx_yy_0, \
                                     ts_xxxx_yz_0, ts_xxxx_zz_0, ts_xxxy_xx_0, ts_xxxy_xy_0, ts_xxxy_xz_0, ts_xxxy_yy_0, \
                                     ts_xxxy_yz_0, ts_xxxy_zz_0, ts_xxxz_xx_0, ts_xxxz_xy_0, ts_xxxz_xz_0, ts_xxxz_yy_0, \
                                     ts_xxxz_yz_0, ts_xxxz_zz_0, ts_xxyy_xx_0, ts_xxyy_xy_0, ts_xxyy_xz_0, ts_xxyy_yy_0, \
                                     ts_xxyy_yz_0, ts_xxyy_zz_0, ts_xxyz_xx_0, ts_xxyz_xy_0, ts_xxyz_xz_0, ts_xxyz_yy_0, \
                                     ts_xxyz_yz_0, ts_xxyz_zz_0, ts_xxzz_xx_0, ts_xxzz_xy_0, ts_xxzz_xz_0, ts_xxzz_yy_0, \
                                     ts_xxzz_yz_0, ts_xxzz_zz_0, ts_xy_xx_0, ts_xy_xy_0, ts_xy_xz_0, ts_xy_yy_0, \
                                     ts_xy_yz_0, ts_xy_zz_0, ts_xyyy_xx_0, ts_xyyy_xy_0, ts_xyyy_xz_0, ts_xyyy_yy_0, \
                                     ts_xyyy_yz_0, ts_xyyy_zz_0, ts_xyyz_xx_0, ts_xyyz_xy_0, ts_xyyz_xz_0, ts_xz_xx_0, \
                                     ts_xz_xy_0, ts_xz_xz_0, ts_xz_yy_0, ts_xz_yz_0, ts_xz_zz_0, ts_yy_xx_0, ts_yy_xy_0, \
                                     ts_yy_xz_0, ts_yy_yy_0, ts_yy_yz_0, ts_yy_zz_0, ts_yz_xx_0, ts_yz_xy_0, ts_yz_xz_0, \
                                     ts_yz_yy_0, ts_yz_yz_0, ts_yz_zz_0, ts_zz_xx_0, ts_zz_xy_0, ts_zz_xz_0, ts_zz_yy_0, \
                                     ts_zz_yz_0, ts_zz_zz_0, tt_xx_xx_0, tt_xx_xy_0, tt_xx_xz_0, tt_xx_yy_0, tt_xx_yz_0, \
                                     tt_xx_zz_0, tt_xxx_x_0, tt_xxx_xx_0, tt_xxx_xy_0, tt_xxx_xz_0, tt_xxx_y_0, \
                                     tt_xxx_yy_0, tt_xxx_yz_0, tt_xxx_z_0, tt_xxx_zz_0, tt_xxxx_xx_0, tt_xxxx_xy_0, \
                                     tt_xxxx_xz_0, tt_xxxx_yy_0, tt_xxxx_yz_0, tt_xxxx_zz_0, tt_xxxy_xx_0, tt_xxxy_xy_0, \
                                     tt_xxxy_xz_0, tt_xxxy_yy_0, tt_xxxy_yz_0, tt_xxxy_zz_0, tt_xxxz_xx_0, tt_xxxz_xy_0, \
                                     tt_xxxz_xz_0, tt_xxxz_yy_0, tt_xxxz_yz_0, tt_xxxz_zz_0, tt_xxy_x_0, tt_xxy_xx_0, \
                                     tt_xxy_xy_0, tt_xxy_xz_0, tt_xxy_y_0, tt_xxy_yy_0, tt_xxy_yz_0, tt_xxy_z_0, \
                                     tt_xxy_zz_0, tt_xxyy_xx_0, tt_xxyy_xy_0, tt_xxyy_xz_0, tt_xxyy_yy_0, tt_xxyy_yz_0, \
                                     tt_xxyy_zz_0, tt_xxyz_xx_0, tt_xxyz_xy_0, tt_xxyz_xz_0, tt_xxyz_yy_0, tt_xxyz_yz_0, \
                                     tt_xxyz_zz_0, tt_xxz_x_0, tt_xxz_xx_0, tt_xxz_xy_0, tt_xxz_xz_0, tt_xxz_y_0, \
                                     tt_xxz_yy_0, tt_xxz_yz_0, tt_xxz_z_0, tt_xxz_zz_0, tt_xxzz_xx_0, tt_xxzz_xy_0, \
                                     tt_xxzz_xz_0, tt_xxzz_yy_0, tt_xxzz_yz_0, tt_xxzz_zz_0, tt_xy_xx_0, tt_xy_xy_0, \
                                     tt_xy_xz_0, tt_xy_yy_0, tt_xy_yz_0, tt_xy_zz_0, tt_xyy_x_0, tt_xyy_xx_0, \
                                     tt_xyy_xy_0, tt_xyy_xz_0, tt_xyy_y_0, tt_xyy_yy_0, tt_xyy_yz_0, tt_xyy_z_0, \
                                     tt_xyy_zz_0, tt_xyyy_xx_0, tt_xyyy_xy_0, tt_xyyy_xz_0, tt_xyyy_yy_0, tt_xyyy_yz_0, \
                                     tt_xyyy_zz_0, tt_xyyz_xx_0, tt_xyyz_xy_0, tt_xyyz_xz_0, tt_xyz_x_0, tt_xyz_xx_0, \
                                     tt_xyz_xy_0, tt_xyz_xz_0, tt_xyz_y_0, tt_xyz_yy_0, tt_xyz_yz_0, tt_xyz_z_0, \
                                     tt_xyz_zz_0, tt_xz_xx_0, tt_xz_xy_0, tt_xz_xz_0, tt_xz_yy_0, tt_xz_yz_0, tt_xz_zz_0, \
                                     tt_xzz_x_0, tt_xzz_xx_0, tt_xzz_xy_0, tt_xzz_xz_0, tt_xzz_y_0, tt_xzz_yy_0, \
                                     tt_xzz_yz_0, tt_xzz_z_0, tt_xzz_zz_0, tt_yy_xx_0, tt_yy_xy_0, tt_yy_xz_0, \
                                     tt_yy_yy_0, tt_yy_yz_0, tt_yy_zz_0, tt_yyy_x_0, tt_yyy_xx_0, tt_yyy_xy_0, \
                                     tt_yyy_xz_0, tt_yyy_y_0, tt_yyy_yy_0, tt_yyy_yz_0, tt_yyy_z_0, tt_yyy_zz_0, \
                                     tt_yyz_x_0, tt_yyz_xx_0, tt_yyz_xy_0, tt_yyz_xz_0, tt_yyz_y_0, tt_yyz_z_0, \
                                     tt_yz_xx_0, tt_yz_xy_0, tt_yz_xz_0, tt_yz_yy_0, tt_yz_yz_0, tt_yz_zz_0, tt_zz_xx_0, \
                                     tt_zz_xy_0, tt_zz_xz_0, tt_zz_yy_0, tt_zz_yz_0, tt_zz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xxxx_xx_0[j] = pa_x[j] * tt_xxx_xx_0[j] + 1.5 * fl1_fx * tt_xx_xx_0[j] + fl1_fx * tt_xxx_x_0[j] + 2.0 * fl1_fz * ts_xxxx_xx_0[j] -
                              3.0 * fl1_fz * fl1_fga * ts_xx_xx_0[j];

            tt_xxxx_xy_0[j] = pa_x[j] * tt_xxx_xy_0[j] + 1.5 * fl1_fx * tt_xx_xy_0[j] + 0.5 * fl1_fx * tt_xxx_y_0[j] +
                              2.0 * fl1_fz * ts_xxxx_xy_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_xy_0[j];

            tt_xxxx_xz_0[j] = pa_x[j] * tt_xxx_xz_0[j] + 1.5 * fl1_fx * tt_xx_xz_0[j] + 0.5 * fl1_fx * tt_xxx_z_0[j] +
                              2.0 * fl1_fz * ts_xxxx_xz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_xz_0[j];

            tt_xxxx_yy_0[j] =
                pa_x[j] * tt_xxx_yy_0[j] + 1.5 * fl1_fx * tt_xx_yy_0[j] + 2.0 * fl1_fz * ts_xxxx_yy_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_yy_0[j];

            tt_xxxx_yz_0[j] =
                pa_x[j] * tt_xxx_yz_0[j] + 1.5 * fl1_fx * tt_xx_yz_0[j] + 2.0 * fl1_fz * ts_xxxx_yz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_yz_0[j];

            tt_xxxx_zz_0[j] =
                pa_x[j] * tt_xxx_zz_0[j] + 1.5 * fl1_fx * tt_xx_zz_0[j] + 2.0 * fl1_fz * ts_xxxx_zz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_xx_zz_0[j];

            tt_xxxy_xx_0[j] = pa_x[j] * tt_xxy_xx_0[j] + fl1_fx * tt_xy_xx_0[j] + fl1_fx * tt_xxy_x_0[j] + 2.0 * fl1_fz * ts_xxxy_xx_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_xy_xx_0[j];

            tt_xxxy_xy_0[j] = pa_x[j] * tt_xxy_xy_0[j] + fl1_fx * tt_xy_xy_0[j] + 0.5 * fl1_fx * tt_xxy_y_0[j] + 2.0 * fl1_fz * ts_xxxy_xy_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_xy_xy_0[j];

            tt_xxxy_xz_0[j] = pa_x[j] * tt_xxy_xz_0[j] + fl1_fx * tt_xy_xz_0[j] + 0.5 * fl1_fx * tt_xxy_z_0[j] + 2.0 * fl1_fz * ts_xxxy_xz_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_xy_xz_0[j];

            tt_xxxy_yy_0[j] =
                pa_x[j] * tt_xxy_yy_0[j] + fl1_fx * tt_xy_yy_0[j] + 2.0 * fl1_fz * ts_xxxy_yy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_yy_0[j];

            tt_xxxy_yz_0[j] =
                pa_x[j] * tt_xxy_yz_0[j] + fl1_fx * tt_xy_yz_0[j] + 2.0 * fl1_fz * ts_xxxy_yz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_yz_0[j];

            tt_xxxy_zz_0[j] =
                pa_x[j] * tt_xxy_zz_0[j] + fl1_fx * tt_xy_zz_0[j] + 2.0 * fl1_fz * ts_xxxy_zz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xy_zz_0[j];

            tt_xxxz_xx_0[j] = pa_x[j] * tt_xxz_xx_0[j] + fl1_fx * tt_xz_xx_0[j] + fl1_fx * tt_xxz_x_0[j] + 2.0 * fl1_fz * ts_xxxz_xx_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_xz_xx_0[j];

            tt_xxxz_xy_0[j] = pa_x[j] * tt_xxz_xy_0[j] + fl1_fx * tt_xz_xy_0[j] + 0.5 * fl1_fx * tt_xxz_y_0[j] + 2.0 * fl1_fz * ts_xxxz_xy_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_xz_xy_0[j];

            tt_xxxz_xz_0[j] = pa_x[j] * tt_xxz_xz_0[j] + fl1_fx * tt_xz_xz_0[j] + 0.5 * fl1_fx * tt_xxz_z_0[j] + 2.0 * fl1_fz * ts_xxxz_xz_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_xz_xz_0[j];

            tt_xxxz_yy_0[j] =
                pa_x[j] * tt_xxz_yy_0[j] + fl1_fx * tt_xz_yy_0[j] + 2.0 * fl1_fz * ts_xxxz_yy_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_yy_0[j];

            tt_xxxz_yz_0[j] =
                pa_x[j] * tt_xxz_yz_0[j] + fl1_fx * tt_xz_yz_0[j] + 2.0 * fl1_fz * ts_xxxz_yz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_yz_0[j];

            tt_xxxz_zz_0[j] =
                pa_x[j] * tt_xxz_zz_0[j] + fl1_fx * tt_xz_zz_0[j] + 2.0 * fl1_fz * ts_xxxz_zz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_xz_zz_0[j];

            tt_xxyy_xx_0[j] = pa_x[j] * tt_xyy_xx_0[j] + 0.5 * fl1_fx * tt_yy_xx_0[j] + fl1_fx * tt_xyy_x_0[j] + 2.0 * fl1_fz * ts_xxyy_xx_0[j] -
                              fl1_fz * fl1_fga * ts_yy_xx_0[j];

            tt_xxyy_xy_0[j] = pa_x[j] * tt_xyy_xy_0[j] + 0.5 * fl1_fx * tt_yy_xy_0[j] + 0.5 * fl1_fx * tt_xyy_y_0[j] +
                              2.0 * fl1_fz * ts_xxyy_xy_0[j] - fl1_fz * fl1_fga * ts_yy_xy_0[j];

            tt_xxyy_xz_0[j] = pa_x[j] * tt_xyy_xz_0[j] + 0.5 * fl1_fx * tt_yy_xz_0[j] + 0.5 * fl1_fx * tt_xyy_z_0[j] +
                              2.0 * fl1_fz * ts_xxyy_xz_0[j] - fl1_fz * fl1_fga * ts_yy_xz_0[j];

            tt_xxyy_yy_0[j] =
                pa_x[j] * tt_xyy_yy_0[j] + 0.5 * fl1_fx * tt_yy_yy_0[j] + 2.0 * fl1_fz * ts_xxyy_yy_0[j] - fl1_fz * fl1_fga * ts_yy_yy_0[j];

            tt_xxyy_yz_0[j] =
                pa_x[j] * tt_xyy_yz_0[j] + 0.5 * fl1_fx * tt_yy_yz_0[j] + 2.0 * fl1_fz * ts_xxyy_yz_0[j] - fl1_fz * fl1_fga * ts_yy_yz_0[j];

            tt_xxyy_zz_0[j] =
                pa_x[j] * tt_xyy_zz_0[j] + 0.5 * fl1_fx * tt_yy_zz_0[j] + 2.0 * fl1_fz * ts_xxyy_zz_0[j] - fl1_fz * fl1_fga * ts_yy_zz_0[j];

            tt_xxyz_xx_0[j] = pa_x[j] * tt_xyz_xx_0[j] + 0.5 * fl1_fx * tt_yz_xx_0[j] + fl1_fx * tt_xyz_x_0[j] + 2.0 * fl1_fz * ts_xxyz_xx_0[j] -
                              fl1_fz * fl1_fga * ts_yz_xx_0[j];

            tt_xxyz_xy_0[j] = pa_x[j] * tt_xyz_xy_0[j] + 0.5 * fl1_fx * tt_yz_xy_0[j] + 0.5 * fl1_fx * tt_xyz_y_0[j] +
                              2.0 * fl1_fz * ts_xxyz_xy_0[j] - fl1_fz * fl1_fga * ts_yz_xy_0[j];

            tt_xxyz_xz_0[j] = pa_x[j] * tt_xyz_xz_0[j] + 0.5 * fl1_fx * tt_yz_xz_0[j] + 0.5 * fl1_fx * tt_xyz_z_0[j] +
                              2.0 * fl1_fz * ts_xxyz_xz_0[j] - fl1_fz * fl1_fga * ts_yz_xz_0[j];

            tt_xxyz_yy_0[j] =
                pa_x[j] * tt_xyz_yy_0[j] + 0.5 * fl1_fx * tt_yz_yy_0[j] + 2.0 * fl1_fz * ts_xxyz_yy_0[j] - fl1_fz * fl1_fga * ts_yz_yy_0[j];

            tt_xxyz_yz_0[j] =
                pa_x[j] * tt_xyz_yz_0[j] + 0.5 * fl1_fx * tt_yz_yz_0[j] + 2.0 * fl1_fz * ts_xxyz_yz_0[j] - fl1_fz * fl1_fga * ts_yz_yz_0[j];

            tt_xxyz_zz_0[j] =
                pa_x[j] * tt_xyz_zz_0[j] + 0.5 * fl1_fx * tt_yz_zz_0[j] + 2.0 * fl1_fz * ts_xxyz_zz_0[j] - fl1_fz * fl1_fga * ts_yz_zz_0[j];

            tt_xxzz_xx_0[j] = pa_x[j] * tt_xzz_xx_0[j] + 0.5 * fl1_fx * tt_zz_xx_0[j] + fl1_fx * tt_xzz_x_0[j] + 2.0 * fl1_fz * ts_xxzz_xx_0[j] -
                              fl1_fz * fl1_fga * ts_zz_xx_0[j];

            tt_xxzz_xy_0[j] = pa_x[j] * tt_xzz_xy_0[j] + 0.5 * fl1_fx * tt_zz_xy_0[j] + 0.5 * fl1_fx * tt_xzz_y_0[j] +
                              2.0 * fl1_fz * ts_xxzz_xy_0[j] - fl1_fz * fl1_fga * ts_zz_xy_0[j];

            tt_xxzz_xz_0[j] = pa_x[j] * tt_xzz_xz_0[j] + 0.5 * fl1_fx * tt_zz_xz_0[j] + 0.5 * fl1_fx * tt_xzz_z_0[j] +
                              2.0 * fl1_fz * ts_xxzz_xz_0[j] - fl1_fz * fl1_fga * ts_zz_xz_0[j];

            tt_xxzz_yy_0[j] =
                pa_x[j] * tt_xzz_yy_0[j] + 0.5 * fl1_fx * tt_zz_yy_0[j] + 2.0 * fl1_fz * ts_xxzz_yy_0[j] - fl1_fz * fl1_fga * ts_zz_yy_0[j];

            tt_xxzz_yz_0[j] =
                pa_x[j] * tt_xzz_yz_0[j] + 0.5 * fl1_fx * tt_zz_yz_0[j] + 2.0 * fl1_fz * ts_xxzz_yz_0[j] - fl1_fz * fl1_fga * ts_zz_yz_0[j];

            tt_xxzz_zz_0[j] =
                pa_x[j] * tt_xzz_zz_0[j] + 0.5 * fl1_fx * tt_zz_zz_0[j] + 2.0 * fl1_fz * ts_xxzz_zz_0[j] - fl1_fz * fl1_fga * ts_zz_zz_0[j];

            tt_xyyy_xx_0[j] = pa_x[j] * tt_yyy_xx_0[j] + fl1_fx * tt_yyy_x_0[j] + 2.0 * fl1_fz * ts_xyyy_xx_0[j];

            tt_xyyy_xy_0[j] = pa_x[j] * tt_yyy_xy_0[j] + 0.5 * fl1_fx * tt_yyy_y_0[j] + 2.0 * fl1_fz * ts_xyyy_xy_0[j];

            tt_xyyy_xz_0[j] = pa_x[j] * tt_yyy_xz_0[j] + 0.5 * fl1_fx * tt_yyy_z_0[j] + 2.0 * fl1_fz * ts_xyyy_xz_0[j];

            tt_xyyy_yy_0[j] = pa_x[j] * tt_yyy_yy_0[j] + 2.0 * fl1_fz * ts_xyyy_yy_0[j];

            tt_xyyy_yz_0[j] = pa_x[j] * tt_yyy_yz_0[j] + 2.0 * fl1_fz * ts_xyyy_yz_0[j];

            tt_xyyy_zz_0[j] = pa_x[j] * tt_yyy_zz_0[j] + 2.0 * fl1_fz * ts_xyyy_zz_0[j];

            tt_xyyz_xx_0[j] = pa_x[j] * tt_yyz_xx_0[j] + fl1_fx * tt_yyz_x_0[j] + 2.0 * fl1_fz * ts_xyyz_xx_0[j];

            tt_xyyz_xy_0[j] = pa_x[j] * tt_yyz_xy_0[j] + 0.5 * fl1_fx * tt_yyz_y_0[j] + 2.0 * fl1_fz * ts_xyyz_xy_0[j];

            tt_xyyz_xz_0[j] = pa_x[j] * tt_yyz_xz_0[j] + 0.5 * fl1_fx * tt_yyz_z_0[j] + 2.0 * fl1_fz * ts_xyyz_xz_0[j];
        }

        idx++;
    }
}

void
compKineticEnergyForGD_45_90(CMemBlock2D<double>&       primBuffer,
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

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    // set up index of integral

    auto pidx_t_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    // check if integral is needed in recursion expansion

    if (pidx_t_4_2_m0 == -1) return;

    // set up indexes of auxilary integral

    auto pidx_t_3_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_4_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {4, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_3_1_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {3, -1, -1, -1}, {1, -1, -1, -1}, 1, 1, 0));

    auto pidx_t_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Kinetic Energy"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

    auto pidx_s_2_2_m0 = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true, {2, -1, -1, -1}, {2, -1, -1, -1}, 1, 1, 0));

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

        auto pa_z = paDistances.data(3 * idx + 2);

        // set up pointers to auxilary integrals

        auto tt_yyy_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 36);

        auto tt_yyy_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 37);

        auto tt_yyy_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 38);

        auto tt_yyy_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 39);

        auto tt_yyy_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 40);

        auto tt_yyy_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 41);

        auto tt_yyz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 42);

        auto tt_yyz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 43);

        auto tt_yyz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 44);

        auto tt_yyz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 45);

        auto tt_yyz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 46);

        auto tt_yyz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 47);

        auto tt_yzz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 48);

        auto tt_yzz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 49);

        auto tt_yzz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 50);

        auto tt_yzz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 51);

        auto tt_yzz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 52);

        auto tt_yzz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 53);

        auto tt_zzz_xx_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 54);

        auto tt_zzz_xy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 55);

        auto tt_zzz_xz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 56);

        auto tt_zzz_yy_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 57);

        auto tt_zzz_yz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 58);

        auto tt_zzz_zz_0 = primBuffer.data(pidx_t_3_2_m0 + 60 * idx + 59);

        auto tt_yy_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 18);

        auto tt_yy_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 19);

        auto tt_yy_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 20);

        auto tt_yy_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 21);

        auto tt_yy_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 22);

        auto tt_yy_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 23);

        auto tt_yz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 24);

        auto tt_yz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 25);

        auto tt_yz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 26);

        auto tt_yz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 27);

        auto tt_yz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 28);

        auto tt_yz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 29);

        auto tt_zz_xx_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 30);

        auto tt_zz_xy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 31);

        auto tt_zz_xz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 32);

        auto tt_zz_yy_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 33);

        auto tt_zz_yz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 34);

        auto tt_zz_zz_0 = primBuffer.data(pidx_t_2_2_m0 + 36 * idx + 35);

        auto tt_yyy_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 18);

        auto tt_yyy_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 19);

        auto tt_yyy_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 20);

        auto tt_yyz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 21);

        auto tt_yyz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 22);

        auto tt_yyz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 23);

        auto tt_yzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 24);

        auto tt_yzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 25);

        auto tt_yzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 26);

        auto tt_zzz_x_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 27);

        auto tt_zzz_y_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 28);

        auto tt_zzz_z_0 = primBuffer.data(pidx_t_3_1_m0 + 30 * idx + 29);

        auto ts_xyyz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 45);

        auto ts_xyyz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 46);

        auto ts_xyyz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 47);

        auto ts_xyzz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 48);

        auto ts_xyzz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 49);

        auto ts_xyzz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 50);

        auto ts_xyzz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 51);

        auto ts_xyzz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 52);

        auto ts_xyzz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 53);

        auto ts_xzzz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 54);

        auto ts_xzzz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 55);

        auto ts_xzzz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 56);

        auto ts_xzzz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 57);

        auto ts_xzzz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 58);

        auto ts_xzzz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 59);

        auto ts_yyyy_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 60);

        auto ts_yyyy_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 61);

        auto ts_yyyy_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 62);

        auto ts_yyyy_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 63);

        auto ts_yyyy_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 64);

        auto ts_yyyy_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 65);

        auto ts_yyyz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 66);

        auto ts_yyyz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 67);

        auto ts_yyyz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 68);

        auto ts_yyyz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 69);

        auto ts_yyyz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 70);

        auto ts_yyyz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 71);

        auto ts_yyzz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 72);

        auto ts_yyzz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 73);

        auto ts_yyzz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 74);

        auto ts_yyzz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 75);

        auto ts_yyzz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 76);

        auto ts_yyzz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 77);

        auto ts_yzzz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 78);

        auto ts_yzzz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 79);

        auto ts_yzzz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 80);

        auto ts_yzzz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 81);

        auto ts_yzzz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 82);

        auto ts_yzzz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 83);

        auto ts_zzzz_xx_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 84);

        auto ts_zzzz_xy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 85);

        auto ts_zzzz_xz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 86);

        auto ts_zzzz_yy_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 87);

        auto ts_zzzz_yz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 88);

        auto ts_zzzz_zz_0 = primBuffer.data(pidx_s_4_2_m0 + 90 * idx + 89);

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

        auto tt_xyyz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 45);

        auto tt_xyyz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 46);

        auto tt_xyyz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 47);

        auto tt_xyzz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 48);

        auto tt_xyzz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 49);

        auto tt_xyzz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 50);

        auto tt_xyzz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 51);

        auto tt_xyzz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 52);

        auto tt_xyzz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 53);

        auto tt_xzzz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 54);

        auto tt_xzzz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 55);

        auto tt_xzzz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 56);

        auto tt_xzzz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 57);

        auto tt_xzzz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 58);

        auto tt_xzzz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 59);

        auto tt_yyyy_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 60);

        auto tt_yyyy_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 61);

        auto tt_yyyy_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 62);

        auto tt_yyyy_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 63);

        auto tt_yyyy_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 64);

        auto tt_yyyy_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 65);

        auto tt_yyyz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 66);

        auto tt_yyyz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 67);

        auto tt_yyyz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 68);

        auto tt_yyyz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 69);

        auto tt_yyyz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 70);

        auto tt_yyyz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 71);

        auto tt_yyzz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 72);

        auto tt_yyzz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 73);

        auto tt_yyzz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 74);

        auto tt_yyzz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 75);

        auto tt_yyzz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 76);

        auto tt_yyzz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 77);

        auto tt_yzzz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 78);

        auto tt_yzzz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 79);

        auto tt_yzzz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 80);

        auto tt_yzzz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 81);

        auto tt_yzzz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 82);

        auto tt_yzzz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 83);

        auto tt_zzzz_xx_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 84);

        auto tt_zzzz_xy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 85);

        auto tt_zzzz_xz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 86);

        auto tt_zzzz_yy_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 87);

        auto tt_zzzz_yz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 88);

        auto tt_zzzz_zz_0 = primBuffer.data(pidx_t_4_2_m0 + 90 * idx + 89);

        // Batch of Integrals (45,90)

        #pragma omp simd aligned(fga, fx, fz, pa_x, pa_y, pa_z, ts_xyyz_yy_0, ts_xyyz_yz_0, ts_xyyz_zz_0, \
                                     ts_xyzz_xx_0, ts_xyzz_xy_0, ts_xyzz_xz_0, ts_xyzz_yy_0, ts_xyzz_yz_0, ts_xyzz_zz_0, \
                                     ts_xzzz_xx_0, ts_xzzz_xy_0, ts_xzzz_xz_0, ts_xzzz_yy_0, ts_xzzz_yz_0, ts_xzzz_zz_0, \
                                     ts_yy_xx_0, ts_yy_xy_0, ts_yy_xz_0, ts_yy_yy_0, ts_yy_yz_0, ts_yy_zz_0, \
                                     ts_yyyy_xx_0, ts_yyyy_xy_0, ts_yyyy_xz_0, ts_yyyy_yy_0, ts_yyyy_yz_0, ts_yyyy_zz_0, \
                                     ts_yyyz_xx_0, ts_yyyz_xy_0, ts_yyyz_xz_0, ts_yyyz_yy_0, ts_yyyz_yz_0, ts_yyyz_zz_0, \
                                     ts_yyzz_xx_0, ts_yyzz_xy_0, ts_yyzz_xz_0, ts_yyzz_yy_0, ts_yyzz_yz_0, ts_yyzz_zz_0, \
                                     ts_yz_xx_0, ts_yz_xy_0, ts_yz_xz_0, ts_yz_yy_0, ts_yz_yz_0, ts_yz_zz_0, \
                                     ts_yzzz_xx_0, ts_yzzz_xy_0, ts_yzzz_xz_0, ts_yzzz_yy_0, ts_yzzz_yz_0, ts_yzzz_zz_0, \
                                     ts_zz_xx_0, ts_zz_xy_0, ts_zz_xz_0, ts_zz_yy_0, ts_zz_yz_0, ts_zz_zz_0, \
                                     ts_zzzz_xx_0, ts_zzzz_xy_0, ts_zzzz_xz_0, ts_zzzz_yy_0, ts_zzzz_yz_0, ts_zzzz_zz_0, \
                                     tt_xyyz_yy_0, tt_xyyz_yz_0, tt_xyyz_zz_0, tt_xyzz_xx_0, tt_xyzz_xy_0, tt_xyzz_xz_0, \
                                     tt_xyzz_yy_0, tt_xyzz_yz_0, tt_xyzz_zz_0, tt_xzzz_xx_0, tt_xzzz_xy_0, tt_xzzz_xz_0, \
                                     tt_xzzz_yy_0, tt_xzzz_yz_0, tt_xzzz_zz_0, tt_yy_xx_0, tt_yy_xy_0, tt_yy_xz_0, \
                                     tt_yy_yy_0, tt_yy_yz_0, tt_yy_zz_0, tt_yyy_x_0, tt_yyy_xx_0, tt_yyy_xy_0, \
                                     tt_yyy_xz_0, tt_yyy_y_0, tt_yyy_yy_0, tt_yyy_yz_0, tt_yyy_z_0, tt_yyy_zz_0, \
                                     tt_yyyy_xx_0, tt_yyyy_xy_0, tt_yyyy_xz_0, tt_yyyy_yy_0, tt_yyyy_yz_0, tt_yyyy_zz_0, \
                                     tt_yyyz_xx_0, tt_yyyz_xy_0, tt_yyyz_xz_0, tt_yyyz_yy_0, tt_yyyz_yz_0, tt_yyyz_zz_0, \
                                     tt_yyz_x_0, tt_yyz_xx_0, tt_yyz_xy_0, tt_yyz_xz_0, tt_yyz_y_0, tt_yyz_yy_0, \
                                     tt_yyz_yz_0, tt_yyz_z_0, tt_yyz_zz_0, tt_yyzz_xx_0, tt_yyzz_xy_0, tt_yyzz_xz_0, \
                                     tt_yyzz_yy_0, tt_yyzz_yz_0, tt_yyzz_zz_0, tt_yz_xx_0, tt_yz_xy_0, tt_yz_xz_0, \
                                     tt_yz_yy_0, tt_yz_yz_0, tt_yz_zz_0, tt_yzz_x_0, tt_yzz_xx_0, tt_yzz_xy_0, \
                                     tt_yzz_xz_0, tt_yzz_y_0, tt_yzz_yy_0, tt_yzz_yz_0, tt_yzz_z_0, tt_yzz_zz_0, \
                                     tt_yzzz_xx_0, tt_yzzz_xy_0, tt_yzzz_xz_0, tt_yzzz_yy_0, tt_yzzz_yz_0, tt_yzzz_zz_0, \
                                     tt_zz_xx_0, tt_zz_xy_0, tt_zz_xz_0, tt_zz_yy_0, tt_zz_yz_0, tt_zz_zz_0, tt_zzz_x_0, \
                                     tt_zzz_xx_0, tt_zzz_xy_0, tt_zzz_xz_0, tt_zzz_y_0, tt_zzz_yy_0, tt_zzz_yz_0, \
                                     tt_zzz_z_0, tt_zzz_zz_0, tt_zzzz_xx_0, tt_zzzz_xy_0, tt_zzzz_xz_0, tt_zzzz_yy_0, \
                                     tt_zzzz_yz_0, tt_zzzz_zz_0: VLX_ALIGN)
        for (int32_t j = 0; j < nprim; j++)
        {
            double fl1_fga = fga[j];

            double fl1_fx = fx[j];

            double fl1_fz = fz[j];

            tt_xyyz_yy_0[j] = pa_x[j] * tt_yyz_yy_0[j] + 2.0 * fl1_fz * ts_xyyz_yy_0[j];

            tt_xyyz_yz_0[j] = pa_x[j] * tt_yyz_yz_0[j] + 2.0 * fl1_fz * ts_xyyz_yz_0[j];

            tt_xyyz_zz_0[j] = pa_x[j] * tt_yyz_zz_0[j] + 2.0 * fl1_fz * ts_xyyz_zz_0[j];

            tt_xyzz_xx_0[j] = pa_x[j] * tt_yzz_xx_0[j] + fl1_fx * tt_yzz_x_0[j] + 2.0 * fl1_fz * ts_xyzz_xx_0[j];

            tt_xyzz_xy_0[j] = pa_x[j] * tt_yzz_xy_0[j] + 0.5 * fl1_fx * tt_yzz_y_0[j] + 2.0 * fl1_fz * ts_xyzz_xy_0[j];

            tt_xyzz_xz_0[j] = pa_x[j] * tt_yzz_xz_0[j] + 0.5 * fl1_fx * tt_yzz_z_0[j] + 2.0 * fl1_fz * ts_xyzz_xz_0[j];

            tt_xyzz_yy_0[j] = pa_x[j] * tt_yzz_yy_0[j] + 2.0 * fl1_fz * ts_xyzz_yy_0[j];

            tt_xyzz_yz_0[j] = pa_x[j] * tt_yzz_yz_0[j] + 2.0 * fl1_fz * ts_xyzz_yz_0[j];

            tt_xyzz_zz_0[j] = pa_x[j] * tt_yzz_zz_0[j] + 2.0 * fl1_fz * ts_xyzz_zz_0[j];

            tt_xzzz_xx_0[j] = pa_x[j] * tt_zzz_xx_0[j] + fl1_fx * tt_zzz_x_0[j] + 2.0 * fl1_fz * ts_xzzz_xx_0[j];

            tt_xzzz_xy_0[j] = pa_x[j] * tt_zzz_xy_0[j] + 0.5 * fl1_fx * tt_zzz_y_0[j] + 2.0 * fl1_fz * ts_xzzz_xy_0[j];

            tt_xzzz_xz_0[j] = pa_x[j] * tt_zzz_xz_0[j] + 0.5 * fl1_fx * tt_zzz_z_0[j] + 2.0 * fl1_fz * ts_xzzz_xz_0[j];

            tt_xzzz_yy_0[j] = pa_x[j] * tt_zzz_yy_0[j] + 2.0 * fl1_fz * ts_xzzz_yy_0[j];

            tt_xzzz_yz_0[j] = pa_x[j] * tt_zzz_yz_0[j] + 2.0 * fl1_fz * ts_xzzz_yz_0[j];

            tt_xzzz_zz_0[j] = pa_x[j] * tt_zzz_zz_0[j] + 2.0 * fl1_fz * ts_xzzz_zz_0[j];

            tt_yyyy_xx_0[j] =
                pa_y[j] * tt_yyy_xx_0[j] + 1.5 * fl1_fx * tt_yy_xx_0[j] + 2.0 * fl1_fz * ts_yyyy_xx_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_xx_0[j];

            tt_yyyy_xy_0[j] = pa_y[j] * tt_yyy_xy_0[j] + 1.5 * fl1_fx * tt_yy_xy_0[j] + 0.5 * fl1_fx * tt_yyy_x_0[j] +
                              2.0 * fl1_fz * ts_yyyy_xy_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_xy_0[j];

            tt_yyyy_xz_0[j] =
                pa_y[j] * tt_yyy_xz_0[j] + 1.5 * fl1_fx * tt_yy_xz_0[j] + 2.0 * fl1_fz * ts_yyyy_xz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_xz_0[j];

            tt_yyyy_yy_0[j] = pa_y[j] * tt_yyy_yy_0[j] + 1.5 * fl1_fx * tt_yy_yy_0[j] + fl1_fx * tt_yyy_y_0[j] + 2.0 * fl1_fz * ts_yyyy_yy_0[j] -
                              3.0 * fl1_fz * fl1_fga * ts_yy_yy_0[j];

            tt_yyyy_yz_0[j] = pa_y[j] * tt_yyy_yz_0[j] + 1.5 * fl1_fx * tt_yy_yz_0[j] + 0.5 * fl1_fx * tt_yyy_z_0[j] +
                              2.0 * fl1_fz * ts_yyyy_yz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_yz_0[j];

            tt_yyyy_zz_0[j] =
                pa_y[j] * tt_yyy_zz_0[j] + 1.5 * fl1_fx * tt_yy_zz_0[j] + 2.0 * fl1_fz * ts_yyyy_zz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_yy_zz_0[j];

            tt_yyyz_xx_0[j] =
                pa_y[j] * tt_yyz_xx_0[j] + fl1_fx * tt_yz_xx_0[j] + 2.0 * fl1_fz * ts_yyyz_xx_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_xx_0[j];

            tt_yyyz_xy_0[j] = pa_y[j] * tt_yyz_xy_0[j] + fl1_fx * tt_yz_xy_0[j] + 0.5 * fl1_fx * tt_yyz_x_0[j] + 2.0 * fl1_fz * ts_yyyz_xy_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_yz_xy_0[j];

            tt_yyyz_xz_0[j] =
                pa_y[j] * tt_yyz_xz_0[j] + fl1_fx * tt_yz_xz_0[j] + 2.0 * fl1_fz * ts_yyyz_xz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_xz_0[j];

            tt_yyyz_yy_0[j] = pa_y[j] * tt_yyz_yy_0[j] + fl1_fx * tt_yz_yy_0[j] + fl1_fx * tt_yyz_y_0[j] + 2.0 * fl1_fz * ts_yyyz_yy_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_yz_yy_0[j];

            tt_yyyz_yz_0[j] = pa_y[j] * tt_yyz_yz_0[j] + fl1_fx * tt_yz_yz_0[j] + 0.5 * fl1_fx * tt_yyz_z_0[j] + 2.0 * fl1_fz * ts_yyyz_yz_0[j] -
                              2.0 * fl1_fz * fl1_fga * ts_yz_yz_0[j];

            tt_yyyz_zz_0[j] =
                pa_y[j] * tt_yyz_zz_0[j] + fl1_fx * tt_yz_zz_0[j] + 2.0 * fl1_fz * ts_yyyz_zz_0[j] - 2.0 * fl1_fz * fl1_fga * ts_yz_zz_0[j];

            tt_yyzz_xx_0[j] =
                pa_y[j] * tt_yzz_xx_0[j] + 0.5 * fl1_fx * tt_zz_xx_0[j] + 2.0 * fl1_fz * ts_yyzz_xx_0[j] - fl1_fz * fl1_fga * ts_zz_xx_0[j];

            tt_yyzz_xy_0[j] = pa_y[j] * tt_yzz_xy_0[j] + 0.5 * fl1_fx * tt_zz_xy_0[j] + 0.5 * fl1_fx * tt_yzz_x_0[j] +
                              2.0 * fl1_fz * ts_yyzz_xy_0[j] - fl1_fz * fl1_fga * ts_zz_xy_0[j];

            tt_yyzz_xz_0[j] =
                pa_y[j] * tt_yzz_xz_0[j] + 0.5 * fl1_fx * tt_zz_xz_0[j] + 2.0 * fl1_fz * ts_yyzz_xz_0[j] - fl1_fz * fl1_fga * ts_zz_xz_0[j];

            tt_yyzz_yy_0[j] = pa_y[j] * tt_yzz_yy_0[j] + 0.5 * fl1_fx * tt_zz_yy_0[j] + fl1_fx * tt_yzz_y_0[j] + 2.0 * fl1_fz * ts_yyzz_yy_0[j] -
                              fl1_fz * fl1_fga * ts_zz_yy_0[j];

            tt_yyzz_yz_0[j] = pa_y[j] * tt_yzz_yz_0[j] + 0.5 * fl1_fx * tt_zz_yz_0[j] + 0.5 * fl1_fx * tt_yzz_z_0[j] +
                              2.0 * fl1_fz * ts_yyzz_yz_0[j] - fl1_fz * fl1_fga * ts_zz_yz_0[j];

            tt_yyzz_zz_0[j] =
                pa_y[j] * tt_yzz_zz_0[j] + 0.5 * fl1_fx * tt_zz_zz_0[j] + 2.0 * fl1_fz * ts_yyzz_zz_0[j] - fl1_fz * fl1_fga * ts_zz_zz_0[j];

            tt_yzzz_xx_0[j] = pa_y[j] * tt_zzz_xx_0[j] + 2.0 * fl1_fz * ts_yzzz_xx_0[j];

            tt_yzzz_xy_0[j] = pa_y[j] * tt_zzz_xy_0[j] + 0.5 * fl1_fx * tt_zzz_x_0[j] + 2.0 * fl1_fz * ts_yzzz_xy_0[j];

            tt_yzzz_xz_0[j] = pa_y[j] * tt_zzz_xz_0[j] + 2.0 * fl1_fz * ts_yzzz_xz_0[j];

            tt_yzzz_yy_0[j] = pa_y[j] * tt_zzz_yy_0[j] + fl1_fx * tt_zzz_y_0[j] + 2.0 * fl1_fz * ts_yzzz_yy_0[j];

            tt_yzzz_yz_0[j] = pa_y[j] * tt_zzz_yz_0[j] + 0.5 * fl1_fx * tt_zzz_z_0[j] + 2.0 * fl1_fz * ts_yzzz_yz_0[j];

            tt_yzzz_zz_0[j] = pa_y[j] * tt_zzz_zz_0[j] + 2.0 * fl1_fz * ts_yzzz_zz_0[j];

            tt_zzzz_xx_0[j] =
                pa_z[j] * tt_zzz_xx_0[j] + 1.5 * fl1_fx * tt_zz_xx_0[j] + 2.0 * fl1_fz * ts_zzzz_xx_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_xx_0[j];

            tt_zzzz_xy_0[j] =
                pa_z[j] * tt_zzz_xy_0[j] + 1.5 * fl1_fx * tt_zz_xy_0[j] + 2.0 * fl1_fz * ts_zzzz_xy_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_xy_0[j];

            tt_zzzz_xz_0[j] = pa_z[j] * tt_zzz_xz_0[j] + 1.5 * fl1_fx * tt_zz_xz_0[j] + 0.5 * fl1_fx * tt_zzz_x_0[j] +
                              2.0 * fl1_fz * ts_zzzz_xz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_xz_0[j];

            tt_zzzz_yy_0[j] =
                pa_z[j] * tt_zzz_yy_0[j] + 1.5 * fl1_fx * tt_zz_yy_0[j] + 2.0 * fl1_fz * ts_zzzz_yy_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_yy_0[j];

            tt_zzzz_yz_0[j] = pa_z[j] * tt_zzz_yz_0[j] + 1.5 * fl1_fx * tt_zz_yz_0[j] + 0.5 * fl1_fx * tt_zzz_y_0[j] +
                              2.0 * fl1_fz * ts_zzzz_yz_0[j] - 3.0 * fl1_fz * fl1_fga * ts_zz_yz_0[j];

            tt_zzzz_zz_0[j] = pa_z[j] * tt_zzz_zz_0[j] + 1.5 * fl1_fx * tt_zz_zz_0[j] + fl1_fx * tt_zzz_z_0[j] + 2.0 * fl1_fz * ts_zzzz_zz_0[j] -
                              3.0 * fl1_fz * fl1_fga * ts_zz_zz_0[j];
        }

        idx++;
    }
}

}  // namespace kinrecfunc
