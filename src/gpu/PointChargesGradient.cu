//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "GpuRuntime.hpp"


#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

#include "BoysFuncGPU.hpp"
#include "PointChargesGradient.hpp"
#include "ErrorHandler.hpp"
#include "GpuConstants.hpp"
#include "GpuSafeChecks.hpp"
#include "GtoFunc.hpp"
#include "GtoInfo.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "MatrixFunc.hpp"
#include "MultiTimer.hpp"
#include "StringFormat.hpp"

namespace gpu {  // gpu namespace

__global__ void
computePointChargesGradientSS(double*         grad_x,
                          const int32_t  grad_cart_ind,
                          const double*   s_prim_info,
                          const int32_t  s_prim_count,
                          const double*   ss_mat_D_local,
                          const int32_t* first_inds_local,
                          const int32_t* second_inds_local,
                          const int32_t  ss_prim_pair_count_local,
                          const int32_t* prim_cart_ao_to_atom_inds,
                          const double*   points_info,
                          const int32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    // each thread computes a primitive V matrix element

    const int32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (ij < ss_prim_pair_count_local)
    {
        const auto i = first_inds_local[ij];
        const auto j = second_inds_local[ij];

        const auto a_i = rawValue(s_prim_info, i + s_prim_count * 0);
        const auto c_i = rawValue(s_prim_info, i + s_prim_count * 1);
        const double r_i[3] = {rawValue(s_prim_info, i + s_prim_count * 2),
                               rawValue(s_prim_info, i + s_prim_count * 3),
                               rawValue(s_prim_info, i + s_prim_count * 4)};

        const auto a_j = rawValue(s_prim_info, j + s_prim_count * 0);
        const auto c_j = rawValue(s_prim_info, j + s_prim_count * 1);
        const double r_j[3] = {rawValue(s_prim_info, j + s_prim_count * 2),
                               rawValue(s_prim_info, j + s_prim_count * 3),
                               rawValue(s_prim_info, j + s_prim_count * 4)};

        const auto r2_ij = (r_j[0] - r_i[0]) * (r_j[0] - r_i[0]) +
                           (r_j[1] - r_i[1]) * (r_j[1] - r_i[1]) +
                           (r_j[2] - r_i[2]) * (r_j[2] - r_i[2]);

        // Electron. J. Theor. Chem., Vol. 2, 66–70 (1997)
        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI * inv_S1, 1.5) * exp(-a_i * a_j * inv_S1 * r2_ij);

        const auto ij_factor = (static_cast<double>(i != j) + 1.0);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        const auto PB_x = (-a_i  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        double grad_i = 0.0, grad_j = 0.0;

        for (int32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {(a_i * r_i[0] + a_j * r_j[0]) * inv_S1 - x_c,
                                  (a_i * r_i[1] + a_j * r_j[1]) * inv_S1 - y_c,
                                  (a_i * r_i[2] + a_j * r_j[2]) * inv_S1 - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F1_t[2];

            gpu::computeBoysFunction(F1_t, (a_i + a_j) * r2_PC, 1, boys_func_table, boys_func_ft);

            // Note: minus sign from electron charge

            const auto V_ij_for_i = (-1.0) * q_c * (

                F1_t[0] * (

                    2.0 * a_i * (
                        PA_x
                    )

                )

                + F1_t[1] * (

                    (-2.0) * a_i * (
                        PC[grad_cart_ind]
                    )

                )

            );

            // Note: minus sign from electron charge

            const auto V_ij_for_j = (-1.0) * q_c * (

                F1_t[0] * (

                    2.0 * a_j * (
                        PB_x
                    )

                )

                + F1_t[1] * (

                    (-2.0) * a_j * (
                        PC[grad_cart_ind]
                    )

                )

            );

            grad_i += theta * S_ij_00 * V_ij_for_i;
            grad_j += theta * S_ij_00 * V_ij_for_j;
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i * ss_mat_D_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[j], grad_j * ss_mat_D_local[ij] * ij_factor);
    }
}

__global__ void
computePointChargesGradientSP(double*         grad_x,
                          const int32_t  grad_cart_ind,
                          const double*   s_prim_info,
                          const int32_t  s_prim_count,
                          const double*   p_prim_info,
                          const int32_t  p_prim_count,
                          const double*   sp_mat_D_local,
                          const int32_t* sp_first_inds_local,
                          const int32_t* sp_second_inds_local,
                          const int32_t  sp_prim_pair_count_local,
                          const int32_t* prim_cart_ao_to_atom_inds,
                          const double*   points_info,
                          const int32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive V matrix element

    const int32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < sp_prim_pair_count_local)
    {
        const auto i = rawValue(sp_first_inds_local, ij);
        const auto j = rawValue(sp_second_inds_local, ij);

        const auto a_i = rawValue(s_prim_info, i + s_prim_count * 0);
        const auto c_i = rawValue(s_prim_info, i + s_prim_count * 1);
        const double r_i[3] = {rawValue(s_prim_info, i + s_prim_count * 2),
                               rawValue(s_prim_info, i + s_prim_count * 3),
                               rawValue(s_prim_info, i + s_prim_count * 4)};

        const auto a_j = rawValue(p_prim_info, j / 3 + p_prim_count * 0);
        const auto c_j = rawValue(p_prim_info, j / 3 + p_prim_count * 1);
        const double r_j[3] = {rawValue(p_prim_info, j / 3 + p_prim_count * 2),
                               rawValue(p_prim_info, j / 3 + p_prim_count * 3),
                               rawValue(p_prim_info, j / 3 + p_prim_count * 4)};

        const auto b0 = j % 3;

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI * inv_S1, 1.5) * exp(-a_i * a_j * inv_S1 * r2_ij);

        const auto ij_factor = 2.0;

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        const auto PB_x = (-a_i  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PB_0 = (-a_i * inv_S1) * rij[b0];

        double grad_i = 0.0, grad_j = 0.0;

        for (int32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {(a_i * r_i[0] + a_j * r_j[0]) * inv_S1 - x_c,
                                  (a_i * r_i[1] + a_j * r_j[1]) * inv_S1 - y_c,
                                  (a_i * r_i[2] + a_j * r_j[2]) * inv_S1 - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F2_t[3];

            gpu::computeBoysFunction(F2_t, (a_i + a_j) * r2_PC, 2, boys_func_table, boys_func_ft);

            // Note: minus sign from electron charge

            const auto V_ij_for_i = (-1.0) * q_c * (

                F2_t[0] * (

                    1.0 * inv_S1 * a_i * (
                        delta[b0][grad_cart_ind]
                    )

                    + 2.0 * a_i * (
                        PA_x * PB_0
                    )

                )

                + F2_t[1] * (

                    (-1.0) * inv_S1 * a_i * (
                        delta[b0][grad_cart_ind]
                    )

                    + (-2.0) * a_i * (
                        PA_x * PC[b0]
                        + PB_0 * PC[grad_cart_ind]
                    )

                )

                + F2_t[2] * (

                    2.0 * a_i * (
                        PC[b0] * PC[grad_cart_ind]
                    )

                )

            );

            // Note: minus sign from electron charge

            const auto V_ij_for_j = (-1.0) * q_c * (

                F2_t[0] * (

                    (-1.0) * (
                        delta[b0][grad_cart_ind]
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b0][grad_cart_ind]
                    )

                    + 2.0 * a_j * (
                        PB_0 * PB_x
                    )

                )

                + F2_t[1] * (

                    (-1.0) / (a_i + a_j) * a_j * (
                        delta[b0][grad_cart_ind]
                    )

                    + (-2.0) * a_j * (
                        PB_0 * PC[grad_cart_ind]
                        + PB_x * PC[b0]
                    )

                )

                + F2_t[2] * (

                    2.0 * a_j * (
                        PC[b0] * PC[grad_cart_ind]
                    )

                )

            );


            grad_i += theta * S_ij_00 * V_ij_for_i;
            grad_j += theta * S_ij_00 * V_ij_for_j;
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i * sp_mat_D_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + j], grad_j * sp_mat_D_local[ij] * ij_factor);
    }
}

__global__ void
computePointChargesGradientSD(double*         grad_x,
                          const int32_t  grad_cart_ind,
                          const double*   s_prim_info,
                          const int32_t  s_prim_count,
                          const double*   d_prim_info,
                          const int32_t  d_prim_count,
                          const double*   sd_mat_D_local,
                          const int32_t* sd_first_inds_local,
                          const int32_t* sd_second_inds_local,
                          const int32_t  sd_prim_pair_count_local,
                          const int32_t* prim_cart_ao_to_atom_inds,
                          const int32_t  p_prim_count,
                          const double*   points_info,
                          const int32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    __shared__ int32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

    const int32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        d_cart_inds[0][0] = 0; d_cart_inds[0][1] = 0;
        d_cart_inds[1][0] = 0; d_cart_inds[1][1] = 1;
        d_cart_inds[2][0] = 0; d_cart_inds[2][1] = 2;
        d_cart_inds[3][0] = 1; d_cart_inds[3][1] = 1;
        d_cart_inds[4][0] = 1; d_cart_inds[4][1] = 2;
        d_cart_inds[5][0] = 2; d_cart_inds[5][1] = 2;

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < sd_prim_pair_count_local)
    {
        const auto i = rawValue(sd_first_inds_local, ij);
        const auto j = rawValue(sd_second_inds_local, ij);

        const auto a_i = rawValue(s_prim_info, i + s_prim_count * 0);
        const auto c_i = rawValue(s_prim_info, i + s_prim_count * 1);
        const double r_i[3] = {rawValue(s_prim_info, i + s_prim_count * 2),
                               rawValue(s_prim_info, i + s_prim_count * 3),
                               rawValue(s_prim_info, i + s_prim_count * 4)};

        const auto a_j = rawValue(d_prim_info, j / 6 + d_prim_count * 0);
        const auto c_j = rawValue(d_prim_info, j / 6 + d_prim_count * 1);
        const double r_j[3] = {rawValue(d_prim_info, j / 6 + d_prim_count * 2),
                               rawValue(d_prim_info, j / 6 + d_prim_count * 3),
                               rawValue(d_prim_info, j / 6 + d_prim_count * 4)};

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI * inv_S1, 1.5) * exp(-a_i * a_j * inv_S1 * r2_ij);

        const auto ij_factor = 2.0;

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        const auto PB_x = (-a_i  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PB_0 = (-a_i * inv_S1) * rij[b0];
        const auto PB_1 = (-a_i * inv_S1) * rij[b1];

        double grad_i = 0.0, grad_j = 0.0;

        for (int32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {(a_i * r_i[0] + a_j * r_j[0]) * inv_S1 - x_c,
                                  (a_i * r_i[1] + a_j * r_j[1]) * inv_S1 - y_c,
                                  (a_i * r_i[2] + a_j * r_j[2]) * inv_S1 - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, (a_i + a_j) * r2_PC, 3, boys_func_table, boys_func_ft);

            // Note: minus sign from electron charge

            const auto V_ij_for_i = (-1.0) * q_c * (

                F3_t[0] * (

                    1.0 * inv_S1 * a_i * (
                        delta[b0][b1] * (PA_x)
                        + delta[b1][grad_cart_ind] * (PB_0)
                        + delta[b0][grad_cart_ind] * (PB_1)
                    )

                    + 2.0 * a_i * (
                        PA_x * PB_0 * PB_1
                    )

                )

                + F3_t[1] * (

                    (-1.0) * inv_S1 * a_i * (
                        delta[b0][b1] * (PA_x + PC[grad_cart_ind])
                        + delta[b1][grad_cart_ind] * (PB_0 + PC[b0])
                        + delta[b0][grad_cart_ind] * (PB_1 + PC[b1])
                    )

                    + (-2.0) * a_i * (
                        PA_x * PB_0 * PC[b1]
                        + PA_x * PB_1 * PC[b0]
                        + PB_0 * PB_1 * PC[grad_cart_ind]
                    )

                )

                + F3_t[2] * (

                    1.0 * inv_S1 * a_i * (
                        delta[b1][grad_cart_ind] * (PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[b1])
                        + delta[b0][b1] * (PC[grad_cart_ind])
                    )

                    + 2.0 * a_i * (
                        PA_x * PC[b0] * PC[b1]
                        + PB_0 * PC[b1] * PC[grad_cart_ind]
                        + PB_1 * PC[b0] * PC[grad_cart_ind]
                    )

                )

                + F3_t[3] * (

                    (-2.0) * a_i * (
                        PC[b0] * PC[b1] * PC[grad_cart_ind]
                    )

                )

            );

            // Note: minus sign from electron charge

            const auto V_ij_for_j = (-1.0) * q_c * (

                F3_t[0] * (

                    (-1.0) * (
                        delta[b1][grad_cart_ind] * (PB_0)
                        + delta[b0][grad_cart_ind] * (PB_1)
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PB_0)
                        + delta[b0][grad_cart_ind] * (PB_1)
                        + delta[b0][b1] * (PB_x)
                    )

                    + 2.0 * a_j * (
                        PB_0 * PB_1 * PB_x
                    )

                )

                + F3_t[1] * (

                    (-1.0) / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PB_0 + PC[b0])
                        + delta[b0][grad_cart_ind] * (PB_1 + PC[b1])
                        + delta[b0][b1] * (PB_x + PC[grad_cart_ind])
                    )

                    + (-2.0) * a_j * (
                        PB_0 * PB_1 * PC[grad_cart_ind]
                        + PB_0 * PB_x * PC[b1]
                        + PB_1 * PB_x * PC[b0]
                    )

                    + (
                        delta[b1][grad_cart_ind] * (PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[b1])
                    )

                )

                + F3_t[2] * (

                    1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[b1])
                        + delta[b0][b1] * (PC[grad_cart_ind])
                    )

                    + 2.0 * a_j * (
                        PB_0 * PC[b1] * PC[grad_cart_ind]
                        + PB_1 * PC[b0] * PC[grad_cart_ind]
                        + PB_x * PC[b0] * PC[b1]
                    )

                )

                + F3_t[3] * (

                    (-2.0) * a_j * (
                        PC[b0] * PC[b1] * PC[grad_cart_ind]
                    )

                )

            );


            grad_i += theta * S_ij_00 * V_ij_for_i;
            grad_j += theta * S_ij_00 * V_ij_for_j;
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i * sd_mat_D_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + j], grad_j * sd_mat_D_local[ij] * ij_factor);
    }
}

__global__ void
computePointChargesGradientPP(double*         grad_x,
                          const int32_t  grad_cart_ind,
                          const double*   p_prim_info,
                          const int32_t  p_prim_count,
                          const double*   pp_mat_D_local,
                          const int32_t* pp_first_inds_local,
                          const int32_t* pp_second_inds_local,
                          const int32_t  pp_prim_pair_count_local,
                          const int32_t* prim_cart_ao_to_atom_inds,
                          const int32_t  s_prim_count,
                          const double*   points_info,
                          const int32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    __shared__ double   delta[3][3];

    // each thread computes a primitive V matrix element

    const int32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < pp_prim_pair_count_local)
    {
        const auto i = rawValue(pp_first_inds_local, ij);
        const auto j = rawValue(pp_second_inds_local, ij);

        const auto a_i = rawValue(p_prim_info, i / 3 + p_prim_count * 0);
        const auto c_i = rawValue(p_prim_info, i / 3 + p_prim_count * 1);
        const double r_i[3] = {rawValue(p_prim_info, i / 3 + p_prim_count * 2),
                               rawValue(p_prim_info, i / 3 + p_prim_count * 3),
                               rawValue(p_prim_info, i / 3 + p_prim_count * 4)};

        const auto a_j = rawValue(p_prim_info, j / 3 + p_prim_count * 0);
        const auto c_j = rawValue(p_prim_info, j / 3 + p_prim_count * 1);
        const double r_j[3] = {rawValue(p_prim_info, j / 3 + p_prim_count * 2),
                               rawValue(p_prim_info, j / 3 + p_prim_count * 3),
                               rawValue(p_prim_info, j / 3 + p_prim_count * 4)};

        const auto a0 = i % 3;
        const auto b0 = j % 3;

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI * inv_S1, 1.5) * exp(-a_i * a_j * inv_S1 * r2_ij);

        const auto ij_factor = (static_cast<double>(i != j) + 1.0);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        const auto PB_x = (-a_i  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PA_0 = (a_j * inv_S1) * rij[a0];

        const auto PB_0 = (-a_i * inv_S1) * rij[b0];

        double grad_i = 0.0, grad_j = 0.0;

        for (int32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {(a_i * r_i[0] + a_j * r_j[0]) * inv_S1 - x_c,
                                  (a_i * r_i[1] + a_j * r_j[1]) * inv_S1 - y_c,
                                  (a_i * r_i[2] + a_j * r_j[2]) * inv_S1 - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F3_t[4];

            gpu::computeBoysFunction(F3_t, (a_i + a_j) * r2_PC, 3, boys_func_table, boys_func_ft);

            // Note: minus sign from electron charge

            const auto V_ij_for_i = (-1.0) * q_c * (

                F3_t[0] * (

                    (-1.0) * (
                        delta[a0][grad_cart_ind] * (PB_0)
                    )

                    + 1.0 * inv_S1 * a_i * (
                        delta[b0][grad_cart_ind] * (PA_0)
                        + delta[a0][b0] * (PA_x)
                        + delta[a0][grad_cart_ind] * (PB_0)
                    )

                    + 2.0 * a_i * (
                        PA_0 * PA_x * PB_0
                    )

                )

                + F3_t[1] * (

                    (-1.0) * inv_S1 * a_i * (
                        delta[b0][grad_cart_ind] * (PA_0 + PC[a0])
                        + delta[a0][b0] * (PA_x + PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PB_0 + PC[b0])
                    )

                    + (-2.0) * a_i * (
                        PA_0 * PA_x * PC[b0]
                        + PA_0 * PB_0 * PC[grad_cart_ind]
                        + PA_x * PB_0 * PC[a0]
                    )

                    + (
                        delta[a0][grad_cart_ind] * (PC[b0])
                    )

                )

                + F3_t[2] * (

                    1.0 * inv_S1 * a_i * (
                        delta[b0][grad_cart_ind] * (PC[a0])
                        + delta[a0][grad_cart_ind] * (PC[b0])
                        + delta[a0][b0] * (PC[grad_cart_ind])
                    )

                    + 2.0 * a_i * (
                        PA_0 * PC[b0] * PC[grad_cart_ind]
                        + PA_x * PC[a0] * PC[b0]
                        + PB_0 * PC[a0] * PC[grad_cart_ind]
                    )

                )

                + F3_t[3] * (

                    (-2.0) * a_i * (
                        PC[a0] * PC[b0] * PC[grad_cart_ind]
                    )

                )

            );

            // Note: minus sign from electron charge

            const auto V_ij_for_j = (-1.0) * q_c * (

                F3_t[0] * (

                    (-1.0) * (
                        delta[b0][grad_cart_ind] * (PA_0)
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b0][grad_cart_ind] * (PA_0)
                        + delta[a0][grad_cart_ind] * (PB_0)
                        + delta[a0][b0] * (PB_x)
                    )

                    + 2.0 * a_j * (
                        PA_0 * PB_0 * PB_x
                    )

                )

                + F3_t[1] * (

                    (-1.0) / (a_i + a_j) * a_j * (
                        delta[b0][grad_cart_ind] * (PA_0 + PC[a0])
                        + delta[a0][grad_cart_ind] * (PB_0 + PC[b0])
                        + delta[a0][b0] * (PB_x + PC[grad_cart_ind])
                    )

                    + (-2.0) * a_j * (
                        PA_0 * PB_0 * PC[grad_cart_ind]
                        + PA_0 * PB_x * PC[b0]
                        + PB_0 * PB_x * PC[a0]
                    )

                    + (
                        delta[b0][grad_cart_ind] * (PC[a0])
                    )

                )

                + F3_t[2] * (

                    1.0 / (a_i + a_j) * a_j * (
                        delta[b0][grad_cart_ind] * (PC[a0])
                        + delta[a0][grad_cart_ind] * (PC[b0])
                        + delta[a0][b0] * (PC[grad_cart_ind])
                    )

                    + 2.0 * a_j * (
                        PA_0 * PC[b0] * PC[grad_cart_ind]
                        + PB_0 * PC[a0] * PC[grad_cart_ind]
                        + PB_x * PC[a0] * PC[b0]
                    )

                )

                + F3_t[3] * (

                    (-2.0) * a_j * (
                        PC[a0] * PC[b0] * PC[grad_cart_ind]
                    )

                )

            );


            grad_i += theta * S_ij_00 * V_ij_for_i;
            grad_j += theta * S_ij_00 * V_ij_for_j;
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + i], grad_i * pp_mat_D_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + j], grad_j * pp_mat_D_local[ij] * ij_factor);
    }
}

__global__ void
computePointChargesGradientPD(double*         grad_x,
                          const int32_t  grad_cart_ind,
                          const double*   p_prim_info,
                          const int32_t  p_prim_count,
                          const double*   d_prim_info,
                          const int32_t  d_prim_count,
                          const double*   pd_mat_D_local,
                          const int32_t* pd_first_inds_local,
                          const int32_t* pd_second_inds_local,
                          const int32_t  pd_prim_pair_count_local,
                          const int32_t* prim_cart_ao_to_atom_inds,
                          const int32_t  s_prim_count,
                          const double*   points_info,
                          const int32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    __shared__ int32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

    const int32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        d_cart_inds[0][0] = 0; d_cart_inds[0][1] = 0;
        d_cart_inds[1][0] = 0; d_cart_inds[1][1] = 1;
        d_cart_inds[2][0] = 0; d_cart_inds[2][1] = 2;
        d_cart_inds[3][0] = 1; d_cart_inds[3][1] = 1;
        d_cart_inds[4][0] = 1; d_cart_inds[4][1] = 2;
        d_cart_inds[5][0] = 2; d_cart_inds[5][1] = 2;

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < pd_prim_pair_count_local)
    {
        const auto i = rawValue(pd_first_inds_local, ij);
        const auto j = rawValue(pd_second_inds_local, ij);

        const auto a_i = rawValue(p_prim_info, i / 3 + p_prim_count * 0);
        const auto c_i = rawValue(p_prim_info, i / 3 + p_prim_count * 1);
        const double r_i[3] = {rawValue(p_prim_info, i / 3 + p_prim_count * 2),
                               rawValue(p_prim_info, i / 3 + p_prim_count * 3),
                               rawValue(p_prim_info, i / 3 + p_prim_count * 4)};

        const auto a_j = rawValue(d_prim_info, j / 6 + d_prim_count * 0);
        const auto c_j = rawValue(d_prim_info, j / 6 + d_prim_count * 1);
        const double r_j[3] = {rawValue(d_prim_info, j / 6 + d_prim_count * 2),
                               rawValue(d_prim_info, j / 6 + d_prim_count * 3),
                               rawValue(d_prim_info, j / 6 + d_prim_count * 4)};

        const auto a0 = i % 3;

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI * inv_S1, 1.5) * exp(-a_i * a_j * inv_S1 * r2_ij);

        const auto ij_factor = 2.0;

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        const auto PB_x = (-a_i  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PA_0 = (a_j * inv_S1) * rij[a0];

        const auto PB_0 = (-a_i * inv_S1) * rij[b0];
        const auto PB_1 = (-a_i * inv_S1) * rij[b1];

        double grad_i = 0.0, grad_j = 0.0;

        for (int32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {(a_i * r_i[0] + a_j * r_j[0]) * inv_S1 - x_c,
                                  (a_i * r_i[1] + a_j * r_j[1]) * inv_S1 - y_c,
                                  (a_i * r_i[2] + a_j * r_j[2]) * inv_S1 - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F4_t[5];

            gpu::computeBoysFunction(F4_t, (a_i + a_j) * r2_PC, 4, boys_func_table, boys_func_ft);

            // Note: minus sign from electron charge

            const auto V_ij_for_i = (-1.0) * q_c * (

                F4_t[0] * (

                    (-0.5) * inv_S1 * (
                        delta[a0][grad_cart_ind] * delta[b0][b1]
                    )

                    + (-1.0) * (
                        delta[a0][grad_cart_ind] * (PB_0 * PB_1)
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                    )

                    + 1.0 * inv_S1 * a_i * (
                        delta[b0][b1] * (PA_0 * PA_x)
                        + delta[b1][grad_cart_ind] * (PA_0 * PB_0)
                        + delta[b0][grad_cart_ind] * (PA_0 * PB_1)
                        + delta[a0][b1] * (PA_x * PB_0)
                        + delta[a0][b0] * (PA_x * PB_1)
                        + delta[a0][grad_cart_ind] * (PB_0 * PB_1)
                    )

                    + 2.0 * a_i * (
                        PA_0 * PA_x * PB_0 * PB_1
                    )

                )

                + F4_t[1] * (

                    (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                    )

                    + (-1.0) * inv_S1 * a_i * (
                        delta[b0][b1] * (PA_0 * PA_x + PA_0 * PC[grad_cart_ind] + PA_x * PC[a0])
                        + delta[b1][grad_cart_ind] * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                        + delta[a0][b1] * (PA_x * PB_0 + PA_x * PC[b0] + PB_0 * PC[grad_cart_ind])
                        + delta[a0][b0] * (PA_x * PB_1 + PA_x * PC[b1] + PB_1 * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                    )

                    + (-2.0) * a_i * (
                        PA_0 * PA_x * PB_0 * PC[b1]
                        + PA_0 * PA_x * PB_1 * PC[b0]
                        + PA_0 * PB_0 * PB_1 * PC[grad_cart_ind]
                        + PA_x * PB_0 * PB_1 * PC[a0]
                    )

                    + 0.5 * inv_S1 * (
                        delta[a0][grad_cart_ind] * delta[b0][b1]
                    )

                    + (
                        delta[a0][grad_cart_ind] * (PB_0 * PC[b1] + PB_1 * PC[b0])
                    )

                )

                + F4_t[2] * (

                    (-1.0) * (
                        delta[a0][grad_cart_ind] * (PC[b0] * PC[b1])
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                    )

                    + 1.0 * inv_S1 * a_i * (
                        delta[b1][grad_cart_ind] * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                        + delta[b0][b1] * (PA_0 * PC[grad_cart_ind] + PA_x * PC[a0] + PC[a0] * PC[grad_cart_ind])
                        + delta[a0][b1] * (PA_x * PC[b0] + PB_0 * PC[grad_cart_ind] + PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PA_x * PC[b1] + PB_1 * PC[grad_cart_ind] + PC[b1] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                    )

                    + 2.0 * a_i * (
                        PA_0 * PA_x * PC[b0] * PC[b1]
                        + PA_0 * PB_0 * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PB_1 * PC[b0] * PC[grad_cart_ind]
                        + PA_x * PB_0 * PC[a0] * PC[b1]
                        + PA_x * PB_1 * PC[a0] * PC[b0]
                        + PB_0 * PB_1 * PC[a0] * PC[grad_cart_ind]
                    )

                )

                + F4_t[3] * (

                    (-1.0) * inv_S1 * a_i * (
                        delta[b1][grad_cart_ind] * (PC[a0] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[a0] * PC[b1])
                        + delta[b0][b1] * (PC[a0] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PC[b0] * PC[b1])
                        + delta[a0][b1] * (PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PC[b1] * PC[grad_cart_ind])
                    )

                    + (-2.0) * a_i * (
                        PA_0 * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PA_x * PC[a0] * PC[b0] * PC[b1]
                        + PB_0 * PC[a0] * PC[b1] * PC[grad_cart_ind]
                        + PB_1 * PC[a0] * PC[b0] * PC[grad_cart_ind]
                    )

                )

                + F4_t[4] * (

                    2.0 * a_i * (
                        PC[a0] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                    )

                )

            );

            // Note: minus sign from electron charge

            const auto V_ij_for_j = (-1.0) * q_c * (

                F4_t[0] * (

                    (-0.5) / (a_i + a_j) * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind])
                    )

                    + (-1.0) * (
                        delta[b1][grad_cart_ind] * (PA_0 * PB_0)
                        + delta[b0][grad_cart_ind] * (PA_0 * PB_1)
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PB_0)
                        + delta[b0][grad_cart_ind] * (PA_0 * PB_1)
                        + delta[b0][b1] * (PA_0 * PB_x)
                        + delta[a0][grad_cart_ind] * (PB_0 * PB_1)
                        + delta[a0][b1] * (PB_0 * PB_x)
                        + delta[a0][b0] * (PB_1 * PB_x)
                    )

                    + 2.0 * a_j * (
                        PA_0 * PB_0 * PB_1 * PB_x
                    )

                )

                + F4_t[1] * (

                    (-1.0) / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                    )

                    + (-1.0) / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PB_0 + PA_0 * PC[b0] + PB_0 * PC[a0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PB_1 + PA_0 * PC[b1] + PB_1 * PC[a0])
                        + delta[b0][b1] * (PA_0 * PB_x + PA_0 * PC[grad_cart_ind] + PB_x * PC[a0])
                        + delta[a0][grad_cart_ind] * (PB_0 * PB_1 + PB_0 * PC[b1] + PB_1 * PC[b0])
                        + delta[a0][b1] * (PB_0 * PB_x + PB_0 * PC[grad_cart_ind] + PB_x * PC[b0])
                        + delta[a0][b0] * (PB_1 * PB_x + PB_1 * PC[grad_cart_ind] + PB_x * PC[b1])
                    )

                    + (-2.0) * a_j * (
                        PA_0 * PB_0 * PB_1 * PC[grad_cart_ind]
                        + PA_0 * PB_0 * PB_x * PC[b1]
                        + PA_0 * PB_1 * PB_x * PC[b0]
                        + PB_0 * PB_1 * PB_x * PC[a0]
                    )

                    + 0.5 / (a_i + a_j) * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind])
                    )

                    + (
                        delta[b1][grad_cart_ind] * (PA_0 * PC[b0] + PB_0 * PC[a0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PC[b1] + PB_1 * PC[a0])
                    )

                )

                + F4_t[2] * (

                    (-1.0) * (
                        delta[b1][grad_cart_ind] * (PC[a0] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[a0] * PC[b1])
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1])
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PC[b0] + PB_0 * PC[a0] + PC[a0] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PC[b1] + PB_1 * PC[a0] + PC[a0] * PC[b1])
                        + delta[b0][b1] * (PA_0 * PC[grad_cart_ind] + PB_x * PC[a0] + PC[a0] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PB_0 * PC[b1] + PB_1 * PC[b0] + PC[b0] * PC[b1])
                        + delta[a0][b1] * (PB_0 * PC[grad_cart_ind] + PB_x * PC[b0] + PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PB_1 * PC[grad_cart_ind] + PB_x * PC[b1] + PC[b1] * PC[grad_cart_ind])
                    )

                    + 2.0 * a_j * (
                        PA_0 * PB_0 * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PB_1 * PC[b0] * PC[grad_cart_ind]
                        + PA_0 * PB_x * PC[b0] * PC[b1]
                        + PB_0 * PB_1 * PC[a0] * PC[grad_cart_ind]
                        + PB_0 * PB_x * PC[a0] * PC[b1]
                        + PB_1 * PB_x * PC[a0] * PC[b0]
                    )

                )

                + F4_t[3] * (

                    (-1.0) / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PC[a0] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[a0] * PC[b1])
                        + delta[b0][b1] * (PC[a0] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PC[b0] * PC[b1])
                        + delta[a0][b1] * (PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PC[b1] * PC[grad_cart_ind])
                    )

                    + (-2.0) * a_j * (
                        PA_0 * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PB_0 * PC[a0] * PC[b1] * PC[grad_cart_ind]
                        + PB_1 * PC[a0] * PC[b0] * PC[grad_cart_ind]
                        + PB_x * PC[a0] * PC[b0] * PC[b1]
                    )

                )

                + F4_t[4] * (

                    2.0 * a_j * (
                        PC[a0] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                    )

                )

            );


            grad_i += theta * S_ij_00 * V_ij_for_i;
            grad_j += theta * S_ij_00 * V_ij_for_j;
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + i], grad_i * pd_mat_D_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + j], grad_j * pd_mat_D_local[ij] * ij_factor);
    }
}

__global__ void
computePointChargesGradientDD(double*         grad_x,
                          const int32_t  grad_cart_ind,
                          const double*   d_prim_info,
                          const int32_t  d_prim_count,
                          const double*   dd_mat_D_local,
                          const int32_t* dd_first_inds_local,
                          const int32_t* dd_second_inds_local,
                          const int32_t  dd_prim_pair_count_local,
                          const int32_t* prim_cart_ao_to_atom_inds,
                          const int32_t  s_prim_count,
                          const int32_t  p_prim_count,
                          const double*   points_info,
                          const int32_t  npoints,
                          const double*   boys_func_table,
                          const double*   boys_func_ft)
{
    __shared__ int32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    // each thread computes a primitive S/T matrix element

    const int32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadIdx.x == 0)
    {
        d_cart_inds[0][0] = 0; d_cart_inds[0][1] = 0;
        d_cart_inds[1][0] = 0; d_cart_inds[1][1] = 1;
        d_cart_inds[2][0] = 0; d_cart_inds[2][1] = 2;
        d_cart_inds[3][0] = 1; d_cart_inds[3][1] = 1;
        d_cart_inds[4][0] = 1; d_cart_inds[4][1] = 2;
        d_cart_inds[5][0] = 2; d_cart_inds[5][1] = 2;

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;
    }

    __syncthreads();

    if (ij < dd_prim_pair_count_local)
    {
        const auto i = rawValue(dd_first_inds_local, ij);
        const auto j = rawValue(dd_second_inds_local, ij);

        // TODO: improve memory access pattern

        const auto a_i = rawValue(d_prim_info, i / 6 + d_prim_count * 0);
        const auto c_i = rawValue(d_prim_info, i / 6 + d_prim_count * 1);
        const double r_i[3] = {rawValue(d_prim_info, i / 6 + d_prim_count * 2),
                               rawValue(d_prim_info, i / 6 + d_prim_count * 3),
                               rawValue(d_prim_info, i / 6 + d_prim_count * 4)};

        const auto a_j = rawValue(d_prim_info, j / 6 + d_prim_count * 0);
        const auto c_j = rawValue(d_prim_info, j / 6 + d_prim_count * 1);
        const double r_j[3] = {rawValue(d_prim_info, j / 6 + d_prim_count * 2),
                               rawValue(d_prim_info, j / 6 + d_prim_count * 3),
                               rawValue(d_prim_info, j / 6 + d_prim_count * 4)};

        const auto a0 = d_cart_inds[i % 6][0];
        const auto a1 = d_cart_inds[i % 6][1];

        const auto b0 = d_cart_inds[j % 6][0];
        const auto b1 = d_cart_inds[j % 6][1];

        const double rij[3] = {r_j[0] - r_i[0], r_j[1] - r_i[1], r_j[2] - r_i[2]};

        const auto r2_ij = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto inv_S1 = 1.0 / (a_i + a_j);

        const auto S_ij_00 = c_i * c_j * pow(MATH_CONST_PI * inv_S1, 1.5) * exp(-a_i * a_j * inv_S1 * r2_ij);

        const auto ij_factor = (static_cast<double>(i != j) + 1.0);

        const auto theta = MATH_CONST_TWO_OVER_SQRT_PI * sqrt(a_i + a_j);

        const auto PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        const auto PB_x = (-a_i  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);

        const auto PA_0 = (a_j * inv_S1) * rij[a0];
        const auto PA_1 = (a_j * inv_S1) * rij[a1];

        const auto PB_0 = (-a_i * inv_S1) * rij[b0];
        const auto PB_1 = (-a_i * inv_S1) * rij[b1];

        double grad_i = 0.0, grad_j = 0.0;

        // TODO: use 2D block to scan the points

        for (int32_t c = 0; c < npoints; c++)
        {
            const auto x_c = points_info[c + npoints * 0];
            const auto y_c = points_info[c + npoints * 1];
            const auto z_c = points_info[c + npoints * 2];
            const auto q_c = points_info[c + npoints * 3];

            const double PC[3] = {(a_i * r_i[0] + a_j * r_j[0]) * inv_S1 - x_c,
                                  (a_i * r_i[1] + a_j * r_j[1]) * inv_S1 - y_c,
                                  (a_i * r_i[2] + a_j * r_j[2]) * inv_S1 - z_c};

            const auto r2_PC = PC[0] * PC[0] + PC[1] * PC[1] + PC[2] * PC[2];

            double F5_t[6];

            gpu::computeBoysFunction(F5_t, (a_i + a_j) * r2_PC, 5, boys_func_table, boys_func_ft);

            // Note: minus sign from electron charge

            const auto V_ij_for_i = (-1.0) * q_c * (

                F5_t[0] * (

                    (-0.5) * inv_S1 * (
                        delta[a1][grad_cart_ind] * delta[b0][b1] * (PA_0)
                        + delta[a0][grad_cart_ind] * delta[b0][b1] * (PA_1)
                        + (delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0)
                        + (delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1)
                    )

                    + (-1.0) * (
                        delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PB_1)
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PB_1)
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PA_0)
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PA_1)
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_x)
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0)
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1)
                    )

                    + 1.0 * inv_S1 * a_i * (
                        delta[b0][b1] * (PA_0 * PA_1 * PA_x)
                        + delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PB_0)
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PB_1)
                        + delta[a1][b1] * (PA_0 * PA_x * PB_0)
                        + delta[a1][b0] * (PA_0 * PA_x * PB_1)
                        + delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PB_1)
                        + delta[a0][b1] * (PA_1 * PA_x * PB_0)
                        + delta[a0][b0] * (PA_1 * PA_x * PB_1)
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PB_1)
                        + delta[a0][a1] * (PA_x * PB_0 * PB_1)
                    )

                    + 2.0 * a_i * (
                        PA_0 * PA_1 * PA_x * PB_0 * PB_1
                    )

                )

                + F5_t[1] * (

                    0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_x * (-2.0) + PC[grad_cart_ind] * (-1.0))
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                    )

                    + (-1.0) * inv_S1 * a_i * (
                        delta[b0][b1] * (PA_0 * PA_1 * PA_x + PA_0 * PA_1 * PC[grad_cart_ind] + PA_0 * PA_x * PC[a1] + PA_1 * PA_x * PC[a0])
                        + delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                        + delta[a1][b1] * (PA_0 * PA_x * PB_0 + PA_0 * PA_x * PC[b0] + PA_0 * PB_0 * PC[grad_cart_ind] + PA_x * PB_0 * PC[a0])
                        + delta[a1][b0] * (PA_0 * PA_x * PB_1 + PA_0 * PA_x * PC[b1] + PA_0 * PB_1 * PC[grad_cart_ind] + PA_x * PB_1 * PC[a0])
                        + delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                        + delta[a0][b1] * (PA_1 * PA_x * PB_0 + PA_1 * PA_x * PC[b0] + PA_1 * PB_0 * PC[grad_cart_ind] + PA_x * PB_0 * PC[a1])
                        + delta[a0][b0] * (PA_1 * PA_x * PB_1 + PA_1 * PA_x * PC[b1] + PA_1 * PB_1 * PC[grad_cart_ind] + PA_x * PB_1 * PC[a1])
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                        + delta[a0][a1] * (PA_x * PB_0 * PB_1 + PA_x * PB_0 * PC[b1] + PA_x * PB_1 * PC[b0] + PB_0 * PB_1 * PC[grad_cart_ind])
                    )

                    + (-2.0) * a_i * (
                        PA_0 * PA_1 * PA_x * PB_0 * PC[b1]
                        + PA_0 * PA_1 * PA_x * PB_1 * PC[b0]
                        + PA_0 * PA_1 * PB_0 * PB_1 * PC[grad_cart_ind]
                        + PA_0 * PA_x * PB_0 * PB_1 * PC[a1]
                        + PA_1 * PA_x * PB_0 * PB_1 * PC[a0]
                    )

                    + 0.5 * inv_S1 * (
                        delta[a1][grad_cart_ind] * delta[b0][b1] * (PA_0 + PC[a0])
                        + delta[a0][grad_cart_ind] * delta[b0][b1] * (PA_1 + PC[a1])
                        + (delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0 + PC[b0])
                        + (delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1 + PC[b1])
                    )

                    + (
                        delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                    )

                )

                + F5_t[2] * (

                    (-0.5) * inv_S1 * (
                        delta[a1][grad_cart_ind] * delta[b0][b1] * (PC[a0])
                        + delta[a0][grad_cart_ind] * delta[b0][b1] * (PC[a1])
                        + (delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PC[b0])
                        + (delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PC[b1])
                    )

                    + (-1.0) * (
                        delta[a1][grad_cart_ind] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                        + delta[a0][grad_cart_ind] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PA_0 + PC[a0] * 2.0)
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PA_1 + PC[a1] * 2.0)
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PA_x + PC[grad_cart_ind] * 2.0)
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0 + PC[b0] * 2.0)
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1 + PC[b1] * 2.0)
                    )

                    + 1.0 * inv_S1 * a_i * (
                        delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                        + delta[b0][b1] * (PA_0 * PA_1 * PC[grad_cart_ind] + PA_0 * PA_x * PC[a1] + PA_0 * PC[a1] * PC[grad_cart_ind] + PA_1 * PA_x * PC[a0] + PA_1 * PC[a0] * PC[grad_cart_ind] + PA_x * PC[a0] * PC[a1])
                        + delta[a1][b1] * (PA_0 * PA_x * PC[b0] + PA_0 * PB_0 * PC[grad_cart_ind] + PA_0 * PC[b0] * PC[grad_cart_ind] + PA_x * PB_0 * PC[a0] + PA_x * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[grad_cart_ind])
                        + delta[a1][b0] * (PA_0 * PA_x * PC[b1] + PA_0 * PB_1 * PC[grad_cart_ind] + PA_0 * PC[b1] * PC[grad_cart_ind] + PA_x * PB_1 * PC[a0] + PA_x * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[grad_cart_ind])
                        + delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                        + delta[a0][b1] * (PA_1 * PA_x * PC[b0] + PA_1 * PB_0 * PC[grad_cart_ind] + PA_1 * PC[b0] * PC[grad_cart_ind] + PA_x * PB_0 * PC[a1] + PA_x * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PA_1 * PA_x * PC[b1] + PA_1 * PB_1 * PC[grad_cart_ind] + PA_1 * PC[b1] * PC[grad_cart_ind] + PA_x * PB_1 * PC[a1] + PA_x * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                        + delta[a0][a1] * (PA_x * PB_0 * PC[b1] + PA_x * PB_1 * PC[b0] + PA_x * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[grad_cart_ind] + PB_0 * PC[b1] * PC[grad_cart_ind] + PB_1 * PC[b0] * PC[grad_cart_ind])
                    )

                    + 2.0 * a_i * (
                        PA_0 * PA_1 * PA_x * PC[b0] * PC[b1]
                        + PA_0 * PA_1 * PB_0 * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PA_1 * PB_1 * PC[b0] * PC[grad_cart_ind]
                        + PA_0 * PA_x * PB_0 * PC[a1] * PC[b1]
                        + PA_0 * PA_x * PB_1 * PC[a1] * PC[b0]
                        + PA_0 * PB_0 * PB_1 * PC[a1] * PC[grad_cart_ind]
                        + PA_1 * PA_x * PB_0 * PC[a0] * PC[b1]
                        + PA_1 * PA_x * PB_1 * PC[a0] * PC[b0]
                        + PA_1 * PB_0 * PB_1 * PC[a0] * PC[grad_cart_ind]
                        + PA_x * PB_0 * PB_1 * PC[a0] * PC[a1]
                    )

                )

                + F5_t[3] * (

                    (-0.5) / ( (a_i + a_j) * (a_i + a_j) ) * a_i * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PC[a0])
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PC[a1])
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PC[b0])
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PC[b1])
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[grad_cart_ind])
                    )

                    + (-1.0) * inv_S1 * a_i * (
                        delta[b1][grad_cart_ind] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                        + delta[b0][b1] * (PA_0 * PC[a1] * PC[grad_cart_ind] + PA_1 * PC[a0] * PC[grad_cart_ind] + PA_x * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[grad_cart_ind])
                        + delta[a1][grad_cart_ind] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                        + delta[a1][b1] * (PA_0 * PC[b0] * PC[grad_cart_ind] + PA_x * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[grad_cart_ind] + PC[a0] * PC[b0] * PC[grad_cart_ind])
                        + delta[a1][b0] * (PA_0 * PC[b1] * PC[grad_cart_ind] + PA_x * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[grad_cart_ind] + PC[a0] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                        + delta[a0][b1] * (PA_1 * PC[b0] * PC[grad_cart_ind] + PA_x * PC[a1] * PC[b0] + PB_0 * PC[a1] * PC[grad_cart_ind] + PC[a1] * PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PA_1 * PC[b1] * PC[grad_cart_ind] + PA_x * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[grad_cart_ind] + PC[a1] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][a1] * (PA_x * PC[b0] * PC[b1] + PB_0 * PC[b1] * PC[grad_cart_ind] + PB_1 * PC[b0] * PC[grad_cart_ind] + PC[b0] * PC[b1] * PC[grad_cart_ind])
                    )

                    + (-2.0) * a_i * (
                        PA_0 * PA_1 * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PA_x * PC[a1] * PC[b0] * PC[b1]
                        + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[grad_cart_ind]
                        + PA_1 * PA_x * PC[a0] * PC[b0] * PC[b1]
                        + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[grad_cart_ind]
                        + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[grad_cart_ind]
                        + PA_x * PB_0 * PC[a0] * PC[a1] * PC[b1]
                        + PA_x * PB_1 * PC[a0] * PC[a1] * PC[b0]
                        + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[grad_cart_ind]
                    )

                    + (
                        delta[a1][grad_cart_ind] * (PC[a0] * PC[b0] * PC[b1])
                        + delta[a0][grad_cart_ind] * (PC[a1] * PC[b0] * PC[b1])
                    )

                )

                + F5_t[4] * (

                    1.0 * inv_S1 * a_i * (
                        delta[b1][grad_cart_ind] * (PC[a0] * PC[a1] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[a0] * PC[a1] * PC[b1])
                        + delta[b0][b1] * (PC[a0] * PC[a1] * PC[grad_cart_ind])
                        + delta[a1][grad_cart_ind] * (PC[a0] * PC[b0] * PC[b1])
                        + delta[a1][b1] * (PC[a0] * PC[b0] * PC[grad_cart_ind])
                        + delta[a1][b0] * (PC[a0] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PC[a1] * PC[b0] * PC[b1])
                        + delta[a0][b1] * (PC[a1] * PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PC[a1] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][a1] * (PC[b0] * PC[b1] * PC[grad_cart_ind])
                    )

                    + 2.0 * a_i * (
                        PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PA_x * PC[a0] * PC[a1] * PC[b0] * PC[b1]
                        + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[grad_cart_ind]
                        + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[grad_cart_ind]
                    )

                )

                + F5_t[5] * (

                    (-2.0) * a_i * (
                        PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                    )

                )

            );

            // Note: minus sign from electron charge

            const auto V_ij_for_j = (-1.0) * q_c * (

                F5_t[0] * (

                    (-0.5) / (a_i + a_j) * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind]) * (PA_0)
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind]) * (PA_1)
                        + delta[a0][a1] * delta[b1][grad_cart_ind] * (PB_0)
                        + delta[a0][a1] * delta[b0][grad_cart_ind] * (PB_1)
                    )

                    + (-1.0) * (
                        delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PB_0)
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PB_1)
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PA_0)
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PA_1)
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0)
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1)
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_x)
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PB_0)
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PB_1)
                        + delta[b0][b1] * (PA_0 * PA_1 * PB_x)
                        + delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PB_1)
                        + delta[a1][b1] * (PA_0 * PB_0 * PB_x)
                        + delta[a1][b0] * (PA_0 * PB_1 * PB_x)
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PB_1)
                        + delta[a0][b1] * (PA_1 * PB_0 * PB_x)
                        + delta[a0][b0] * (PA_1 * PB_1 * PB_x)
                        + delta[a0][a1] * (PB_0 * PB_1 * PB_x)
                    )

                    + 2.0 * a_j * (
                        PA_0 * PA_1 * PB_0 * PB_1 * PB_x
                    )

                )

                + F5_t[1] * (

                    0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PA_0 * (-2.0) + PC[a0] * (-1.0))
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PA_1 * (-2.0) + PC[a1] * (-1.0))
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0 * (-2.0) + PC[b0] * (-1.0))
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1 * (-2.0) + PC[b1] * (-1.0))
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_x * (-2.0) + PC[grad_cart_ind] * (-1.0))
                    )

                    + (-1.0) / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PB_0 + PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PB_1 + PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                        + delta[b0][b1] * (PA_0 * PA_1 * PB_x + PA_0 * PA_1 * PC[grad_cart_ind] + PA_0 * PB_x * PC[a1] + PA_1 * PB_x * PC[a0])
                        + delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PB_1 + PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a0])
                        + delta[a1][b1] * (PA_0 * PB_0 * PB_x + PA_0 * PB_0 * PC[grad_cart_ind] + PA_0 * PB_x * PC[b0] + PB_0 * PB_x * PC[a0])
                        + delta[a1][b0] * (PA_0 * PB_1 * PB_x + PA_0 * PB_1 * PC[grad_cart_ind] + PA_0 * PB_x * PC[b1] + PB_1 * PB_x * PC[a0])
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PB_1 + PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PB_0 * PB_1 * PC[a1])
                        + delta[a0][b1] * (PA_1 * PB_0 * PB_x + PA_1 * PB_0 * PC[grad_cart_ind] + PA_1 * PB_x * PC[b0] + PB_0 * PB_x * PC[a1])
                        + delta[a0][b0] * (PA_1 * PB_1 * PB_x + PA_1 * PB_1 * PC[grad_cart_ind] + PA_1 * PB_x * PC[b1] + PB_1 * PB_x * PC[a1])
                        + delta[a0][a1] * (PB_0 * PB_1 * PB_x + PB_0 * PB_1 * PC[grad_cart_ind] + PB_0 * PB_x * PC[b1] + PB_1 * PB_x * PC[b0])
                    )

                    + (-2.0) * a_j * (
                        PA_0 * PA_1 * PB_0 * PB_1 * PC[grad_cart_ind]
                        + PA_0 * PA_1 * PB_0 * PB_x * PC[b1]
                        + PA_0 * PA_1 * PB_1 * PB_x * PC[b0]
                        + PA_0 * PB_0 * PB_1 * PB_x * PC[a1]
                        + PA_1 * PB_0 * PB_1 * PB_x * PC[a0]
                    )

                    + 0.5 / (a_i + a_j) * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind]) * (PA_0 + PC[a0])
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind]) * (PA_1 + PC[a1])
                        + delta[a0][a1] * delta[b1][grad_cart_ind] * (PB_0 + PC[b0])
                        + delta[a0][a1] * delta[b0][grad_cart_ind] * (PB_1 + PC[b1])
                    )

                    + (
                        delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_1 * PB_0 * PC[a0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_1 * PB_1 * PC[a0])
                    )

                )

                + F5_t[2] * (

                    (-0.5) / (a_i + a_j) * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind]) * (PC[a0])
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind]) * (PC[a1])
                        + delta[a0][a1] * delta[b1][grad_cart_ind] * (PC[b0])
                        + delta[a0][a1] * delta[b0][grad_cart_ind] * (PC[b1])
                    )

                    + (-1.0) * (
                        delta[b1][grad_cart_ind] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                        + delta[b0][grad_cart_ind] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                    )

                    + 0.5 / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PA_0 + PC[a0] * 2.0)
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PA_1 + PC[a1] * 2.0)
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PB_0 + PC[b0] * 2.0)
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PB_1 + PC[b1] * 2.0)
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PB_x + PC[grad_cart_ind] * 2.0)
                    )

                    + 1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PA_1 * PC[b0] + PA_0 * PB_0 * PC[a1] + PA_0 * PC[a1] * PC[b0] + PA_1 * PB_0 * PC[a0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1])
                        + delta[b0][grad_cart_ind] * (PA_0 * PA_1 * PC[b1] + PA_0 * PB_1 * PC[a1] + PA_0 * PC[a1] * PC[b1] + PA_1 * PB_1 * PC[a0] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1])
                        + delta[b0][b1] * (PA_0 * PA_1 * PC[grad_cart_ind] + PA_0 * PB_x * PC[a1] + PA_0 * PC[a1] * PC[grad_cart_ind] + PA_1 * PB_x * PC[a0] + PA_1 * PC[a0] * PC[grad_cart_ind] + PB_x * PC[a0] * PC[a1])
                        + delta[a1][grad_cart_ind] * (PA_0 * PB_0 * PC[b1] + PA_0 * PB_1 * PC[b0] + PA_0 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a0] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0])
                        + delta[a1][b1] * (PA_0 * PB_0 * PC[grad_cart_ind] + PA_0 * PB_x * PC[b0] + PA_0 * PC[b0] * PC[grad_cart_ind] + PB_0 * PB_x * PC[a0] + PB_0 * PC[a0] * PC[grad_cart_ind] + PB_x * PC[a0] * PC[b0])
                        + delta[a1][b0] * (PA_0 * PB_1 * PC[grad_cart_ind] + PA_0 * PB_x * PC[b1] + PA_0 * PC[b1] * PC[grad_cart_ind] + PB_1 * PB_x * PC[a0] + PB_1 * PC[a0] * PC[grad_cart_ind] + PB_x * PC[a0] * PC[b1])
                        + delta[a0][grad_cart_ind] * (PA_1 * PB_0 * PC[b1] + PA_1 * PB_1 * PC[b0] + PA_1 * PC[b0] * PC[b1] + PB_0 * PB_1 * PC[a1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0])
                        + delta[a0][b1] * (PA_1 * PB_0 * PC[grad_cart_ind] + PA_1 * PB_x * PC[b0] + PA_1 * PC[b0] * PC[grad_cart_ind] + PB_0 * PB_x * PC[a1] + PB_0 * PC[a1] * PC[grad_cart_ind] + PB_x * PC[a1] * PC[b0])
                        + delta[a0][b0] * (PA_1 * PB_1 * PC[grad_cart_ind] + PA_1 * PB_x * PC[b1] + PA_1 * PC[b1] * PC[grad_cart_ind] + PB_1 * PB_x * PC[a1] + PB_1 * PC[a1] * PC[grad_cart_ind] + PB_x * PC[a1] * PC[b1])
                        + delta[a0][a1] * (PB_0 * PB_1 * PC[grad_cart_ind] + PB_0 * PB_x * PC[b1] + PB_0 * PC[b1] * PC[grad_cart_ind] + PB_1 * PB_x * PC[b0] + PB_1 * PC[b0] * PC[grad_cart_ind] + PB_x * PC[b0] * PC[b1])
                    )

                    + 2.0 * a_j * (
                        PA_0 * PA_1 * PB_0 * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PA_1 * PB_1 * PC[b0] * PC[grad_cart_ind]
                        + PA_0 * PA_1 * PB_x * PC[b0] * PC[b1]
                        + PA_0 * PB_0 * PB_1 * PC[a1] * PC[grad_cart_ind]
                        + PA_0 * PB_0 * PB_x * PC[a1] * PC[b1]
                        + PA_0 * PB_1 * PB_x * PC[a1] * PC[b0]
                        + PA_1 * PB_0 * PB_1 * PC[a0] * PC[grad_cart_ind]
                        + PA_1 * PB_0 * PB_x * PC[a0] * PC[b1]
                        + PA_1 * PB_1 * PB_x * PC[a0] * PC[b0]
                        + PB_0 * PB_1 * PB_x * PC[a0] * PC[a1]
                    )

                )

                + F5_t[3] * (

                    (-0.5) / ( (a_i + a_j) * (a_i + a_j) ) * a_j * (
                        (delta[a1][b0] * delta[b1][grad_cart_ind] + delta[a1][b1] * delta[b0][grad_cart_ind] + delta[a1][grad_cart_ind] * delta[b0][b1]) * (PC[a0])
                        + (delta[a0][b0] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[b0][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[b0][b1]) * (PC[a1])
                        + (delta[a0][a1] * delta[b1][grad_cart_ind] + delta[a0][b1] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b1]) * (PC[b0])
                        + (delta[a0][a1] * delta[b0][grad_cart_ind] + delta[a0][b0] * delta[a1][grad_cart_ind] + delta[a0][grad_cart_ind] * delta[a1][b0]) * (PC[b1])
                        + (delta[a0][a1] * delta[b0][b1] + delta[a0][b0] * delta[a1][b1] + delta[a0][b1] * delta[a1][b0]) * (PC[grad_cart_ind])
                    )

                    + (-1.0) / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PA_0 * PC[a1] * PC[b0] + PA_1 * PC[a0] * PC[b0] + PB_0 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PA_0 * PC[a1] * PC[b1] + PA_1 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[b1])
                        + delta[b0][b1] * (PA_0 * PC[a1] * PC[grad_cart_ind] + PA_1 * PC[a0] * PC[grad_cart_ind] + PB_x * PC[a0] * PC[a1] + PC[a0] * PC[a1] * PC[grad_cart_ind])
                        + delta[a1][grad_cart_ind] * (PA_0 * PC[b0] * PC[b1] + PB_0 * PC[a0] * PC[b1] + PB_1 * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[b1])
                        + delta[a1][b1] * (PA_0 * PC[b0] * PC[grad_cart_ind] + PB_0 * PC[a0] * PC[grad_cart_ind] + PB_x * PC[a0] * PC[b0] + PC[a0] * PC[b0] * PC[grad_cart_ind])
                        + delta[a1][b0] * (PA_0 * PC[b1] * PC[grad_cart_ind] + PB_1 * PC[a0] * PC[grad_cart_ind] + PB_x * PC[a0] * PC[b1] + PC[a0] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PA_1 * PC[b0] * PC[b1] + PB_0 * PC[a1] * PC[b1] + PB_1 * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[b1])
                        + delta[a0][b1] * (PA_1 * PC[b0] * PC[grad_cart_ind] + PB_0 * PC[a1] * PC[grad_cart_ind] + PB_x * PC[a1] * PC[b0] + PC[a1] * PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PA_1 * PC[b1] * PC[grad_cart_ind] + PB_1 * PC[a1] * PC[grad_cart_ind] + PB_x * PC[a1] * PC[b1] + PC[a1] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][a1] * (PB_0 * PC[b1] * PC[grad_cart_ind] + PB_1 * PC[b0] * PC[grad_cart_ind] + PB_x * PC[b0] * PC[b1] + PC[b0] * PC[b1] * PC[grad_cart_ind])
                    )

                    + (-2.0) * a_j * (
                        PA_0 * PA_1 * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PB_0 * PC[a1] * PC[b1] * PC[grad_cart_ind]
                        + PA_0 * PB_1 * PC[a1] * PC[b0] * PC[grad_cart_ind]
                        + PA_0 * PB_x * PC[a1] * PC[b0] * PC[b1]
                        + PA_1 * PB_0 * PC[a0] * PC[b1] * PC[grad_cart_ind]
                        + PA_1 * PB_1 * PC[a0] * PC[b0] * PC[grad_cart_ind]
                        + PA_1 * PB_x * PC[a0] * PC[b0] * PC[b1]
                        + PB_0 * PB_1 * PC[a0] * PC[a1] * PC[grad_cart_ind]
                        + PB_0 * PB_x * PC[a0] * PC[a1] * PC[b1]
                        + PB_1 * PB_x * PC[a0] * PC[a1] * PC[b0]
                    )

                    + (
                        delta[b1][grad_cart_ind] * (PC[a0] * PC[a1] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[a0] * PC[a1] * PC[b1])
                    )

                )

                + F5_t[4] * (

                    1.0 / (a_i + a_j) * a_j * (
                        delta[b1][grad_cart_ind] * (PC[a0] * PC[a1] * PC[b0])
                        + delta[b0][grad_cart_ind] * (PC[a0] * PC[a1] * PC[b1])
                        + delta[b0][b1] * (PC[a0] * PC[a1] * PC[grad_cart_ind])
                        + delta[a1][grad_cart_ind] * (PC[a0] * PC[b0] * PC[b1])
                        + delta[a1][b1] * (PC[a0] * PC[b0] * PC[grad_cart_ind])
                        + delta[a1][b0] * (PC[a0] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][grad_cart_ind] * (PC[a1] * PC[b0] * PC[b1])
                        + delta[a0][b1] * (PC[a1] * PC[b0] * PC[grad_cart_ind])
                        + delta[a0][b0] * (PC[a1] * PC[b1] * PC[grad_cart_ind])
                        + delta[a0][a1] * (PC[b0] * PC[b1] * PC[grad_cart_ind])
                    )

                    + 2.0 * a_j * (
                        PA_0 * PC[a1] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PA_1 * PC[a0] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                        + PB_0 * PC[a0] * PC[a1] * PC[b1] * PC[grad_cart_ind]
                        + PB_1 * PC[a0] * PC[a1] * PC[b0] * PC[grad_cart_ind]
                        + PB_x * PC[a0] * PC[a1] * PC[b0] * PC[b1]
                    )

                )

                + F5_t[5] * (

                    (-2.0) * a_j * (
                        PC[a0] * PC[a1] * PC[b0] * PC[b1] * PC[grad_cart_ind]
                    )

                )

            );


            grad_i += theta * S_ij_00 * V_ij_for_i;
            grad_j += theta * S_ij_00 * V_ij_for_j;
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + i], grad_i * dd_mat_D_local[ij] * ij_factor);
        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + j], grad_j * dd_mat_D_local[ij] * ij_factor);
    }
}

}  // namespace gpu
