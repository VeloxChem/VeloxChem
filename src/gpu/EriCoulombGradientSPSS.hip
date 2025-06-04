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

#include <hip/hip_runtime.h>

#include "BoysFuncGPU.hpp"
#include "EriCoulombGradientSPSS.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombGradientSPSS_I_0(double*         grad_x,
                       const uint32_t  grad_cart_ind,
                       const double    prefac_coulomb,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double*   sp_mat_D_local,
                       const double*   ss_mat_D,
                       const double*   sp_mat_Q_local,
                       const double*   ss_mat_Q,
                       const uint32_t* sp_first_inds_local,
                       const uint32_t* sp_second_inds_local,
                       const double*   sp_pair_data_local,
                       const uint32_t  sp_prim_pair_count_local,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const double*   ss_pair_data,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* prim_cart_ao_to_atom_inds,
                       const double*   boys_func_table,
                       const double*   boys_func_ft,
                       const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_0, PA_x, PB_x;
    uint32_t i, j, b0;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

    }

    __syncthreads();

    if (ij < sp_prim_pair_count_local)
    {
        i = sp_first_inds_local[ij];
        j = sp_second_inds_local[ij];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_j = p_prim_info[j / 3 + p_prim_count * 0];

        r_j[0] = p_prim_info[j / 3 + p_prim_count * 2];
        r_j[1] = p_prim_info[j / 3 + p_prim_count * 3];
        r_j[2] = p_prim_info[j / 3 + p_prim_count * 4];

        S1 = a_i + a_j;
        inv_S1 = 1.0 / S1;

        S_ij_00 = sp_pair_data_local[ij];

        ij_factor_D = 2.0 * sp_mat_D_local[ij];

        PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        PB_x = (-a_i * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);


        b0 = j % 3;

        PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= sp_prim_pair_count_local) || (fabs(sp_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = ss_first_inds[kl];
        const auto l = ss_second_inds[kl];

        const auto a_k = s_prim_info[k + s_prim_count * 0];

        const double r_k[3] = {s_prim_info[k + s_prim_count * 2],
                               s_prim_info[k + s_prim_count * 3],
                               s_prim_info[k + s_prim_count * 4]};

        const auto a_l = s_prim_info[l + s_prim_count * 0];

        const double r_l[3] = {s_prim_info[l + s_prim_count * 2],
                               s_prim_info[l + s_prim_count * 3],
                               s_prim_info[l + s_prim_count * 4]};

        const auto S_kl_00 = ss_pair_data[kl];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F2_t[3];

        gpu::computeBoysFunction(F2_t, S1 * S2 * inv_S4 * r2_PQ, 2, boys_func_table, boys_func_ft);



        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // mu grad

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    F2_t[0] * 2.0 * a_i * (
                        + PA_x * PB_0
                    )

                    + F2_t[0] * inv_S1 * a_i * (
                        delta[b0][grad_cart_ind]
                    )

                    + F2_t[1] * (-1.0) * S2 * inv_S1 * inv_S4 * a_i * (
                        delta[b0][grad_cart_ind]
                    )

                    + F2_t[1] * 2.0 * S2 * inv_S4 * a_i * (
                        + PA_x * PQ[b0]
                        + PB_0 * PQ[grad_cart_ind]
                    )

                    + F2_t[2] * 2.0 * S2 * S2 * inv_S4 * inv_S4 * a_i * (
                        PQ[b0] * PQ[grad_cart_ind]
                    )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;

    }


    __syncthreads();

    if ((threadIdx.y == 0) && (ij < sp_prim_pair_count_local))
    {
        double grad_i_x = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            grad_i_x += ERIs[n][threadIdx.x];
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i_x * ij_factor_D * 2.0 * prefac_coulomb);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombGradientSPSS_J_0(double*         grad_x,
                       const uint32_t  grad_cart_ind,
                       const double    prefac_coulomb,
                       const double*   s_prim_info,
                       const uint32_t  s_prim_count,
                       const double*   p_prim_info,
                       const uint32_t  p_prim_count,
                       const double*   sp_mat_D_local,
                       const double*   ss_mat_D,
                       const double*   sp_mat_Q_local,
                       const double*   ss_mat_Q,
                       const uint32_t* sp_first_inds_local,
                       const uint32_t* sp_second_inds_local,
                       const double*   sp_pair_data_local,
                       const uint32_t  sp_prim_pair_count_local,
                       const uint32_t* ss_first_inds,
                       const uint32_t* ss_second_inds,
                       const double*   ss_pair_data,
                       const uint32_t  ss_prim_pair_count,
                       const uint32_t* prim_cart_ao_to_atom_inds,
                       const double*   boys_func_table,
                       const double*   boys_func_ft,
                       const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ double   delta[3][3];

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PB_0, PA_x, PB_x;
    uint32_t i, j, b0;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {

        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

    }

    __syncthreads();

    if (ij < sp_prim_pair_count_local)
    {
        i = sp_first_inds_local[ij];
        j = sp_second_inds_local[ij];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_j = p_prim_info[j / 3 + p_prim_count * 0];

        r_j[0] = p_prim_info[j / 3 + p_prim_count * 2];
        r_j[1] = p_prim_info[j / 3 + p_prim_count * 3];
        r_j[2] = p_prim_info[j / 3 + p_prim_count * 4];

        S1 = a_i + a_j;
        inv_S1 = 1.0 / S1;

        S_ij_00 = sp_pair_data_local[ij];

        ij_factor_D = 2.0 * sp_mat_D_local[ij];

        PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        PB_x = (-a_i * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);


        b0 = j % 3;

        PB_0 = (-a_i * inv_S1) * (r_j[b0] - r_i[b0]);

    }

    for (uint32_t m = 0; m < (ss_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= ss_prim_pair_count) || (ij >= sp_prim_pair_count_local) || (fabs(sp_mat_Q_local[ij] * ss_mat_Q[kl] * ss_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = ss_first_inds[kl];
        const auto l = ss_second_inds[kl];

        const auto a_k = s_prim_info[k + s_prim_count * 0];

        const double r_k[3] = {s_prim_info[k + s_prim_count * 2],
                               s_prim_info[k + s_prim_count * 3],
                               s_prim_info[k + s_prim_count * 4]};

        const auto a_l = s_prim_info[l + s_prim_count * 0];

        const double r_l[3] = {s_prim_info[l + s_prim_count * 2],
                               s_prim_info[l + s_prim_count * 3],
                               s_prim_info[l + s_prim_count * 4]};

        const auto S_kl_00 = ss_pair_data[kl];


        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F2_t[3];

        gpu::computeBoysFunction(F2_t, S1 * S2 * inv_S4 * r2_PQ, 2, boys_func_table, boys_func_ft);



        double kl_factor = (static_cast<double>(k != l) + 1.0);

        // nu grad

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    F2_t[0] * 2.0 * a_j * (
                        + PB_0 * PB_x
                    )

                    + F2_t[0] * (-1.0) * (
                        delta[b0][grad_cart_ind]
                    )

                    + F2_t[0] * inv_S1 * a_j * (
                        delta[b0][grad_cart_ind]
                    )

                    + F2_t[1] * (-1.0) * S2 * inv_S1 * inv_S4 * a_j * (
                        delta[b0][grad_cart_ind]
                    )

                    + F2_t[1] * 2.0 * S2 * inv_S4 * a_j * (
                        + PB_0 * PQ[grad_cart_ind]
                        + PB_x * PQ[b0]
                    )

                    + F2_t[2] * 2.0 * S2 * S2 * inv_S4 * inv_S4 * a_j * (
                        PQ[b0] * PQ[grad_cart_ind]
                    )


                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * ss_mat_D[kl] * kl_factor;
    }


    __syncthreads();

    if ((threadIdx.y == 0) && (ij < sp_prim_pair_count_local))
    {
        double grad_j_x = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            grad_j_x += ERIs[n][threadIdx.x];
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + j], grad_j_x * ij_factor_D * 2.0 * prefac_coulomb);
    }
}


}  // namespace gpu
