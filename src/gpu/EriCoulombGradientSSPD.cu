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



#include "BoysFuncGPU.hpp"
#include "EriCoulombGradientSSPD.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_J)
computeCoulombGradientSSPD_I_0(double*         grad_x,
                               const uint32_t  grad_cart_ind,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   d_prim_info,
                               const uint32_t  d_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   pd_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   pd_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* pd_first_inds,
                               const uint32_t* pd_second_inds,
                               const double*   pd_pair_data,
                               const uint32_t  pd_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PB_x;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
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

    if (ij < ss_prim_pair_count_local)
    {
        i = ss_first_inds_local[ij];
        j = ss_second_inds_local[ij];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_j = s_prim_info[j + s_prim_count * 0];

        r_j[0] = s_prim_info[j + s_prim_count * 2];
        r_j[1] = s_prim_info[j + s_prim_count * 3];
        r_j[2] = s_prim_info[j + s_prim_count * 4];

        S1 = a_i + a_j;
        inv_S1 = 1.0 / S1;

        S_ij_00 = ss_pair_data_local[ij];

        ij_factor_D = (static_cast<double>(i != j) + 1.0) * ss_mat_D_local[ij];

        PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        PB_x = (-a_i * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);




    }

    for (uint32_t m = 0; m < (pd_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= pd_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * pd_mat_Q[kl] * pd_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = pd_first_inds[kl];
        const auto l = pd_second_inds[kl];

        const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];

        const double r_k[3] = {p_prim_info[k / 3 + p_prim_count * 2],
                               p_prim_info[k / 3 + p_prim_count * 3],
                               p_prim_info[k / 3 + p_prim_count * 4]};

        const auto a_l = d_prim_info[l / 6 + d_prim_count * 0];

        const double r_l[3] = {d_prim_info[l / 6 + d_prim_count * 2],
                               d_prim_info[l / 6 + d_prim_count * 3],
                               d_prim_info[l / 6 + d_prim_count * 4]};

        const auto S_kl_00 = pd_pair_data[kl];

        const auto c0 = k % 3;
        const auto d0 = d_cart_inds[l % 6][0];
        const auto d1 = d_cart_inds[l % 6][1];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F4_t[5];

        gpu::computeBoysFunction(F4_t, S1 * S2 * inv_S4 * r2_PQ, 4, boys_func_table, boys_func_ft);

        const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
        const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);
        const auto QD_1 = (-a_k * inv_S2) * (r_l[d1] - r_k[d1]);


        double kl_factor = 2.0;

        // mu grad

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    F4_t[0] * inv_S2 * a_i * (
                        delta[c0][d1] * (PA_x * QD_0)
                        + delta[c0][d0] * (PA_x * QD_1)
                        + delta[d0][d1] * (PA_x * QC_0)
                    )

                    + F4_t[0] * 2.0 * a_i * (
                        + PA_x * QC_0 * QD_0 * QD_1
                    )

                    + F4_t[1] * (-1.0) * S1 * inv_S2 * inv_S4 * a_i * (
                        delta[d0][d1] * (PA_x * PQ[c0] + PA_x * QC_0)
                        + delta[c0][d1] * (PA_x * PQ[d0] + PA_x * QD_0)
                        + delta[c0][d0] * (PA_x * PQ[d1] + PA_x * QD_1)
                    )

                    + F4_t[1] * inv_S4 * a_i * (
                        delta[c0][d1] * (PQ[grad_cart_ind] * QD_0)
                        + delta[c0][d0] * (PQ[grad_cart_ind] * QD_1)
                        + delta[d1][grad_cart_ind] * (QC_0 * QD_0)
                        + delta[d0][grad_cart_ind] * (QC_0 * QD_1)
                        + delta[c0][grad_cart_ind] * (QD_0 * QD_1)
                        + delta[d0][d1] * (PQ[grad_cart_ind] * QC_0)
                    )

                    + F4_t[1] * 2.0 * S1 * inv_S4 * a_i * (
                        + PA_x * PQ[c0] * QD_0 * QD_1 * (-1.0)
                        + PA_x * PQ[d0] * QC_0 * QD_1 * (-1.0)
                        + PA_x * PQ[d1] * QC_0 * QD_0 * (-1.0)
                    )

                    + F4_t[1] * 2.0 * S2 * inv_S4 * a_i * (
                        + PQ[grad_cart_ind] * QC_0 * QD_0 * QD_1
                    )

                    + F4_t[1] * 0.5 * inv_S2 * inv_S4 * a_i * (
                        (delta[c0][d0] * delta[d1][grad_cart_ind] + delta[c0][d1] * delta[d0][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[d0][d1])
                    )

                    + F4_t[2] * (-0.5) * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * (
                        (delta[c0][d0] * delta[d1][grad_cart_ind] + delta[c0][d1] * delta[d0][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[d0][d1])
                    )

                    + F4_t[2] * (-1.0) * S1 * inv_S4 * inv_S4 * a_i * (
                        delta[d0][d1] * (PQ[c0] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QC_0)
                        + delta[d1][grad_cart_ind] * (PQ[c0] * QD_0 + PQ[d0] * QC_0)
                        + delta[d0][grad_cart_ind] * (PQ[c0] * QD_1 + PQ[d1] * QC_0)
                        + delta[c0][d1] * (PQ[d0] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QD_0)
                        + delta[c0][grad_cart_ind] * (PQ[d0] * QD_1 + PQ[d1] * QD_0)
                        + delta[c0][d0] * (PQ[d1] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QD_1)
                    )

                    + F4_t[2] * 2.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * (
                        + PA_x * PQ[c0] * PQ[d0] * QD_1
                        + PA_x * PQ[c0] * PQ[d1] * QD_0
                        + PA_x * PQ[d0] * PQ[d1] * QC_0
                    )

                    + F4_t[2] * 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                        + PQ[c0] * PQ[grad_cart_ind] * QD_0 * QD_1 * (-1.0)
                        + PQ[d0] * PQ[grad_cart_ind] * QC_0 * QD_1 * (-1.0)
                        + PQ[d1] * PQ[grad_cart_ind] * QC_0 * QD_0 * (-1.0)
                    )

                    + F4_t[2] * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * (
                        delta[d0][d1] * (PA_x * PQ[c0])
                        + delta[c0][d1] * (PA_x * PQ[d0])
                        + delta[c0][d0] * (PA_x * PQ[d1])
                    )

                    + F4_t[3] * (-2.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                        PA_x * PQ[c0] * PQ[d0] * PQ[d1]
                    )

                    + F4_t[3] * 2.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                        + PQ[c0] * PQ[d0] * PQ[grad_cart_ind] * QD_1
                        + PQ[c0] * PQ[d1] * PQ[grad_cart_ind] * QD_0
                        + PQ[d0] * PQ[d1] * PQ[grad_cart_ind] * QC_0
                    )

                    + F4_t[3] * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                        delta[d1][grad_cart_ind] * (PQ[c0] * PQ[d0])
                        + delta[d0][grad_cart_ind] * (PQ[c0] * PQ[d1])
                        + delta[d0][d1] * (PQ[c0] * PQ[grad_cart_ind])
                        + delta[c0][grad_cart_ind] * (PQ[d0] * PQ[d1])
                        + delta[c0][d1] * (PQ[d0] * PQ[grad_cart_ind])
                        + delta[c0][d0] * (PQ[d1] * PQ[grad_cart_ind])
                    )

                    + F4_t[4] * (-2.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                        PQ[c0] * PQ[d0] * PQ[d1] * PQ[grad_cart_ind]
                    )

                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * pd_mat_D[kl] * kl_factor;

    }


    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
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
computeCoulombGradientSSPD_J_0(double*         grad_x,
                               const uint32_t  grad_cart_ind,
                               const double    prefac_coulomb,
                               const double*   s_prim_info,
                               const uint32_t  s_prim_count,
                               const double*   p_prim_info,
                               const uint32_t  p_prim_count,
                               const double*   d_prim_info,
                               const uint32_t  d_prim_count,
                               const double*   ss_mat_D_local,
                               const double*   pd_mat_D,
                               const double*   ss_mat_Q_local,
                               const double*   pd_mat_Q,
                               const uint32_t* ss_first_inds_local,
                               const uint32_t* ss_second_inds_local,
                               const double*   ss_pair_data_local,
                               const uint32_t  ss_prim_pair_count_local,
                               const uint32_t* pd_first_inds,
                               const uint32_t* pd_second_inds,
                               const double*   pd_pair_data,
                               const uint32_t  pd_prim_pair_count,
                               const uint32_t* prim_cart_ao_to_atom_inds,
                               const double*   boys_func_table,
                               const double*   boys_func_ft,
                               const double    eri_threshold)
{
    // each thread row scans over [ij|??] and sum up to a primitive J matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM][TILE_DIM + 1];
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    const uint32_t ij = blockDim.x * blockIdx.x + threadIdx.x;

    double a_i, a_j, r_i[3], r_j[3], S_ij_00, S1, inv_S1, ij_factor_D;
    double PA_x, PB_x;
    uint32_t i, j;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
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

    if (ij < ss_prim_pair_count_local)
    {
        i = ss_first_inds_local[ij];
        j = ss_second_inds_local[ij];

        a_i = s_prim_info[i + s_prim_count * 0];

        r_i[0] = s_prim_info[i + s_prim_count * 2];
        r_i[1] = s_prim_info[i + s_prim_count * 3];
        r_i[2] = s_prim_info[i + s_prim_count * 4];

        a_j = s_prim_info[j + s_prim_count * 0];

        r_j[0] = s_prim_info[j + s_prim_count * 2];
        r_j[1] = s_prim_info[j + s_prim_count * 3];
        r_j[2] = s_prim_info[j + s_prim_count * 4];

        S1 = a_i + a_j;
        inv_S1 = 1.0 / S1;

        S_ij_00 = ss_pair_data_local[ij];

        ij_factor_D = (static_cast<double>(i != j) + 1.0) * ss_mat_D_local[ij];

        PA_x = (a_j  * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);
        PB_x = (-a_i * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);




    }

    for (uint32_t m = 0; m < (pd_prim_pair_count + TILE_DIM - 1) / TILE_DIM; m++)
    {
        const uint32_t kl = m * TILE_DIM + threadIdx.y;

        if ((kl >= pd_prim_pair_count) || (ij >= ss_prim_pair_count_local) || (fabs(ss_mat_Q_local[ij] * pd_mat_Q[kl] * pd_mat_D[kl]) <= eri_threshold))
        {
            break;
        }

        const auto k = pd_first_inds[kl];
        const auto l = pd_second_inds[kl];

        const auto a_k = p_prim_info[k / 3 + p_prim_count * 0];

        const double r_k[3] = {p_prim_info[k / 3 + p_prim_count * 2],
                               p_prim_info[k / 3 + p_prim_count * 3],
                               p_prim_info[k / 3 + p_prim_count * 4]};

        const auto a_l = d_prim_info[l / 6 + d_prim_count * 0];

        const double r_l[3] = {d_prim_info[l / 6 + d_prim_count * 2],
                               d_prim_info[l / 6 + d_prim_count * 3],
                               d_prim_info[l / 6 + d_prim_count * 4]};

        const auto S_kl_00 = pd_pair_data[kl];

        const auto c0 = k % 3;
        const auto d0 = d_cart_inds[l % 6][0];
        const auto d1 = d_cart_inds[l % 6][1];

        // J. Chem. Phys. 84, 3963-3974 (1986)

        const auto S2 = a_k + a_l;

        const auto inv_S2 = 1.0 / S2;
        const auto inv_S4 = 1.0 / (S1 + S2);

        const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                              (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                              (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

        const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

        const auto Lambda = sqrt(4.0 * S1 * S2 * MATH_CONST_INV_PI * inv_S4);

        double F4_t[5];

        gpu::computeBoysFunction(F4_t, S1 * S2 * inv_S4 * r2_PQ, 4, boys_func_table, boys_func_ft);

        const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
        const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);
        const auto QD_1 = (-a_k * inv_S2) * (r_l[d1] - r_k[d1]);


        double kl_factor = 2.0;

        // nu grad

        const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                    F4_t[0] * inv_S2 * a_j * (
                        delta[d0][d1] * (PB_x * QC_0)
                        + delta[c0][d1] * (PB_x * QD_0)
                        + delta[c0][d0] * (PB_x * QD_1)
                    )

                    + F4_t[0] * 2.0 * a_j * (
                        + PB_x * QC_0 * QD_0 * QD_1
                    )

                    + F4_t[1] * (-1.0) * S1 * inv_S2 * inv_S4 * a_j * (
                        delta[d0][d1] * (PB_x * PQ[c0] + PB_x * QC_0)
                        + delta[c0][d1] * (PB_x * PQ[d0] + PB_x * QD_0)
                        + delta[c0][d0] * (PB_x * PQ[d1] + PB_x * QD_1)
                    )

                    + F4_t[1] * inv_S4 * a_j * (
                        delta[c0][d1] * (PQ[grad_cart_ind] * QD_0)
                        + delta[c0][d0] * (PQ[grad_cart_ind] * QD_1)
                        + delta[d1][grad_cart_ind] * (QC_0 * QD_0)
                        + delta[d0][grad_cart_ind] * (QC_0 * QD_1)
                        + delta[c0][grad_cart_ind] * (QD_0 * QD_1)
                        + delta[d0][d1] * (PQ[grad_cart_ind] * QC_0)
                    )

                    + F4_t[1] * 2.0 * S1 * inv_S4 * a_j * (
                        + PB_x * PQ[c0] * QD_0 * QD_1 * (-1.0)
                        + PB_x * PQ[d0] * QC_0 * QD_1 * (-1.0)
                        + PB_x * PQ[d1] * QC_0 * QD_0 * (-1.0)
                    )

                    + F4_t[1] * 2.0 * S2 * inv_S4 * a_j * (
                        + PQ[grad_cart_ind] * QC_0 * QD_0 * QD_1
                    )

                    + F4_t[1] * 0.5 * inv_S2 * inv_S4 * a_j * (
                        (delta[c0][d0] * delta[d1][grad_cart_ind] + delta[c0][d1] * delta[d0][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[d0][d1])
                    )

                    + F4_t[2] * (-0.5) * S1 * inv_S2 * inv_S4 * inv_S4 * a_j * (
                        (delta[c0][d0] * delta[d1][grad_cart_ind] + delta[c0][d1] * delta[d0][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[d0][d1])
                    )

                    + F4_t[2] * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_j * (
                        delta[d0][d1] * (PB_x * PQ[c0])
                        + delta[c0][d1] * (PB_x * PQ[d0])
                        + delta[c0][d0] * (PB_x * PQ[d1])
                    )

                    + F4_t[2] * (-1.0) * S1 * inv_S4 * inv_S4 * a_j * (
                        delta[d0][d1] * (PQ[c0] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QC_0)
                        + delta[d1][grad_cart_ind] * (PQ[c0] * QD_0 + PQ[d0] * QC_0)
                        + delta[d0][grad_cart_ind] * (PQ[c0] * QD_1 + PQ[d1] * QC_0)
                        + delta[c0][d1] * (PQ[d0] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QD_0)
                        + delta[c0][grad_cart_ind] * (PQ[d0] * QD_1 + PQ[d1] * QD_0)
                        + delta[c0][d0] * (PQ[d1] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QD_1)
                    )

                    + F4_t[2] * 2.0 * S1 * S1 * inv_S4 * inv_S4 * a_j * (
                        + PB_x * PQ[c0] * PQ[d0] * QD_1
                        + PB_x * PQ[c0] * PQ[d1] * QD_0
                        + PB_x * PQ[d0] * PQ[d1] * QC_0
                    )

                    + F4_t[2] * 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_j * (
                        + PQ[c0] * PQ[grad_cart_ind] * QD_0 * QD_1 * (-1.0)
                        + PQ[d0] * PQ[grad_cart_ind] * QC_0 * QD_1 * (-1.0)
                        + PQ[d1] * PQ[grad_cart_ind] * QC_0 * QD_0 * (-1.0)
                    )

                    + F4_t[3] * 2.0 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * (
                        + PB_x * PQ[c0] * PQ[d0] * PQ[d1] * (-1.0)
                    )

                    + F4_t[3] * 2.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_j * (
                        + PQ[c0] * PQ[d0] * PQ[grad_cart_ind] * QD_1
                        + PQ[c0] * PQ[d1] * PQ[grad_cart_ind] * QD_0
                        + PQ[d0] * PQ[d1] * PQ[grad_cart_ind] * QC_0
                    )

                    + F4_t[3] * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_j * (
                        delta[d1][grad_cart_ind] * (PQ[c0] * PQ[d0])
                        + delta[d0][grad_cart_ind] * (PQ[c0] * PQ[d1])
                        + delta[d0][d1] * (PQ[c0] * PQ[grad_cart_ind])
                        + delta[c0][grad_cart_ind] * (PQ[d0] * PQ[d1])
                        + delta[c0][d1] * (PQ[d0] * PQ[grad_cart_ind])
                        + delta[c0][d0] * (PQ[d1] * PQ[grad_cart_ind])
                    )

                    + F4_t[4] * (-2.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_j * (
                        PQ[c0] * PQ[d0] * PQ[d1] * PQ[grad_cart_ind]
                    )


                );

        ERIs[threadIdx.y][threadIdx.x] += eri_ijkl * pd_mat_D[kl] * kl_factor;
    }


    __syncthreads();

    if ((threadIdx.y == 0) && (ij < ss_prim_pair_count_local))
    {
        double grad_j_x = 0.0;

        for (uint32_t n = 0; n < TILE_DIM; n++)
        {
            grad_j_x += ERIs[n][threadIdx.x];
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[j], grad_j_x * ij_factor_D * 2.0 * prefac_coulomb);
    }
}


}  // namespace gpu
