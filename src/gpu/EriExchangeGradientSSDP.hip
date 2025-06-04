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
#include "EriExchangeGradientSSDP.hpp"

namespace gpu {  // gpu namespace

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeGradientSSDP_I_0(double*         grad_x,
                        const uint32_t  grad_cart_ind,
                        const double    frac_exact_exchange,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const double*   D_ik_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_dp,
                        const uint32_t* prim_cart_ao_to_atom_inds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold)
{
    // each thread block scans over [i?|k?] and sum up to a primitive K matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM_Y_K][TILE_DIM_X_K + 1];
    __shared__ uint32_t i, k, count_i, count_k, displ_i, displ_k;
    __shared__ double   a_i, r_i[3], a_k, r_k[3], ik_factor_D;
    __shared__ uint32_t c0, c1;
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    const uint32_t ik = blockIdx.x;

    double d2 = 1.0;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        count_i = 0;
        count_k = 0;

        d_cart_inds[0][0] = 0; d_cart_inds[0][1] = 0;
        d_cart_inds[1][0] = 0; d_cart_inds[1][1] = 1;
        d_cart_inds[2][0] = 0; d_cart_inds[2][1] = 2;
        d_cart_inds[3][0] = 1; d_cart_inds[3][1] = 1;
        d_cart_inds[4][0] = 1; d_cart_inds[4][1] = 2;
        d_cart_inds[5][0] = 2; d_cart_inds[5][1] = 2;
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        if (ik < pair_inds_count_for_K_sd)
        {
            i = pair_inds_i_for_K_sd[ik];
            k = pair_inds_k_for_K_sd[ik];

            count_i = pair_counts_K_ss[i];
            count_k = pair_counts_K_dp[k];

            displ_i = pair_displs_K_ss[i];
            displ_k = pair_displs_K_dp[k];

            a_i = s_prim_info[i + s_prim_count * 0];

            r_i[0] = s_prim_info[i + s_prim_count * 2];
            r_i[1] = s_prim_info[i + s_prim_count * 3];
            r_i[2] = s_prim_info[i + s_prim_count * 4];

            a_k = d_prim_info[k / 6 + d_prim_count * 0];

            r_k[0] = d_prim_info[k / 6 + d_prim_count * 2];
            r_k[1] = d_prim_info[k / 6 + d_prim_count * 3];
            r_k[2] = d_prim_info[k / 6 + d_prim_count * 4];

            ik_factor_D = 2.0 * D_ik_for_K_sd[ik];

            c0 = d_cart_inds[k % 6][0];
            c1 = d_cart_inds[k % 6][1];

        }
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x;
        uint32_t j_prim, j_cgto;

        if ((ik < pair_inds_count_for_K_sd) && (j < count_i))
        {
            Q_ij   = Q_K_ss[displ_i + j];

            j_prim = D_inds_K_ss[displ_i + j];

            j_cgto = s_prim_aoinds[j_prim];

            a_j = s_prim_info[j_prim + s_prim_count * 0];

            r_j[0] = s_prim_info[j_prim + s_prim_count * 2];
            r_j[1] = s_prim_info[j_prim + s_prim_count * 3];
            r_j[2] = s_prim_info[j_prim + s_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_ss[displ_i + j];

            PA_x = (a_j * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);



        }


        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            if ((ik >= pair_inds_count_for_K_sd) || (j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_dp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            const auto l_prim = D_inds_K_dp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_dp[displ_k + l];

            const auto d0 = l_prim % 3;

            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F4_t[5];

            gpu::computeBoysFunction(F4_t, rho * d2 * r2_PQ, 4, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F4_t[1] *= d2;
                F4_t[2] *= d2 * d2;
                F4_t[3] *= d2 * d2 * d2;
                F4_t[4] *= d2 * d2 * d2 * d2;
            }

            const double QC_x = (a_l * inv_S2) * (r_l[grad_cart_ind] - r_k[grad_cart_ind]);

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QC_1 = (a_l * inv_S2) * (r_l[c1] - r_k[c1]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);


            // mu grad

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                            F4_t[0] * inv_S2 * a_i * (
                                delta[c0][c1] * (PA_x * QD_0)
                                + delta[c1][d0] * (PA_x * QC_0)
                                + delta[c0][d0] * (PA_x * QC_1)
                            )

                            + F4_t[0] * 2.0 * a_i * (
                                + PA_x * QC_0 * QC_1 * QD_0
                            )

                            + F4_t[1] * (-1.0) * S1 * inv_S2 * inv_S4 * a_i * (
                                delta[c1][d0] * (PA_x * PQ[c0] + PA_x * QC_0)
                                + delta[c0][d0] * (PA_x * PQ[c1] + PA_x * QC_1)
                                + delta[c0][c1] * (PA_x * PQ[d0] + PA_x * QD_0)
                            )

                            + F4_t[1] * inv_S4 * a_i * (
                                delta[c0][c1] * (PQ[grad_cart_ind] * QD_0)
                                + delta[c1][grad_cart_ind] * (QC_0 * QD_0)
                                + delta[c0][grad_cart_ind] * (QC_1 * QD_0)
                                + delta[c1][d0] * (PQ[grad_cart_ind] * QC_0)
                                + delta[c0][d0] * (PQ[grad_cart_ind] * QC_1)
                                + delta[d0][grad_cart_ind] * (QC_0 * QC_1)
                            )

                            + F4_t[1] * 2.0 * S1 * inv_S4 * a_i * (
                                + PA_x * PQ[c0] * QC_1 * QD_0 * (-1.0)
                                + PA_x * PQ[c1] * QC_0 * QD_0 * (-1.0)
                                + PA_x * PQ[d0] * QC_0 * QC_1 * (-1.0)
                            )

                            + F4_t[1] * 2.0 * S2 * inv_S4 * a_i * (
                                + PQ[grad_cart_ind] * QC_0 * QC_1 * QD_0
                            )

                            + F4_t[1] * 0.5 * inv_S2 * inv_S4 * a_i * (
                                (delta[c0][c1] * delta[d0][grad_cart_ind] + delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[2] * (-0.5) * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * (
                                (delta[c0][c1] * delta[d0][grad_cart_ind] + delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[2] * (-1.0) * S1 * inv_S4 * inv_S4 * a_i * (
                                delta[c1][d0] * (PQ[c0] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QC_0)
                                + delta[d0][grad_cart_ind] * (PQ[c0] * QC_1 + PQ[c1] * QC_0)
                                + delta[c1][grad_cart_ind] * (PQ[c0] * QD_0 + PQ[d0] * QC_0)
                                + delta[c0][d0] * (PQ[c1] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QC_1)
                                + delta[c0][grad_cart_ind] * (PQ[c1] * QD_0 + PQ[d0] * QC_1)
                                + delta[c0][c1] * (PQ[d0] * PQ[grad_cart_ind] + PQ[grad_cart_ind] * QD_0)
                            )

                            + F4_t[2] * 2.0 * S1 * S1 * inv_S4 * inv_S4 * a_i * (
                                + PA_x * PQ[c0] * PQ[c1] * QD_0
                                + PA_x * PQ[c0] * PQ[d0] * QC_1
                                + PA_x * PQ[c1] * PQ[d0] * QC_0
                            )

                            + F4_t[2] * 2.0 * S1 * S2 * inv_S4 * inv_S4 * a_i * (
                                + PQ[c0] * PQ[grad_cart_ind] * QC_1 * QD_0 * (-1.0)
                                + PQ[c1] * PQ[grad_cart_ind] * QC_0 * QD_0 * (-1.0)
                                + PQ[d0] * PQ[grad_cart_ind] * QC_0 * QC_1 * (-1.0)
                            )

                            + F4_t[2] * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_i * (
                                delta[c1][d0] * (PA_x * PQ[c0])
                                + delta[c0][d0] * (PA_x * PQ[c1])
                                + delta[c0][c1] * (PA_x * PQ[d0])
                            )

                            + F4_t[3] * (-2.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                                PA_x * PQ[c0] * PQ[c1] * PQ[d0]
                            )

                            + F4_t[3] * 2.0 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                                + PQ[c0] * PQ[c1] * PQ[grad_cart_ind] * QD_0
                                + PQ[c0] * PQ[d0] * PQ[grad_cart_ind] * QC_1
                                + PQ[c1] * PQ[d0] * PQ[grad_cart_ind] * QC_0
                            )

                            + F4_t[3] * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                                delta[d0][grad_cart_ind] * (PQ[c0] * PQ[c1])
                                + delta[c1][grad_cart_ind] * (PQ[c0] * PQ[d0])
                                + delta[c1][d0] * (PQ[c0] * PQ[grad_cart_ind])
                                + delta[c0][grad_cart_ind] * (PQ[c1] * PQ[d0])
                                + delta[c0][d0] * (PQ[c1] * PQ[grad_cart_ind])
                                + delta[c0][c1] * (PQ[d0] * PQ[grad_cart_ind])
                            )

                            + F4_t[4] * (-2.0) * S1 * S1 * S1 * S2 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_i * (
                                PQ[c0] * PQ[c1] * PQ[d0] * PQ[grad_cart_ind]
                            )

                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[j_cgto * naos + l_cgto];

        }
    }


    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sd))
    {
        double grad_i_x = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                grad_i_x += ERIs[y][x];
            }
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[i], grad_i_x * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}

__global__ void __launch_bounds__(TILE_SIZE_K)
computeExchangeGradientSSDP_K_0(double*         grad_x,
                        const uint32_t  grad_cart_ind,
                        const double    frac_exact_exchange,
                        const uint32_t* pair_inds_i_for_K_sd,
                        const uint32_t* pair_inds_k_for_K_sd,
                        const double*   D_ik_for_K_sd,
                        const uint32_t  pair_inds_count_for_K_sd,
                        const double*   s_prim_info,
                        const uint32_t* s_prim_aoinds,
                        const uint32_t  s_prim_count,
                        const double*   p_prim_info,
                        const uint32_t* p_prim_aoinds,
                        const uint32_t  p_prim_count,
                        const double*   d_prim_info,
                        const uint32_t* d_prim_aoinds,
                        const uint32_t  d_prim_count,
                        const double    sp_max_D,
                        const double*   mat_D_full_AO,
                        const uint32_t  naos,
                        const double*   Q_K_ss,
                        const double*   Q_K_dp,
                        const uint32_t* D_inds_K_ss,
                        const uint32_t* D_inds_K_dp,
                        const uint32_t* pair_displs_K_ss,
                        const uint32_t* pair_displs_K_dp,
                        const uint32_t* pair_counts_K_ss,
                        const uint32_t* pair_counts_K_dp,
                        const double*   pair_data_K_ss,
                        const double*   pair_data_K_dp,
                        const uint32_t* prim_cart_ao_to_atom_inds,
                        const double*   boys_func_table,
                        const double*   boys_func_ft,
                        const double    omega,
                        const double    eri_threshold)
{
    // each thread block scans over [i?|k?] and sum up to a primitive K matrix element
    // J. Chem. Theory Comput. 2009, 5, 4, 1004-1015

    __shared__ double   ERIs[TILE_DIM_Y_K][TILE_DIM_X_K + 1];
    __shared__ uint32_t i, k, count_i, count_k, displ_i, displ_k;
    __shared__ double   a_i, r_i[3], a_k, r_k[3], ik_factor_D;
    __shared__ uint32_t c0, c1;
    __shared__ uint32_t d_cart_inds[6][2];
    __shared__ double   delta[3][3];

    const uint32_t ik = blockIdx.x;

    double d2 = 1.0;

    ERIs[threadIdx.y][threadIdx.x] = 0.0;

    if ((threadIdx.y == 0) && (threadIdx.x == 0))
    {
        count_i = 0;
        count_k = 0;

        d_cart_inds[0][0] = 0; d_cart_inds[0][1] = 0;
        d_cart_inds[1][0] = 0; d_cart_inds[1][1] = 1;
        d_cart_inds[2][0] = 0; d_cart_inds[2][1] = 2;
        d_cart_inds[3][0] = 1; d_cart_inds[3][1] = 1;
        d_cart_inds[4][0] = 1; d_cart_inds[4][1] = 2;
        d_cart_inds[5][0] = 2; d_cart_inds[5][1] = 2;
        delta[0][0] = 1.0; delta[0][1] = 0.0; delta[0][2] = 0.0;
        delta[1][0] = 0.0; delta[1][1] = 1.0; delta[1][2] = 0.0;
        delta[2][0] = 0.0; delta[2][1] = 0.0; delta[2][2] = 1.0;

        if (ik < pair_inds_count_for_K_sd)
        {
            i = pair_inds_i_for_K_sd[ik];
            k = pair_inds_k_for_K_sd[ik];

            count_i = pair_counts_K_ss[i];
            count_k = pair_counts_K_dp[k];

            displ_i = pair_displs_K_ss[i];
            displ_k = pair_displs_K_dp[k];

            a_i = s_prim_info[i + s_prim_count * 0];

            r_i[0] = s_prim_info[i + s_prim_count * 2];
            r_i[1] = s_prim_info[i + s_prim_count * 3];
            r_i[2] = s_prim_info[i + s_prim_count * 4];

            a_k = d_prim_info[k / 6 + d_prim_count * 0];

            r_k[0] = d_prim_info[k / 6 + d_prim_count * 2];
            r_k[1] = d_prim_info[k / 6 + d_prim_count * 3];
            r_k[2] = d_prim_info[k / 6 + d_prim_count * 4];

            ik_factor_D = 2.0 * D_ik_for_K_sd[ik];

            c0 = d_cart_inds[k % 6][0];
            c1 = d_cart_inds[k % 6][1];

        }
    }

    __syncthreads();

    for (uint32_t m = 0; m < (count_i + TILE_DIM_Y_K - 1) / TILE_DIM_Y_K; m++)
    {
        const uint32_t j = m * TILE_DIM_Y_K + threadIdx.y;

        // sync threads before starting a new scan
        __syncthreads();

        double Q_ij, a_j, r_j[3], S_ij_00, S1, inv_S1;
        double PA_x;
        uint32_t j_prim, j_cgto;

        if ((ik < pair_inds_count_for_K_sd) && (j < count_i))
        {
            Q_ij   = Q_K_ss[displ_i + j];

            j_prim = D_inds_K_ss[displ_i + j];

            j_cgto = s_prim_aoinds[j_prim];

            a_j = s_prim_info[j_prim + s_prim_count * 0];

            r_j[0] = s_prim_info[j_prim + s_prim_count * 2];
            r_j[1] = s_prim_info[j_prim + s_prim_count * 3];
            r_j[2] = s_prim_info[j_prim + s_prim_count * 4];

            S1 = a_i + a_j;
            inv_S1 = 1.0 / S1;

            S_ij_00 = pair_data_K_ss[displ_i + j];

            PA_x = (a_j * inv_S1) * (r_j[grad_cart_ind] - r_i[grad_cart_ind]);



        }


        for (uint32_t n = 0; n < (count_k + TILE_DIM_X_K - 1) / TILE_DIM_X_K; n++)
        {
            const uint32_t l = n * TILE_DIM_X_K + threadIdx.x;

            if ((ik >= pair_inds_count_for_K_sd) || (j >= count_i) || (l >= count_k) || (fabs(Q_ij * Q_K_dp[displ_k + l] * sp_max_D) <= eri_threshold))
            {
                break;
            }

            const auto l_prim = D_inds_K_dp[displ_k + l];

            const auto l_cgto = p_prim_aoinds[(l_prim / 3) + p_prim_count * (l_prim % 3)];

            const auto a_l = p_prim_info[l_prim / 3 + p_prim_count * 0];

            const double r_l[3] = {p_prim_info[l_prim / 3 + p_prim_count * 2],
                                   p_prim_info[l_prim / 3 + p_prim_count * 3],
                                   p_prim_info[l_prim / 3 + p_prim_count * 4]};

            const auto S_kl_00 = pair_data_K_dp[displ_k + l];

            const auto d0 = l_prim % 3;

            // J. Chem. Phys. 84, 3963-3974 (1986)

            const auto S2 = a_k + a_l;

            const auto inv_S2 = 1.0 / S2;
            const auto inv_S4 = 1.0 / (S1 + S2);

            const double PQ[3] = {(a_k * r_k[0] + a_l * r_l[0]) * inv_S2 - (a_i * r_i[0] + a_j * r_j[0]) * inv_S1,
                                  (a_k * r_k[1] + a_l * r_l[1]) * inv_S2 - (a_i * r_i[1] + a_j * r_j[1]) * inv_S1,
                                  (a_k * r_k[2] + a_l * r_l[2]) * inv_S2 - (a_i * r_i[2] + a_j * r_j[2]) * inv_S1};

            const auto r2_PQ = PQ[0] * PQ[0] + PQ[1] * PQ[1] + PQ[2] * PQ[2];

            const auto rho = S1 * S2 * inv_S4;

            if (omega != 0.0) d2 = omega * omega / (rho + omega * omega);

            const auto Lambda = sqrt(4.0 * rho * d2 * MATH_CONST_INV_PI);

            double F4_t[5];

            gpu::computeBoysFunction(F4_t, rho * d2 * r2_PQ, 4, boys_func_table, boys_func_ft);

            if (omega != 0.0)
            {
                F4_t[1] *= d2;
                F4_t[2] *= d2 * d2;
                F4_t[3] *= d2 * d2 * d2;
                F4_t[4] *= d2 * d2 * d2 * d2;
            }

            const double QC_x = (a_l * inv_S2) * (r_l[grad_cart_ind] - r_k[grad_cart_ind]);

            const auto QC_0 = (a_l * inv_S2) * (r_l[c0] - r_k[c0]);
            const auto QC_1 = (a_l * inv_S2) * (r_l[c1] - r_k[c1]);
            const auto QD_0 = (-a_k * inv_S2) * (r_l[d0] - r_k[d0]);


            // lambda grad

            const double eri_ijkl = Lambda * S_ij_00 * S_kl_00 * (

                            F4_t[0] * inv_S2 * a_k * (
                                delta[c1][grad_cart_ind] * (QC_0 * QD_0)
                                + delta[c0][grad_cart_ind] * (QC_1 * QD_0)
                                + delta[c0][c1] * (QC_x * QD_0)
                                + delta[d0][grad_cart_ind] * (QC_0 * QC_1)
                                + delta[c1][d0] * (QC_0 * QC_x)
                                + delta[c0][d0] * (QC_1 * QC_x)
                            )

                            + F4_t[0] * (-0.5) * inv_S2 * (
                                (delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[0] * 2.0 * a_k * (
                                + QC_0 * QC_1 * QC_x * QD_0
                            )

                            + F4_t[0] * (-1.0) * (
                                delta[c1][grad_cart_ind] * (QC_0 * QD_0)
                                + delta[c0][grad_cart_ind] * (QC_1 * QD_0)
                            )

                            + F4_t[0] * 0.5 * inv_S2 * inv_S2 * a_k * (
                                (delta[c0][c1] * delta[d0][grad_cart_ind] + delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[1] * (-1.0) * S1 * inv_S2 * inv_S2 * inv_S4 * a_k * (
                                (delta[c0][c1] * delta[d0][grad_cart_ind] + delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[1] * (-1.0) * S1 * inv_S2 * inv_S4 * a_k * (
                                delta[d0][grad_cart_ind] * (PQ[c0] * QC_1 + PQ[c1] * QC_0 + QC_0 * QC_1)
                                + delta[c1][grad_cart_ind] * (PQ[c0] * QD_0 + PQ[d0] * QC_0 + QC_0 * QD_0)
                                + delta[c1][d0] * (PQ[c0] * QC_x + PQ[grad_cart_ind] * QC_0 + QC_0 * QC_x)
                                + delta[c0][grad_cart_ind] * (PQ[c1] * QD_0 + PQ[d0] * QC_1 + QC_1 * QD_0)
                                + delta[c0][d0] * (PQ[c1] * QC_x + PQ[grad_cart_ind] * QC_1 + QC_1 * QC_x)
                                + delta[c0][c1] * (PQ[d0] * QC_x + PQ[grad_cart_ind] * QD_0 + QC_x * QD_0)
                            )

                            + F4_t[1] * 2.0 * S1 * inv_S4 * a_k * (
                                + PQ[c0] * QC_1 * QC_x * QD_0 * (-1.0)
                                + PQ[c1] * QC_0 * QC_x * QD_0 * (-1.0)
                                + PQ[d0] * QC_0 * QC_1 * QC_x * (-1.0)
                                + PQ[grad_cart_ind] * QC_0 * QC_1 * QD_0 * (-1.0)
                            )

                            + F4_t[1] * S1 * inv_S4 * (
                                delta[c1][grad_cart_ind] * (PQ[c0] * QD_0 + PQ[d0] * QC_0)
                                + delta[c0][grad_cart_ind] * (PQ[c1] * QD_0 + PQ[d0] * QC_1)
                            )

                            + F4_t[1] * 0.5 * S1 * inv_S2 * inv_S4 * (
                                (delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[2] * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * a_k * (
                                delta[c1][grad_cart_ind] * (PQ[c0] * PQ[d0] + PQ[c0] * QD_0 + PQ[d0] * QC_0)
                                + delta[c0][grad_cart_ind] * (PQ[c1] * PQ[d0] + PQ[c1] * QD_0 + PQ[d0] * QC_1)
                                + delta[c0][c1] * (PQ[d0] * PQ[grad_cart_ind] + PQ[d0] * QC_x + PQ[grad_cart_ind] * QD_0)
                                + delta[d0][grad_cart_ind] * (PQ[c0] * PQ[c1] + PQ[c0] * QC_1 + PQ[c1] * QC_0)
                                + delta[c1][d0] * (PQ[c0] * PQ[grad_cart_ind] + PQ[c0] * QC_x + PQ[grad_cart_ind] * QC_0)
                                + delta[c0][d0] * (PQ[c1] * PQ[grad_cart_ind] + PQ[c1] * QC_x + PQ[grad_cart_ind] * QC_1)
                            )

                            + F4_t[2] * 2.0 * S1 * S1 * inv_S4 * inv_S4 * a_k * (
                                + PQ[c0] * PQ[c1] * QC_x * QD_0
                                + PQ[c0] * PQ[d0] * QC_1 * QC_x
                                + PQ[c0] * PQ[grad_cart_ind] * QC_1 * QD_0
                                + PQ[c1] * PQ[d0] * QC_0 * QC_x
                                + PQ[c1] * PQ[grad_cart_ind] * QC_0 * QD_0
                                + PQ[d0] * PQ[grad_cart_ind] * QC_0 * QC_1
                            )

                            + F4_t[2] * (-1.0) * S1 * S1 * inv_S4 * inv_S4 * (
                                delta[c1][grad_cart_ind] * (PQ[c0] * PQ[d0])
                                + delta[c0][grad_cart_ind] * (PQ[c1] * PQ[d0])
                            )

                            + F4_t[2] * 0.5 * S1 * S1 * inv_S2 * inv_S2 * inv_S4 * inv_S4 * a_k * (
                                (delta[c0][c1] * delta[d0][grad_cart_ind] + delta[c0][d0] * delta[c1][grad_cart_ind] + delta[c0][grad_cart_ind] * delta[c1][d0])
                            )

                            + F4_t[3] * (-1.0) * S1 * S1 * S1 * inv_S2 * inv_S4 * inv_S4 * inv_S4 * a_k * (
                                delta[d0][grad_cart_ind] * (PQ[c0] * PQ[c1])
                                + delta[c1][grad_cart_ind] * (PQ[c0] * PQ[d0])
                                + delta[c1][d0] * (PQ[c0] * PQ[grad_cart_ind])
                                + delta[c0][grad_cart_ind] * (PQ[c1] * PQ[d0])
                                + delta[c0][d0] * (PQ[c1] * PQ[grad_cart_ind])
                                + delta[c0][c1] * (PQ[d0] * PQ[grad_cart_ind])
                            )

                            + F4_t[3] * (-2.0) * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * a_k * (
                                PQ[c0] * PQ[c1] * PQ[d0] * QC_x
                                + PQ[c0] * PQ[c1] * PQ[grad_cart_ind] * QD_0
                                + PQ[c0] * PQ[d0] * PQ[grad_cart_ind] * QC_1
                                + PQ[c1] * PQ[d0] * PQ[grad_cart_ind] * QC_0
                            )

                            + F4_t[4] * 2.0 * S1 * S1 * S1 * S1 * inv_S4 * inv_S4 * inv_S4 * inv_S4 * a_k * (
                                PQ[c0] * PQ[c1] * PQ[d0] * PQ[grad_cart_ind]
                            )


                    );

            ERIs[threadIdx.y][threadIdx.x] -= eri_ijkl * mat_D_full_AO[l_cgto * naos + j_cgto];
        }
    }


    __syncthreads();

    if ((threadIdx.y == 0) && (threadIdx.x == 0) && (ik < pair_inds_count_for_K_sd))
    {
        double grad_k_x = 0.0;

        for (uint32_t y = 0; y < TILE_DIM_Y_K; y++)
        {
            for (uint32_t x = 0; x < TILE_DIM_X_K; x++)
            {
                grad_k_x += ERIs[y][x];
            }
        }

        atomicAdd(grad_x + prim_cart_ao_to_atom_inds[s_prim_count + p_prim_count * 3 + k], grad_k_x * ik_factor_D * 2.0 * frac_exact_exchange);
    }
}


}  // namespace gpu
